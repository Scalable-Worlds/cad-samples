// Observe the files opened by the simulation and its children, including
// data opened by VPI after elaboration. Use the kernel's resolved file path.
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ptrace.h>
#include <sys/syscall.h>
#include <sys/user.h>
#include <sys/wait.h>
#include <unistd.h>

static bool beneath(const char *path, const char *directory) {
    size_t n = strlen(directory);
    return strncmp(path, directory, n) == 0 && (path[n] == '/' || path[n] == 0);
}

static bool system_resource(const char *path) {
    if (beneath(path, "/usr/lib") || beneath(path, "/usr/lib64"))
        return strstr(path, ".so") != NULL || strstr(path, ".vpi") != NULL ||
               beneath(path, "/usr/lib/locale");
    return strcmp(path, "/etc/ld.so.cache") == 0 || strcmp(path, "/etc/localtime") == 0 ||
           beneath(path, "/usr/share/locale") || beneath(path, "/usr/share/zoneinfo") ||
           beneath(path, "/dev");
}

static int check_open(pid_t pid, const char *scratch, const char *delivery) {
    struct user_regs_struct regs;
    if (ptrace(PTRACE_GETREGS, pid, NULL, &regs) == -1) return -1;
    if ((long)regs.rax < 0) return 0;
    int flags;
    if (regs.orig_rax == SYS_openat) flags = regs.rdx;
    else if (regs.orig_rax == SYS_open) flags = regs.rsi;
    else return 0;
    if ((flags & O_ACCMODE) == O_WRONLY || (flags & O_DIRECTORY)) return 0;

    char descriptor[96], path[PATH_MAX + 1];
    snprintf(descriptor, sizeof(descriptor), "/proc/%ld/fd/%ld", (long)pid, (long)regs.rax);
    ssize_t length = readlink(descriptor, path, sizeof(path)-1);
    if (length == -1) return -1;
    path[length] = 0;
    if (beneath(path, scratch) || beneath(path, delivery) || system_resource(path)) return 0;

    // Present an unshipped file as absent. Optional diagnostics or a fallback
    // resource that the library can supply locally can still complete.
    fprintf(stderr, "External simulation resource unavailable: %s\n", path);
    regs.rax = -ENOENT;
    return ptrace(PTRACE_SETREGS, pid, NULL, &regs);
}

int main(int argc, char **argv) {
    if (argc < 5) return 2;
    pid_t root = fork();
    if (root == -1) return 2;
    if (root == 0) {
        if (ptrace(PTRACE_TRACEME, 0, NULL, NULL) < 0) _exit(2);
        raise(SIGSTOP);
        execv(argv[3], argv + 3);
        _exit(126);
    }
    int status, result = 2;
    if (waitpid(root, &status, 0) < 0 || !WIFSTOPPED(status)) return 2;
    long options = PTRACE_O_TRACESYSGOOD | PTRACE_O_TRACEFORK |
                   PTRACE_O_TRACEVFORK | PTRACE_O_TRACECLONE;
    if (ptrace(PTRACE_SETOPTIONS, root, NULL, options) < 0 ||
        ptrace(PTRACE_SYSCALL, root, NULL, NULL) < 0) return 2;
    for (;;) {
        pid_t pid = waitpid(-1, &status, __WALL);
        if (pid == -1) {
            if (errno == EINTR) continue;
            return errno == ECHILD ? result : 2;
        }
        if (WIFEXITED(status) || WIFSIGNALED(status)) {
            if (pid == root) {
                result = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
                // The calling supervisor reaps any remaining descendants.
                return result;
            }
            continue;
        }
        int sig = WSTOPSIG(status);
        if (sig == (SIGTRAP | 0x80)) {
            if (check_open(pid, argv[1], argv[2]) < 0) {
                perror("simulation file trace");
                return 2;
            }
            sig = 0;
        } else if (sig == SIGTRAP || sig == SIGSTOP) {
            sig = 0;
        }
        if (ptrace(PTRACE_SYSCALL, pid, NULL, (long)sig) < 0 && errno != ESRCH) return 2;
    }
}

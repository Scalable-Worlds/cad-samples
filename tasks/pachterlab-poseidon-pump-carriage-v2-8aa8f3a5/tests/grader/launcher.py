"""Launch an inspection command and reap the process tree it leaves behind."""

import ctypes
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile


def descendants():
    """Enumerate descendants through every thread's child list."""
    seen = set()
    pending = [os.getpid()]
    while pending:
        parent = pending.pop()
        for child_list in Path(f"/proc/{parent}/task").glob("*/children"):
            try:
                children = child_list.read_bytes().split()
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                continue
            for child in map(int, children):
                if child not in seen:
                    seen.add(child)
                    pending.append(child)
                    yield child


def reap_descendants():
    while True:
        signalled = False
        for pid in descendants():
            try:
                # Only signal an adopted child whose PID we ourselves hold
                # until waitpid. Deeper descendants are adopted in turn.
                os.waitid(os.P_PID, pid, os.WEXITED | os.WNOHANG | os.WNOWAIT)
            except ChildProcessError:
                continue
            os.kill(pid, signal.SIGKILL)
            signalled = True
        try:
            os.waitpid(-1, 0 if signalled else os.WNOHANG)
        except ChildProcessError:
            return


def main():
    libc = ctypes.CDLL(None, use_errno=True)
    pr_set_child_subreaper = 36
    if libc.prctl(pr_set_child_subreaper, 1, 0, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), "Cannot adopt orphaned workers")
    command = sys.argv[2:]
    env = dict(os.environ)
    with tempfile.TemporaryDirectory(dir="/tmp") as private:
        if command[0] == "--model-only":
            command = command[1:]
            subprocess.run([sys.executable, "-I", Path(__file__).with_name("file_access.py"),
                            private], check=True, cwd="/")
            library = Path(private) / "file_access.so"
            # Only the delivered process and its descendants see the simulated read errors.
            os.chmod(private, 0o755)
            env["LD_PRELOAD"] = str(library) + " " + env.get("LD_PRELOAD", "")
        try:
            child = subprocess.run(
                command, cwd=sys.argv[1], stdin=subprocess.DEVNULL, env=env,
                user=65534, group=65534, extra_groups=[], start_new_session=True,
            )
            return child.returncode if child.returncode >= 0 else 128 - child.returncode
        finally:
            reap_descendants()


if __name__ == "__main__":
    sys.exit(main())

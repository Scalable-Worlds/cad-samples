"""Run one unprivileged call and reap its descendants, including detached children."""
import ctypes
import os
from pathlib import Path
import signal
import subprocess
import sys


def reap_descendants():
    children = Path(f"/proc/self/task/{os.getpid()}/children")
    while True:
        try:
            pid, _ = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if pid:
            continue
        # This supervisor has one invocation. A child which detached from the
        # process group, and its own children, reparent here as parents exit.
        for child in children.read_text().split():
            try:
                os.kill(int(child), signal.SIGKILL)
            except ProcessLookupError:
                pass


def main():
    # An intervening subreaper could adopt this process before Python starts.
    parent = int(os.environ.pop("INVOKER_PARENT_PID", os.getppid()))
    stop_signals = (signal.SIGTERM, signal.SIGINT, signal.SIGHUP)

    def stop(signum, frame):
        # Finish reaping even if cancellation is signalled more than once.
        for sig in stop_signals:
            signal.signal(sig, signal.SIG_IGN)
        raise SystemExit(128 + signum)

    for sig in stop_signals:
        signal.signal(sig, stop)
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Cannot supervise descendant processes")
    if libc.prctl(1, signal.SIGTERM, 0, 0, 0) != 0:  # PR_SET_PDEATHSIG
        raise OSError(ctypes.get_errno(), "Cannot follow the calling process")
    if os.getppid() != parent:
        return 128 + signal.SIGTERM
    proc = None
    try:
        proc = subprocess.Popen(sys.argv[1:], user=65534, group=65534,
                                extra_groups=[], start_new_session=True)
        return proc.wait()
    finally:
        for sig in stop_signals:
            signal.signal(sig, signal.SIG_IGN)
        if proc is not None:
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        reap_descendants()


if __name__ == "__main__":
    sys.exit(main())

"""Reap the descendants of a single unprivileged command on Linux."""

import ctypes
import os
from pathlib import Path
import signal
import subprocess
import sys


def reap_descendants():
    # No other work is launched by this single-threaded subreaper. Even a
    # worker that created its own session is adopted when its parent exits.
    children = Path("/proc/thread-self/children")
    while True:
        try:
            finished, _ = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if finished:
            continue
        for child in map(int, children.read_text().split()):
            # The PID belongs to this supervisor until waitpid reaps it.
            try:
                os.kill(child, signal.SIGKILL)
            except ProcessLookupError:
                pass
            os.waitpid(child, 0)


def run(command):
    libc = ctypes.CDLL(None, use_errno=True)
    pr_set_child_subreaper = 36
    if libc.prctl(pr_set_child_subreaper, 1, 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    try:
        return subprocess.run(command, user=65534, group=65534, extra_groups=[],
                              start_new_session=True).returncode
    finally:
        reap_descendants()


if __name__ == "__main__":
    status = run(sys.argv[1:])
    sys.exit(128 - status if status < 0 else status)

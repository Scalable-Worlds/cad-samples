"""Run one editor and reap its complete process tree before inspecting its files."""

import ctypes
import os
from pathlib import Path
import signal
import subprocess
import sys


def run(command):
    # Orphaned grandchildren, including helpers that start a new session, become
    # children of this single-invocation supervisor. Other processes stay alone.
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Unable to supervise the editor")
    process = subprocess.Popen(
        command, stdin=subprocess.DEVNULL, start_new_session=True,
        user=65534, group=65534, extra_groups=[],
    )
    try:
        # Retain the leader PID until its process group has been signalled.
        os.waitid(os.P_PID, process.pid, os.WEXITED | os.WNOWAIT)
    finally:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
        child_list = Path(f"/proc/self/task/{os.getpid()}/children")
        while children := child_list.read_text().split():
            for value in children:
                pid = int(value)
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                os.waitpid(pid, 0)
    return process.returncode


if __name__ == "__main__":
    sys.exit(run(sys.argv[1:]))

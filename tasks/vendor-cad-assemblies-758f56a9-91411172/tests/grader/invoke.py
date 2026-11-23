"""Launch one study and reap all of its descendants, including detached workers."""

import ctypes
import os
import signal
import subprocess
import sys
from pathlib import Path


def reap_descendants():
    """As the dedicated subreaper, adopt and drain each level of the process tree."""
    children = Path(f"/proc/self/task/{os.getpid()}/children")
    while True:
        try:
            pid, _ = os.waitpid(-1, os.WNOHANG)
        except ChildProcessError:
            return
        if pid:
            continue
        for child in children.read_text().split():
            try:
                os.kill(int(child), signal.SIGKILL)
            except ProcessLookupError:
                pass


def main():
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Cannot adopt study workers")
    scratch, *command = sys.argv[1:]
    env = {**os.environ, "HOME": scratch, "TMPDIR": scratch, "MPLCONFIGDIR": scratch}
    env.pop("PYTEST_CURRENT_TEST", None)
    with open(f"{scratch}/stdout", "wb") as out, open(f"{scratch}/stderr", "wb") as err:
        try:
            return subprocess.run(
                command, cwd=scratch, env=env, stdin=subprocess.DEVNULL, stdout=out,
                stderr=err, user=65534, group=65534, extra_groups=[], start_new_session=True,
            ).returncode
        finally:
            reap_descendants()


if __name__ == "__main__":
    sys.exit(main())

"""Generate fabrication-admissible stock maps for independent proof ratings."""

import numpy as np


def challenge_layout(joint_flags, case):
    """Pack intact strips using the measured joint centers in each orientation."""
    strips = np.asarray(joint_flags).reshape(-1, 54)
    oriented = np.stack((strips, strips[:, ::-1]), axis=1).reshape(-1, 54)
    # A joint rules out five stations on the next layer: differences of 0, 1 or 2.
    blocked = oriented.copy()
    for offset in (1, 2):
        blocked[:, offset:] |= oriented[:, :-offset]
        blocked[:, :-offset] |= oriented[:, offset:]

    # The last two maps expose ordinary-stock tension and competing width-interface traces.
    rng = np.random.default_rng((691, 921, 8, 0)[case])
    available = np.ones_like(oriented[:, 0], dtype=bool)
    layout = np.empty((3, 54, 2, 12, 2), dtype=int)
    stations = np.arange(54)
    halves = np.arange(2)
    for member in range(3):
        previous = None
        for layer in range(12):
            eligible = available.copy()
            if previous is not None:
                eligible &= ~np.any(oriented & blocked[previous], axis=1)
            chosen = int(rng.choice(np.flatnonzero(eligible)))
            strip, turned = divmod(chosen, 2)
            available[2 * strip : 2 * strip + 2] = False
            rows = 54 * strip + (stations[::-1] if turned else stations)
            layout[member, :, :, layer, 0] = rows[:, None]
            layout[member, :, :, layer, 1] = halves[::-1] if turned else halves
            previous = chosen
    return layout

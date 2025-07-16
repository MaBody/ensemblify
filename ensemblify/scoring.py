import numpy as np


def descriptive_stats(scores: np.ndarray):
    stats = {}
    funcs = {
        "mean": np.mean,
        "std": lambda x: np.std(x, ddof=1),
        "min": np.min,
        "p25": lambda x: np.percentile(x, 25),
        "p50": lambda x: np.percentile(x, 50),
        "p75": lambda x: np.percentile(x, 75),
        "max": np.max,
    }
    for name, func in funcs.items():
        stats[name] = func(scores)
    return stats

"""
Evolution of resistance to mono or dual targeted therapy (Bozic I et al., eLife 2013;
2:e00747, doi:10.7554/eLife.00747).

Stochastic simulation: src/resistance/resistance_core.h (C++, all cores), built as
build/libresistance.so by `make build`. Analytic results of the paper (equation 1 and the formulas before it) are used to
validate the simulation and are shown next to it.
"""

import ctypes
import math
import os
import time
from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np

_LIB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build", "libresistance.so")
OUTCOMES = ("eradicated", "failed", "undecided")
MAX_REPLICATES = 1_000_000


class _Params(ctypes.Structure):
    _fields_ = [(name, ctypes.c_double) for name in (
        "M", "u", "n1", "n2", "n12", "b", "d", "bp", "dp", "dt", "k_exact", "fail_size", "t_max_treat")] + [
        ("dual", ctypes.c_int), ("exact_events", ctypes.c_int), ("cost", ctypes.c_double),
        ("sequential", ctypes.c_int)]


def _load():
    try:
        lib = ctypes.CDLL(_LIB_PATH)
    except OSError:  # pragma: no cover - depends on the build
        return None
    lib.resistance_run.restype = ctypes.c_int
    lib.resistance_run.argtypes = [
        ctypes.POINTER(_Params), ctypes.c_double, ctypes.c_double, ctypes.c_uint64, ctypes.c_int64,
        ctypes.c_int64] + [np.ctypeslib.ndpointer(dtype=t, flags="C_CONTIGUOUS") for t in (
            np.int32, np.int32, np.float64, np.float64, np.float64, np.float64)]
    lib.resistance_threads.restype = ctypes.c_int
    return lib


_LIB = _load()


@dataclass
class Scenario:
    M: float = 1e9
    u: float = 1e-9
    n1: float = 50
    n2: float = 50
    n12: float = 0
    b: float = 0.14
    d: float = 0.13
    bp: float = 0.14
    dp: float = 0.17
    dual: bool = True
    dt: float = 0.1
    k_exact: float = 1000.0
    fail_size: float = 1e5
    t_max_treat: float = 20000.0
    exact_events: bool = False   # small populations event by event (as in the paper) instead of tau-leaping
    cost: float = 0.0            # cost of resistance per mutation: birth (b - d)(1 - cost)^k + d (Table 2)
    sequential: bool = False     # dual: drug 1 until relapse (back at M cells), then drug 2 (Figure 5)


def simulate(sc: Scenario, n_replicates: int, seed: int, first_replicate: int = 0) -> Dict[str, np.ndarray]:
    """Replicates first_replicate .. first_replicate + n - 1 on the CPU (all cores)."""
    if _LIB is None:
        raise RuntimeError("build/libresistance.so is not built (run: make build)")
    p = _Params(sc.M, sc.u, sc.n1, sc.n2, sc.n12, sc.b, sc.d, sc.bp, sc.dp, sc.dt, sc.k_exact,
                sc.fail_size, sc.t_max_treat, 1 if sc.dual else 0, 1 if sc.exact_events else 0,
                sc.cost, 1 if (sc.dual and sc.sequential) else 0)
    n = int(n_replicates)
    out = {"outcome": np.zeros(n, np.int32), "restarts": np.zeros(n, np.int32),
           "t_detect": np.zeros(n), "t_end": np.zeros(n), "at_detect": np.zeros(4 * n), "t_switch": np.zeros(n)}
    rc = _LIB.resistance_run(ctypes.byref(p), math.exp((sc.b - sc.d) * sc.dt), math.exp((sc.bp - sc.dp) * sc.dt),
                             int(seed), int(first_replicate), n, out["outcome"], out["restarts"],
                             out["t_detect"], out["t_end"], out["at_detect"], out["t_switch"])
    if rc != 0:
        raise ValueError("invalid parameters")
    out["at_detect"] = out["at_detect"].reshape(n, 4)
    return out


def threads() -> int:
    return int(_LIB.resistance_threads()) if _LIB is not None else 0


# --------------------------------------------------------------------------- analytic (paper)
def analytic(sc: Scenario) -> Dict[str, Optional[float]]:
    """Bozic 2013: expected cells resistant to the therapy at detection (X), probability
    that such resistance is present at detection (p_res), and probability of eradication
    (equation 1). Dual therapy only (the monotherapy formulas are the n1 = n2 = 0 limit
    with n12 = n); s = 1 - d/b, s' = 1 - d'/b'."""
    M, u, b, d = sc.M, sc.u, sc.b, sc.d
    s = 1.0 - d / b
    sp = 1.0 - sc.dp / sc.bp
    if sc.dual:
        n1, n2, n12 = sc.n1, sc.n2, sc.n12
    else:   # one drug: resistance in one step via any of the n1 (+ n12) mutations
        n1, n2, n12 = 0.0, 0.0, sc.n1 + sc.n12
    mu = u / s * math.log(M * s)
    X = M * (n12 * mu + (n1 * n2 + n12 / 2.0 * (n1 + n2 - n12)) * mu * mu)
    p1 = math.exp(M * u * n12 * math.log(s) / (1.0 - s))
    p2 = math.exp(M * u * u * (2 * n1 * n2 + n12 * (n1 + n2)) * math.log(s) * math.log(M * s) / (s * (1.0 - s)))
    result = {"expected_resistant_at_detection": X, "p_resistance_at_detection": 1.0 - p1 * p2,
              "s": s, "s_treatment": sp}
    if sp >= 0:   # equation (1) assumes treatment shrinks sensitive cells (s' < 0)
        result["p_eradication"] = None
        return result
    p1_up = math.exp(-M * u * n12)
    p1_down = math.exp(M * u * n12 * s / sp)
    c = (sp - s) / (s * sp)
    p2_up = math.exp(M * u * u * c * (
        n1 * (n2 + n12) * math.log(1.0 / (s * M) + u * (n2 + n12) * c)
        + n2 * (n1 + n12) * math.log(1.0 / (s * M) + u * (n1 + n12) * c)))
    p2_down = math.exp(-M * u * u * (2 * n1 * n2 + n12 * (n1 + n2)) * s / (sp * sp))
    result["p_eradication"] = p1_up * p1_down * p2_up * p2_down
    return result


def wilson_interval(k: int, n: int, z: float = 1.96):
    """95 % confidence interval of a proportion k/n (Wilson score)."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def summarize(sc: Scenario, out: Dict[str, np.ndarray]) -> Dict[str, Any]:
    n = len(out["outcome"])
    counts = {name: int(np.sum(out["outcome"] == i)) for i, name in enumerate(OUTCOMES)}
    resistant_type = 3 if sc.dual else None
    at = out["at_detect"]
    resistant = at[:, 3] if sc.dual else at[:, 1] + at[:, 3]
    decided = counts["eradicated"] + counts["failed"]
    return {
        "replicates": n,
        "outcomes": counts,
        "p_eradication": counts["eradicated"] / n if n else None,
        "p_eradication_ci95": wilson_interval(counts["eradicated"], n),
        "p_resistance_at_detection": float(np.mean(resistant > 0)) if n else None,
        "mean_resistant_at_detection": float(np.mean(resistant)) if n else None,
        "mean_days_to_detection": float(np.mean(out["t_detect"])) if n else None,
        "decided": decided,
        "resistant_type": resistant_type,
    }

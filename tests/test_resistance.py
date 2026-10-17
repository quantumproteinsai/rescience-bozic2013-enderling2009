"""Resistance evolution (Bozic et al. eLife 2013): quick checks against the paper's analytic
results (equation 1) at its Figure 3 parameters (the full sweep is scripts/bozic_data.py),
determinism and independence from the number of threads. Run: make test"""

import os
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
import resistance_model as model  # noqa: E402

pytestmark = pytest.mark.skipif(model._LIB is None, reason="build/libresistance.so not built (make build)")


@pytest.mark.parametrize("dual,n12,M", [(True, 1, 1e8), (True, 2, 1e7), (False, 0, 1e7)])
def test_eradication_probability_matches_equation_1(dual, n12, M):
    sc = model.Scenario(M=M, n12=n12, dual=dual)
    out = model.simulate(sc, 3000, seed=5)
    s = model.summarize(sc, out)
    low, high = s["p_eradication_ci95"]
    expected = model.analytic(sc)["p_eradication"]
    # 99.9 % interval to keep the test from failing by chance (3.3 SD)
    se = np.sqrt(max(s["p_eradication"] * (1 - s["p_eradication"]), 1e-4) / 3000)
    assert abs(s["p_eradication"] - expected) < 3.3 * se + 1e-3, (s["p_eradication"], expected, (low, high))
    assert s["outcomes"]["undecided"] == 0


def test_replicates_are_reproducible_and_independent_of_batching():
    sc = model.Scenario(M=1e7, n12=1)
    a = model.simulate(sc, 600, seed=9)
    b = model.simulate(sc, 200, seed=9, first_replicate=400)
    for key in ("outcome", "restarts", "t_detect", "t_end", "at_detect"):
        np.testing.assert_array_equal(a[key][400:], b[key])
    c = model.simulate(sc, 600, seed=10)
    assert not np.array_equal(a["at_detect"], c["at_detect"])


def test_results_do_not_depend_on_the_number_of_threads():
    code = ("import numpy as np; import resistance_model as model; "
            "o = model.simulate(model.Scenario(M=1e7, n12=1), 300, seed=4); "
            "print(int(o['outcome'].sum()), repr(float(o['at_detect'].sum())), repr(float(o['t_end'].sum())))")
    runs = [subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True,
                           env=dict(os.environ, OMP_NUM_THREADS=str(t)), cwd=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python")).stdout
            for t in (1, 8)]
    assert runs[0] == runs[1]


def test_no_cross_resistance_means_near_certain_eradication_at_1e8():
    sc = model.Scenario(M=1e8, n12=0)
    s = model.summarize(sc, model.simulate(sc, 1000, seed=2))
    assert s["p_eradication"] > 0.99 and model.analytic(sc)["p_eradication"] > 0.99


def test_sequential_therapy_fails_where_simultaneous_cures():
    # Figure 5C/D (n12 = 0): sequential fails in all lesions, simultaneous cures > 99 %
    # (M = 1e8 to keep the test fast: a lesion without drug-1 resistant cells can still be
    # cured by drug 1 alone, which is rare)
    seq = model.Scenario(M=1e8, n12=0, sequential=True)
    out = model.simulate(seq, 200, seed=6)
    s = model.summarize(seq, out)
    assert s["outcomes"]["failed"] >= 195
    failed = out["outcome"] == 1
    assert (out["t_switch"][failed] > 0).all()
    sim = model.Scenario(M=1e8, n12=0)
    assert model.summarize(sim, model.simulate(sim, 200, seed=6))["p_eradication"] > 0.98


def test_cost_of_resistance_zero_is_the_default_model():
    a = model.simulate(model.Scenario(M=1e7, n12=1), 300, seed=4)
    b = model.simulate(model.Scenario(M=1e7, n12=1, cost=0.0), 300, seed=4)
    assert all(np.array_equal(a[k], b[k]) for k in a)
    c = model.simulate(model.Scenario(M=1e7, n12=1, cost=0.1), 300, seed=4)
    assert c["at_detect"][:, 1:].sum() < a["at_detect"][:, 1:].sum()

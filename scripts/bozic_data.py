"""Data for the Bozic et al. (2013) replication: probability of eradication vs lesion size
(their Figure 3 conditions) and the two-step resistance check. Run: make data-bozic"""
import csv, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import resistance_model as m

OUT = os.path.join(HERE, "..", "data")
LESIONS, SEED = 4000, 11
M_VALUES = [1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10]
CASES = [("dual", 0), ("dual", 1), ("dual", 2), ("mono", 0)]

with open(os.path.join(OUT, "bozic_eradication.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["therapy", "n12", "M", "lesions", "eradicated", "p_erad", "ci_low", "ci_high", "p_erad_eq1",
                "p_res_sim", "p_res_formula", "seconds"])
    for therapy, n12 in CASES:
        for M in M_VALUES:
            sc = m.Scenario(M=M, n12=n12, dual=therapy == "dual")
            t = time.perf_counter()
            s = m.summarize(sc, m.simulate(sc, LESIONS, seed=SEED))
            a = m.analytic(sc)
            row = [therapy, n12, M, LESIONS, s["outcomes"]["eradicated"], s["p_eradication"], *s["p_eradication_ci95"],
                   a["p_eradication"], s["p_resistance_at_detection"], a["p_resistance_at_detection"],
                   round(time.perf_counter() - t, 2)]
            w.writerow(row)
            print(therapy, n12, "%.0e" % M, "sim %.4f eq1 %.4f" % (s["p_eradication"], a["p_eradication"]), flush=True)

# Two-step resistance present at detection (dual, n12 = 0, M = 1e9): time step and threshold
with open(os.path.join(OUT, "bozic_two_step.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["method", "dt", "k_exact", "lesions", "with_resistance", "p_res_sim", "ci_low", "ci_high", "p_res_formula"])
    for exact, dt, k in ((False, 0.1, 1000.0), (False, 0.1, 1e4), (False, 0.02, 1000.0), (True, 0.1, 1000.0)):
        sc = m.Scenario(M=1e9, n12=0, dt=dt, k_exact=k, exact_events=exact)
        out = m.simulate(sc, 20000, seed=33)
        k_res = int((out["at_detect"][:, 3] > 0).sum())
        lo, hi = m.wilson_interval(k_res, 20000)
        w.writerow(["exact events" if exact else "tau-leaping", dt, k, 20000, k_res, k_res / 20000, lo, hi,
                    m.analytic(sc)["p_resistance_at_detection"]])
        print("two-step dt %.2f k %g: %.5f [%.5f, %.5f]" % (dt, k, k_res / 20000, lo, hi), flush=True)

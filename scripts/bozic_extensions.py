"""Bozic et al. (2013) beyond equation 1: cost of resistance (their Table 2), sequential vs
simultaneous dual therapy (Figure 5) and a therapy that suppresses division instead of
killing (text: 26 % -> 29 %). Simulation only, as in the paper. Run: make data-extensions"""
import csv, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import resistance_model as m

OUT = os.path.join(HERE, "..", "data")
LESIONS = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
SEED = 21

# Table 2: n1 = n2 = 50, d = b - 0.01, d' = b + 0.03, u = 1e-9; probability of FAILURE.
# Paper values (10^6 runs each) for c = 0, 1, 5, 10 %.
TABLE2 = [(0, 1e9, 0.14, (0.0, 0.0, 0.0, 0.0)), (0, 1e9, 1.0, (0.01, 0.01, 0.01, 0.0)),
          (1, 1e9, 0.14, (0.74, 0.73, 0.72, 0.70)), (1, 1e9, 1.0, (0.74, 0.74, 0.72, 0.70)),
          (0, 1e11, 0.14, (0.12, 0.11, 0.08, 0.06)), (0, 1e11, 1.0, (0.53, 0.51, 0.42, 0.32)),
          (1, 1e11, 0.14, (1.0, 1.0, 1.0, 1.0)), (1, 1e11, 1.0, (1.0, 1.0, 1.0, 1.0))]
COSTS = (0.0, 0.01, 0.05, 0.10)

with open(os.path.join(OUT, "bozic_cost.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["n12", "M", "b", "cost", "lesions", "failed", "undecided", "p_fail", "ci_low", "ci_high",
                "p_fail_paper", "seconds"])
    for n12, M, b, paper in TABLE2:
        for c, ref in zip(COSTS, paper):
            # time step: b dt = 0.014 as in the default (b = 0.14, dt = 0.1)
            sc = m.Scenario(M=M, n12=n12, b=b, d=b - 0.01, bp=b, dp=b + 0.03, cost=c, dt=0.014 / b)
            t = time.perf_counter()
            s = m.summarize(sc, m.simulate(sc, LESIONS, seed=SEED))
            k = s["outcomes"]["failed"]
            lo, hi = m.wilson_interval(k, LESIONS)
            w.writerow([n12, M, b, c, LESIONS, k, s["outcomes"]["undecided"], k / LESIONS, lo, hi, ref,
                        round(time.perf_counter() - t, 1)])
            f.flush()
            print("n12 %d M %.0e b %.2f c %.2f: fail %.3f [%.3f, %.3f] paper %.2f" % (n12, M, b, c, k / LESIONS, lo, hi, ref),
                  flush=True)

# Figure 5 (M = 1e9, b 0.14, d 0.13, d' 0.17) and division-suppressing therapy
with open(os.path.join(OUT, "bozic_sequential.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["case", "n12", "lesions", "eradicated", "p_erad", "ci_low", "ci_high", "p_11_at_start", "paper"])
    cases = [("simultaneous", 1, {}, "26.3 % eradicated"), ("sequential", 1, {"sequential": True}, "0 % (73.7 % 11 at start)"),
             ("simultaneous", 0, {}, "> 99 % eradicated"), ("sequential", 0, {"sequential": True}, "0 %"),
             ("suppress division b' 0.10 d' 0.13", 1, {"bp": 0.10, "dp": 0.13}, "29 % (26 % by killing)")]
    for name, n12, kw, paper in cases:
        sc = m.Scenario(M=1e9, n12=n12, **kw)
        out = m.simulate(sc, LESIONS, seed=SEED)
        s = m.summarize(sc, out)
        k = s["outcomes"]["eradicated"]
        lo, hi = m.wilson_interval(k, LESIONS)
        p11 = float((out["at_detect"][:, 3] > 0).mean())
        w.writerow([name, n12, LESIONS, k, k / LESIONS, lo, hi, p11, paper])
        print("%s n12 %d: erad %.4f [%.4f, %.4f], 11 at start %.3f; paper %s" % (name, n12, k / LESIONS, lo, hi, p11, paper),
              flush=True)
print("DONE")

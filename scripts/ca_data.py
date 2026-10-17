"""Data for the cancer stem cell cellular automaton replication (Enderling et al. 2009;
Poleszczuk & Enderling 2014), with the asynchronous reference implementation
(src/ca/ca_reference.cpp, built in build/). Run: make data-ca"""
import csv, os, re, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "..", "build", "ca_reference")
OUT = os.path.join(HERE, "..", "data")
PM_CELL_WIDTHS = 10.0 / 24.0   # 100 um/day, 10 um sites, 1 h steps

def run(tumours, days, seed, alpha, pm, L=1000, dump=None):
    env = dict(os.environ)
    if dump:
        env["CA_DUMP"] = dump
    out = subprocess.run([BIN, str(tumours), str(days), str(seed), str(alpha), "0.1", "10", str(L), "0", repr(pm)],
                         capture_output=True, text=True, check=True, env=env).stdout
    rows = [dict(zip(("tumour", "cells", "stem", "quiescent", "edge"), map(int, m.groups())))
            for m in re.finditer(r"TUMOUR (\d+) cells (\d+) stem (\d+) quiescent (\d+) edge (\d+)", out)]
    assert len(rows) == tumours
    return rows

with open(os.path.join(OUT, "ca_paradox.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["alpha", "tumour", "cells", "stem", "quiescent", "edge"])
    for alpha in (0.0, 0.01, 0.05, 0.1, 0.2, 0.3):
        for r in run(32, 180, 17, alpha, PM_CELL_WIDTHS):
            w.writerow([alpha, r["tumour"], r["cells"], r["stem"], r["quiescent"], r["edge"]])
        print("paradox alpha", alpha, "done", flush=True)

with open(os.path.join(OUT, "ca_migration.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pm", "tumour", "cells", "stem", "quiescent", "edge"])
    for pm in (PM_CELL_WIDTHS, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        for r in run(16, 180, 23, 0.01, pm, L=1400):
            w.writerow([pm, r["tumour"], r["cells"], r["stem"], r["quiescent"], r["edge"]])
        print("migration pm", round(pm, 4), "done", flush=True)

for alpha in (0.01, 0.2):
    run(1, 180, 5, alpha, PM_CELL_WIDTHS, L=700, dump=os.path.join(OUT, "ca_tumour_alpha%g.pgm" % alpha))
print("images done")

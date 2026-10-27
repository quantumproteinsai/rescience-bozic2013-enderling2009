"""Quick reproduction check for reviewers (make quick-check, a few minutes on a laptop).

Recomputes a representative subset of the published results from scratch and compares
them EXACTLY with the data files in data/ (every data point has its own fixed seed, so a
subset regenerates the same numbers as the full run). Nothing in data/ is modified.

  - Bozic, Figure 3: three points of the eradication curves (4,000 lesions each)
  - Bozic, Table 1: patients N20 and N22 (1,000 replicates per lesion)
  - Bozic, Table 2: one cost-of-resistance entry (10,000 lesions)
  - Bozic, Figure 5: simultaneous dual therapy with one cross-resistance mutation
  - Enderling: the first 4 of the 32 tumours at alpha = 0 and alpha = 0.2 (180 days)

The automaton results are bit-identical only with GCC/libstdc++ (see README).
"""
import csv
import os
import re
import subprocess
import sys
import time
import warnings
import zlib

import numpy as np
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import resistance_model as m   # noqa: E402

DATA = os.path.join(HERE, "..", "data")
CA_BIN = os.path.join(HERE, "..", "build", "ca_reference")
results = []


def read(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def check(label, got, expected):
    ok = got == expected
    results.append(ok)
    print("%-62s %-24s %s" % (label, "got %s" % got, "PASS" if ok else "FAIL (data: %s)" % expected), flush=True)


t0 = time.perf_counter()

# ---------------------------------------------------------------- Figure 3 (bozic_data.py: seed 11)
erad = read("bozic_eradication.csv")
for therapy, n12, M in (("dual", 1, 1e8), ("dual", 1, 1e9), ("mono", 0, 1e7)):
    row = next(r for r in erad if r["therapy"] == therapy and int(r["n12"]) == n12 and float(r["M"]) == M)
    sc = m.Scenario(M=M, n12=n12, dual=therapy == "dual")
    s = m.summarize(sc, m.simulate(sc, int(row["lesions"]), seed=11))
    check("Fig. 3: %s, n12=%d, M=%.0e: lesions eradicated" % (therapy, n12, M),
          s["outcomes"]["eradicated"], int(row["eradicated"]))

# ---------------------------------------------------------------- Table 1 (bozic_patients.py)
warnings.filterwarnings("ignore", module="openpyxl")   # harmless notes about the spreadsheet
wb = openpyxl.load_workbook(os.path.join(DATA, "bozic2013_supplementary_file2.xlsx"), data_only=True)
rows = list(wb["Summary"].iter_rows(values_only=True))
patients = [p for p in rows[1][1:] if p]
lesions = {p: [] for p in patients}
for r in rows[2:]:
    if r[0] and str(r[0]).startswith("Met."):
        for p, v in zip(patients, r[1:]):
            if isinstance(v, (int, float)) and v > 0:
                lesions[p].append(float(v))
pat = {r["patient"]: r for r in read("bozic_patients.csv")}
for p in ("N20", "N22"):
    reps = int(pat[p]["replicates"])
    for n12, col in ((1, "sim_n12_1"), (0, "sim_n12_0")):
        cured = np.ones(reps, bool)
        for i, M in enumerate(lesions[p]):
            seed = zlib.crc32(("%s-%d-%d" % (p, i, n12)).encode())
            cured &= m.simulate(m.Scenario(M=M, n12=n12, dual=True), reps, seed=seed)["outcome"] == 0
        check("Table 1: patient %s, n12=%d: failure probability" % (p, n12),
              int((~cured).sum()) / reps, float(pat[p][col]))

# ---------------------------------------------------------------- Table 2 (bozic_extensions.py: seed 21)
row = next(r for r in read("bozic_cost.csv")
           if r["n12"] == "1" and float(r["M"]) == 1e9 and float(r["b"]) == 0.14 and float(r["cost"]) == 0.05)
b = 0.14
sc = m.Scenario(M=1e9, n12=1, b=b, d=b - 0.01, bp=b, dp=b + 0.03, cost=0.05, dt=0.014 / b)
s = m.summarize(sc, m.simulate(sc, int(row["lesions"]), seed=21))
check("Table 2: n12=1, M=1e9, b=0.14, cost 5 %: lesions failed", s["outcomes"]["failed"], int(row["failed"]))

# ---------------------------------------------------------------- Figure 5 (bozic_extensions.py: seed 21)
row = next(r for r in read("bozic_sequential.csv") if r["case"] == "simultaneous" and r["n12"] == "1")
sc = m.Scenario(M=1e9, n12=1)
s = m.summarize(sc, m.simulate(sc, int(row["lesions"]), seed=21))
check("Fig. 5: simultaneous, n12=1, M=1e9: lesions eradicated", s["outcomes"]["eradicated"], int(row["eradicated"]))

# ---------------------------------------------------------------- Enderling (ca_data.py: seed 17, L 1000)
par = read("ca_paradox.csv")
for alpha in (0.0, 0.2):
    out = subprocess.run([CA_BIN, "4", "180", "17", repr(alpha), "0.1", "10", "1000", "0", repr(10.0 / 24.0)],
                         capture_output=True, text=True, check=True).stdout
    got = [(int(a), int(c), int(s)) for a, c, s in re.findall(r"TUMOUR (\d+) cells (\d+) stem (\d+)", out)]
    want = [(int(r["tumour"]), int(r["cells"]), int(r["stem"])) for r in par
            if float(r["alpha"]) == alpha and int(r["tumour"]) < 4]
    check("Enderling: alpha=%g, tumours 0-3: (cells, stem cells)" % alpha,
          [g[1:] for g in got], [w[1:] for w in want])

print("\n%d of %d checks reproduce the published data exactly (%.0f s)."
      % (sum(results), len(results), time.perf_counter() - t0))
sys.exit(0 if all(results) else 1)

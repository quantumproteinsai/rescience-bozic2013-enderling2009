"""Patient-level probability of treatment failure (Bozic et al. 2013, Table 1): a patient is
cured only if every lesion is eradicated. Lesion sizes: their Supplementary file 2
(data/bozic2013_supplementary_file2.xlsx, eLife, CC-BY 3.0). Computed (a) as the paper does,
1 - product of equation (1) over the lesions, and (b) by simulating every lesion.
Run: make data-patients"""
import csv, os, sys, time, zlib
import numpy as np
import openpyxl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
import resistance_model as m

DATA = os.path.join(HERE, "..", "data")
REPS = int(os.environ.get("BOZIC_PATIENT_REPS", "1000"))

# Table 1 of the paper: failure probability (monotherapy, dual n12 = 1, dual n12 = 0)
TABLE1 = {"N1": (1, 1, 0.283), "N2": (1, 1, 0.26), "N3": (1, 1, 0.203), "N4": (1, 1, 0.172), "N5": (1, 1, 0.128),
          "N6": (1, 1, 0.12), "N7": (1, 1, 0.112), "N8": (1, 1, 0.092), "N9": (1, 1, 0.08), "N10": (1, 1, 0.069),
          "N11": (1, 1, 0.068), "N12": (1, 1, 0.061), "N13": (1, 1, 0.052), "N14": (1, 1, 0.051), "N15": (1, 1, 0.042),
          "N16": (1, 1, 0.028), "N17": (1, 1, 0.027), "N18": (1, 1, 0.026), "N19": (1, 1, 0.023),
          "N20": (1, 0.881, 0.002), "N21": (1, 0.828, 0.002), "N22": (1, 0.677, 0.001)}

wb = openpyxl.load_workbook(os.path.join(DATA, "bozic2013_supplementary_file2.xlsx"), data_only=True)
ws = wb["Summary"]
rows = list(ws.iter_rows(values_only=True))
patients = [p for p in rows[1][1:] if p]
lesions = {p: [] for p in patients}
for row in rows[2:]:
    if not (row[0] and str(row[0]).startswith("Met.")):
        continue
    for p, v in zip(patients, row[1:]):
        if isinstance(v, (int, float)) and v > 0:
            lesions[p].append(float(v))

def failure_analytic(cells, dual, n12):
    p = 1.0
    for M in cells:
        p *= m.analytic(m.Scenario(M=M, n12=n12, dual=dual))["p_eradication"]
    return 1.0 - p

def failure_simulated(patient, cells, n12):
    cured = np.ones(REPS, bool)
    for i, M in enumerate(cells):
        seed = zlib.crc32(("%s-%d-%d" % (patient, i, n12)).encode())
        out = m.simulate(m.Scenario(M=M, n12=n12, dual=True), REPS, seed=seed)
        cured &= out["outcome"] == 0
    k = int((~cured).sum())
    return k / REPS, m.wilson_interval(k, REPS)

with open(os.path.join(DATA, "bozic_patients.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["patient", "lesions", "burden", "table1_mono", "table1_n12_1", "table1_n12_0", "eq1_mono", "eq1_n12_1",
                "eq1_n12_0", "sim_n12_1", "sim_n12_1_lo", "sim_n12_1_hi", "sim_n12_0", "sim_n12_0_lo", "sim_n12_0_hi",
                "replicates", "seconds"])
    for p in patients:
        t = time.perf_counter()
        cells = lesions[p]
        eq = (failure_analytic(cells, False, 0), failure_analytic(cells, True, 1), failure_analytic(cells, True, 0))
        s1, ci1 = failure_simulated(p, cells, 1)
        s0, ci0 = failure_simulated(p, cells, 0)
        w.writerow([p, len(cells), sum(cells), *TABLE1[p], *eq, s1, *ci1, s0, *ci0, REPS, round(time.perf_counter() - t, 1)])
        f.flush()
        print("%s: %d lesions, burden %.2g | Table 1 %s | eq.1 %.3f %.3f %.3f | sim n12=1 %.3f, n12=0 %.3f"
              % (p, len(cells), sum(cells), TABLE1[p], *eq, s1, s0), flush=True)

"""Figures, tables and the numbers cited in main.tex (results.tex), from data/.
Run: make figures"""
import csv
import math
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
FIG = os.path.join(HERE, "..", "figures")
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
COLOURS = {("dual", 0): "#1f77b4", ("dual", 1): "#d62728", ("dual", 2): "#ff7f0e", ("mono", 0): "#555555"}
LABELS = {("dual", 0): "dual, $n_{12}=0$", ("dual", 1): "dual, $n_{12}=1$", ("dual", 2): "dual, $n_{12}=2$",
          ("mono", 0): "monotherapy"}


def read(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def fmt_int(x):
    return "{:,}".format(int(round(x))).replace(",", "{,}")


macros = {}

# ----------------------------------------------------------------------------- Bozic
rows = read("bozic_eradication.csv")
groups = defaultdict(list)
for r in rows:
    groups[(r["therapy"], int(r["n12"]))].append(r)
inside, points, max_diff = 0, 0, 0.0
fig, ax = plt.subplots(figsize=(5.2, 3.4))
for key, rs in groups.items():
    M = np.array([float(r["M"]) for r in rs])
    p = np.array([float(r["p_erad"]) for r in rs])
    lo = np.array([float(r["ci_low"]) for r in rs])
    hi = np.array([float(r["ci_high"]) for r in rs])
    eq = np.array([float(r["p_erad_eq1"]) for r in rs])
    for a, l, h in zip(eq, lo, hi):
        points += 1
        inside += l <= a <= h
    max_diff = max(max_diff, float(np.max(np.abs(p - eq))))
    Mf = np.logspace(np.log10(M.min()), np.log10(M.max()), 200)
    ax.errorbar(M, p, yerr=[np.maximum(p - lo, 0), np.maximum(hi - p, 0)], fmt="o", ms=3.5, color=COLOURS[key], capsize=2, lw=1)
    ax.plot(M, eq, "-", color=COLOURS[key], lw=1.2, label=LABELS[key])
ax.set_xscale("log")
ax.set_xlabel("cells at the start of treatment, $M$")
ax.set_ylabel("probability of eradication")
ax.set_ylim(-0.03, 1.03)
ax.legend(frameon=False, fontsize=8)
fig.savefig(os.path.join(FIG, "bozic_eradication.pdf"), metadata={"CreationDate": None})
plt.close(fig)

at = {(r["therapy"], int(r["n12"]), float(r["M"])): r for r in rows}
macros.update({
    "BozicInside": str(inside), "BozicPoints": str(points), "BozicMaxDiff": "{:.3f}".format(max_diff),
    "BozicLesions": fmt_int(float(rows[0]["lesions"])),
    "BozicNoCross": "{:.3f}".format(float(at[("dual", 0, 1e9)]["p_erad"])),
    "BozicOneCross": "{:.2f}".format(float(at[("dual", 1, 1e9)]["p_erad"])),
})

two = read("bozic_two_step.csv")
base = two[0]
macros["BozicTwoStepSim"] = "{:.4f} [{:.4f}, {:.4f}]".format(float(base["p_res_sim"]), float(base["ci_low"]),
                                                              float(base["ci_high"]))
macros["BozicTwoStepFormula"] = "{:.4f}".format(float(base["p_res_formula"]))
devs = [1 - float(r["p_res_sim"]) / float(r["p_res_formula"]) for r in two]
macros["BozicTwoStepDev"] = "{:.0f}--{:.0f}\\,\\%".format(100 * min(devs), 100 * max(devs))
with open(os.path.join(FIG, "two_step_table.tex"), "w") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\toprule\nsmall populations & $\\Delta t$ (day) & threshold $K$ & simulated (95\\,\\% CI) & formula $1-p_1p_2$ \\\\\n\\midrule\n")
    for r in two:
        f.write("{} & {} & {} & {:.4f} [{:.4f}, {:.4f}] & {:.4f} \\\\\n".format(
            r["method"], r["dt"], fmt_int(float(r["k_exact"])), float(r["p_res_sim"]), float(r["ci_low"]), float(r["ci_high"]),
            float(r["p_res_formula"])))
    f.write("\\bottomrule\n\\end{tabular}\n")

# ----------------------------------------------------------------------------- patients (Table 1)
pts = read("bozic_patients.csv")
match = [r for r in pts if max(abs(float(r["eq1_" + k]) - float(r["table1_" + k])) for k in ("mono", "n12_1", "n12_0")) <= 0.0015]
with open(os.path.join(FIG, "patients_table.tex"), "w") as f:
    f.write("\\begin{tabular}{lrr|rrr|rrr}\n\\toprule\n & & & \\multicolumn{3}{c|}{dual, $n_{12}=1$} & \\multicolumn{3}{c}{dual, $n_{12}=0$} \\\\\n"
            "patient & lesions & cells & Table 1 & eq.~(1) & simulated & Table 1 & eq.~(1) & simulated \\\\\n\\midrule\n")
    for r in pts:
        f.write("{} & {} & \\num{{{:.2g}}} & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} & {:.3f} \\\\\n".format(
            r["patient"], r["lesions"], float(r["burden"]), float(r["table1_n12_1"]), float(r["eq1_n12_1"]),
            float(r["sim_n12_1"]), float(r["table1_n12_0"]), float(r["eq1_n12_0"]), float(r["sim_n12_0"])))
    f.write("\\bottomrule\n\\end{tabular}\n")
n20 = next(r for r in pts if r["patient"] == "N20")
inside_sim = sum(1 for r in pts for k in ("n12_1", "n12_0")
                 if float(r["sim_%s_lo" % k]) - 1e-12 <= float(r["eq1_" + k]) <= float(r["sim_%s_hi" % k]) + 1e-12)
macros.update({
    "BozicPatients": str(len(pts)), "BozicPatientsMatch": str(len(match)),
    "BozicPatientReps": fmt_int(float(pts[0]["replicates"])),
    "BozicPatientsSimInside": str(inside_sim), "BozicPatientsSimTotal": str(2 * len(pts)),
    "BozicNTwentyBurden": "\\num{%.3g}" % float(n20["burden"]),
    "BozicNTwentyEq": "%.3f" % float(n20["eq1_n12_1"]),
})

# ----------------------------------------------------------------------------- cost of resistance, sequential therapy
cost = read("bozic_cost.csv")
consistent = 0
cells = {}
for r in cost:
    p, lo, hi, ref = (float(r[k]) for k in ("p_fail", "ci_low", "ci_high", "p_fail_paper"))
    # the paper's values are rounded to two decimals
    consistent += lo <= ref + 0.005 and hi >= ref - 0.005
    cells[(r["n12"], float(r["M"]), float(r["b"]), float(r["cost"]))] = (p, ref)
macros["BozicCostConsistent"] = "{} of {}".format(consistent, len(cost))
macros["BozicCostLesions"] = fmt_int(int(cost[0]["lesions"]))
conds = []
for k in cells:
    if k[:3] not in conds:
        conds.append(k[:3])
costs = sorted({k[3] for k in cells})
with open(os.path.join(FIG, "cost_table.tex"), "w") as f:
    f.write("\\begin{tabular}{rrr" + "c" * len(costs) + "}\n\\toprule\n$n_{12}$ & $M$ & $b$ & "
            + " & ".join("$c={:g}\\,\\%$".format(100 * c) for c in costs) + " \\\\\n\\midrule\n")
    for n12, M, b in conds:
        f.write("{} & $10^{{{:.0f}}}$ & {:g} & ".format(n12, math.log10(M), b)
                + " & ".join("{:.3f} ({:.2f})".format(*cells[(n12, M, b, c)]) for c in costs) + " \\\\\n")
    f.write("\\bottomrule\n\\end{tabular}\n")
seq = read("bozic_sequential.csv")
names = {"simultaneous": "simultaneous", "sequential": "sequential"}
with open(os.path.join(FIG, "sequential_table.tex"), "w") as f:
    f.write("\\begin{tabular}{lrll}\n\\toprule\nTherapy & $n_{12}$ & eradicated (95\\,\\% CI) & paper \\\\\n\\midrule\n")
    for r in seq:
        label = "division suppressed ($b'=0.10$, $d'=0.13$)" if r["case"].startswith("suppress") else r["case"]
        hi = float(r["ci_high"])
        f.write("{} & {} & {:.3f} [{:.3f}, {} ] & {} \\\\\n".format(
            label, r["n12"], float(r["p_erad"]), float(r["ci_low"]), "{:.3f}".format(hi) if hi >= 0.001 else "{:.4f}".format(hi),
            r["paper"].replace("%", "\\%").replace(">", "$>$")))
    f.write("\\bottomrule\n\\end{tabular}\n")
seqd = {(r["case"], r["n12"]): r for r in seq}
macros["BozicSimulOneCross"] = "{:.1f}\\,\\%".format(100 * float(seqd[("simultaneous", "1")]["p_erad"]))
macros["BozicSuppress"] = "{:.1f}\\,\\%".format(100 * float(seqd[("suppress division b' 0.10 d' 0.13", "1")]["p_erad"]))
macros["BozicSeqLesions"] = fmt_int(int(seq[0]["lesions"]))

# ----------------------------------------------------------------------------- automaton
def summary(rows, key):
    out = defaultdict(lambda: defaultdict(list))
    for r in rows:
        k = float(r[key])
        for col in ("cells", "stem", "quiescent", "edge"):
            out[k][col].append(float(r[col]))
    return {k: {c: np.array(v) for c, v in cols.items()} for k, cols in out.items()}


par = summary(read("ca_paradox.csv"), "alpha")
alphas = sorted(par)
n_t = len(par[alphas[0]]["cells"])
mean = lambda a: float(a.mean())
se = lambda a: float(a.std(ddof=1) / math.sqrt(len(a)))

with open(os.path.join(FIG, "paradox_table.tex"), "w") as f:
    f.write("\\begin{tabular}{rrrr}\n\\toprule\n$\\alpha$ & cells & stem cells & without free neighbour \\\\\n\\midrule\n")
    for a in alphas:
        d = par[a]
        q = 100 * d["quiescent"] / d["cells"]
        f.write("{:g} & {} $\\pm$ {} & {:.0f} $\\pm$ {:.0f} & {:.1f}\\,\\% \\\\\n".format(
            a, fmt_int(mean(d["cells"])), fmt_int(se(d["cells"])), mean(d["stem"]), se(d["stem"]), q.mean()))
    f.write("\\bottomrule\n\\end{tabular}\n")


def pgm(path):
    with open(path, "rb") as fh:
        raw = fh.read()
    parts = raw.split(b"\n", 3)
    w, h = map(int, parts[1].split())
    return np.frombuffer(parts[3], np.uint8).reshape(h, w)


fig = plt.figure(figsize=(7.2, 2.6))
for i, a in enumerate((0.01, 0.2)):
    ax = fig.add_axes([0.0 + 0.27 * i, 0.08, 0.26, 0.84])
    img = pgm(os.path.join(DATA, "ca_tumour_alpha{:g}.pgm".format(a)))
    rgb = np.ones(img.shape + (3,))
    rgb[img == 128] = (0.62, 0.10, 0.12)
    rgb[img == 255] = (0.98, 0.80, 0.08)
    ax.imshow(rgb, interpolation="nearest")
    ax.set_axis_off()
    ax.set_title("({}) $\\alpha = {:g}$, {} cells".format("ab"[i], a, fmt_int((img > 0).sum()).replace("{,}", ",")),
                 fontsize=8)
ax = fig.add_axes([0.62, 0.18, 0.33, 0.72])
cells = [mean(par[a]["cells"]) for a in alphas]
ax.errorbar(alphas, cells, yerr=[se(par[a]["cells"]) for a in alphas], fmt="o-", color="#b22222", ms=3.5,
            capsize=2, label="all cells")
ax.set_xlabel("death probability at division, $\\alpha$")
ax.set_ylabel("cells after 180 days")
ax2 = ax.twinx()
ax2.errorbar(alphas, [mean(par[a]["stem"]) for a in alphas], yerr=[se(par[a]["stem"]) for a in alphas],
             fmt="s--", color="#c79a00", ms=3.5, capsize=2, label="stem cells")
ax2.set_ylabel("stem cells")
ax2.spines["right"].set_visible(True)
ax.set_title("(c)", fontsize=8)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=7, loc="upper left")
fig.savefig(os.path.join(FIG, "ca_paradox.pdf"), metadata={"CreationDate": None})
plt.close(fig)

mig = summary(read("ca_migration.csv"), "pm")
pms = sorted(mig)
mcells = np.array([mean(mig[p]["cells"]) for p in pms])
fig, ax = plt.subplots(figsize=(4.0, 2.8))
ax.errorbar(pms, mcells, yerr=[se(mig[p]["cells"]) for p in pms], fmt="o-", color="#333333", ms=3.5, capsize=2)
ax.axhline(140000, ls="--", color="#d62728", lw=1)
ax.text(pms[0], 150000, "reported: ~140,000", color="#d62728", fontsize=7)
ax.set_xlabel("migration probability per step, $p_m$")
ax.set_ylabel("cells after 180 days")
fig.savefig(os.path.join(FIG, "ca_migration.pdf"), metadata={"CreationDate": None})
plt.close(fig)

# p_m at which the mean size crosses 140,000 (linear interpolation)
match = None
for (p0, c0), (p1, c1) in zip(zip(pms, mcells), zip(pms[1:], mcells[1:])):
    if c0 <= 140000 <= c1:
        match = p0 + (140000 - c0) * (p1 - p0) / (c1 - c0)
        break
q01 = par[0.01]
macros.update({
    "CATumours": str(n_t),
    "CAParadoxLow": fmt_int(mean(par[0.0]["cells"])),
    "CAParadoxHigh": fmt_int(mean(par[0.2]["cells"])),
    "CAStemRatio": "{:.0f}".format(mean(par[0.2]["stem"]) / mean(par[0.0]["stem"])),
    "CAQuiescent": "{:.0f}\\,\\%".format(float((100 * q01["quiescent"] / q01["cells"]).mean())),
    "CASizeCellWidths": fmt_int(mean(q01["cells"])),
    "CASizeLiteral": fmt_int(mean(mig[1.0]["cells"])),
    "CAPmMatch": "{:.2f}".format(match) if match is not None else "n/a",
    "CAEdgeTouches": str(int(sum(par[a]["edge"].sum() for a in alphas))),
})

with open(os.path.join(HERE, "..", "paper", "results.tex"), "w") as f:
    f.write("% Written by scripts/make_figures.py from data/ - do not edit by hand\n")
    for k, v in sorted(macros.items()):
        f.write("\\newcommand{\\%s}{%s}\n" % (k, v))
print("\n".join("%s = %s" % kv for kv in sorted(macros.items())))

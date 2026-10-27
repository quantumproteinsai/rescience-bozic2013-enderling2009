# Independent replication of two stochastic models of cancer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23437598.svg)](https://doi.org/10.5281/zenodo.23437598)

Code, data and manuscript for the preprint

> *Independent replication of two stochastic models of cancer: resistance to combination
> therapy (Bozic et al., 2013) and cancer stem cell driven tumour growth (Enderling et al., 2009).*
> Bertrand Mercier des Rochettes, 2026. [Preprint DOI, to be added]
>
> Code and data archived at Zenodo: [doi:10.5281/zenodo.23437598](https://doi.org/10.5281/zenodo.23437598) (all versions)

Both models were re-implemented from the published descriptions only:

- **Resistance to targeted combination therapy:** the multitype branching process of
  Bozic et al., *eLife* 2013;2:e00747. Its probability of eradication (equation 1) is
  compared with simulations, together with the patient-level failure probabilities
  (Table 1).
- **Cancer stem cell cellular automaton:** Enderling et al., *Cancer Res* 2009;69:8814,
  as specified by Poleszczuk & Enderling, *Applied Mathematics* 2014;5:144. Covers the
  death-rate "paradox", the fraction of quiescent cells, and tumour size against the
  migration rule.

## Contents

| Path | What it is |
|---|---|
| `src/resistance/` | Resistance simulation (C++, OpenMP): hybrid tau-leaping, plus an exact event-by-event option |
| `src/ca/ca_reference.cpp` | Cellular automaton, asynchronous updates (C++, OpenMP across tumours) |
| `python/resistance_model.py` | Python interface to the resistance simulation; the paper's analytic formulas |
| `scripts/` | Generate the data (`bozic_data.py`, `bozic_patients.py`, `ca_data.py`); figures, tables and the numbers quoted in the text (`make_figures.py`) |
| `data/` | All simulation output used in the paper, with run logs. Also `bozic2013_supplementary_file2.xlsx`, the lesion sizes from Bozic et al. (eLife, CC BY 4.0), unmodified |
| `figures/` | Figures (PDF) and tables (LaTeX) written by `make_figures.py` |
| `rescience/` | The same work as two ReScience C articles: `bozic2013/` and `enderling2009/` (template in `rescience/template/`, GPL-3+ and font licences) |
| `paper/` | Manuscript. `results.tex` holds every simulated number as a macro, written by `make_figures.py` |
| `tests/` | Quick checks of the resistance simulation |

## Reproducing

Requirements: Linux or macOS, `g++` with OpenMP, Python ≥ 3.10, and a LaTeX distribution
with `latexmk` and these packages: `lmodern`, `siunitx`, `xurl`, `booktabs`, `caption`,
`placeins`, `hyperref` (on Debian/Ubuntu, `texlive-latex-extra` and `lmodern`).

```bash
make paper     # manuscript from the data shipped in data/      (minutes)
make test      # quick checks of the resistance simulation       (about 1 min)
make quick-ca  # quick check of the automaton's death-rate paradox (about 1 min on 2 cores)
make all       # rebuild programs, regenerate ALL data, figures and manuscript
```

`make quick-ca` runs 4 tumours at α = 0 and α = 0.2 and writes only to `data/quick/`; the
data used in the paper are not touched. It checks the central qualitative result (a higher
death rate gives larger tumours), not the numbers of the paper, which come from the full
ensembles (32 tumours per death rate, `make data-ca`).

`make all` takes several hours on an 8-thread laptop CPU, almost all of it in the
automaton ensembles (`make data-ca`). Each data set can be regenerated on its own with
`make data-bozic`, `make data-patients` or `make data-ca`. The Python packages are
installed into `.venv/` on first use.

### Exact reproducibility

Every simulation is seeded, and the seeds are in the scripts.

- **Resistance model:** every lesion has its own random stream (xoshiro256**), so the
  results do not depend on the number of threads (checked in `tests/`).
- **Automaton:** uses `std::mt19937_64` with `std::shuffle` and
  `std::uniform_real_distribution`. Their output depends on the C++ standard library,
  so the published numbers are reproduced bit for bit with GCC/libstdc++. Other
  compilers give statistically equivalent, not identical, results.

The data in this repository were produced with GCC 15.2 on Linux (Intel i5-8279U).

## Licence

Code: MIT (see `LICENSE`). Manuscript text, figures and generated data: CC BY 4.0. The
supplementary file of Bozic et al. is redistributed under its original CC BY 4.0
licence.

## Citation

See `CITATION.cff`. Please cite the original papers of the replicated models as well.

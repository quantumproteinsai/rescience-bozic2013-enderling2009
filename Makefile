# Replication of Bozic et al. (2013) and Enderling et al. (2009) / Poleszczuk & Enderling (2014).
#   make paper     compile the manuscript from the data in data/ (minutes)
#   make all       rebuild the programs, regenerate ALL data (hours), figures and manuscript
#   make test      quick checks of the resistance simulation
#   make quick-check  reproduce a subset of the published data exactly (a few minutes)
#   make rescience the two ReScience C articles (rescience/*/article.pdf; needs XeLaTeX)
#   make quick-ca  quick check of the automaton's death-rate paradox (4 tumours; writes data/quick/ only)
# Requirements: g++ with OpenMP, Python 3.10+, a LaTeX distribution with latexmk
# (LaTeX packages: lmodern, siunitx, xurl, booktabs, caption, placeins, hyperref).

CXX     ?= g++
PYTHON  := .venv/bin/python
# -ffp-contract=off, no -ffast-math: reproducible floating point
CXXFLAGS = -O2 -fopenmp -ffp-contract=off -Wall -Wextra

.PHONY: all build venv rescience quick-check data data-bozic data-patients data-ca quick-ca figures paper test clean

all: build data figures paper

build: build/libresistance.so build/ca_reference

build/libresistance.so: src/resistance/resistance_cpu.cpp src/resistance/resistance_core.h
	mkdir -p build
	$(CXX) $(CXXFLAGS) -fPIC -shared -o $@ src/resistance/resistance_cpu.cpp

build/ca_reference: src/ca/ca_reference.cpp
	mkdir -p build
	$(CXX) -O2 -fopenmp -o $@ $<

venv: .venv/.done
.venv/.done: requirements.txt
	python3 -m venv .venv
	$(PYTHON) -m pip install -q -r requirements.txt
	touch $@

data: data-bozic data-patients data-extensions data-ca

data-bozic: build venv
	$(PYTHON) scripts/bozic_data.py | tee data/bozic_log.txt

data-patients: build venv
	$(PYTHON) scripts/bozic_patients.py | tee data/bozic_patients_log.txt

data-extensions: build venv
	$(PYTHON) scripts/bozic_extensions.py 10000 | tee data/bozic_extensions_log.txt

data-ca: build venv
	$(PYTHON) scripts/ca_data.py | tee data/ca_log.txt

quick-ca: build venv
	CA_QUICK=1 $(PYTHON) scripts/ca_data.py

figures: venv
	$(PYTHON) scripts/make_figures.py

paper: figures
	cd paper && latexmk -pdf -interaction=nonstopmode -outdir=build main.tex
	@echo "Manuscript: paper/build/main.pdf"

quick-check: build venv
	$(PYTHON) scripts/quick_check.py

rescience: figures
	cd rescience/bozic2013 && $(MAKE) article.pdf
	cd rescience/enderling2009 && $(MAKE) article.pdf

test: build venv
	$(PYTHON) -m pytest -q tests

clean:
	rm -rf build paper/build

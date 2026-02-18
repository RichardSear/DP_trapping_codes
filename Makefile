# make -n DEST_DIR=..  overrides the setting below
# Warren and Sear 2025/2026

DEST_DIR = .

PYTHON3 = /usr/bin/env python3

MAIN_FIGS = fig1 fig2 fig3 fig4 fig5
SUPP_FIGS = figS1 figS2a figS2b figS3 figS4ab figS4cd

default: main_figs supp_figs

main_figs: $(MAIN_FIGS)
supp_figs: $(SUPP_FIGS)


vardp100k.ods: raw_analyse.py data/vardp100k.dat.gz
	$(PYTHON3) $^ -o $@

varrc100k.ods: raw_analyse.py data/varrc100k.dat.gz
	$(PYTHON3) $^ --col rc -o $@

poredp10k.ods: raw_analyse.py data/poredp10k.dat.gz
	$(PYTHON3) $^ -o $@

poredp100k.ods: raw_analyse.py data/poredp100k.dat.gz
	$(PYTHON3) $^ -o $@

porerc10k.ods: raw_analyse.py data/porerc10k.dat.gz
	$(PYTHON3) $^ --col rc -o $@

poremsteps.ods: raw_analyse.py data/poremsteps.dat.gz
	$(PYTHON3) $^ --col code -o $@


fig1: fig1.py
	$(PYTHON3) $^ -o $(DEST_DIR)/fig1.pdf

fig2: fig2.py
	$(PYTHON3) $^ -o $(DEST_DIR)/fig2.pdf

fig3: fig3.py vardp100k.ods
	$(PYTHON3) $^ -j -o $(DEST_DIR)/fig3.pdf

fig4: fig4.py
	$(PYTHON3) $^ -o $(DEST_DIR)/fig4.pdf

fig5: fig5.py poredp100k.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/fig5.pdf


figS1: figS1.py data/vardp100k.dat.gz
	$(PYTHON3) $^ -j --dpi 300 -o $(DEST_DIR)/figS1.png

figS2a: figS2a.py varrc100k.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/figS2a.pdf

figS2b: figS2b.py porerc10k.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/figS2b.pdf


figS3: figS3.py poremsteps.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/figS3.pdf

figS4ab: figS4ab.py vardp100k.ods
	$(PYTHON3) $^ -j -o $(DEST_DIR)/figS4ab.pdf

figS4cd: figS4cd.py poredp10k.ods
	$(PYTHON3) $^ -j -o $(DEST_DIR)/figS4cd.pdf

clean:
	rm -f *~

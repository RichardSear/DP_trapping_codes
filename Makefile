# make -n DEST_DIR=..  overrides the setting below
# Warren and Sear 2025/2026

DEST_DIR = .
DATA_DIR = data

PYTHON = /usr/bin/env python3

MAIN_FIGS = fig1 fig2 fig3 fig4 fig5
SUPP_FIGS = figS1 figS2 figS3 figS4 figS5 figS6

default: main_figs supp_figs

main_figs: $(MAIN_FIGS)
supp_figs: $(SUPP_FIGS)


vardp100k.ods: raw_analyse.py $(DATA_DIR)/vardp100k.dat.gz
	$(PYTHON) $^ -o $@

varrc100k.ods: raw_analyse.py $(DATA_DIR)/varrc100k.dat.gz
	$(PYTHON) $^ --col rc -o $@

poredp10k.ods: raw_analyse.py $(DATA_DIR)/poredp10k.dat.gz
	$(PYTHON) $^ -o $@

poredp100k.ods: raw_analyse.py $(DATA_DIR)/poredp100k.dat.gz
	$(PYTHON) $^ -o $@

porerc10k.ods: raw_analyse.py $(DATA_DIR)/porerc10k.dat.gz
	$(PYTHON) $^ --col rc -o $@

poremsteps.ods: raw_analyse.py $(DATA_DIR)/poremsteps.dat.gz
	$(PYTHON) $^ --col code -o $@


fig1: fig1.py
	$(PYTHON) $^ -o $(DEST_DIR)/fig1.pdf

fig2: fig2.py
	$(PYTHON) $^ -o $(DEST_DIR)/fig2.pdf

fig3: fig3.py vardp100k.ods
	$(PYTHON) $^ -o $(DEST_DIR)/fig3.pdf

fig4: fig4.py
	$(PYTHON) $^ -o $(DEST_DIR)/fig4.pdf

fig5: fig5.py poredp100k.ods
	$(PYTHON) $^ -o $(DEST_DIR)/fig5.pdf


figS1: figS1.py $(DATA_DIR)/vardp100k.dat.gz
	$(PYTHON) $^ -j --dpi 300 -o $(DEST_DIR)/figS1.png

figS2: figS2.py varrc100k.ods
	$(PYTHON) $^ -o $(DEST_DIR)/figS2.pdf

figS3: figS3.py poremsteps.ods
	$(PYTHON) $^ -o $(DEST_DIR)/figS3.pdf

figS4: figS4.py porerc10k.ods 
	$(PYTHON) $^ -o $(DEST_DIR)/figS4.pdf

figS5: figS5.py vardp100k.ods
	$(PYTHON) $^ -j -o $(DEST_DIR)/figS5.pdf

figS6: figS6.py poredp10k.ods
	$(PYTHON) $^ -j -o $(DEST_DIR)/figS6.pdf

clean:
	rm -f *~

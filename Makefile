# make -n DEST_DIR=..  overrides the setting below

DEST_DIR = .
DATA_DIR = data
ODS_DIR = .

PYTHON = /usr/bin/env python3

MAIN_FIGS = fig1 fig2 fig3 fig4
SUPP_FIGS = figS1 figS2 figS3 figS4 figS5 figS6

PIPETTE_ODS = vardp100k.ods varrc100k.ods
PORE_ODS = poredp100k.ods porerc10k.ods poremsteps.ods

default: spreadsheets main_figs supp_figs

spreadsheets: $(PIPETTE_ODS) $(PORE_ODS)

main_figs: $(MAIN_FIGS)
supp_figs: $(SUPP_FIGS)

vardp100k.ods: raw_analyse.py $(DATA_DIR)/vardp100k.dat.gz
	$(PYTHON) $^ -o $(ODS_DIR)/$@

varrc100k.ods: raw_analyse.py $(DATA_DIR)/varrc100k.dat.gz
	$(PYTHON) $^ --col rc -o $(ODS_DIR)/$@

poredp100k.ods: raw_analyse.py $(DATA_DIR)/poredp100k.dat.gz
	$(PYTHON) $^ -o $(ODS_DIR)/$@

porerc10k.ods: raw_analyse.py $(DATA_DIR)/porerc10k.dat.gz
	$(PYTHON) $^ --col rc -o $(ODS_DIR)/$@

poremsteps.ods: raw_analyse.py $(DATA_DIR)/poremsteps.dat.gz
	$(PYTHON) $^ --col code -o $(ODS_DIR)/$@


fig1: fig1.py
	$(PYTHON) $^ -o $(DEST_DIR)/fig1.pdf

fig2: fig2.py
	$(PYTHON) $^ -o $(DEST_DIR)/fig2.pdf

fig3: fig3.py $(ODS_DIR)/vardp100k.ods
	$(PYTHON) $^ -o $(DEST_DIR)/fig3.pdf

fig4: fig4.py
	$(PYTHON) $^ -o $(DEST_DIR)/fig4.pdf

fig5: fig5.py $(ODS_DIR)/poredp100k.ods
	$(PYTHON) $^ -o $(DEST_DIR)/fig5.pdf


figS1: figS1.py $(DATA_DIR)/vardp100k.dat.gz
	$(PYTHON) $^ -j --dpi 300 -o $(DEST_DIR)/figS1.png

figS2: figS2.py $(ODS_DIR)/varrc100k.ods
	$(PYTHON) $^ -o $(DEST_DIR)/figS2.pdf

figS3: figS3.py $(ODS_DIR)/poremsteps.ods
	$(PYTHON) $^ -o $(DEST_DIR)/figS3.pdf

figS4: figS4.py $(ODS_DIR)/porerc10k.ods 
	$(PYTHON) $^ -o $(DEST_DIR)/figS4.pdf

figS5: figS5.py $(ODS_DIR)/vardp100k.ods
	$(PYTHON) $^ -j -o $(DEST_DIR)/figS5.pdf

figS6: figS6.py $(ODS_DIR)/poredp10k.ods
	$(PYTHON) $^ -j -o $(DEST_DIR)/figS6.pdf

clean:
	rm -f *~

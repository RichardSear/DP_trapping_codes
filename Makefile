# make -n DEST_DIR=..  overrides the setting below

DATA_DIR = data
DEST_DIR = .

MAIN_FIGS = fig1 fig2 fig3 fig4
SUPP_FIGS = figS1 figS2 figS3 figS4 figS5 figS6

ODS_PIPETTE = vardp100k.ods varrc100k.ods
ODS_PORE = poredp100k.ods porerc10k.ods poremsteps.ods

default: spreadsheets figs

spreadsheets: $(ODS_PIPETTE) $(ODS_PORE)

figs: $(MAIN_FIGS) $(SUPP_FIGS)


vardp100k.ods: raw_analyse.py
	./raw_analyse.py $(DATA_DIR)/$(subst .ods,.dat.gz,$@) -o $@

varrc100k.ods: raw_analyse.py
	./raw_analyse.py --col rc $(DATA_DIR)/$(subst .ods,.dat.gz,$@) -o $@

poredp100k.ods: raw_analyse.py
	./raw_analyse.py $(DATA_DIR)/$(subst .ods,.dat.gz,$@) -o $@

porerc10k.ods: raw_analyse.py
	./raw_analyse.py --col rc $(DATA_DIR)/$(subst .ods,.dat.gz,$@) -o $@

poremsteps.ods: raw_analyse.py
	./raw_analyse.py --col  $(DATA_DIR)/$(subst .ods,.dat.gz,$@) -o $@


fig1: fig1_make.py
	./fig1_make.py -o $(DEST_DIR)/fig1.pdf

fig2: fig2_make.py
	./fig2_make.py -o $(DEST_DIR)/fig2.pdf

fig3: fig3_make.py vardp100k.ods
	./fig3_make.py vardp100k.ods -o $(DEST_DIR)/fig3.pdf

fig4: fig4_make.py
	./fig4_make.py -o $(DEST_DIR)/fig4.pdf

fig4: fig5_make.py poredp100k.ods
	./fig5_make.py poredp100k.ods -o $(DEST_DIR)/fig5.pdf


figS1: figS1_make.py
	./figS1_make.py $(DATA_DIR)/vardp100k.dat.gz -j --dpi 300 -o $(DEST_DIR)/figS1.png

figS2: figS2_make.py varrc100k.ods 
	./figS2_make.py varrc100k.ods -o $(DEST_DIR)/figS2.pdf

figS3: figS3_make.py poremsteps.ods
	./figS3_make.py poremsteps.ods -o $(DEST_DIR)/figS3.pdf

figS4: figS4_make.py porerc10k.ods 
	./figS4_make.py porerc10k.ods -o $(DEST_DIR)/figS4.pdf

figS5: figS5_make.py vardp100k.ods
	./figS5_make.py vardp100k.ods -j -o $(DEST_DIR)/figS5.pdf

figS6: figS6_make.py poredp10k.ods
	./figS6_make.py poredp10k.ods -j -o $(DEST_DIR)/figS6.pdf

clean:
	rm -f *~

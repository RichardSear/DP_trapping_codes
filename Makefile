PAPER = trappingms
BIBLIO = trapping
FIGS =

PDFL = pdflatex -synctex=1
BTEX = bibtex

#PDFL = pdflatex -synctex=1 -interaction=batchmode
#BTEX = bibtex -terse

MAIN_FIGS = fig1 fig2 fig3 fig4
SUPP_FIGS = figS1 figS2 figS3 figS4 figS5 figS6

default: spreadsheets figs

spreadsheets:
	./raw_analyse.py --col Dp data/vardp100k.dat.gz -o vardp100k.ods
	./raw_analyse.py --col Dp data/poredp100k.dat.gz -o poredp100k.ods
	./raw_analyse.py --col rc data/varrc100k.dat.gz -o varrc100k.ods
	./raw_analyse.py --col rc data/porerc10k.dat.gz -o porerc10k.ods
	./raw_analyse.py --col code data/poremsteps.dat.gz -o poremsteps.ods

figs: $(MAIN_FIGS) $(SUPP_FIGS)

fig1: fig1_make.py
	./fig1_make.py -o ../fig2.pdf

fig2: fig2_make.py
	./fig2_make.py -o ../fig2.pdf

fig3: fig3_make.py vardp100k.ods
	./fig3_make.py vardp100k.ods -o ../fig3.pdf

fig4: fig4_make.py
	./fig4_make.py -o ../fig4.pdf

figS1: figS1_make.py
	./figS1_make.py data/vardp100k.dat.gz -j --dpi 300 -o ../figS1.png

figS2: figS2_make.py varrc100k.ods 
	./figS2_make.py varrc100k.ods -o ../figS2.pdf

figS3: figS3_make.py poremsteps.ods
	./figS3_make.py poremsteps.ods -o ../figS3.pdf

figS4: figS4_make.py porerc10k.ods 
	./figS4_make.py porerc10k.ods -o ../figS4.pdf

figS5: figS5_make.py vardp100k.ods
	./figS5_make.py vardp100k.ods -j -o ../figS5.pdf

figS6: figS6_make.py poredp10k.ods
	./figS6_make.py poredp10k.ods -j -o ../figS6.pdf

clean:
	rm -f *~

# make -n DEST_DIR=..  overrides the setting below
# Warren and Sear 2025/2026

DEST_DIR = .

PYTHON3 = /usr/bin/python3

MAIN_FIGS = fig_schem fig_pip fig_pip_bd fig_pore fig_pore_bd
SUPP_FIGS_1 = fig_bimodal fig_pip_rc fig_pore_rc fig_pip_extra
SUPP_FIGS_2 = fig_pore_msteps fig_pore_extra fig_pore_bif_fixedp
SUPP_FIGS_EXTRA = fig_pore_bif fig_pore_bif_bd

default: main_figs supp_figs

main_figs: $(MAIN_FIGS)
supp_figs: $(SUPP_FIGS_1) $(SUPP_FIGS_2)
extra_figs: $(SUPP_FIGS_EXTRA)

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

porebif10k25.ods: raw_analyse.py data/porebif10k.dat.gz
	$(PYTHON3) $^ -f k=25 -o $@

porebif10k30.ods: raw_analyse.py data/porebif10k.dat.gz
	$(PYTHON3) $^ -f k=30 -o $@

porebif10k35.ods: raw_analyse.py data/porebiff10k.dat.gz
	$(PYTHON3) $^ -o $@

fig_schem: fig_schem.py
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_schem.pdf

fig_pip: fig_pip.py
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_pip.pdf

fig_pip_bd: fig_pip_bd.py vardp100k.ods
	$(PYTHON3) $^ -j -o $(DEST_DIR)/fig_pip_bd.pdf

fig_pore: fig_pore.py
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_pore.pdf

fig_pore_bd: fig_pore_bd.py poredp100k.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_pore_bd.pdf

fig_bimodal: fig_bimodal.py data/vardp100k.dat.gz
	$(PYTHON3) $^ -j --dpi 300 -o $(DEST_DIR)/fig_bimodal.png

fig_pip_rc: fig_pip_rc.py varrc100k.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_pip_rc.pdf

fig_pore_rc: fig_pore_rc.py porerc10k.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_pore_rc.pdf

fig_pore_msteps: fig_pore_msteps.py poremsteps.ods
	$(PYTHON3) $^ -o $(DEST_DIR)/fig_pore_msteps.pdf

fig_pip_extra: fig_pip_extra.py vardp100k.ods
	$(PYTHON3) $^ -j -o $(DEST_DIR)/fig_pip_extra.pdf

fig_pore_extra: fig_pore_extra.py poredp10k.ods
	$(PYTHON3) $^ -j -o $(DEST_DIR)/fig_pore_extra.pdf

fig_pore_bif_fixedp: fig_pore_bif_fixedp.py
	$(PYTHON3) $< -o $(DEST_DIR)/fig_pore_bif_fixedp.pdf

fig_pore_bif: fig_pore_bif.py
	$(PYTHON3) $< -o $(DEST_DIR)/fig_pore_bif.pdf

fig_pore_bif_bd: fig_pore_bif_bd.py porebif10k25.ods porebif10k30.ods porebif10k35.ods
	$(PYTHON3) $< porebif10k{k}.ods -o $(DEST_DIR)/fig_pore_bif_bd.pdf

clean:
	rm -f *~

#!/usr/bin/bash

echo ./fig1_make.py -o ../fig1.pdf
echo ./fig2_make.py -o ../fig2.pdf
echo ./fig3_make.py vardp100k.ods -o ../fig3.pdf
echo ./fig4_make.py -o ../fig4.pdf

echo ./figS1_make.py data/vardp100k.dat.gz -j --dpi 300 -o ../figS1.png
echo ./figS2_make.py varrc100k.ods -o ../figS2.pdf
echo ./figS3_make.py poremsteps.ods -o ../figS3.pdf
echo ./figS4_make.py porerc10k.ods -o ../figS4.pdf
echo ./figS5_make.py vardp100k.ods -j -o ../figS5.pdf
echo ./figS6_make.py poredp10k.ods -j -o ../figS6.pdf

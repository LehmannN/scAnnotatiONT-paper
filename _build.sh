#!/bin/sh

# Remove previously generated output
rm -ri docs/

# Uncomment next line to activate the conda environment for this study (see README.md file to create the conda environment)
#conda activate scAnnotatiONT

# Generate a HTML version of the report
set -ev
Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')" >& _build.log

# Generate a PDF version of the report (not by default)
#Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::pdf_book')"

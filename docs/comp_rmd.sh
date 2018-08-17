#!/bin/bash

# Set the first argument as variable file
file=$1.Rmd

filename=$(basename "$file")
extension="${filename##*.}"
filename="${filename%.*}"

Rscript -e "rmarkdown::render('$file')"


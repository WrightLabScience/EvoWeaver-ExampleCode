#!/bin/bash

# Rscript JobScript.R ${1} ${2}
Rscript JobScript.R ${1} ${2}

if [ -e *.fasta ]
then
  exit 0
else
  exit 1
fi
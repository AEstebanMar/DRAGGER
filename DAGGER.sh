#! /usr/bin/env bash

###########################################################################################
######## DAGGER: Drug repositioning by Analysis of GWAS and Gene Expression in R ##########
###########################################################################################

# sudo apt install unzip
# sudo apt-get install r-base

### Analysis parameters.

echo -e "\nLaunching DAGGER 1.2\n============================================================================================"

export pval=0.001 ### The p-value cutoff for the GWAS data. Default is 0.001
Qcutoff=5 ### The top percentaje Q-value that will be selected from the GTEx data. Default is 5.

cd ./src

#./preprocessing.sh
#Rscript ./Analysis.R $Qcutoff
./stats_plots.R

echo -e "\n============================================================================================\nDAGGER has finished execution. See Output folder for results.\n============================================================================================\n"


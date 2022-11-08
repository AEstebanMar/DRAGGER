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

./preprocessing.sh
Rscript ./Analysis.R $Qcutoff
./stats_plots.R

echo -e "\n============================================================================================\nDAGGER has finished execution. See Output folder for results.\n============================================================================================\n"
### POR HACER: NO VALE EL REPS PORQUE A PARTIR DE CIERTO P-VALOR DE GWAS APARECE UN RS QUE NO ESTÁ EN GTEX.
### ESTUDIAR SI CON UN MERGE CONSEGUIMOS RESULTADOS EQUIVALENTES, Y SI NO, ELIMINAR EL RS A MANO (COMPARACIÓN
### VECTORIAL DE LOS NOMBRES DE LOS RS QUE ELIMINE DE LA TABLA GWAS LOS QUE NO APAREZCAN EN AMBOS)

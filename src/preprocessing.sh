#! /usr/bin/env bash

### Procesamiento de GWAS. Necesario tener el paquete unzip instalado. Si no está instalado:

# Extraer archivo de GWAS.

echo -e "\nStarting preprocessing...\n"

cd ../data/raw

unzip metaGWAS_repli16dbs_20190930.1tbl.zip

# Contar las líneas del archivo bruto de GWAS. 

totalpols=$(wc -l < metaGWAS_repli16dbs_20190930.1tbl.rs)

# Head del archivo para ver su estructura.

nawk -F, '{if(NR<11) {print}}' metaGWAS_repli16dbs_20190930.1tbl.rs > ../processed/head_GWAS.txt

# Filtrado para tomar los polimorfismos con p-valor menor a (pval).

awk -v pval=$pval '{ if ($9 <= pval) {print} }' metaGWAS_repli16dbs_20190930.1tbl.rs > ../processed/GWAS_filtered.txt

# Recuento de RS tras el filtrado.

filteredpols=$(wc -l < ../processed/GWAS_filtered.txt)

echo -e "\n" $filteredpols "polymorphisms passed filter\n"

# Selección aleatoria de 41919 polimorfismos. Se usarán para buscar falsos positivos más adelante en el análisis.

awk -v filteredpols=$filteredpols 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=filteredpols; i++){x=int(rand()*NR) + 1; print a[x];}}' metaGWAS_repli16dbs_20190930.1tbl.rs > ../processed/GWAS_random.txt
tar -xvf GTEx_Analysis_v8_eQTL.tar -C ../temp
cd ../temp/GTEx_Analysis_v8_eQTL
gunzip *egenes.txt.gz
awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' OFS=, *egenes.txt > ../../processed/Merged_eQTL.txt
TotalGTEx=$(wc -l ../../processed/Merged_eQTL.txt)

rm -r ../GTEx_Analysis_v8_eQTL

cd ../../../Output/

echo "Total polymorphisms from GWAS:" > Summary.txt
echo $totalpols >> Summary.txt ### Un unique a esto en R ha hecho ver que realmente son menos. Hay unas diez mil repeticiones.
echo "Selected p-value:" >> Summary.txt
echo $pval >> Summary.txt
echo "Total filtered polymorphisms:" >> Summary.txt
echo $filteredpols >> Summary.txt
echo -e "\nPreprocessing done!\n\nStarting R script..."

cd ../src/
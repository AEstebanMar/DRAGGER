#! /usr/bin/env bash

### GWAS preprocessing ###

echo -e "\nStarting preprocessing...\n"

cd ../data/raw

unzip metaGWAS_repli16dbs_20190930.1tbl.zip

totalpols=$(wc -l < metaGWAS_repli16dbs_20190930.1tbl.rs)


nawk -F, '{if(NR<11) {print}}' metaGWAS_repli16dbs_20190930.1tbl.rs > ../processed/head_GWAS.txt

awk -v pval=$pval '{ if ($9 <= pval) {print} }' metaGWAS_repli16dbs_20190930.1tbl.rs > ../processed/GWAS_filtered.txt

filteredpols=$(wc -l < ../processed/GWAS_filtered.txt) # Comprobar si esto cuenta los únicos

uniqpols=$(uniq -c ../processed/GWAS_filtered.txt | wc -l)

duplicates=$((filteredpols-uniqpols))

echo -e "\n" $filteredpols "polymorphisms passed filter\n"

# Random selection of 41,919 polymorphisms. Will be used to detect false positives further in the analysis.

awk -v filteredpols=$filteredpols 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=filteredpols; i++){x=int(rand()*NR) + 1; print a[x];}}' metaGWAS_repli16dbs_20190930.1tbl.rs > ../processed/GWAS_random.txt
tar -xvf GTEx_Analysis_v8_eQTL.tar -C ../temp
cd ../temp/GTEx_Analysis_v8_eQTL
gunzip *egenes.txt.gz
awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' OFS=, *egenes.txt > ../../processed/Merged_eQTL.txt
TotalGTEx=$(wc -l ../../processed/Merged_eQTL.txt)

rm -r ../GTEx_Analysis_v8_eQTL

cd ../../../Output/

echo "Total polymorphisms from GWAS:" > Summary.txt
echo $totalpols >> Summary.txt
echo "Selected p-value:" >> Summary.txt
echo $pval >> Summary.txt
echo "Total filtered polymorphisms:" >> Summary.txt
echo $filteredpols >> Summary.txt

cd ../src/
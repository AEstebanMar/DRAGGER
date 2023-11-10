## bash code to prepare `GTEx` dataset. Result is used by GTEx.R
## File downloaded from https://www.gtexportal.org/home/downloads/adult-gtex#qtl

tar -xvf GTEx_Analysis_v8_eQTL.tar -C ../temp
cd ../temp/GTEx_Analysis_v8_eQTL
gunzip *egenes.txt.gz
awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' OFS=, *egenes.txt > Merged_eQTL.txt
rm -r ../GTEx_Analysis_v8_eQTL
#!/bin/bash

### Procesamiento de GWAS. Necesario tener el paquete unzip instalado. Si no está instalado:

# sudo apt install unzip

# Extraer archivo de GWAS.

unzip metaGWAS_repli16dbs_20190930.1tbl.zip

# Contar las líneas del archivo bruto de GWAS. 

wc -l metaGWAS_repli16dbs_20190930.1tbl.rs

# Head del archivo para ver su estructura.

nawk -F, '{if(NR<11) {print}}' metaGWAS_repli16dbs_20190930.1tbl.rs > ../head_GWAS.txt

# Imprimir el head en la consola de bash.

awk '{print}' ../head_GWAS.txt

# Filtrado para tomar los polimorfismos con p-valor menor a 0.001.

awk '{ if ($9 <= 0.001) {print} }' metaGWAS_repli16dbs_20190930.1tbl.rs > ../GWAS_0.001.txt

# Recuento de RS tras el filtrado. Resultado: 41919 polimorfismos.

wc -l ../GWAS_0.001.txt

# Selección aleatoria de 41919 polimorfismos. Se usarán para buscar falsos positivos más adelante en el análisis.

awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=41919; i++){x=int(rand()*NR) + 1; print a[x];}}' metaGWAS_repli16dbs_20190930.1tbl.rs > ../GWAS_random.txt

# Recuento de RS tras el filtrado para ver que todo va bien.

wc -l ../GWAS_random.txt

### Procesamiento de GTEx

# Se extrae la carpeta

tar -xvf GTEx_Analysis_v8_eQTL.tar -C ../

# Cambiar de directorio a la carpeta con los archivos de texto del GTEx.

cd ../GTEx_Analysis_v8_eQTL

# Primero se extraen los archivos .gz

gunzip *egenes.txt.gz

# Ahora se unen los archivos

awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' OFS=, *egenes.txt > ../eQTL_Unidos.txt

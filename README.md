# Análisis de GWAS para el reposicionamiento de fármacos en la Enfermedad de Alzheimer
Script de análisis de datos de GWAS, expresión génica (GTEx) y farmacológica (DGIdb) para el reposicionamiento de fármacos para tratar la enfermedad de Alzheimer. Elaborado para mi TFM en la Universidad de Málaga (Máster en Biología Molecular y Celular).

Archivos necesarios para ejecutar el script:

Los tres primeros archivos se obtienen mediante un corto procesamiento en AWK. Explicado en el archivo preprocessing.sh.

Archivos necesarios:

Datos genómicos: metaGWAS_repli16dbs_20190930.1tbl.rs (Bellenguez *et al*, 2022). Extraer el archivo comprimido.

Datos de expresión génica: GTEx_Analysis_v8_eQTL.tar (https://www.gtexportal.org/home/datasets). Extraer archivos, se trabajará sobre los terminados en egenes.gz. (Información de GTEx).

Datos farmacológicos: interactions.tsv (datos de fármacos descargados de https://www.dgidb.org/downloads)

Datos de Alzforum descargados desde: https://www.alzforum.org/therapeutics/search?fda_statuses=&target_types=&therapy_types=&conditions=&keywords-entry=&keywords=#results

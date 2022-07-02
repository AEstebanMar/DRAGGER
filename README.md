# Análisis de GWAS para el reposicionamiento de fármacos en la Enfermedad de Alzheimer
Script de análisis de datos de GWAS, expresión génica (GTEx) y farmacológica (DGIdb) para el reposicionamiento de fármacos para tratar la enfermedad de Alzheimer. Elaborado para mi TFM en la Universidad de Málaga (Máster en Biología Molecular y Celular).

Archivos necesarios para ejecutar el script:

Datos genómicos: metaGWAS_repli16dbs_20190930.1tbl.zip (Bellenguez *et al*, 2022).

Datos de expresión génica: GTEx_Analysis_v8_eQTL.tar (https://www.gtexportal.org/home/datasets).

Datos farmacológicos: interactions.tsv (datos de fármacos descargados de https://www.dgidb.org/downloads)

Datos de Alzforum descargados desde: https://www.alzforum.org/therapeutics/search?fda_statuses=&target_types=&therapy_types=&conditions=&keywords-entry=&keywords=#results. Descargar en formato texto. Archivo incluido en este repositorio.

Para ejecutar el programa, crear una carpeta en la que guardar el script de AlzGWAS.R. En esa carpeta, crear una subcarpeta para los archivos metaGWAS_repli16dbs_20190930.1tbl.rs y GTEx_Analysis_v8_eQTL.tar. Guardar en ella el script preprocessing.sh y ejecutarlo. Una vez ejecutado, salir de la subcarpeta y ejecutar el script AlzGWAS.R.

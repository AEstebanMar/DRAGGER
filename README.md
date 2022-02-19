# AlzGWAS
Script de comparación de datos de GWAS, expresión génica (GTEx) y farmacológica para el reposicionamiento de fármacos para tratar la enfermedad de Alzheimer. Elaborado para mi TFM en la Universidad de Málaga (Máster en Biología Molecular y Celular)

Archivos necesarios para ejecutar el script:

Los tres primeros archivos se obtienen mediante un corto procesamiento en AWK. Explicado en el archivo bashscript.txt

eQTL_Unidos.txt (Información de GTEx). Datos original descargados desde https://www.gtexportal.org/home/datasets (GTEx_Analysis_v8_eQTL.tar)

GWAS_0.001.txt (RS asociados significativamente según el GWAS). Datos obtenidos de la publicación 

GWAS_random.txt (RS escogidos aleatoriamente del GWAS).

head_GWAS.txt (head del GWAS completo, necesario para asignar nombres a las columnas en R)

interactions.tsv (datos de fármacos descargados de https://www.dgidb.org/downloads)

Datos de Alzforum descargados desde: https://www.alzforum.org/therapeutics/search?fda_statuses=&target_types=&therapy_types=&conditions=&keywords-entry=&keywords=#results

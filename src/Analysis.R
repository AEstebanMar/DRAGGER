#! /usr/bin/env Rscript

### AlzGWAS: An R Tool for Drug Repurposing by Functional Annotation of GWAS in Alzheimer's Disease

Qcutoff = as.numeric(tail(commandArgs(),1))
setwd('../data/raw')
message('\n============================================================================================\nLoading input data')
message('\nLoading Alzforum table')
Alzforum = read.table ("Alzforum Therapeutics.txt", sep = "\t", header=TRUE, dec=".", fill = TRUE, quote = "")
message('\nLoading DGIdb table')
Interactions = read.csv("interactions.tsv", sep = "\t")
setwd('../processed')
# There are some duplicates in the significant RS table. This will be corrected later, as it will be easier
# after some downscaling.
message('\nLoading GWAS data')
Sign = read.table("GWAS_filtered.txt", header=FALSE, sep="", dec=".")
Rand = read.table("GWAS_random.txt", header=FALSE, sep="", dec=".")
message('\nLoading GTEx data (this might take a bit)')
GTEx = read.table ("Merged_eQTL.txt", header=TRUE, sep="", dec=".", fill = TRUE)
TotalGenes = length(unique(GTEx$gene_name))
GWASheader = read.table ("head_GWAS.txt", header=FALSE, sep="", dec=".") [1,]	### Tal y como he escrito los códigos anteriores la tabla del GWAS se ha generado sin nombres de columnas. Esta línea lo corrige.
message('\nInput data loaded!')
colnames(Sign) = colnames(Rand) = GWASheader
saveRDS(Sign, file = "Sign.rds")
saveRDS(Rand, file = "Rand.rds")

######################################### GWAS-GTEx merging ###############################################

message('\nApplying Q-value filter...')

RSSign = Sign [,"RS"]					

SortedGTEx = GTEx[order(GTEx$qval),]			

# Se escoge el 5 % más significativo de los valores de GTEx. Considerar cambio a criterio simple de p-valor, no tan restrictivo.

TopQval = SortedGTEx[seq(nrow(SortedGTEx)*(Qcutoff/100)),]

# Se ajusta el nombre de la columna de los RS, que por cómo se ha generado el archivo de GTEx no es adecuado.

colnames(TopQval)[19] = "rs_id"
TopGenes = length(unique(TopQval$gene_name))
saveRDS(TopQval, "TopQval.rds")

# Borro estos marcos de datos del área de trabajo, ya que no vuelven a hacer falta y ocupan mucho espacio en la RAM.

rm(ls=GTEx,ls=SortedGTEx)

message('Done!')
message('\nMatching GWAS and GTEx Data...')
# Almaceno un vector con los RS del top 5 por ciento del GTEx.
									
RSTop <- TopQval[,"rs_id"]

MatchGTEx <- TopQval[RSTop %in% RSSign,]	

# Almacena las filas de RSSign que se han encontrado.

MatchGWAS <- Sign[RSSign %in% RSTop,]		

# Eliminación de RS repetidos en tabla de GWAS.

MatchGWAS <- MatchGWAS[!duplicated(MatchGWAS$RS),]

# Ahora las ordeno en función del RS, es necesario que estén en el mismo orden porque se van a unir en función del RS.

MatchGTEx <- MatchGTEx[order(MatchGTEx$rs_id),]

MatchGWAS <- MatchGWAS[order(MatchGWAS$RS),]

# Es posible que ambos marcos de datos generados tengan un número distinto de elementos, ya que hay RS que aparecen varias veces en el GTEx (están relacionados con más de un gen) mientras
# que en el GWAS deberían aparecer una sola vez. Si es el caso, repetimos las filas necesarias hasta igualarlas en número. Para arreglarlo voy a anotar los RS de cada una, a ordenar
# los vectores resultantes y usar los datos de tabulación de los RS de GTEx para generar un vector que me diga cuántas veces hay que repetir cada fila del GWAS.

MatchRSGTEx = MatchGTEx[,"rs_id"]	### Anoto los RS que aparecen en los marcos de datos.

MatchRSGWAS = MatchGWAS[,"RS"]

MatchRSGTEx = sort(MatchRSGTEx)	### Las ordeno numéricamente por el RS

MatchRSGWAS = sort(MatchRSGWAS)

# Recuentos de cuántas veces aparece cada RS en la tabla GTEx. Necesario para unir ambas tablas.

TablaGTEx = table(MatchRSGTEx)

# Este resultado es un vector con cuántas veces se repite cada RS en la tabla de GTEx. Como en la tabla de GWAS cada RS aparece una sola vez, tiene tantas filas como
# elementos tiene este vector. Por tanto, podemos utilizar el vector para indicar cuántas veces debe repetirse cada elemento de la tabla. Esta parte del código se basa
# en la respuesta de Andrew a la pregunta planteada en https://stackoverflow.com/questions/29743691/duplicate-rows-in-a-data-frame-in-r

reps = rep(1:nrow(MatchGWAS), TablaGTEx)

MatchGWASajustado = MatchGWAS[reps,]

# Se crea la matriz con las coincidencias, uniendo la información de ambos archivos.

result = cbind(MatchGTEx, MatchGWASajustado)			

# Ahora mismo la primera columna es algo fea, así que la retoco con ayuda del paquete "stringr"

colnames(result)[1]="Tissue"

# Considerar sustituir esto por  sub(" .v8.egenes.", "", test) o algo así para eliminar la necesidad del paquete stringr.

Unoydos = stringr::str_split_fixed(result$Tissue, ".v8.egenes.txt,", 2)

result = result[-1]

colnames(Unoydos) = c("Tissue","gene_id")

finalresult = cbind(Unoydos,result)

SigMerged = finalresult[order(finalresult$Tissue),]

matchedSigPols <- length(SigMerged$RS)

saveRDS(SigMerged, file='SigMerged.rds')

### Repetimos para los aleatorios.

# Creo dos vectores con los RS de cada grupo. Importante eliminar los NA, son RS que no están genotipados.

RSRand = Rand [,"RS"]				

MatchGTEx = TopQval [RSTop %in% RSRand,]	

# Almacena las filas de RSSign que se han encontrado.	

MatchRand = Rand [RSRand %in% RSTop,]

# Eliminación de RS repetidos en tabla de GWAS.

MatchRand = MatchRand[!duplicated(MatchRand),]

# Ahora las ordeno en función del RS, es necesario que estén en el mismo orden porque se van a unir en función del RS.

MatchGTEx = MatchGTEx [order(MatchGTEx$rs_id),]		

MatchRand = MatchRand [order(MatchRand$RS),]

# Es posible que ambos marcos de datos generados tengan un número distinto de elementos, ya que hay RS que aparecen varias veces en el GTEx (están relacionados con más de un gen) mientras
# que en el GWAS deberían aparecer una sola vez. Si es el caso, repetimos las filas necesarias hasta igualarlas en número. Para arreglarlo voy a anotar los RS de cada una, a ordenar
# los vectores resultantes y usar los datos de tabulación de los RS de GTEx para generar un vector que me diga cuántas veces hay que repetir cada fila del GWAS.

# Anoto los RS que aparecen en los marcos de datos.

MatchRSGTEx = MatchGTEx[,"rs_id"]	

MatchRSRand = MatchRand[,"RS"]

# Las ordeno numéricamente por el RS

MatchRSGTEx = sort(MatchRSGTEx)	

MatchRSRand = sort(MatchRSRand)

# Recuentos de cuántas veces aparece cada RS en la tabla GTEx. Necesario para unir ambas tablas.

TablaGTEx = table(MatchRSGTEx)

# Este resultado es un vector con cuántas veces se repite cada RS en la tabla de GTEx. Como en la tabla de GWAS cada RS aparece una sola vez, tiene tantas filas como

# elementos tiene este vector. Por tanto, podemos utilizar el vector para indicar cuántas veces debe repetirse cada elemento de la tabla. Esta parte del código se basa

# en la respuesta de Andrew a la pregunta planteada en https://stackoverflow.com/questions/29743691/duplicate-rows-in-a-data-frame-in-r

reps = rep(1:nrow(MatchRand), TablaGTEx)

MatchRandajustado = MatchRand[reps,]

# Se crea la matriz con las coincidencias, uniendo la información de ambos archivos.

result = cbind(MatchGTEx,MatchRandajustado)			

# Ahora mismo la primera columna es algo fea, así que la retoco con ayuda del paquete "stringr"

colnames(result)[1]="Tissue"

Unoydos = stringr::str_split_fixed(result$Tissue, ".v8.egenes.txt,", 2)

result = result[-1]

colnames(Unoydos) = c("Tissue","gene_id")

finalresult = cbind(Unoydos,result)

RandMerged = finalresult[order(finalresult$Tissue),]

matchedRandPols <- length(RandMerged$RS)

saveRDS(RandMerged, file = "RandMerged.rds")

#eQTLs <- data.frame("eQTL"=c(unique(SigMerged$rs_id),unique(RandMerged$rs_id)), "non-eQTL"=unique())

#eQTLs = data.frame("eQTL" = c(439, 107),"no eQTL" = c(41480, 41812))
#rownames(eQTLs) = c("Significant SNPs", "Random SNPs")

# Se genera un archivo de texto con la lista de dianas aleatorias. Esto servirá para construir la red en STRING.

setwd('../../output/tables')

write.table(RandMerged[,"gene_name"], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE, "RandomTargets.txt")

### Depuración para eliminar los probables falsos positivos.

# Se anotan las ID en Ensembl de los genes encontrados.

SigGenes = SigMerged[,"gene_id"]

RandGenes = RandMerged[,"gene_id"]

message('Done!')
message('\nRemoving false positives...')

# Se anotan las posiciones en la tabla de significativos de los falsos positivos, que son los genes del archivo de resultados significativos que se encuentran también en el archivo de aleatorios.

FPs = which (SigGenes %in% RandGenes)

# Se eliminan dichos genes de SigGenes, almacenando el resultado en una nueva variable.
{
	if(length(FPs)>0)
	{
		DepMerge = SigMerged[-FPs,]
	}
	else
	{
		DepMerge = SigMerged
	}

}
# Se guarda el resultado en un archivo. El condicional se asegura de que ha salido bien, de ser así vale TRUE y el archivo se guarda.

{
	if (!any(DepMerge[,"gene_id"] %in% RandGenes))

	{

	message('Depuration complete!')

	}
}

# Se genera un archivo de texto con la lista de dianas. Esto servirá para construir la red en STRING.

PotentialTargets = length(unique(DepMerge$gene_name))

DepEQTL = length(unique(DepMerge$rs_id))

write.table(DepMerge[,"gene_name"], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE, "Targets.txt")

### Análisis de fármacos. Nota: aquí también hago análisis aleatorio, pero realmente no debería ser aplicable. Los genes falsos positivos ya se han retirado,
### así que aquí simplemente eliminaré los fármacos que TAMBIÉN tengan un efecto sobre un gen no implicado en la enfermedad. Esto sólo tiene sentido a nivel de
### efectos secundarios, algo que no es cribable a este nivel.

message('\nStarting drug interaction analysis...')

# Se sustituyen las casillas vacías para evitar errores al construir el archivo final.

Interactions [Interactions==""] = "Missing"

# Anotación de filas que coinciden.

# Primero las filas del GWAS-GTEx para las que se ha encontrado fármaco en la base de datos.

FoundGenes = DepMerge$gene_name %in% Interactions$gene_name

Druggable = DepMerge [FoundGenes,]

# Después las filas con los fármacos encontrados. Se elimina la columna gene_name del marco de datos de fármacos, que ya está en el GWAS-GTEx.

FoundDrugs = Interactions$gene_name %in% DepMerge$gene_name

Drugs = Interactions [FoundDrugs,]

# Como hay varios fármacos por gen, hay más fármacos que genes. Para unir las dos matrices tengo que ajustarlo de la misma forma que ajusté la matriz del GWAS.

# Primero las ordeno según el nombre del gen.

Druggable = Druggable [order(Druggable$gene_name),]

Drugs = Drugs [order(Drugs$gene_name),]

# Anoto dichos nombres.

DruggableGenes = Druggable [,"gene_name"]	

DrugGenes = Drugs [,"gene_name"]

# Por la forma en que se unen las matrices, el nombre del gen aparecerá sólo una vez. Para comprobar más adelante que todo ha salido bien, añado otra columna con esta información

# que no se eliminará.

Drugs = cbind (DrugGenes, Drugs)

# Se crea la matriz con las coincidencias, uniendo la información de ambos archivos.

Pharma = merge (Druggable, Drugs, by.x = "gene_name", by.y = "gene_name")

# Comprobación:

if (all (Pharma$gene_name == Pharma$DrugGenes) )

{

	if (any (colnames (Pharma) == "DrugGenes"))

		{

		Pharma = Pharma [-which (colnames (Pharma) == "DrugGenes") ]

		}

	message('Drugs identified successfully!')

}

### Escoger fármacos con efecto deseado

message('\nFiltering for desired interaction...')

# Beta positiva o negativa. Almacenado de esta forma, TRUE significa protector y FALSE significa riesgo.

Betas = Pharma$beta < 0

# Pendiente positiva o negativa. Almacenado de esta forma, TRUE significa aumento en la expresión y FALSE significa disminución.

Slopes = Pharma$slope > 0

# Se crea un vector para decidir si se necesita un inhibidor o un activador.

Interaction = Betas == Slopes

# Si un valor de Betas coincide con el de Slopes, significa que es un RS protector que activa expresión génica o un RS de riesgo que la disminuye. En cualquier caso, hará falta un activador.

Interaction [Interaction == TRUE] = "activator"

# Si un valor de Betas NO coincide con el de Slopes, significa que es un RS protector que inhibe expresión génica o un RS de riesgo que la aumenta. En cualquier caso, hará falta un inhibidor.

Interaction [Interaction == FALSE] = "inhibitor"

# Se crea el vector que decide si se recomienda o no el fármaco. Un valor de TRUE significa que se recomienda, los valores de FALSE se deben procesar un poco más.

Recommended = Interaction == Pharma$interaction_types

# Las posiciones que valen FALSE porque el fármaco detectado no es activador ni inhibidor (o falta información) se cambian por una interrogación.

Recommended [!Pharma$interaction_types %in% c("inhibitor","activator")] = "?"

RecPharma = cbind (Pharma, Interaction, Recommended)

DruggableTargets = length(unique(RecPharma$gene_name))

Candidates = RecPharma[Recommended==TRUE,]

TreatableTargets = length(unique(Candidates$gene_name))

message('Filtering complete!')

### Comparo estos resultados con la tabla de Alzforum para saber cuáles se han probado y cuáles no.
	
message('\nRemoving Alzforum matches...')

# Ahora comparo con el archivo de candidatos. Probando con la metformina (350 en Candidates, aparece como 351 porque tiene header)
# Para ello hay que asegurarse de que los nombres están escritos en el mismo formato, lo que se consigue con los siguientes comandos.

# Los nombres de los candidatos no están exclusivamente con el fármaco, sino con la forma concreta (por ejemplo, no aparece "Metformin" sino "Metformin hydrochloride").
# En este primer comando se elimina este problema dividiendo el nombre y quedándonos con la primera palabra. Además, se pasa todo a mayúsculas.
	
CandidatosNombre = toupper(stringr::str_split_fixed(Candidates$drug_claim_primary_name, " ", 2)[,1])

AlzforumNombre = toupper(stringr::str_split_fixed(Alzforum$Name, " ", 2)[,1])

# A continuación se eliminan los caracteres especiales (como el símbolo de marca registrada) y los espacios, que por el método en el que se han descargado
# los datos de Alzforum pueden aparecer. Para ello se sustituyen todos por apóstrofos (').

CandidatosNombre = gsub ( "[^0-9A-Za-z///' ]" , "'" , CandidatosNombre)

AlzforumNombre = gsub ( "[^0-9A-Za-z///' ]" , "'" , AlzforumNombre)

# Y para finalizar se eliminan los apóstrofos.

CandidatosNombre = gsub ( "'" , "" , CandidatosNombre)

AlzforumNombre = gsub ( "'" , "" , AlzforumNombre)

# Nota importante: este método no es perfecto. Hay alguna excepción en la lista de candidatos que no se procesa bien de esta manera. Sin embargo, en este caso no es relevante.
# Los errores son en fármacos cuyo nombre es del tipo "Compound 31" o "EGFR inhibitor", pero se trata de moléculas que todavía no están del todo
# desarrolladas y por tanto no se han probado en nada. Por si acaso, he comprobado estos casos a mano, y no aparecen en la lista de Alzforum.

# Ahora se hace la comparación. De la forma en la que se ha implementado, se creará un vector en el que aparecerá TRUE si el fármaco en cuestión no se ha probado y por tanto vale para nuestro ensayo.

Nuevo = !CandidatosNombre%in%AlzforumNombre

# Se añade este vector al archivo de candidatos.

NewCandidates = Candidates[Nuevo,]

# Se genera el archivo de candidatos.

write.table(NewCandidates, row.names = FALSE, sep = "\t", "Results.tsv")

message('Done!')

message('\nAnalysis complete!\n')

setwd('../')
cat("Total GTEx Genes:\n", file="Summary.txt", append = TRUE)
cat(TotalGenes, file="Summary.txt", append = TRUE)
cat("\nQ-value cutoff:\n", file="Summary.txt", append = TRUE)
cat(Qcutoff, file="Summary.txt", append = TRUE)
cat(" %\n", file="Summary.txt", append = TRUE)
cat("Filtered genes:\n", file="Summary.txt", append = TRUE)
cat(TopGenes, file="Summary.txt", append = TRUE)
cat("\nMax Q-value:\n", file="Summary.txt", append = TRUE)
cat(TopQval[nrow(TopQval),"qval"], file="Summary.txt", append = TRUE)
cat("\nSignificant eQTLs:\n", file="Summary.txt", append = TRUE)
cat(matchedSigPols, file="Summary.txt", append = TRUE)
cat("\nRandom eQTLs:\n", file="Summary.txt", append = TRUE)
cat(matchedRandPols, file="Summary.txt", append = TRUE)
cat("\nFinal eQTLs:\n", file="Summary.txt", append = TRUE)
cat(DepEQTL, file="Summary.txt", append = TRUE)
cat("\nTotal potential targets:\n", file="Summary.txt", append = TRUE)
cat(PotentialTargets, file="Summary.txt", append = TRUE)
cat("\nTotal druggable targets:\n", file="Summary.txt", append = TRUE)
cat(DruggableTargets, file="Summary.txt", append = TRUE)
cat("\nTotal treatable targets:\n", file="Summary.txt", append = TRUE)
cat(TreatableTargets, file="Summary.txt", append = TRUE)


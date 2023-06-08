#! /usr/bin/env Rscript

### AlzGWAS: An R Tool for Drug Repurposing by Functional Annotation of GWAS in Alzheimer's Disease

Qcutoff <- as.numeric(tail(commandArgs(),1))
setwd('../data/raw')
message('\n============================================================================================\nLoading input data')
message('\nLoading Alzforum table')
Alzforum <- read.table ("Alzforum Therapeutics.txt", sep = "\t", header=TRUE, dec=".", fill = TRUE, quote = "")
message('\nLoading DGIdb table')
Interactions <- read.csv("interactions.tsv", sep = "\t")
setwd('../processed')
# There are some duplicates in the significant RS table. This will be corrected later, as it will be easier
# after some downscaling.
message('\nLoading GWAS data')
Sign <- read.table("GWAS_filtered.txt", header=FALSE, sep="", dec=".")
Rand <- read.table("GWAS_random.txt", header=FALSE, sep="", dec=".")
message('\nLoading GTEx data (this might take a bit)')
GTEx <- read.table ("Merged_eQTL.txt", header=TRUE, sep="", dec=".", fill = TRUE)
TotalGenes <- length(unique(GTEx$gene_name))
TotalEQTL <- length(GTEx$rs_id)
GWASheader <- read.table("header_GWAS.txt", header=FALSE, sep="", dec=".")	### Tal y como he escrito los códigos anteriores la tabla del GWAS se ha generado sin nombres de columnas. Esta línea lo corrige.
message('\nInput data loaded!')
colnames(Sign) = colnames(Rand) <- GWASheader
saveRDS(Sign, file <- "Sign.rds")
saveRDS(Rand, file <- "Rand.rds")

######################################### GWAS-GTEx merging ###############################################

message('\nApplying Q-value filter...')

RSSign <- Sign [,"RS"]					

SortedGTEx <- GTEx[order(GTEx$qval),]			

# Top 5 % GTEx selection

TopQval <- SortedGTEx[seq(nrow(SortedGTEx)*(Qcutoff/100)),]

colnames(TopQval)[19] = "rs_id"
TopGenes <- length(unique(TopQval$gene_name))

saveRDS(TopQval, "TopQval.rds")

# Cleanup

rm(ls=GTEx,ls=SortedGTEx)

message('Done!')
message('\nMatching GWAS and GTEx Data...')
									
RSTop <- TopQval[,"rs_id"]
TopEQTL <- length(RSTop)

MatchGTEx <- TopQval[RSTop %in% RSSign,]	

# Saving of RSSign rows with matches in RSTop.

MatchGWAS <- Sign[RSSign %in% RSTop,]		

# Eliminación de RS repetidos en tabla de GWAS.

MatchGWAS <- MatchGWAS[!duplicated(MatchGWAS$RS),]

# Sorted by RS. Necessary for merging.

MatchGTEx <- MatchGTEx[order(MatchGTEx$rs_id),]

MatchGWAS <- MatchGWAS[order(MatchGWAS$RS),]

# Correction due to SNPs appearing more than once in the GTEx dataset, as they can be linked to the expression of several genes.

MatchRSGTEx = MatchGTEx[,"rs_id"]

MatchRSGWAS = MatchGWAS[,"RS"]

MatchRSGTEx = sort(MatchRSGTEx)	# Sorting both with the same criterion to allow merging.

MatchRSGWAS = sort(MatchRSGWAS)


TablaGTEx = table(MatchRSGTEx)

# The result is a vector of the times an RS is repeated in the GTEx table. The GWAS RS are sorted in the same way, and they both
# have the same list of RS. This allows to combine the repetitions with the GWAS SNPs to duplicate the necessary rows.
# Approach based on Andrew's response to the following post:
# https://stackoverflow.com/questions/29743691/duplicate-rows-in-a-data-frame-in-r

reps = rep(1:nrow(MatchGWAS), TablaGTEx)

MatchGWASajustado = MatchGWAS[reps,]

result = cbind(MatchGTEx, MatchGWASajustado)			

colnames(result)[1]="Tissue"

Unoydos = stringr::str_split_fixed(result$Tissue, ".v8.egenes.txt,", 2)

result = result[-1]

colnames(Unoydos) = c("Tissue","gene_id")

finalresult = cbind(Unoydos,result)

SigMerged = finalresult[order(finalresult$Tissue),]

matchedSigPols <- length(SigMerged$RS)

saveRDS(SigMerged, file='SigMerged.rds')

### Same process applied to the random RS.

RSRand = Rand [,"RS"]				

MatchGTEx = TopQval [RSTop %in% RSRand,]	

MatchRand = Rand [RSRand %in% RSTop,]

MatchRand = MatchRand[!duplicated(MatchRand),]

MatchGTEx = MatchGTEx [order(MatchGTEx$rs_id),]		

MatchRand = MatchRand [order(MatchRand$RS),]

MatchRSGTEx = MatchGTEx[,"rs_id"]	

MatchRSRand = MatchRand[,"RS"]

MatchRSGTEx = sort(MatchRSGTEx)	

MatchRSRand = sort(MatchRSRand)

TablaGTEx = table(MatchRSGTEx)

reps = rep(1:nrow(MatchRand), TablaGTEx)

MatchRandajustado = MatchRand[reps,]

result = cbind(MatchGTEx,MatchRandajustado)			

colnames(result)[1]="Tissue"

Unoydos = stringr::str_split_fixed(result$Tissue, ".v8.egenes.txt,", 2)

result = result[-1]

colnames(Unoydos) = c("Tissue","gene_id")

finalresult = cbind(Unoydos,result)

RandMerged = finalresult[order(finalresult$Tissue),]

matchedRandPols <- length(RandMerged$RS)

saveRDS(RandMerged, file = "RandMerged.rds")

# A text file is generated for the random targets. Will be used for STRING network analysis.

setwd('../../output/tables')

write.table(RandMerged[,"gene_name"], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE, "RandomTargets.txt")

### Removal of false positives.

SigGenes = SigMerged[,"gene_id"]

RandGenes = RandMerged[,"gene_id"]

message('Done!')
message('\nRemoving false positives...')

# Significant genes found in the random dataset are saved as FP (False Positives).

FPs = which (SigGenes %in% RandGenes)

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


{
	if (!any(DepMerge[,"gene_id"] %in% RandGenes))

	{

	message('Depuration complete!')

	}
}

# Final targets list saved in a text file for STRING network analysis.

PotentialTargets = length(unique(DepMerge$gene_name))

DepEQTL = length(unique(DepMerge$rs_id))

write.table(DepMerge[,"gene_name"], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE, "Targets.txt")

### Drug analysis. A random analysis is also performed, removing coincidences for the results. This reduces the likelihood
### of the found drugs having side effects, an important step on the priorization. Future versions of DAGGER will include
### an option to disable this step if a larger drugs list is desired.

message('\nStarting drug interaction analysis...')

# Empty fields are switched for "Missing" string.

Interactions [Interactions==""] <- "Missing"

# Store GWAS-GTEx rows for which a drug has been found.

FoundGenes <- DepMerge$gene_name %in% Interactions$gene_name

Druggable <- DepMerge [FoundGenes,]

# Store DGIdb rows with mach in GWAS-GTEx table.

FoundDrugs <- Interactions$gene_name %in% DepMerge$gene_name

Drugs <- Interactions [FoundDrugs,]

# Adjustments neccesary due to there being more drugs than genes.

Druggable <- Druggable[order(Druggable$gene_name),]

Drugs <- Drugs[order(Drugs$gene_name),]

DruggableGenes <- Druggable [,"gene_name"]	

DrugGenes <- Drugs [,"gene_name"]

# Due to the way the tables are merged, the gene name column used for merging will be deleted.
# To prevent this information from being lost, the column is duplicated first.

Drugs <- cbind(DrugGenes, Drugs)


Pharma <- merge(Druggable, Drugs, by.x="gene_name", by.y="gene_name")

# Final check:

if(all(Pharma$gene_name==Pharma$DrugGenes))

{

	if(any(colnames(Pharma)=="DrugGenes"))

		{

		Pharma <- Pharma[-which(colnames(Pharma)=="DrugGenes")]

		}

	message('Drugs identified successfully!')

}

### Filtering for drugs with desired effect

message('\nFiltering for desired interaction...')

Betas <- Pharma$beta < 0


Slopes <- Pharma$slope > 0

Interaction <- Betas == Slopes

# If the "Betas" value matches the "Slopes" value, it means a protector SNP increasing expression
# levels or a risk SNP decreasing them has been found. In both cases, an activator is desired.

Interaction [Interaction == TRUE] = "activator"

# If they do NOT match, the inverse situations are true. Thus, an inhibitor is desired.

Interaction [Interaction == FALSE] = "inhibitor"

Recommended = Interaction == Pharma$interaction_types

# Cases where the known interaction is not "activator" or "inhibitor", a question mark is introduced.
# Deciding if said drugs are appropriate requires further analysis, currently not included in DAGGER.

Recommended [!Pharma$interaction_types %in% c("inhibitor","activator")] = "?"

RecPharma = cbind (Pharma, Interaction, Recommended)

DruggableTargets = length(unique(RecPharma$gene_name))

Candidates = RecPharma[Recommended==TRUE,]

TreatableTargets = length(unique(Candidates$gene_name))

message('Filtering complete!')

### Found drugs are compared with Alzforum Therapeutics Alzheimer Disease table to check for already tested drugs.
	
message('\nRemoving Alzforum matches...')

# Unification of name formats. First word is selected, as we will consider entries such as "Metformin" and
# "Metformin Hydroclorate" to be the same drug.

CandidatesName = toupper(stringr::str_split_fixed(Candidates$drug_claim_primary_name, " ", 2)[,1])

AlzforumName = toupper(stringr::str_split_fixed(Alzforum$Name, " ", 2)[,1])

# Elimination of special characters, which are problematic.

CandidatesName = gsub ( "[^0-9A-Za-z///' ]" , "'" , CandidatesName)

AlzforumName = gsub ( "[^0-9A-Za-z///' ]" , "'" , AlzforumName)

# Apostrophes are removed.

CandidatesName = gsub ( "'" , "" , CandidatesName)

AlzforumName = gsub ( "'" , "" , AlzforumName)

# This approach has a limitation. Drugs with names such as "Compund 31" or "EGFR inhibitor" are not properly considered
# in Candidates-Alzforum matching. However, such names are assigned to little-known drugs, which are yet to be tested
# for any disease and therefore will not have a match. In any case, they are in very small number, making manual
# confirmation very easy.

New = !CandidatesName%in%AlzforumName

NewCandidates = Candidates[New,]

write.table(NewCandidates, row.names = FALSE, sep = "\t", "Results.tsv")

message('Done!')

message('\nAnalysis complete!\n')

setwd('../')

cat("\nTotal GTEx Genes:\n", file="Summary.txt", append = TRUE)
cat(TotalGenes, file="Summary.txt", append = TRUE)
cat("\nTotal GTEx eQTL:\n", file="Summary.txt", append = TRUE)
cat(TotalEQTL, file="Summary.txt", append = TRUE)
cat("\nQ-value cutoff:\n", file="Summary.txt", append = TRUE)
cat(Qcutoff, file="Summary.txt", append = TRUE)
cat(" %\n", file="Summary.txt", append = TRUE)
cat("Filtered genes:\n", file="Summary.txt", append = TRUE)
cat(TopGenes, file="Summary.txt", append = TRUE)
cat("\nFiltered GTEx eQTL:\n", file="Summary.txt", append = TRUE)
cat(TopEQTL, file="Summary.txt", append = TRUE)
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


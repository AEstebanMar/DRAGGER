
### FUNCTIONS

remove_duplicate_rs <- function(df) {

	if(!is.null(df$p_value)) {
		message("Sorting input by statistical significance")
		df <- df[order(df$p_value), ]
	}
	if(is.null(df$rs_id)) {
		stop('RS ID column not found or not properly parsed. Please run input
			through DAGGER::parse_column_names()')
	}

	if(!any(duplicated(df$rs_id))) {
		return(df)
	} else {
		message("Duplicates found in dataset, choosing first occurrence")
		res <- df[!duplicated(df$rs_id), ]
		return(res)
	}		
}

filter_significance <- function(df, value = 0.05) {
	if(is.null(df$p_value)) {
		message("No statistical significance column found in input or improperly
			parsed. You might want to run it through DAGGER::parse_column_names
			first. Returning it as-is")
		return(df)
	}
	res <- df[df$p_value <= value, ]
	return(res)
}

DAGGER <- function(GWAS, GTEx, DGI) {
	GWAS <- remove_duplicate_rs(parse_column_names(GWAS))
	GTEx <- filter_significance(parse_column_names(GTEx), 0.05)
	DGI <- parse_column_names(DGI)

	GWAS-GTEx <- merge(GWAS, GTEx, by = "gene_id")
	return(GWAS-GTEx)
}

### Drug analysis. A random analysis is also performed, removing coincidences for the results. This reduces the likelihood
### of the found drugs having side effects, an important step on the priorization. Future versions of DAGGER will include
### an option to disable this step if a larger drugs list is desired.

# message('\nStarting drug interaction analysis...')

# # Empty fields are switched for "Missing" string.

# Interactions [Interactions==""] <- "Missing"

# # Store GWAS-GTEx rows for which a drug has been found.

# FoundGenes <- DepMerge$gene_name %in% Interactions$gene_name

# Druggable <- DepMerge [FoundGenes,]

# # Store DGIdb rows with mach in GWAS-GTEx table.

# FoundDrugs <- Interactions$gene_name %in% DepMerge$gene_name

# Drugs <- Interactions [FoundDrugs,]

# # Adjustments neccesary due to there being more drugs than genes.

# Druggable <- Druggable[order(Druggable$gene_name),]

# Drugs <- Drugs[order(Drugs$gene_name),]

# DruggableGenes <- Druggable [,"gene_name"]	

# DrugGenes <- Drugs [,"gene_name"]

# # Due to the way the tables are merged, the gene name column used for merging will be deleted.
# # To prevent this information from being lost, the column is duplicated first.

# Drugs <- cbind(DrugGenes, Drugs)


# Pharma <- merge(Druggable, Drugs, by.x="gene_name", by.y="gene_name")

# # Final check:

# if(all(Pharma$gene_name==Pharma$DrugGenes))

# {

# 	if(any(colnames(Pharma)=="DrugGenes"))

# 		{

# 		Pharma <- Pharma[-which(colnames(Pharma)=="DrugGenes")]

# 		}

# 	message('Drugs identified successfully!')

# }

# ### Filtering for drugs with desired effect

# message('\nFiltering for desired interaction...')

# Betas <- Pharma$beta < 0


# Slopes <- Pharma$slope > 0

# Interaction <- Betas == Slopes

# # If the "Betas" value matches the "Slopes" value, it means a protector SNP increasing expression
# # levels or a risk SNP decreasing them has been found. In both cases, an activator is desired.

# Interaction [Interaction == TRUE] = "activator"

# # If they do NOT match, the inverse situations are true. Thus, an inhibitor is desired.

# Interaction [Interaction == FALSE] = "inhibitor"

# Recommended = Interaction == Pharma$interaction_types

# # Cases where the known interaction is not "activator" or "inhibitor", a question mark is introduced.
# # Deciding if said drugs are appropriate requires further analysis, currently not included in DAGGER.

# Recommended [!Pharma$interaction_types %in% c("inhibitor","activator")] = "?"

# RecPharma = cbind (Pharma, Interaction, Recommended)

# DruggableTargets = length(unique(RecPharma$gene_name))

# Candidates = RecPharma[Recommended==TRUE,]

# TreatableTargets = length(unique(Candidates$gene_name))

# message('Filtering complete!')

# ### Found drugs are compared with Alzforum Therapeutics Alzheimer Disease table to check for already tested drugs.
	
# message('\nRemoving Alzforum matches...')

# # Unification of name formats. First word is selected, as we will consider entries such as "Metformin" and
# # "Metformin Hydroclorate" to be the same drug.

# CandidatesName = toupper(stringr::str_split_fixed(Candidates$drug_claim_primary_name, " ", 2)[,1])

# AlzforumName = toupper(stringr::str_split_fixed(Alzforum$Name, " ", 2)[,1])

# # Elimination of special characters, which are problematic.

# CandidatesName = gsub ( "[^0-9A-Za-z///' ]" , "'" , CandidatesName)

# AlzforumName = gsub ( "[^0-9A-Za-z///' ]" , "'" , AlzforumName)

# # Apostrophes are removed.

# CandidatesName = gsub ( "'" , "" , CandidatesName)

# AlzforumName = gsub ( "'" , "" , AlzforumName)

# # This approach has a limitation. Drugs with names such as "Compund 31" or "EGFR inhibitor" are not properly considered
# # in Candidates-Alzforum matching. However, such names are assigned to little-known drugs, which are yet to be tested
# # for any disease and therefore will not have a match. In any case, they are in very small number, making manual
# # confirmation very easy.

# New = !CandidatesName%in%AlzforumName

# NewCandidates = Candidates[New,]

# write.table(NewCandidates, row.names = FALSE, sep = "\t", "Results.tsv")

# message('Done!')

# message('\nAnalysis complete!\n')

# setwd('../')

# cat("\nTotal GTEx Genes:\n", file="Summary.txt", append = TRUE)
# cat(TotalGenes, file="Summary.txt", append = TRUE)
# cat("\nTotal GTEx eQTL:\n", file="Summary.txt", append = TRUE)
# cat(TotalEQTL, file="Summary.txt", append = TRUE)
# cat("\nQ-value cutoff:\n", file="Summary.txt", append = TRUE)
# cat(Qcutoff, file="Summary.txt", append = TRUE)
# cat(" %\n", file="Summary.txt", append = TRUE)
# cat("Filtered genes:\n", file="Summary.txt", append = TRUE)
# cat(TopGenes, file="Summary.txt", append = TRUE)
# cat("\nFiltered GTEx eQTL:\n", file="Summary.txt", append = TRUE)
# cat(TopEQTL, file="Summary.txt", append = TRUE)
# cat("\nMax Q-value:\n", file="Summary.txt", append = TRUE)
# cat(TopQval[nrow(TopQval),"qval"], file="Summary.txt", append = TRUE)
# cat("\nSignificant eQTLs:\n", file="Summary.txt", append = TRUE)
# cat(matchedSigPols, file="Summary.txt", append = TRUE)
# cat("\nRandom eQTLs:\n", file="Summary.txt", append = TRUE)
# cat(matchedRandPols, file="Summary.txt", append = TRUE)
# cat("\nFinal eQTLs:\n", file="Summary.txt", append = TRUE)
# cat(DepEQTL, file="Summary.txt", append = TRUE)
# cat("\nTotal potential targets:\n", file="Summary.txt", append = TRUE)
# cat(PotentialTargets, file="Summary.txt", append = TRUE)
# cat("\nTotal druggable targets:\n", file="Summary.txt", append = TRUE)
# cat(DruggableTargets, file="Summary.txt", append = TRUE)
# cat("\nTotal treatable targets:\n", file="Summary.txt", append = TRUE)
# cat(TreatableTargets, file="Summary.txt", append = TRUE)


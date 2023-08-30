#! /usr/bin/env Rscript


### Script to test gprofiler2 on our GWAS, reuse after testing to absorb into DAGGER

RSSign <- readRDS('../data/processed/Sign.rds')$RS
RSRand <- readRDS('../data/processed/Rand.rds')$RS

# Remove NAs

RSSign <- RSSign[!is.na(RSSign)]
RSRand <- RSRand[!is.na(RSRand)]

# Enrichment analysis

SignRes <- gprofiler2::gost(query=RSSign)$result
RandRes <- gprofiler2::gost(query=RSRand)$result

# Structural analysis

SignStruct <- gprofiler2::gsnpense(query=(RSSign))
RandStruct <- gprofiler2::gsnpense(query=(RSRand))

# Convert to character data frame
save.image('Testing.RData')
SignRes <- apply(SignRes,2,as.character)
RandRes <- apply(RandRes,2,as.character)
SignStruct <- apply(SignStruct,2,as.character)
RandStruct <- apply(RandStruct,2,as.character)

# Save output

write.table(SignRes, file = '../output/tables/SignRes.tsv', row.names=FALSE, sep="\t")
write.table(RandRes, file = '../output/tables/RandRes.tsv', row.names=FALSE, sep="\t")
write.table(SignStruct, file = '../output/tables/SignStruct.tsv', row.names=FALSE, sep="\t")
write.table(RandStruct, file = '../output/tables/RandStruct.tsv', row.names=FALSE, sep="\t")

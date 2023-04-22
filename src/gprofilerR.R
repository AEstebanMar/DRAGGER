#! /usr/bin/env Rscript


### Script to test gprofiler2 on our GWAS, reuse after testing to absorb it into DAGGER


RSSign <- readRDS('../data/processed/Sign.rds')$RS
RSRand <- readRDS('../data/processed/Rand.rds')$RS

SignRes <- gprofiler2::gost(query=RSSign)
RandRes <- gprofiler2::gost(query=RSRand)

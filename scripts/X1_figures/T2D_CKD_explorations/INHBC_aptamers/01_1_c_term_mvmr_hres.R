#' MVMR using c-terminus INHBC levels
library(tidyverse)
library(ggplot2)
library(data.table)
library(glue)
library(MVMR)
library(yaml)
source("custom_mv_harmonize.R")

yml_file  <-  "configs/proj_config.yml"
args  <- read_yaml(yml_file)[['T2D_CKD_results']]

outdir  <- "INHBC_aptamers_followup"

if (!file.exists(file.path(outdir,'mvmr_hres_on_inhbc_cterm.Rds'))){
  bmi_exposure  <- readRDS("bmi_processed.exposure_gwas.Rds")
  bmi_outcome  <- readRDS("bmi_processed.outcome_gwas.Rds")

  WHR_exposure  <- readRDS("WHR_processed.exposure_gwas.Rds")
  WHR_outcome  <- readRDS("WHR_processed.outcome_gwas.Rds")


  t2d_exposure  <- readRDS("T2D_processed.exposure_gwas.Rds")
  t2d_outcome  <- readRDS("T2D_processed.outcome_gwas.Rds")


  inhbc_cterm_outcome  <- read.table(
    "deCODE_protein_MR_inputs/OL_NA_SL_INHBC@15686_49_UP_P55103_MR_outcome.tsv.gz",sep="\t",quote="",header=T,colClasses=c("effect_allele.outcome"='character','other_allele.outcome'='character')
  )


  mvmr_hres  <- custom_mv_harmonize_data(
    all_exposures_exposure_gwas = list(
      'bmi' = bmi_exposure,
      'T2D'=t2d_exposure,
      'WHR'=WHR_exposure
    ),
    all_exposures_outcome_gwas = list(
      'bmi'=bmi_outcome,
      'T2D'=t2d_outcome,
      'WHR'=WHR_outcome
    ),
    outcome_gwas = inhbc_cterm_outcome,
    stringency=2,
    refpan = args$ukb_chrpos_refpan,
    mhc_chr=6,
    mhc_pos=c(28477797,33448354),
    plink=args$plink,
    ancestry='EUR',
    tag_kb=1000,tag_r2=0.8
  )
  saveRDS(
    mvmr_hres,
    file=file.path(outdir,'mvmr_hres_on_inhbc_cterm.Rds')
  )
}else{
  mvmr_hres = readRDS(file.path(outdir,'mvmr_hres_on_inhbc_cterm.Rds'))
}

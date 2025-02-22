#' MR using biometric
library(tidyverse)
library(ggplot2)
library(data.table)
library(glue)
library(yaml)
source("harmonize.R")

yml_file  <-  "configs/proj_config.yml"
args  <- read_yaml(yml_file)[['T2D_CKD_results']]


olink_out <- "UKB_replication"
decode_out  <- "INHBC_aptamers_followup"


whr_exposure  <- readRDS("WHR_processed.exposure_gwas.Rds")
bmi_exposure  <- readRDS("bmi_processed.exposure_gwas.Rds")

cterm_out <- read.table("OL_NA_SL_INHBC@15686_49_UP_P55103_MR_outcome.tsv.gz",sep="\t",header=T,quote="")
full_out <- read.table("OL_NA_SL_INHBC@6408_2_UP_P55103_MR_outcome.tsv.gz",sep="\t",header=T,quote="")
olink_outcome  <- readRDS("OL_INHBC_SL_INHBC@15686_49_UP_P55103_MR.outcome_gwas.Rds")



whr_inhbc_olink_hres  <- custom_harmonize_data(
  exposure_gwas = whr_exposure,
  outcome_gwas = olink_outcome,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)
bmi_inhbc_olink_hres  <- custom_harmonize_data(
  exposure_gwas = bmi_exposure,
  outcome_gwas = olink_outcome,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)

whr_inhbc_olink_mr  <- mr(whr_inhbc_olink_hres)
bmi_inhbc_olink_mr  <- mr(bmi_inhbc_olink_hres)




whr_inhbc_cterm_hres  <- custom_harmonize_data(
  exposure_gwas = whr_exposure,
  outcome_gwas = cterm_out,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)
bmi_inhbc_cterm_hres  <- custom_harmonize_data(
  exposure_gwas = bmi_exposure,
  outcome_gwas = cterm_out,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)

whr_inhbc_cterm_mr  <- mr(whr_inhbc_cterm_hres)
bmi_inhbc_cterm_mr  <- mr(bmi_inhbc_cterm_hres)



whr_inhbc_full_hres  <- custom_harmonize_data(
  exposure_gwas = whr_exposure,
  outcome_gwas = full_out,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)
bmi_inhbc_full_hres  <- custom_harmonize_data(
  exposure_gwas = bmi_exposure,
  outcome_gwas = full_out,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)

whr_inhbc_full_mr  <- mr(whr_inhbc_full_hres)
bmi_inhbc_full_mr  <- mr(bmi_inhbc_full_hres)


saveRDS(
  whr_inhbc_olink_hres,
  file.path(olink_out,'whr_inhbc_olink_hres.Rds')  
)
saveRDS(
  whr_inhbc_olink_mr,
  file.path(olink_out,'whr_inhbc_olink_mr.Rds')  
)
saveRDS(
  bmi_inhbc_olink_hres,
  file.path(olink_out,'bmi_inhbc_olink_hres.Rds')  
)
saveRDS(
  bmi_inhbc_olink_mr,
  file.path(olink_out,'bmi_inhbc_olink_mr.Rds')  
)


saveRDS(
  whr_inhbc_cterm_hres,
  file.path(decode_out,'whr_inhbc_cterm_hres.Rds')  
)
saveRDS(
  whr_inhbc_cterm_mr,
  file.path(decode_out,'whr_inhbc_cterm_mr.Rds')  
)
saveRDS(
  bmi_inhbc_cterm_hres,
  file.path(decode_out,'bmi_inhbc_cterm_hres.Rds')  
)
saveRDS(
  bmi_inhbc_cterm_mr,
  file.path(decode_out,'bmi_inhbc_cterm_mr.Rds')  
)

saveRDS(
  whr_inhbc_full_hres,
  file.path(decode_out,'whr_inhbc_full_hres.Rds')  
)
saveRDS(
  whr_inhbc_full_mr,
  file.path(decode_out,'whr_inhbc_full_mr.Rds')  
)
saveRDS(
  bmi_inhbc_full_hres,
  file.path(decode_out,'bmi_inhbc_full_hres.Rds')  
)
saveRDS(
  bmi_inhbc_full_mr,
  file.path(decode_out,'bmi_inhbc_full_mr.Rds')  
)
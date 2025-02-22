#' run mr_steigger_filter first before runing MR.
#' 
library(tidyverse)
library(glue)
library(data.table)
library(TwoSampleMR)
library(yaml)

source("harmonize.R")
source("compute_F_stat.R")
args  <- read_yaml("configs/proj_config.yml")[['T2D_CKD_results']]

t2d_exp <- readRDS("T2D_processed.exposure_gwas.Rds")

inhbc_out <- read.table("OL_NA_SL_INHBC@15686_49_UP_P55103_MR_outcome.tsv.gz",sep="\t",header=T)

#%% with no proxy search
t2d_on_inhbc_hres <- TwoSampleMR::harmonise_data(
    t2d_exp,
    inhbc_out,
    2
)
t2d_on_inhbc_hres_steiger <- TwoSampleMR::steiger_filtering(t2d_on_inhbc_hres)
t2d_on_inhbc_hres_steiger <- t2d_on_inhbc_hres_steiger %>% filter(steiger_dir)
t2d_on_inhbc_steiger_mr <- TwoSampleMR::mr(t2d_on_inhbc_hres_steiger)
#' MR IVW (nsnp, b, se, pval):
#' Inverse variance weighted  193 0.048528737 0.009878637 8.992546e-07
#%% with proxy
t2d_on_inhbc_hres_proxy <- custom_harmonize_data(
  exposure_gwas = t2d_exp,
  outcome_gwas = inhbc_out,
  stringency=2,
  refpan=args$ukb_chrpos_refpan,
  mhc_chr = 6,
  mhc_pos = c(28477797,33448354),
  plink = args$plink,
  ancestry='EUR',
  tag_kb=1000,tag_r2=0.8
)
t2d_on_inhbc_hres_proxy_steiger <- TwoSampleMR::steiger_filtering(t2d_on_inhbc_hres_proxy)
t2d_on_inhbc_hres_proxy_steiger <- t2d_on_inhbc_hres_proxy_steiger %>% filter(steiger_dir)
t2d_on_inhbc_hres_proxy_steiger_mr <- TwoSampleMR::mr(t2d_on_inhbc_hres_proxy_steiger)
#' MR IVW (nsnp, b, se, pval):
#' Inverse variance weighted 201 0.044496039 0.00985282 6.299741e-06
saveRDS(
    t2d_on_inhbc_hres_proxy_steiger,
    "t2d_on_inhbc_hres_proxy_steiger.Rds"
)
saveRDS(
    t2d_on_inhbc_hres_proxy_steiger_mr,
    "t2d_on_inhbc_hres_proxy_steiger_mr.Rds"
)
#%% proxy search F-stat
t2d_on_inhbc_hres_proxy_steiger <- t2d_on_inhbc_hres_proxy_steiger %>% filter(mr_keep)
t2d_on_inhbc_steiger_r2 <- compute_var_explained(
    beta = t2d_on_inhbc_hres_proxy_steiger$beta.exposure,
    maf = t2d_on_inhbc_hres_proxy_steiger$eaf.exposure,
    se = t2d_on_inhbc_hres_proxy_steiger$se.exposure,
    N = t2d_on_inhbc_hres_proxy_steiger$samplesize.exposure
)
t2d_on_inhbc_steiger_fstat <- compute_fstat(
    r2 = t2d_on_inhbc_steiger_r2,
    n = 146550,
    k = 201
)
#' F-stat for T2D is 79
#%%
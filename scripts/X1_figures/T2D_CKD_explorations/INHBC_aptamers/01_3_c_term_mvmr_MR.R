#' MVMR using c-terminus of INHBC
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
mvmr_hres = readRDS(file.path(outdir,'mvmr_hres_on_inhbc_cterm.Rds'))
cterm_phenocov  <- readRDS(file.path(outdir,'t2d_bmi_whr_phenocov_cterm.Rds'))

mvmr_hres_mvmr_package  <- get_mvmr_inputs_from_twosampleMR(mvmr_hres)
mvmr_res  <- TwoSampleMR::mv_multiple(mvmr_hres)
exposure_beta_cols  <- colnames(mvmr_hres_mvmr_package)[!grepl("outcome",colnames(mvmr_hres_mvmr_package)) & grepl("beta",colnames(mvmr_hres_mvmr_package))]
exposure_se_cols  <- colnames(mvmr_hres_mvmr_package)[!grepl("outcome",colnames(mvmr_hres_mvmr_package)) & grepl("se",colnames(mvmr_hres_mvmr_package))]
mvmr_package_hres_formatted  <- MVMR::format_mvmr(
  BXGs = mvmr_hres_mvmr_package[exposure_beta_cols],
  BYG = mvmr_hres_mvmr_package['outcome_beta'],
  seBXGs = mvmr_hres_mvmr_package[exposure_se_cols],
  seBYG = mvmr_hres_mvmr_package['outcome_se'],
  RSID = mvmr_hres_mvmr_package['SNP']
)


phenocov  <- phenocov_mvmr(cterm_phenocov,mvmr_hres_mvmr_package[exposure_se_cols])
res  <- pleiotropy_mvmr(r_input = mvmr_package_hres_formatted, gencov = phenocov)
compute_isq  <- function(Q,Q_df){
  I_sq  <- 100*(Q - Q_df)/Q
  return(I_sq)
}
compute_isq(res$Qstat,1287)
saveRDS(
  mvmr_res,
  file.path(
    outdir,'T2D_cterm_BMI_WHR_mvmr_res.Rds'
  )
)

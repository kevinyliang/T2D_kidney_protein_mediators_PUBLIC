#' Obtain the covariance matrix of BMI, T2D and WHR
#' 
#' Estimate pheno covariance using estimate_Syy which is correlation of the betas
#' 
#' done for INHCB cterm and INHCB full protein.
#' Results:
#'  Both are the same as expected since it is the same outcome
#' 
#' The allele assignment does not matter
library(tidyverse)
library(ggplot2)
library(data.table)
library(glue)
library(MVMR)
library(yaml)
library(metaCCA)

source("custom_mv_harmonize.R")

yml_file  <-  "configs/proj_config.yml"
args  <- read_yaml(yml_file)[['T2D_CKD_results']]

get_phenocov  <- function(mvmr_hres,outname){
  mvmr_package_hres  <- get_mvmr_inputs_from_twosampleMR(mvmr_hres)
  beta_cols  <- colnames(mvmr_package_hres)[
    grepl("_beta",colnames(mvmr_package_hres)) & colnames(mvmr_package_hres) != 'outcome_beta'
  ]
  se_cols  <- colnames(mvmr_package_hres)[
    grepl("_se",colnames(mvmr_package_hres)) & colnames(mvmr_package_hres) != 'outcome_se'
  ]

  S_XY  <- mvmr_package_hres[c("SNP",beta_cols,se_cols)]
  S_XY$allele_0  <- strsplit(
    S_XY$SNP,":"
  ) %>% lapply(
    X = .,
    FUN = function(x){toupper(x[length(x)-1])}
  ) %>% unlist()
  S_XY$allele_1  <- strsplit(
    S_XY$SNP,":"
  ) %>% lapply(
    X = .,
    FUN = function(x){toupper(x[length(x)])}
  ) %>% unlist()
  rownames(S_XY)  <- S_XY$SNP
  S_xy_trait_cols  <- c()
  for (i in seq_along(beta_cols)){
    S_xy_trait_cols <- c(
      S_xy_trait_cols,
      beta_cols[i],se_cols[i]
    )
  }
  S_XY  <- S_XY %>% dplyr::select(
    allele_0,allele_1,all_of(S_xy_trait_cols)
  )
  colnames(S_XY)  <- gsub("_beta","_b",colnames(S_XY))
  S_XY$allele_0  <- as_factor(S_XY$allele_0)
  S_XY$allele_1  <- as_factor(S_XY$allele_1)

  # estimate phenotype covar
  S_YY  <- estimateSyy(S_XY=S_XY)
  saveRDS(
    S_YY,
    outname
  )
}


bmi_t2d_whr_on_inhbc_mvmr_hres_cterm  <- readRDS("mvmr_hres_on_inhbc_cterm.Rds")
get_phenocov(bmi_t2d_whr_on_inhbc_mvmr_hres_cterm,'t2d_bmi_whr_phenocov_cterm.Rds')


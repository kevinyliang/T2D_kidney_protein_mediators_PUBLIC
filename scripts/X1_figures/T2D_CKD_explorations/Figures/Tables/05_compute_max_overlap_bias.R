#' Compute the bias for each MR result assuming maximum overlap of samples.
#' 
#' Uses the same formula for Yoshiji et al.
#' 
#' 

library(tidyverse)
library(ggplot2)
library(glue)
library(data.table)


#%% all f-stat res
all_sen_res  <- rbind(
  read.table(
    "all_sensitivity_res_T2D_and_prot.tsv",sep="\t",header=T
  ),read.table(
    "all_sensitivity_res_T2D_prot_and_disease.tsv",sep="\t",header=T
  )
)
#%%
#%%
deCODE_sample_size  <- 35559
# min is 21.5
min_decode_f_stat  <- all_sen_res %>% dplyr::filter(
  grepl("OL_",exposure)
) %>% pull(fstat) %>% na.omit() %>% min()
INHBC_fstat  <- 1373.795

#%%
#%%
sample_sizes  <- c(
  "T2D"=80154 +853816,
  "CKD"=41395 +439303,
  'eGFR'=1004040,
  "BUN"=852678
)

F_stat <-   c(
  "T2D"=all_sen_res %>% filter(exposure =='T2D') %>% pull(fstat) %>% min(),
  "Kidney_stones"=all_sen_res %>% filter(exposure =='kidney_stones') %>% pull(fstat) %>% min(),
  'eGFR_decline'=all_sen_res %>% filter(exposure =='eGFR_decline_unadj') %>% pull(fstat) %>% min(),
  'urate'=all_sen_res %>% filter(exposure =='urate') %>% pull(fstat) %>% min(),
  'eGFR'=all_sen_res %>% filter(exposure =='eGFR') %>% pull(fstat) %>% min(),
  "BUN"=all_sen_res %>% filter(exposure =='BUN') %>% pull(fstat) %>% min(),
  'kidney_cancer'=all_sen_res %>% filter(exposure =='kidney_cancer') %>% pull(fstat) %>% min(),
  'CKD'=all_sen_res %>% filter(exposure =='CKD') %>% pull(fstat) %>% min()
)
#%%


#%% compute max bias Disease on protein

decode_overlaps_proportion  <- ifelse(
  deCODE_sample_size <= sample_sizes,deCODE_sample_size/sample_sizes,sample_sizes/deCODE_sample_size
)

one_over_fstat  <- 1/F_stat

max_bias  <- decode_overlaps_proportion * one_over_fstat
max_bias*100

# 0.05%
((deCODE_sample_size/sample_sizes['T2D']) * 1/F_stat['T2D'])*100 
#%%

#%% compute max bias INHBC on disease

decode_overlaps_proportion  <- ifelse(
  deCODE_sample_size <= sample_sizes,deCODE_sample_size/sample_sizes,sample_sizes/deCODE_sample_size
)

one_over_fstat  <- 1/INHBC_fstat

# Max 0.005%
((deCODE_sample_size/sample_sizes['CKD']) * 1/INHBC_fstat)*100 
max_bias  <- decode_overlaps_proportion * one_over_fstat

#%%
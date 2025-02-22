#' Steigger test 
#' 
#' 
library(TwoSampleMR)
library(data.table)
library(glue)
library(yaml)

args  <- read_yaml("configs/proj_config.yml")[['common']]
protein  <- c(
  "OL_NA_SL_AGRN@15483_377_UP_O00468",
  "OL_NA_SL_CTSD@5508_62_UP_P07339",
  "OL_NA_SL_GNPTG@10666_7_UP_Q9UJJ9",
  "OL_NA_SL_INHBC@15686_49_UP_P55103",
  "OL_NA_SL_LPO@4801_13_UP_P22079"
)
protein_trait_pairs  <- list(
  "OL_NA_SL_AGRN@15483_377_UP_O00468" = c("eGFR"),
  "OL_NA_SL_CTSD@5508_62_UP_P07339" = c("eGFR"),
  "OL_NA_SL_INHBC@15686_49_UP_P55103" = c("eGFR","BUN"),
  "OL_NA_SL_LPO@4801_13_UP_P22079" = c("eGFR")
)


all_protein_on_trait_df  <- data.frame()
for (each_exposure in names(protein_trait_pairs)){
  prot_name  <- gsub("OL_NA_SL_|@.*$","",each_exposure)
  outcomes  <- protein_trait_pairs[[each_exposure]]
  for (each_outcome in outcomes){
    hres  <- readRDS(
        glue("{each_exposure}_MR_exposure.csv.gz_ON_{each_outcome}_processed_harmonized.Rds"
      )
    )
    
    steigger_res <- mr_steiger(
      p_exp = hres$pval.exposure,
      p_out = hres$pval.outcome,
      n_exp = hres$samplesize.exposure,
      n_out = hres$samplesize.outcome,
      r_exp = get_r_from_bsen(
        b = hres$beta.exposure,
        se = hres$se.exposure,
        n = hres$samplesize.exposure
      ),
      r_out = get_r_from_bsen(
        b = hres$beta.outcome,
        se = hres$se.outcome,
        n = hres$samplesize.outcome
      )
    )
    all_protein_on_trait_df  <- rbind(
      all_protein_on_trait_df,
      data.frame(
        exposure = prot_name,
        outcome = args$disease_name_plots[[each_outcome]],
        r2_exp = steigger_res$r2_exp,
        r2_out  = steigger_res$r2_out,
        r2_exp_adj = steigger_res$r2_exp_adj,
        r2_out_adj = steigger_res$r2_out_adj,
        correct_causal_direction = steigger_res$correct_causal_direction
      )
    )
  }
}




all_trait_on_protein_df <- data.frame()
for (each_outcome in protein){
  hres  <- readRDS(glue(
    "harmonized_disease_on_prot/T2D_processed.exposure_gwas.Rds_ON_{each_outcome}_MR_outcome.tsv.gz_harmonized.Rds"
  ))
  outcome  <- gsub("OL_NA_SL_|@.*$","",each_outcome)
  
  steigger_res <- mr_steiger(
    p_exp = hres$pval.exposure,
    p_out = hres$pval.outcome,
    n_exp = hres$samplesize.exposure,
    n_out = hres$samplesize.outcome,
    r_exp = get_r_from_lor(
      hres$beta.exposure,
      hres$eaf.exposure,
      ncase = 80154,
      ncontrol =  853816,
      prevalence = 0.092
    ),
    r_out = get_r_from_bsen(
      b = hres$beta.outcome,
      se = hres$se.outcome,
      n = hres$samplesize.outcome
    )
  )
  all_trait_on_protein_df  <- rbind(
    all_trait_on_protein_df,
     data.frame(
      exposure = 'T2D',
      outcome = outcome,
      r2_exp = steigger_res$r2_exp,
      r2_out  = steigger_res$r2_out,
      r2_exp_adj = steigger_res$r2_exp_adj,
      r2_out_adj = steigger_res$r2_out_adj,
      correct_causal_direction = steigger_res$correct_causal_direction
    )
  )
}

#%% GPNTG on CKD
#' CKD prevalence: https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(23)00570-3/fulltext#seccestitle140
#'  12.8% in west/east Europe
hres  <- readRDS(
  glue(
    "OL_NA_SL_GNPTG@10666_7_UP_Q9UJJ9_MR_exposure.csv.gz_ON_CKD_processed_harmonized.Rds"
  )
)


steigger_res <- mr_steiger(
  p_exp = hres$pval.exposure,
  p_out = hres$pval.outcome,
  n_exp = hres$samplesize.exposure,
  n_out = hres$samplesize.outcome,
  r_exp = get_r_from_bsen(
    b = hres$beta.exposure,
    se = hres$se.exposure,
    n = hres$samplesize.exposure
  ),
  r_out = get_r_from_lor(
    hres$beta.outcome,
    hres$eaf.outcome,
    ncase = 41395,
    ncontrol =  439303,
    prevalence = 0.128
  )
)
all_trait_on_protein_df  <- rbind(
  all_trait_on_protein_df,
    data.frame(
    exposure = 'GNPTG',
    outcome = unique(hres$outcome),
    r2_exp = steigger_res$r2_exp,
    r2_out  = steigger_res$r2_out,
    r2_exp_adj = steigger_res$r2_exp_adj,
    r2_out_adj = steigger_res$r2_out_adj,
    correct_causal_direction = steigger_res$correct_causal_direction
  )
)


all_steigger_res  <- rbind(
  all_trait_on_protein_df,
  all_protein_on_trait_df
)
numeric_cols <- c(
  "r2_exp","r2_out","r2_exp_adj","r2_out_adj"
)
for (each_col in numeric_cols){
  all_steigger_res[abs(all_steigger_res[,each_col])>= 0.0001,each_col] <- round(all_steigger_res[abs(all_steigger_res[,each_col])>= 0.0001,each_col],4)
  all_steigger_res[abs(all_steigger_res[,each_col]) < 0.0001 & abs(all_steigger_res[,each_col]) > 0,each_col] <- formatC(all_steigger_res[abs(all_steigger_res[,each_col]) < 0.0001 & abs(all_steigger_res[,each_col]) > 0,each_col],format='e',digits=3)
}



fwrite(
  all_steigger_res,
  "Tab_INHBC_Steigger_res.tsv",sep="\t",col.names=T,quote=F
)


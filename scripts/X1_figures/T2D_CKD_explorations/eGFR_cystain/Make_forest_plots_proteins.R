#' Forest plot of cystaatin based eGFR
library(ggforestplot)
library(yaml)
library(data.table)
library(glue)
library(tidyverse)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,"FigS4_eGFR_cystatin_creatinine.pdf")

args  <- read_yaml("configs/proj_config.yml")[['T2D_CKD_results']]
 

cystatin_out  <- "eGFR_cystatin_rep"

complications   <- c("eGFR_cystatin","eGFR")
trait_order  <- c("eGFR_cystatin","eGFR")
disease_name_plots  <- list(
  "eGFR_cystatin" = 'Cystatin based eGFR',
  "eGFR"='Creatinine based eGFR'
)
proteins  <- "OL_NA_SL_INHBC@15686_49_UP_P55103"
name  <- "OL_NA_SL_INHBC@15686_49_UP_P55103"


disease_on_prot_dir  <- "MR_results_disease_on_prot"
prot_on_disease_dir  <- "MR_results_prot_on_disease"

trait_prot_plot_data  <- rbind(
  readRDS(
    file.path(
      disease_on_prot_dir,
      glue(
        "eGFR_cystatin_processed.exposure_gwas.Rds_ON_{proteins}_MR_outcome.tsv.gz_MR_RES.Rds"
      )
    )
  ) %>% 
  dplyr::filter(
    method %in% c('Inverse variance weighted','Wald ratio')
  ) %>%
  dplyr::select(exposure,b,se,pval) %>% 
  mutate(
    direction = glue("eGFR on INHBC")
  ) %>% rename(
    disease = exposure
  ),
  readRDS(
    file.path(
      disease_on_prot_dir,
      glue(
        "eGFR_processed.exposure_gwas.Rds_ON_{proteins}_MR_outcome.tsv.gz_MR_RES.Rds"
      )
    )
  ) %>% 
  dplyr::filter(
    method %in% c('Inverse variance weighted','Wald ratio')
  ) %>%
  dplyr::select(exposure,b,se,pval) %>% 
  mutate(
    direction = glue("eGFR on INHBC")
  ) %>% rename(
    disease = exposure
  ),
  readRDS(
    file.path(
      prot_on_disease_dir,
      glue(
        "{proteins}_MR_exposure.csv.gz_ON_eGFR_cystatin_processed_MR_RES.Rds"
      )
    )
  ) %>% 
  dplyr::filter(
    method %in% c('Inverse variance weighted','Wald ratio')
  ) %>%
  dplyr::select(outcome,b,se,pval) %>% 
  mutate(
    direction = glue("INHBC on eGFR")
  ) %>% rename(
    disease = outcome
  ),
  readRDS(
    file.path(
      prot_on_disease_dir,
      glue(
        "{proteins}_MR_exposure.csv.gz_ON_eGFR_processed_MR_RES.Rds"
      )
    )
  ) %>% 
  dplyr::filter(
    method %in% c('Inverse variance weighted','Wald ratio')
  ) %>%
  dplyr::select(outcome,b,se,pval) %>% 
  mutate(
    direction = glue("INHBC on eGFR")
  ) %>% rename(
    disease = outcome
  )
)



trait_prot_plot_data$disease  <- unlist(disease_name_plots[trait_prot_plot_data$disease])

trait_prot_plot_data  <- trait_prot_plot_data[order(
  match(
    trait_prot_plot_data$direction,c("INHBC on eGFR","eGFR on INHBC")
  )
),]
trait_prot_plot_data$direction  <- as_factor(trait_prot_plot_data$direction)



plot_title  <- unlist(strsplit(unlist(strsplit(unlist(strsplit(name,"_SL_"))[2],"_UP_"))[1],"@"))


forestplot  <- ggforestplot::forestplot(
  df = trait_prot_plot_data,
  name = disease,
  estimate = b,
  pvalue = pval
) +
  ggforce::facet_col(
    facets = ~direction,
    scales = "free",
    space = "free"
  ) +
  theme_classic() +
  labs(
    x='Estimated causal effects of INHBC on eGFR'
  ) +
  theme(text=element_text(size=20),legend.position='none',plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'))

ggsave(
  filename = tfile,
  device='png',width=9,height=5,units='in',
  plot = forestplot
)
ggsave(
  filename = final_file,
  device='pdf',width=9,height=5,units='in',
  plot = forestplot
)



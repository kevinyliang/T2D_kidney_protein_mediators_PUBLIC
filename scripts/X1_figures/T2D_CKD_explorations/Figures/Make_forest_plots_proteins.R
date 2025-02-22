#' Protein MR foorest plots
library(ggforestplot)
library(yaml)
library(data.table)
library(glue)
library(tidyverse)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')


args  <- read_yaml("configs/proj_config.yml")[['T2D_CKD_results']]
 
outdir_base  <- "figures"

all_disease_on_prot_results  <- read.table(
  "all_disease_on_T2D_prot_mr_results.tsv",sep="\t",header=T
) %>% filter(exposure %in% names(args$disease_name_plots))


all_prot_on_disease_results  <- read.table(
  "all_T2D_prot_on_disease_mr_results.tsv",sep="\t",header=T
) %>% filter(outcome %in% names(args$disease_name_plots))


all_t2d_on_prot_results  <- read.table(
  "all_T2D_on_prot_mr_results.tsv",sep="\t",header=T
)
all_prot_on_t2d_results  <- read.table(
  "all_prot_on_T2D_mr_results.tsv",sep="\t",header=T
)

all_disease_on_prot_results$bonferroni_p  <- p.adjust(all_disease_on_prot_results$pval,method='bonferroni')
all_prot_on_disease_results$bonferroni_p  <- p.adjust(all_prot_on_disease_results$pval,method='bonferroni')
all_t2d_on_prot_results$bonferroni_p  <- p.adjust(all_t2d_on_prot_results$pval,method='bonferroni')
all_prot_on_t2d_results$bonferroni_p  <- p.adjust(all_prot_on_t2d_results$pval,method='bonferroni')



complications   <- c('eGFR','BUN','urate')
disease  <- "T2D"
order  <- c(disease,complications)

proteins  <- "OL_NA_SL_INHBC@15686_49_UP_P55103"
name  <- "OL_NA_SL_INHBC@15686_49_UP_P55103"


trait_prot_plot_data  <- rbind(
  all_disease_on_prot_results %>% filter(
    exposure %in% complications &
    outcome == proteins
  ) %>% 
  dplyr::select(exposure,b,se,bonferroni_p) %>% 
  mutate(
    direction = glue("disease on Protein")
  ) %>% rename(
    disease = exposure
  ),
  all_t2d_on_prot_results %>% filter(
    exposure == disease &
    outcome == proteins
  ) %>% 
  dplyr::select(exposure,b,se,bonferroni_p) %>% 
  mutate(
    direction = glue("disease on Protein")
  ) %>% rename(
    disease = exposure
  ),
  all_prot_on_disease_results %>% filter(
    exposure == proteins,
    outcome %in% complications
  ) %>% 
  dplyr::select(outcome,b,se,bonferroni_p) %>% 
  mutate(
    direction = glue("Protein on disease")
  ) %>% rename(
    disease = outcome
  ),
  all_prot_on_t2d_results %>% filter(
    exposure == proteins,
    outcome == disease
  ) %>% 
  dplyr::select(outcome,b,se,bonferroni_p) %>% 
  mutate(
    direction = glue("Protein on disease")
  ) %>% rename(
    disease = outcome
  )
)

trait_prot_plot_data  <- rbind(
  trait_prot_plot_data[trait_prot_plot_data$disease == disease,],
  trait_prot_plot_data %>% filter(disease %in% complications)
)

trait_prot_plot_data  <- trait_prot_plot_data[
  order(
    match(
      trait_prot_plot_data$disease,
      order
    )
  ),]

trait_prot_plot_data$disease  <- as_factor(unlist(args$disease_name_plots[trait_prot_plot_data$disease]))



plot_title  <- unlist(strsplit(unlist(strsplit(unlist(strsplit(name,"_SL_"))[2],"_UP_"))[1],"@"))
plot_title  <- glue("{plot_title[1]}\nSomalogic seqID: {plot_title[2]}")

forestplot  <- ggforestplot::forestplot(
  df = trait_prot_plot_data,
  name = disease,
  estimate = b,
  pvalue = bonferroni_p,
  colour= direction
) +
  ggforce::facet_col(
    facets = ~direction,
    scales = "free",
    space = "free"
  ) + 
  theme_classic() +
  labs(
    x='Effect size',
    title=plot_title
  ) +
  theme(text=element_text(size=20))

ggsave(
  filename = tfile,
  device='png',width=13,height=8,units='in',
  plot = forestplot
)


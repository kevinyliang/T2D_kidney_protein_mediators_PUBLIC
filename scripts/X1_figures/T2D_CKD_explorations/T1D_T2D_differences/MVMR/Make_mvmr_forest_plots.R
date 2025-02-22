library(ggforestplot)
library(yaml)
library(data.table)
library(glue)
library(tidyverse)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,'Fig5b_MVMR_t2d_t1d_bmi_whr.pdf')

args  <- read_yaml("configs/proj_config.yml")[['T2D_CKD_results']]



mvmr_res  <- readRDS(
  "t1d_T2d_bmi_whr_on_inhbc_mvmr_mr.Rds"
)$result
forest_plot_data  <- data.frame(
  b = mvmr_res$b,
  exposure = unlist(args$all_disease_name_map[mvmr_res$exposure]),
  se = mvmr_res$se,
  p = mvmr_res$pval
)

trait_order  <- c("T2D","T1D","WHR",'bmi')
order <- unlist(args$all_disease_name_map[trait_order])
forest_plot_data  <- forest_plot_data[
  order(
    match(
      forest_plot_data$exposure,
      order
    )
  ),]

p  <- ggforestplot::forestplot(
  forest_plot_data,
  estimate = b,
  se = se,
  pvalue = p,
  name = exposure,
  colour = exposure
) + theme_classic() +
labs(
  x='MVMR causal estimates'
) +
scale_x_continuous(expand = expansion(mult=c(0.1,0.2)))+
theme(
  text = element_text(size=20),
  legend.position='none'
) 

ggsave(
  tfile,
  p,
  device='png',
  width=6,height=3,units='in'
)
ggsave(
  final_file,
  p,
  device='pdf',
  width=6,height=3,units='in'
)






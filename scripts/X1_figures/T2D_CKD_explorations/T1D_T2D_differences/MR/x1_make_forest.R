library(ggforestplot)
library(yaml)
library(data.table)
library(glue)
library(tidyverse)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,"Fig5a_uniMR_t2d_t1d_bmi_whr.pdf")


args  <- read_yaml("proj_config.yml")[['T2D_CKD_results']]


t2d_res  <- readRDS("t2d_inhbc_mr.Rds") %>% filter(method=="Inverse variance weighted")
t1d_res  <- readRDS("t1d_inhbc_mr.Rds") %>% filter(method=="Inverse variance weighted")
bmi_res  <- readRDS("bmi_inhbc_mr.Rds") %>% filter(method=="Inverse variance weighted")
whr_res  <- readRDS("whr_inhbc_mr.Rds") %>% filter(method=="Inverse variance weighted")


name_map  <- c(
  "T2D"="T2D",
  "T1D"="T1D",
  "bmi"="BMI",
  "WHR"="Wasit to hip ratio"
)
forest_plot_data <- reduce(
  .x = list(t2d_res,t1d_res,bmi_res,whr_res),
  .f = rbind
)
forest_plot_data$exposure  <- name_map[forest_plot_data$exposure]



trait_order  <- c("T2D","T1D","WHR",'bmi')
order <- unlist(name_map[trait_order])
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
  pvalue = pval,
  name = exposure,
  colour = exposure
) + theme_classic() +
labs(
  x='Univariable MR causal estimate'
) +
theme(text = element_text(size=20),legend.position='none') 

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
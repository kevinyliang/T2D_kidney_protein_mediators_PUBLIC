# look at protein associations
#' Filters for:
#'  invalid batches
#'  Europeans only
#'  batches 1-6 only

library(tidyverse)
library(glue)
library(ggplot2)
library(data.table)
library(ggforestplot)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')


library(tidyverse)
library(glue)
library(data.table)
library(yaml)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')

yml_file  <-  "configs/proj_config.yml"

args  <- read_yaml(yml_file)[['T2D_CKD_results']]
protein_name  <- 'inhbc'
outdir = file.path(
  args$rt_dir,'results','UKB_disease_observations','INHBC_T1D_T2D'
)
km_info <- read.table(
  file.path(args$rt_dir,'results','UKB_disease_observations','KM_information.tsv'),
  sep="\t",header=T,quote="",colClasses=c("eid"="character")
)
T2D_counts  <- sum(km_info$T2D == 1)
T1D_counts  <- sum(km_info$T1D == 1)
non_diabetics  <- sum(km_info$T1D == 0 & km_info$T2D == 0)
km_info  <- km_info %>% dplyr::mutate(
  label = ifelse(
    T2D == 1, glue("T2D (N = {T2D_counts})"),
    ifelse(
      T1D == 1, glue("T1D (N = {T1D_counts})"),glue("Non diabetics (N = {non_diabetics})")
    )
  )
)
scaled_km_info  <- km_info %>% dplyr::mutate(
  inhbc = scale(inhbc),
  BMI = scale(BMI)
)


#%%

r_t2d_t1d  <- t.test(
  km_info %>% filter(T1D==1) %>% pull(inhbc),
  km_info %>% filter(T2D==1) %>% pull(inhbc)
)
r_t2d_non  <- t.test(
  km_info %>% filter(T1D==0 & T2D == 0) %>% pull(inhbc),
  km_info %>% filter(T2D==1) %>% pull(inhbc)
)
r_t1d_non  <- t.test(
  km_info %>% filter(T1D==0 & T2D == 0) %>% pull(inhbc),
  km_info %>% filter(T1D==1) %>% pull(inhbc)
)

p  <- ggplot(
    km_info
  ) + 
  geom_boxplot(
    aes(x=label,y=inhbc,fill=label)
  )+
  theme_classic() +
  theme(
    text = element_text(size=25),legend.position='none'
  )+
  xlab("") +
  labs(
    y = 'INHBC (NPX levels)'
  )
ggsave(
  tfile,
  plot = p,
  device='png',width=10,height=10,units='in'
)


#%%


# look at protein associations
#' Filters for:
#'  invalid batches
#'  Europeans only
#'  batches 1-6 only

library(tidyverse)
library(yaml)
library(glue)
library(ggplot2)
library(data.table)
library(ggforestplot)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file_a <- file.path(final_dir,"FigS5a_observational_t2d_inhbc.pdf")
final_file_b <- file.path(final_dir,"FigS5b_observational_t2d_inhbc.pdf")

yml_file  <-  "proj_config.yml"

args  <- read_yaml(yml_file)[['T2D_CKD_results']]
protein_name  <- 'inhbc'
outdir = file.path(
  args$rt_dir,'results','UKB_disease_observations','INHBC'
)
km_info <- read.table(
  file.path(args$rt_dir,'results','UKB_disease_observations',glue('KM_information_{toupper(protein_name)}_ALL_Ancestries.tsv')),
  sep="\t",header=T,quote="",colClasses=c("eid"="character")
)
dia  <- sum(km_info$T2D==1)
non_dia  <- sum(km_info$T2D == 0)

km_info$label  <- ifelse(
  km_info$T2D==1,glue('T2D (N = {dia})'),glue('Non-T2D (N = {non_dia})')
)
scaled_km_info  <- km_info %>% dplyr::mutate(
  inhbc = scale(inhbc),
  BMI = scale(BMI)
)
scaled_km_info_w_pc <- scaled_km_info


#%%

r  <- t.test(
  km_info %>% filter(T2D==1) %>% pull(inhbc),
  km_info %>% filter(T2D==0) %>% pull(inhbc)
)

p  <- ggplot(
    km_info
  ) + 
  geom_violin(
    aes(x=label,y=inhbc,fill=label)
  )+
  scale_fill_manual(
    values = setNames(
      c("#f35f5b","#17b4b7"),
      nm=c(
        glue('T2D (N = {dia})'),
        glue('Non-T2D (N = {non_dia})')
      )
    )
  )+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  theme_classic() +
  theme(
    text = element_text(size=25),
    legend.position='none'
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
ggsave(
  final_file_a,
  plot = p,
  device='pdf',width=10,height=10,units='in'
)


#%% linear regression
# pc_var  <- paste0("p22009_a",seq(1,10),collapse="+")
t2d_on_inhbc  <- lm(
  formula = as.formula(glue("inhbc ~ T2D + Age_at_recruitment + Genetic_sex + ancestry")),
  data = scaled_km_info_w_pc
) %>% summary()
forest_plot_data  <- data.frame(
  b = c(
    t2d_on_inhbc$coefficients['T2D','Estimate'],
    t2d_on_inhbc$coefficients['Age_at_recruitment','Estimate'],
    t2d_on_inhbc$coefficients['Genetic_sex','Estimate']
  ),
  se = c(
    t2d_on_inhbc$coefficients['T2D','Std. Error'],
    t2d_on_inhbc$coefficients['Age_at_recruitment','Std. Error'],
    t2d_on_inhbc$coefficients['Genetic_sex','Std. Error']
  ),
  label = c(
    "T2D",
    "Age at recruitment",
    "Genetically determined sex (male)"
  )
)
p  <- ggforestplot::forestplot(
  df = forest_plot_data,
  name = label,
  se = se,
  estimate = b,
  colour = label
) + theme_classic() + theme(text = element_text(size=20),legend.position='none') +
labs(
  x='change in INHBC per\nchange in predictor'
)
ggsave(
  tfile,
  p,
  device='png',width=8,height=2,units='in'
)
ggsave(
  final_file_b,
  p,
  device='pdf',width=8,height=2,units='in'
)

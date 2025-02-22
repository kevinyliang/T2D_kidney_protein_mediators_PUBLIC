
library(tidyverse)
library(glue)
library(data.table)
library(yaml)
library(lubridate)
library(ggsurvfit)
library(survival)
library(gtsummary)
library(tidycmprsk)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,"Fig4b_cox_regression.pdf")

yml_file  <-  "configs/proj_config.yml"

args  <- read_yaml(yml_file)[['T2D_CKD_results']]
protein_name  <- 'inhbc'
outdir = file.path(
  args$rt_dir,'results','UKB_disease_observations',toupper(protein_name)
)
km_info <- read.table(
  file.path(args$rt_dir,'results','UKB_disease_observations',glue('KM_information_{toupper(protein_name)}_ALL_Ancestries.tsv')),
  sep="\t",header=T,quote="",colClasses=c("eid"="character")
)


WHR_info  <- read.table(
  "WHR_information_Mar9_2024.tsv",sep="\t",quote="",header=T,colClasses = c("FID"='character')
)
km_info  <- inner_join(
  km_info,
  WHR_info,
  # by=c("eid"="FID")
  by=c('eid'="FID")
)

km_info  <- km_info %>% dplyr::mutate(
  sex_plot = ifelse(Genetic_sex == 1,"Male","Female")
) %>% mutate(
  WHR = scale(WHR),
  HbA1c = scale(HbA1c),
  BMI = scale(BMI),
  !!sym(protein_name) := scale(!!sym(protein_name)),
  Genetic_sex = as_factor(Genetic_sex)
)
km_info_w_PC <- km_info


#%% fit cox regression

cox_model  <- coxph(
  as.formula(glue("Surv(time_of_followup, disease_status) ~ WHR + BMI + T2D + {protein_name} + sex_plot + Age_at_recruitment + T2D*{protein_name} + ancestry")), data = km_info_w_PC) %>% summary()
forestplot_df  <- data.frame(
  name = rownames(cox_model$coefficients),
  estimate = cox_model$coefficients[,'coef'],
  se = cox_model$coefficients[,'se(coef)'],
  p = cox_model$coefficients[,'Pr(>|z|)']
)
forestplot_df <- forestplot_df %>% filter(
  name %in% c("WHR","BMI","T2D","inhbc","sex_plotMale","Age_at_recruitment","T2D:inhbc")
)

name_map  <- setNames(
  c('WHR','BMI','T2D',toupper(protein_name),'Male','Age at recruitment',glue('T2D * {toupper(protein_name)}')),
  nm = c(
    'WHR','BMI',"T2D",protein_name,'sex_plotMale','Age_at_recruitment',glue('T2D:{protein_name}')
  )
)
forestplot_df$name  <- name_map[forestplot_df$name]
p  <- ggforestplot::forestplot(
  forestplot_df,
  name = name,
  estimate=estimate,
  se=se,
  pvalue=p,
  colour=name,logodds=T
) + 
scale_color_manual(
  values=c(
    "WHR" = "#8d0922",
    "BMI" = "#4178bc",
    "T2D" = "#f76f00",
    "INHBC" = "#7b6dff",
    "Male" = "#17b5b3",
    "Age at recruitment" = "#ff0011",
    "T2D * INHBC" = "#252e48"
  )
)+
theme_classic() + theme(
  text =element_text(size=20),legend.position='none'
) + labs(
  x='Hazard ratio of standardized risk factors',
  y='Risk factors'
) +xlim(c(NA,4.2))

ggsave(
  tfile,p,
  device='png',height=5,width=8,units='in'
)
ggsave(
  final_file,
  p,
  device='pdf',height=5,width=8,units='in'
)
#%% Create table
forestplot_df  <- forestplot_df %>% 
  mutate(
    lower_ci  = estimate - 1.96 * se,
    upper_ci  = estimate + 1.96 * se
  )

nmap <- c(
  "name"="Predictor",
  "estimate"="ln(Hazard ratio)",
  "se"="SE",
  "p"="P",
  "lower_ci"="Lower CI",
  "upper_ci"="Upper CI"
)

colnames(forestplot_df) <- nmap[colnames(forestplot_df)]


numeric_cols <- c(
  "ln(Hazard ratio)","SE","P","Lower CI","Upper CI"
)

for (each_col in numeric_cols){
  forestplot_df[abs(forestplot_df[,each_col])>= 0.0001 & !is.na(forestplot_df[,each_col]),each_col] <- round(forestplot_df[abs(forestplot_df[,each_col])>= 0.0001 & !is.na(forestplot_df[,each_col]),each_col],4)
  forestplot_df[abs(forestplot_df[,each_col]) < 0.0001 & abs(forestplot_df[,each_col]) > 0 & !is.na(forestplot_df[,each_col]),each_col] <- formatC(forestplot_df[abs(forestplot_df[,each_col]) < 0.0001 & abs(forestplot_df[,each_col]) > 0 & !is.na(forestplot_df[,each_col]),each_col],format='e',digits=3)
}







fwrite(
  forestplot_df,
  "cox_regression_results.tsv",
  sep="\t",col.names=T,row.names=F,quote=F
)


#%%
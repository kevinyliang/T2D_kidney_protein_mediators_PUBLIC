library(tidyverse)
library(glue)
library(yaml)
library(ggforestplot)
library(data.table)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,"Fig5c_observational_t2d_t1d_bmi_whr.pdf")

yml_file <- "configs/proj_config.yml"

args  <- read_yaml(yml_file)[['T2D_CKD_results']]

protein <- 'inhbc'
outdir  <-  glue('UKB_disease_observations')

km_information  <- read.table(
  glue('UKB_disease_observations/KM_information_{toupper(protein)}.tsv'),header=T,sep="\t",quote="",colClasses=c("eid"="character")
)%>% rename(FID=eid)

whr_info  <- read.table(
  "UKB_meta_information/WHR_information_Mar9_2024.tsv",sep="\t",quote="",header=T,colClasses=c("FID"="character")
)
prs_info  <- read.table(
  "UKB_diabetes_PRS_Feb27_2024_participant.tsv",sep="\t",quote="",header=T,
  colClasses= c("Participant.ID"='character')
) %>% rename(
  FID = Participant.ID,
  Standard_T2D_PRS = Standard.PRS.for.type.2.diabetes..T2D.,
  Standard_T1D_PRS = Standard.PRS.for.type.1.diabetes..T1D.
)
km_information  <- inner_join(
  km_information,
  whr_info,by=c("FID")
) %>% inner_join(
  .,prs_info,by=c("FID")
)%>%dplyr::select(
  Standard_T2D_PRS,Standard_T1D_PRS,BMI,WHR,!!sym(protein)
)  %>% na.omit()
km_information_scale  <- scale(km_information) %>% as.data.frame()

#%%
adj  <- lm(
  formula = as.formula(
    glue("{protein} ~ Standard_T2D_PRS + Standard_T1D_PRS + BMI + WHR")
  ),
  data = km_information_scale
) %>% summary()

forestplot_data  <- data.frame(
  b = c(
    adj$coefficients['Standard_T2D_PRS','Estimate'],
    adj$coefficients['Standard_T1D_PRS','Estimate'],
    adj$coefficients['WHR','Estimate'],
    adj$coefficients['BMI','Estimate']
  ),
  se = c(
    adj$coefficients['Standard_T2D_PRS','Std. Error'],
    adj$coefficients['Standard_T1D_PRS','Std. Error'],
    adj$coefficients['WHR','Std. Error'],
    adj$coefficients['BMI','Std. Error']
  ),
  p = c(
    adj$coefficients['Standard_T2D_PRS','Pr(>|t|)'],
    adj$coefficients['Standard_T1D_PRS','Pr(>|t|)'],
    adj$coefficients['WHR','Pr(>|t|)'],
    adj$coefficients['BMI','Pr(>|t|)']
  ),
  model = c(
    "T2D PRS",
    "T1D PRS",
    "Waist to hip ratio",
    "BMI"
  )
)

p  <- ggforestplot::forestplot(
  forestplot_data,
  estimate=b,
  se = se,
  pvalue=p,
  name=model,
  colour=model
) + theme_classic() +
theme(text = element_text(size=20),legend.position='none') +
labs(
  x='Multivariable \n linear regression estimate'
) +
scale_color_manual(
  values = c(
    "T2D PRS" = "#7b6dff",
    "T1D PRS" = "#29bab8",
    "Waist to hip ratio" = "#fb303f",
    "BMI" = "#40475e"
  )
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
#%%





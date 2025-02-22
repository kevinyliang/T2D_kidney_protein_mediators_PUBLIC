
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
final_dir <- "figures"
final_file <- file.path(final_dir,'Fig4a_km_curve.pdf')

yml_file  <-  "configs/proj_config.yml"

args  <- read_yaml(yml_file)[['T2D_CKD_results']]
protein_name  <- 'inhbc'


km_info <- read.table(
  file.path(args$rt_dir,'results','UKB_disease_observations',glue('KM_information_{toupper(protein_name)}_ALL_Ancestries.tsv')),
  sep="\t",header=T,quote="",colClasses=c("eid"="character")
)
#%% create groups
protein_quantiles  <- 0.1
top_boundary  <- quantile(km_info%>%pull(!!sym(protein_name)),1-protein_quantiles)
bottom_boundary  <- quantile(km_info%>%pull(!!sym(protein_name)),protein_quantiles)
assertthat::assert_that(
  top_boundary > bottom_boundary
)
km_info  <- km_info %>% mutate(
  groups = ifelse(
    T2D == 1 & !!sym(protein_name) >= top_boundary,glue('T2D, high {toupper(protein_name)}'), ifelse(
      T2D == 1 & !!(sym(protein_name)) <= bottom_boundary,glue("T2D, low {toupper(protein_name)}"), ifelse(
        T2D == 0 & !!sym(protein_name) >= top_boundary, glue("Non-T2D, high {toupper(protein_name)}"), ifelse(
          T2D == 0 & !!sym(protein_name) <= bottom_boundary,glue('Non-T2D, low {toupper(protein_name)}'),NA
        )
      )
    )
  )
)


#%%
p_data  <- km_info %>% filter(!is.na(groups))
p_data  <- p_data[
  order(match(
    p_data$groups,c(
      glue("Non-T2D, low {toupper(protein_name)}"),
      glue("Non-T2D, high {toupper(protein_name)}"),
      glue("T2D, low {toupper(protein_name)}"),
      glue("T2D, high {toupper(protein_name)}")
    )
  )),
]
p_data$groups  <- as_factor(p_data$groups)

p  <- survfit2(
  Surv(time_of_followup, disease_status) ~ groups, data = p_data) %>% 
  ggsurvfit()+
  labs(
    x = "Days of followup",
    y = "Proportion of individuals without kidney diseases"
  ) + 
  add_confidence_interval() +
  add_censor_mark(size = 2, alpha = 0.2) +
  
  theme_classic() + 
  theme(
    plot.margin=unit(c(0.5,0.5,0.5,0.5),'cm'),
    text=element_text(size=20),
    legend.position = 'inside',
    legend.position.inside=c(0.2,0.3)
  ) 

ggsave(
  tfile,
  p,
  device='png',height=7,width=8,units='in'
)
ggsave(
  final_file,
  p,
  device='pdf',height=7,width=8,units='in'
)
#%%





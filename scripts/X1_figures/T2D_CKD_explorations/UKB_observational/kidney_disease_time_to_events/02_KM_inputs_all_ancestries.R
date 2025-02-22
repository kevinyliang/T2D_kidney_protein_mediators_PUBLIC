#' Create input for all observational studies (mainly survival analysis, cox regression, linear regression)
#' 
#' Inclusion criteria:
#'  Prevalent T2D, T1D cases only (i.e., incident T2D or T1D individuals are removed).
#'  Incident Kidney disease only (i.e., prevalent kidney disease removed).
#' 
#' Exclusion criteria:
#'  Diagnosed with T2D and T1D
#'  Have liver cirrhosis
#'  Without INHBC measurements
#'  Belong in batch 0 or 7 of UKB-PPP proteomics release
#'  Withdrawn individuals
#'  Without blood draw dates
#'  Individuals with T1D, T2D, or kidney diseases (i.e., the diseases considered), but without valid diagnosis dates
#'  Individuals without T1D or T2D but have HbA1c levels above 48 mmol/mol (6.5%)
#'
#' Kidney disease defintion prevalent:
#'  Main: ICD-10/OPCS-4
#'  secondary: eGFR < 60, uACR >= 3 mg/mmol
#' 
#' Kidney disease definition incident:
#'  ICD-10/OPCS-4 codes
#' 
#' T2D definition (Incident/Prevalent):
#'  ICD-10 codes.
#' 
#' 


library(glue)
library(data.table)
library(yaml)
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(tidyverse)


yml_file  <-  "configs/proj_config.yml"

args  <- read_yaml(yml_file)[['T2D_CKD_results']]

outdir = file.path(
  args$rt_dir,'results','UKB_disease_observations'
)
system(command = glue("mkdir -p {outdir}"))
protein_name = 'inhbc'

#%% read in covar
ukb_covar_info  <- read.table(
  args$ukb_covars,quote="",header=T,sep="\t",colClasses=c(
    'FID' = 'character'
  )
)
all_blood_draw_dates  <- read.table(
  "all_blood_draw_dates_formatted.tsv",sep="\t",header=T,quote="",colClasses=c("eid"="character")
)
olink_data <- read.table(
  "data_olink_results_i0_Dec2023.csv",sep=",",quote="",header=T,colClasses=setNames(c("character","numeric"),c("eid",protein_name))
) %>% dplyr::select(eid,!!sym(protein_name)) %>% na.omit()

batch_info = read.table(
  "olink_batch_number.dat",sep="\t",quote="",header=T,colClasses=c("Batch"="numeric")
)
plate_info = read.table(
  "UKBPPP_participant_info_participant.tsv",sep="\t",quote="",header=T,
  colClasses=c("eid" = 'character')
)
sample_plate_info = dplyr::inner_join(
  plate_info,batch_info,by=c("p30901_i0"="PlateID")
)

europeans_ind  <- read.table(
  "ukb.eurIDsPCA.plink.txt",sep=" ",header=F,quote="",colClasses=c("V1"="character","V2" = 'character')
) %>%rename(FID = V1, IID = V2) %>% filter(!grepl("^-",FID))%>% mutate(
  ancestry = "EUR"
)
afr_ind  <- read.table(
  "ukb.afrIDsPCA.plink.txt",sep=" ",header=F,quote="",colClasses=c("V1"="character","V2" = 'character')
) %>%rename(FID = V1, IID = V2) %>% filter(!grepl("^-",FID))%>% mutate(
  ancestry = "AFR"
)
amr_ind  <- read.table(
  "ukb.amrIDsPCA.plink.txt",sep=" ",header=F,quote="",colClasses=c("V1"="character","V2" = 'character')
) %>%rename(FID = V1, IID = V2) %>% filter(!grepl("^-",FID))%>% mutate(
  ancestry = "AMR"
)
eas_ind  <- read.table(
  "ukb.easIDsPCA.plink.txt",sep=" ",header=F,quote="",colClasses=c("V1"="character","V2" = 'character')
) %>%rename(FID = V1, IID = V2) %>% filter(!grepl("^-",FID))%>% mutate(
  ancestry = "EAS"
)
sas_ind  <- read.table(
  "ukb.sasIDsPCA.plink.txt",sep=" ",header=F,quote="",colClasses=c("V1"="character","V2" = 'character')
) %>%rename(FID = V1, IID = V2) %>% filter(!grepl("^-",FID)) %>% mutate(
  ancestry = "SAS"
)


all_sample_ancestry <- rbind(
  europeans_ind,
  sas_ind,
  eas_ind,
  afr_ind,
  amr_ind
)


#%% filter samples based on batch, withdrawn and european
keep_ind = sample_plate_info %>% dplyr::filter(
  Batch != 0 & Batch != 7
)
withdrawn  <- read.table(
      "UKB_withdrawal.csv",sep="\t",header=F,quote="",colClasses=c("V1"="character")
    ) %>% pull(V1)
keep_ind  <- keep_ind %>% filter(!eid %in% withdrawn)


#%% read in kidney disease data and format date
#'
all_kidney_disease_incidence = read.table(
  "kidney_diseases_refined_definitions_status.txt.gz",sep="\t",header=T,quote="",colClasses=c("FID"='character')
)
all_kidney_disease_incidence$earliest_date <- 
  as_datetime(
    all_kidney_disease_incidence$earliest_date,format="%Y-%m-%d"
  )
all_kidney_disease_incidence$blood_draw_date <- 
  as_datetime(
    all_kidney_disease_incidence$blood_draw_date,format="%Y-%m-%d"
  )

#%% remove prevalent CKD based on eGFR data
#' CKD threshold: https://onlinelibrary.wiley.com/doi/epdf/10.1002/jcsm.12705
#' https://www.medrxiv.org/content/10.1101/2023.12.13.23299901v1.full.pdf

# identify additional prevalent CKD cases
eGFR_threshold  <- 60 # less than 60 is case
baseline_eGFR_data  <- read.table("all_UKB_eGFRcr_results.tsv",sep="\t",header=T,quote="",colClasses=c("FID"='character'))
eGFR_based_baseline_CKD  <- baseline_eGFR_data %>% filter(eGFR < eGFR_threshold)
assertthat::assert_that(
  max(eGFR_based_baseline_CKD$eGFR) < 60
)

# all samples that were 'no kidney disease', but have eGFR less than 60 at baseline (rounded to the nearest whole number) are considered prevalent caes and are removed
eGFR_based_prevalent_to_removed  <- all_kidney_disease_incidence %>% filter(disease_status == 0 & FID %in% eGFR_based_baseline_CKD$FID)
all_kidney_disease_incidence  <- all_kidney_disease_incidence %>% filter(
  ! FID %in% eGFR_based_prevalent_to_removed$FID
)
check_df  <- inner_join(
  all_kidney_disease_incidence %>% filter(disease_status == 0),
  baseline_eGFR_data,
  by = c("FID")
)
assertthat::assert_that(
  min(check_df$eGFR) >= 60
)
rm(check_df)

#%% remove prevalent cases based on urinary albumin-creatinine ratio
#' https://www.agappe.com/swiss_en/blog-details/albumin-microalbumin-prealbumin.html
#' Notes:
#'  microalbumin units: mg/L
#'  creatinine units:umol/L
#' 
#' Therefore: uACR (mg/mmol) = ([microalbumin mg/L] /[creatinine umol/L]) *(1umol/0.001mmol)
#' 
uACR_threshold  <- 3 # units mg/mmol # greater than 3 is case
UKB_uACR_data  <- read.table("UKB_urinary_ACR_data_participant_June21_2024.tsv",sep="\t",header=T,colClasses=c('Participant.ID'='character')) %>% dplyr::rename(
  microalbumin_urine_i0 = Microalbumin.in.urine...Instance.0,
  creatinine_urine_i0 = Creatinine..enzymatic..in.urine...Instance.0
) %>% dplyr::select(Participant.ID,microalbumin_urine_i0,creatinine_urine_i0)
UKB_uACR_data  <- UKB_uACR_data %>% filter(!is.na(microalbumin_urine_i0) & !is.na(creatinine_urine_i0))

UKB_uACR_data$uACR  <- UKB_uACR_data$microalbumin_urine_i0/UKB_uACR_data$creatinine_urine_i0 / 0.001

uACR_based_CKD  <- UKB_uACR_data %>% filter(uACR >= 3)

uACR_based_CKD_to_remove  <- all_kidney_disease_incidence %>% filter(
  disease_status == 0 & FID %in% uACR_based_CKD$Participant.ID
)
# UKB_uACR_data %>% filter(Participant.ID %in% uACR_based_CKD_to_remove$FID) %>% View()
all_kidney_disease_incidence  <- all_kidney_disease_incidence %>% filter(
  ! FID %in% uACR_based_CKD_to_remove$FID
)
check_df  <- inner_join(
  all_kidney_disease_incidence %>% filter(disease_status == 0),
  UKB_uACR_data,
  by = c("FID"='Participant.ID')
)
assertthat::assert_that(
  max(check_df$uACR) < 3
)
rm(check_df)
#%% re classify those who died with kidney related diagnoses as incident
death_w_kidney_related_diseases  <- read.table("deaths_with_kidney_related_diagnoses.txt",sep="\t",header=T,quote="",colClasses=c("FID"='character'))


all_kidney_disease_incidence  <- all_kidney_disease_incidence %>% 
  mutate(
    disease_status_new = ifelse(
      FID %in% death_w_kidney_related_diseases$FID,1,disease_status
    )
  )
# make sure the only individual changed are those we expect
changed_df  <- all_kidney_disease_incidence %>% filter(
  disease_status != disease_status_new
)
assertthat::assert_that(
  all(changed_df$FID %in% death_w_kidney_related_diseases$FID)
)
all_kidney_disease_incidence <- all_kidney_disease_incidence %>% 
  dplyr::select(-disease_status) %>% rename(
    disease_status = disease_status_new
  )
assertthat::assert_that(
  all(all_kidney_disease_incidence %>% filter(FID %in% death_w_kidney_related_diseases$FID) %>% pull(disease_status) == 1)
)

#%% remove individuals with any liver diseases
liver_disease_incidence  <- read.table(
  "liver_disease_status.txt.gz",sep="\t",header=T,quote="",colClasses = c("FID"="character")
) %>% filter(disease_status == 1)
liver_disease_prevalence  <- read.table(
  "liver_disease_status.txt.gz",sep="\t",header=T,quote="",colClasses=c("FID"='character')
) %>% filter(disease_status == 1)

deaths_w_liver_disease  <- read.table("deaths_with_liver_related_diagnoses.txt",header=T,sep="\t",quote="",colClasses=c("FID"="character"))
liver_disease_to_remove  <- unique(
  c(liver_disease_incidence$FID,liver_disease_prevalence$FID,deaths_w_liver_disease$FID)
)


all_kidney_disease_incidence <- all_kidney_disease_incidence %>% filter(
  ! FID %in% liver_disease_to_remove
)


#%% read in prevalence T2D data and format dates
all_T2D_prevalence  <- read.table(
  "T2D_status.txt.gz",
  sep="\t",header=T,quote="",colClasses=c("FID"='character')
)
all_T2D_individuals  <- all_T2D_prevalence %>% dplyr::select(FID,disease_status) %>% dplyr::rename(T2D = disease_status)

all_T1D_prevalence  <- read.table(
  "T1D_status.txt.gz",
  sep="\t",header=T,quote="",colClasses=c("FID"='character')
)
all_T1D_individuals <- all_T1D_prevalence %>% dplyr::select(FID,disease_status) %>% dplyr::rename(T1D = disease_status)


# inner join will remove any incident T2D or T1D individuals or those we can't figure out T2D or T1D incident/prevalent status for.
# if these folks are control -> kept
all_diabetic_individuals  <- inner_join(
  all_T2D_individuals %>% dplyr::select(FID,T2D),
  all_T1D_individuals  %>% dplyr::select(FID,T1D),
  by = c("FID")
)

#%% format the death date
#' The most recent censor date is earliest of death or follow update
#' 
#' Results:
#'  There seems to be people who have records of death and diagnosis after lost to followup.
#'  Mentioned this field was not updated since 2017, perhaps some updated again?
#'  Only death date is considered (https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=191)
death_and_lost_data  <- read.table("UKB_health_records_meta_data_participant_June20_2024.tsv",sep="\t",header=T,quote="",colClasses = c("Participant.ID"="character")) %>% rename(
  FID = Participant.ID,
  date_lost_to_followup = Date.lost.to.follow.up,
  date_of_death = Date.of.death...Instance.0
) %>% dplyr::select(FID,date_of_death)

death_and_lost_data  <- death_and_lost_data %>% filter(date_of_death != "")

death_and_lost_data$date_of_death <- as_datetime(death_and_lost_data$date_of_death,format=c("%Y-%m-%d"))

censor_date_map  <- setNames(
  death_and_lost_data$date_of_death,nm=death_and_lost_data$FID
)




#%%
#' Merge the data we want

km_info  <- inner_join(
  olink_data %>% filter(eid %in% keep_ind$eid),
  ukb_covar_info,
  by=c("eid"="FID") 
) %>% inner_join(
 all_blood_draw_dates %>% dplyr::select(eid,blood_draw_date),by=c("eid")
)
km_info$blood_draw_date  <- as_datetime(
  km_info$blood_draw_date,format=c("%Y-%m-%d")
)
# remove all prevalent KD cases
km_info  <- inner_join(
  km_info,
  all_kidney_disease_incidence,by=c("eid"="FID")
)
# the death recorded added cases had no blood draw date as it was not originally considered as cases.
#' the blood draw date should be added in later 
check_df  <- km_info %>% filter(disease_status == 1 & !is.na(blood_draw_date.y))
# check consistent blood draw date (i.e., no issue with merging)
assertthat::assert_that(
  all(check_df$blood_draw_date.x == check_df$blood_draw_date.y))
rm(check_df)

km_info  <- km_info %>% dplyr::select(-blood_draw_date.y) %>% rename(blood_draw_date = blood_draw_date.x)
assertthat::assert_that(
  sum(is.na(km_info$blood_draw_date)) == 0
)

km_info  <- inner_join(
  km_info,
  all_diabetic_individuals %>% dplyr::select(FID,T1D,T2D),
  by = c("eid"="FID")
)
km_info$earliest_death_lost_date  <- censor_date_map[km_info$eid]

#%% remove individuals 
# remove double diabetics
km_info  <- km_info %>% filter(
 !(T1D == 1 & T2D == 1) 
)
# remove non-diabetics but abnormally high HbA1c
to_remove  <- km_info %>% filter(
  T2D == 0 & T1D == 0 & HbA1c > 48
)
km_info <- km_info %>% filter(! eid %in% to_remove$eid)
#%% merge with ancestry information

km_info  <- inner_join(
  km_info,
  all_sample_ancestry %>% dplyr::select(-IID),
  by= c("eid" = "FID")
)

#%%
#' Compute the time of followup
# Last follow-up date: Sep 2023 (https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41270)
last_followup_date = as_datetime("2023-09-01",format="%Y-%m-%d")


#' Right censor date is the earliest of occurence of one of the following:
#'  - kidney disease diagnosis dates
#'  - death date
#' If neither happened, then the last followup date was used.
km_info$right_censor_date <- apply(
  km_info,
  MARGIN = 1,
  FUN = function(x){
    if (is.na(x['earliest_death_lost_date']) & is.na(x['earliest_date'])){
      date = NA
    }else{
    date = min(x['earliest_death_lost_date'],x['earliest_date'],na.rm=T)
    }
    date
  }
)
km_info$right_censor_date  <- as_datetime(km_info$right_censor_date,format=c("%Y-%m-%d"))

# Since the right censor dates are either death or diagnoses, all diagnoses should preceed death
# this means that if there is a diagnosis, it should be the right censor date
check_df  <- km_info %>% filter(!is.na(earliest_date))
assertthat::assert_that(
  all(check_df$right_censor_date == check_df$earliest_date)
)
rm(check_df)

# add in the last update for the remaining
km_info$right_censor_date  <- replace(km_info$right_censor_date,is.na(km_info$right_censor_date),last_followup_date)




# all right censor date must be at most the last follow update
assertthat::assert_that(
  all(km_info$right_censor_date <= last_followup_date)
)
# all disease status is 0 or one and no missing
assertthat::assert_that(
  all(sort(unique(all_diabetic_individuals$T1D)) == c(0,1)) &
  sum(is.na(all_diabetic_individuals$T1D)) == 0
)
assertthat::assert_that(
  all(sort(unique(all_diabetic_individuals$T2D)) == c(0,1)) &
  sum(is.na(all_diabetic_individuals$T2D)) == 0
)
assertthat::assert_that(
  all(sort(unique(all_diabetic_individuals$T2D)) == c(0,1)) &
  sum(is.na(km_info$disease_status)) == 0
)
# There were outliers with HbA1c with high HbA1c, but low glucose.
# make sure were not included
assertthat::assert_that(
  nrow(km_info %>% filter(HbA1c > 200)) == 0
)

# make sure no more double diabetics
assertthat::assert_that(
  nrow(
    km_info %>% filter(T2D == 1 & T1D == 1)
  ) == 0
)

# compute follow up time and save
km_info$time_of_followup  <- km_info$right_censor_date - km_info$blood_draw_date



fwrite(
  km_info,
  file.path(outdir,glue('KM_information_{toupper(protein_name)}_ALL_Ancestries_update.tsv')),sep="\t",col.names=T,quote=F
)

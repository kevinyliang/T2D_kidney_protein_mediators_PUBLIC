import os,sys,glob,re,yaml
from collections import namedtuple
import pandas as pd
import numpy as np
import pyreadr

from multiprocessing import Pool
from functools import partial


def get_prevalence(row,ICD_10,OPCS_4,opcs4_date_parsed,icd10_date_parsed,blood_draw_fu0_date_uniqued):
  # print(f"..{row.FID}..",end="")
  
  if ICD_10 != "":
    if isinstance(row.ICD10_summary_diagnoses,str):
      # does this person have the disease?
      # i.e., any ICD-10 code matches
      case = bool(re.search(ICD_10,row.ICD10_summary_diagnoses))
      if not case:
        ICD10_case_control_stat = 0
      else:
        # if this person has disease, which ICD-10 are the relevant ones
        # what is the earliest dates?

        # find the relevant ICD 10 codes
        individual_ICD_10 = [
          x.split(" ")[0].strip() for x in 
          row.ICD10_summary_diagnoses.split("|")
        ]
        # find the coresdponding dates
        icd10_date = icd10_date_parsed.query(f"FID == '{row.FID}'").drop(['FID'],axis=1)
        assert(len(icd10_date) == 1),f"more than 1 ICD10 {row.FID}"
        
        # only keep diagnoses with dates
        earliest_icd10_date = pd.to_datetime(icd10_date.iloc[0,:],format = "%Y-%m-%d").values[0:len(individual_ICD_10)]
        
        # matching each one ID, so can strip the space
        # the space sometime use to find exact ICD code matches
        relevant_ICD10 = [
          bool(re.search(re.sub(" ","",ICD_10),x[0])) and not np.isnat(x[1]) for x in zip(individual_ICD_10,earliest_icd10_date)
        ]
        # this person is a case, but there isn't a date for it
        # can't distinguish incidence/prevalence
        # dropped
        if np.sum(relevant_ICD10) == 0:
          ICD10_case_control_stat = np.nan
        else:
          earliest_icd10_date = earliest_icd10_date[relevant_ICD10]
          earliest_icd10_codes = np.asarray(individual_ICD_10)[relevant_ICD10]
          # among the selected codes and dates, pick the earliest
          earliest_icd10_codes = earliest_icd10_codes[earliest_icd10_date == earliest_icd10_date.min()]
          earliest_icd10_date = earliest_icd10_date[earliest_icd10_date == earliest_icd10_date.min()]
          # even if there are 2 relevant ICD10 code with the same date, it doesn't matter
          # pick 1 
          earliest_icd10_codes = earliest_icd10_codes[0]
          earliest_icd10_date = earliest_icd10_date[0]
          
          # if this person has disease, is it incidence?
          blood_draw_date = blood_draw_fu0_date_uniqued.query(
            f"FID == '{row.FID}'"
          )['date']
          assert(len(blood_draw_date) < 2),f"{row.FID}"
          if len(blood_draw_date) < 1:
            # this is unknown -> dropped
            ICD10_case_control_stat = np.nan
          else:
            blood_draw_date = pd.to_datetime(blood_draw_date.values[0],format="%Y-%m-%d")
            if blood_draw_date > earliest_icd10_date:
              ICD10_case_control_stat = 1
            else:
              # this is not incidence -> dropped
              ICD10_case_control_stat = np.nan
    else:
      ICD10_case_control_stat = 0
  else:
    # not defined by ICD-10
    ICD10_case_control_stat = 0
  
  if OPCS_4 != "":
    if isinstance(row.OPCS4_summary_operations,str):
      # does this person have the disease?
      case = bool(re.search(OPCS_4,row.OPCS4_summary_operations))
      if not case:
        OPCS4_case_control_stat = 0
      else:
        individual_OPCS_4 = [
          x.split(" ")[0].strip() for x in 
          row.OPCS4_summary_operations.split("|")
        ]

        
        OPCS4_date = opcs4_date_parsed.query(f"FID == '{row.FID}'").drop(['FID'],axis=1)
        assert(len(OPCS4_date) == 1),f"more than 1 OPCS4 {row.FID}"
        # find the date of the relevant OPCS_4
        
        earliest_OPCS4_date = pd.to_datetime(OPCS4_date.iloc[0,:],format = "%Y-%m-%d").values[0:len(individual_OPCS_4)]
        
        relevant_OPCS4 = [
          bool(re.search(re.sub(" ","",OPCS_4),x[0])) and not np.isnat(x[1]) for x in zip(individual_OPCS_4,earliest_OPCS4_date)
        ]
        if np.sum(relevant_OPCS4) == 0:
          OPCS4_case_control_stat = np.nan
        else:
          earliest_OPCS4_date = earliest_OPCS4_date[relevant_OPCS4]
          earliest_OPCS4_codes = np.asarray(individual_OPCS_4)[relevant_OPCS4]
          
          # among the selected codes and dates, pick the earliest
          earliest_OPCS4_codes = earliest_OPCS4_codes[earliest_OPCS4_date == earliest_OPCS4_date.min()]
          earliest_OPCS4_date = earliest_OPCS4_date[earliest_OPCS4_date == earliest_OPCS4_date.min()]
          # even if there are 2 relevant OPCS4 code with the same date, it doesn't matter
          # pick 1 
          earliest_OPCS4_codes = earliest_OPCS4_codes[0]
          earliest_OPCS4_date = earliest_OPCS4_date[0]

          
          blood_draw_date = blood_draw_fu0_date_uniqued.query(
            f"FID == '{row.FID}'"
          )['date']
          assert(len(blood_draw_date) < 2),f"{row.FID}"
          if len(blood_draw_date) < 1:
            # this is unknown -> dropped
            OPCS4_case_control_stat = np.nan
          else:
            blood_draw_date = pd.to_datetime(blood_draw_date.values[0],format="%Y-%m-%d")
            if blood_draw_date > earliest_OPCS4_date:
              OPCS4_case_control_stat = 1
            else:
              # this is not incidence -> dropped
              OPCS4_case_control_stat = np.nan
    else:
      OPCS4_case_control_stat = 0
  else:
    # not defined by ICD-10
    OPCS4_case_control_stat = 0

  # if undefined by both, it is dropped
  if np.isnan(OPCS4_case_control_stat) or np.isnan(ICD10_case_control_stat):
    final_status = np.nan
    earliest_date = np.nan
    earliest_code = np.nan
    earliest_src = np.nan
    blood_draw_date = np.nan
  # if either is case, then we figure out which one and get the dates
  elif OPCS4_case_control_stat == 1 or ICD10_case_control_stat == 1:
    final_status = 1
    if OPCS4_case_control_stat == 1 and ICD10_case_control_stat == 1:
      earliest_date = min(earliest_OPCS4_date,earliest_icd10_date)
      if earliest_date == earliest_OPCS4_date:
        earliest_code = earliest_OPCS4_codes
        earliest_src = 'OPCS4'
      elif earliest_date == earliest_icd10_date:
        earliest_code = earliest_icd10_codes
        earliest_src = 'ICD10'
    elif OPCS4_case_control_stat == 1:
      final_status = 1
      earliest_date = earliest_OPCS4_date
      earliest_code = earliest_OPCS4_codes
      earliest_src = 'OPCS4'
    elif ICD10_case_control_stat == 1:
      final_status = 1
      earliest_date = earliest_icd10_date
      earliest_code = earliest_icd10_codes
      earliest_src = 'ICD10'
    else:
      assert(False),f'error with the case control {row}'
  # if both are control then it is settled
  elif OPCS4_case_control_stat == 0 or ICD10_case_control_stat == 0:
    final_status = 0
    earliest_date = np.nan
    earliest_code = np.nan
    earliest_src = np.nan
    blood_draw_date = np.nan
  # unconsidered conditions
  else:
    assert(False),f"What status {row}"
  return final_status,earliest_date,earliest_code,earliest_src,blood_draw_date


def get_incidence(row,ICD_10,OPCS_4,opcs4_date_parsed,icd10_date_parsed,blood_draw_fu0_date_uniqued):
  # print(f"..{row.FID}..",end="")
  if ICD_10 != "":
    if isinstance(row.ICD10_summary_diagnoses,str):
      # does this person have the disease?
      # i.e., any ICD-10 code matches
      case = bool(re.search(ICD_10,row.ICD10_summary_diagnoses))
      if not case:
        ICD10_case_control_stat = 0
      else:
        # if this person has disease, which ICD-10 are the relevant ones
        # what is the earliest dates?
        # find the relevant ICD 10 codes
        individual_ICD_10 = [
          x.split(" ")[0].strip() for x in 
          row.ICD10_summary_diagnoses.split("|")
        ]
        # find the coresdponding dates
        icd10_date = icd10_date_parsed.query(f"FID == '{row.FID}'").drop(['FID'],axis=1)
        assert(len(icd10_date) == 1),f"more than 1 ICD10 {row}"
        
        # only keep diagnoses with dates
        # obtain the vector of dates in the same order
        # since we can match by array structure
        # https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41282
        # https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41280
        earliest_icd10_date = pd.to_datetime(icd10_date.iloc[0,:],format = "%Y-%m-%d").dt.floor('D').values[0:len(individual_ICD_10)]
        
        relevant_ICD10 = [
          bool(re.search(re.sub(" ","",ICD_10),x[0])) and not np.isnat(x[1]) for x in zip(individual_ICD_10,earliest_icd10_date)
        ]
        if np.sum(relevant_ICD10) == 0:
          ICD10_case_control_stat = np.nan
        else:
          earliest_icd10_date = earliest_icd10_date[relevant_ICD10]
          earliest_icd10_codes = np.asarray(individual_ICD_10)[relevant_ICD10]
          # among the selected codes and dates, pick the earliest
          earliest_icd10_codes = earliest_icd10_codes[earliest_icd10_date == earliest_icd10_date.min()]
          earliest_icd10_date = earliest_icd10_date[earliest_icd10_date == earliest_icd10_date.min()]
          # even if there are 2 relevant ICD10 code with the same date, it doesn't matter
          # pick 1 
          earliest_icd10_codes = earliest_icd10_codes[0]
          earliest_icd10_date = earliest_icd10_date[0]
          
          # if this person has disease, is it incidence?
          blood_draw_date = blood_draw_fu0_date_uniqued.query(
            f"FID == '{row.FID}'"
          )['date']
          assert(len(blood_draw_date) < 2),f"{row.FID}"
          if len(blood_draw_date) < 1:
            # this is unknown -> dropped
            ICD10_case_control_stat = np.nan
          else:
            blood_draw_date = pd.to_datetime(blood_draw_date.values[0],format="%Y-%m-%d")
            if blood_draw_date < earliest_icd10_date:
              ICD10_case_control_stat = 1
            else:
              # this is not incidence -> dropped
              ICD10_case_control_stat = np.nan
    else:
      ICD10_case_control_stat = 0
  else:
    # not defined by ICD-10
    ICD10_case_control_stat = 0
  
  if OPCS_4 != "":
    if isinstance(row.OPCS4_summary_operations,str):
      # does this person have the disease?
      case = bool(re.search(OPCS_4,row.OPCS4_summary_operations))
      if not case:
        OPCS4_case_control_stat = 0
      else:
        individual_OPCS_4 = [
          x.split(" ")[0].strip() for x in 
          row.OPCS4_summary_operations.split("|")
        ]
        
        OPCS4_date = opcs4_date_parsed.query(f"FID == '{row.FID}'").drop(['FID'],axis=1)
        assert(len(OPCS4_date) == 1),f"more than 1 OPCS4 {row.FID}"
        # find the date of the relevant OPCS_4
        earliest_OPCS4_date = pd.to_datetime(OPCS4_date.iloc[0,:],format = "%Y-%m-%d").values[0:len(individual_OPCS_4)]
        
        relevant_OPCS4 = [
          bool(re.search(re.sub(" ","",OPCS_4),x[0])) and not np.isnat(x[1]) for x in zip(individual_OPCS_4,earliest_OPCS4_date)
        ]
        if np.sum(relevant_OPCS4) == 0:
          OPCS4_case_control_stat = np.nan
        else:
          earliest_OPCS4_date = earliest_OPCS4_date[relevant_OPCS4]
          earliest_OPCS4_codes = np.asarray(individual_OPCS_4)[relevant_OPCS4]
          # among the selected codes and dates, pick the earliest
          earliest_OPCS4_codes = earliest_OPCS4_codes[earliest_OPCS4_date == earliest_OPCS4_date.min()]
          earliest_OPCS4_date = earliest_OPCS4_date[earliest_OPCS4_date == earliest_OPCS4_date.min()]
          # even if there are 2 relevant OPCS4 code with the same date, it doesn't matter
          # pick 1 
          earliest_OPCS4_codes = earliest_OPCS4_codes[0]
          earliest_OPCS4_date = earliest_OPCS4_date[0]

          
          blood_draw_date = blood_draw_fu0_date_uniqued.query(
            f"FID == '{row.FID}'"
          )['date']
          assert(len(blood_draw_date) < 2),f"{row.FID}"
          if len(blood_draw_date) < 1:
            # this is unknown -> dropped
            OPCS4_case_control_stat = np.nan
          else:
            blood_draw_date = pd.to_datetime(blood_draw_date.values[0],format="%Y-%m-%d")
            if blood_draw_date < earliest_OPCS4_date:
              OPCS4_case_control_stat = 1
            else:
              # this is not incidence -> dropped
              OPCS4_case_control_stat = np.nan
    else:
      OPCS4_case_control_stat = 0
  else:
    # not defined by OPCS-4
    OPCS4_case_control_stat = 0

  # if undefined by both, it is dropped
  if np.isnan(OPCS4_case_control_stat) or np.isnan(ICD10_case_control_stat):
    final_status = np.nan
    earliest_date = np.nan
    earliest_code = np.nan
    earliest_src = np.nan
    blood_draw_date = np.nan
  # if either is case, then we figure out which one and get the dates
  elif OPCS4_case_control_stat == 1 or ICD10_case_control_stat == 1:
    final_status = 1
    if OPCS4_case_control_stat == 1 and ICD10_case_control_stat == 1:
      earliest_date = min(earliest_OPCS4_date,earliest_icd10_date)
      if earliest_date == earliest_OPCS4_date:
        earliest_code = earliest_OPCS4_codes
        earliest_src = 'OPCS4'
      elif earliest_date == earliest_icd10_date:
        earliest_code = earliest_icd10_codes
        earliest_src = 'ICD10'
    elif OPCS4_case_control_stat == 1:
      final_status = 1
      earliest_date = earliest_OPCS4_date
      earliest_code = earliest_OPCS4_codes
      earliest_src = 'OPCS4'
    elif ICD10_case_control_stat == 1:
      final_status = 1
      earliest_date = earliest_icd10_date
      earliest_code = earliest_icd10_codes
      earliest_src = 'ICD10'
    else:
      assert(False),f'error with the case control {row}'
  # if both are control then it is settled
  elif OPCS4_case_control_stat == 0 or ICD10_case_control_stat == 0:
    final_status = 0
    earliest_date = np.nan
    earliest_code = np.nan
    earliest_src = np.nan
    blood_draw_date = np.nan
  # unconsidered conditions
  else:
    assert(False),f"What status {row}"
  return final_status,earliest_date,earliest_code,earliest_src,blood_draw_date

def get_case(row,ICD_10,OPCS_4):
  if ICD_10 != "":
    if isinstance(row.ICD10_summary_diagnoses,str):
      ICD10_case_control_stat = bool(re.search(ICD_10,row.ICD10_summary_diagnoses))
    else:
      ICD10_case_control_stat = False
  else:
    ICD10_case_control_stat = False
    
  if OPCS_4 != "":
    if isinstance(row.OPCS4_summary_operations,str):
      OPCS4_case_control_stat = bool(re.search(OPCS_4,row.OPCS4_summary_operations))
    else:
      OPCS4_case_control_stat = False
  else:
    OPCS4_case_control_stat = False

  case_control_stat = ICD10_case_control_stat or OPCS4_case_control_stat
  return case_control_stat

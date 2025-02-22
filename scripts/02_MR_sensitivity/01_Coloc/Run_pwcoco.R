#' run_pwcoco
#' Prep and run pwCoCo
#' 

source("get_p_from_bse.R")
run_pwcoco <- function(harmonized_res,outdir,protein_gwas_file,disease_gwas_file,refpan_dir){
  protein_name  <- harmonized_res$exposure %>% unique()
  disease_name  <- harmonized_res$outcome %>% unique()
  outdir  <- file.path(outdir,glue("pwcoco_Prot_{protein_name}_Disease_{disease_name}"))
  system(command = glue("mkdir -p {outdir}"))

  
  soma_info  <- unlist(strsplit(unlist(
    strsplit((strsplit(protein_name,"_UP_") %>% unlist())[1],"_SL_")
  )[2],"@"))



  protein_gwas  <- read.table(protein_gwas_file,sep="\t",header=T,colClasses = c("effect_allele.outcome"="character",'other_allele.outcome'='character'))

  protein_gwas$pos.outcome  <- as.numeric(protein_gwas$pos.outcome)
  
  disease_gwas  <- readRDS(disease_gwas_file)

  for (each_mr_snp in 1:nrow(harmonized_res)){
    mr_snp_info  <- harmonized_res[each_mr_snp,]
    mr_snp_info$pos  <- lapply(
      strsplit(mr_snp_info$SNP,":"),FUN=function(x){x[2]}
    ) %>% unlist()
    chr  <- mr_snp_info$chr.exposure
    pos  <- as.numeric(mr_snp_info$pos)
    snp  <- mr_snp_info$SNP %>% toupper()
    pwcoco_out_file  <- file.path(outdir,glue('coco_res_{snp}'))
    if (file.exists(glue("{pwcoco_out_file}.coloc"))){
      next
    }
    tdir  <- file.path(outdir,glue('COCO_inputs{snp}'))
    system(command=glue("mkdir -p {tdir}"))
    min_pos = max(pos - 500000,1)
    max_pos = pos + 500000
    rel_protein_region  <- protein_gwas %>% 
      dplyr::filter(pos.outcome >= min_pos & pos.outcome <= max_pos & chr.outcome == chr)
    colnames(rel_protein_region)  <- gsub("outcome","exposure",colnames(rel_protein_region))
    
    disease_gwas  <- disease_gwas %>% dplyr::filter(chr.outcome == chr)
    disease_gwas$pos  <- lapply(
      strsplit(disease_gwas$SNP,":"),FUN=function(x){x[2]}
    ) %>% unlist()
    disease_gwas$pos  <- as.numeric(disease_gwas$pos)
    rel_disease_region  <- disease_gwas %>% 
      dplyr::filter(pos >= min_pos & pos <= max_pos)

    # prepare the data to harmonize so the beta and everything are mapped
    hres  <- TwoSampleMR::harmonise_data(
      rel_protein_region,
      rel_disease_region,
      action=2
    )
    rel_exp_col  <- colnames(hres)[grepl("exposure",colnames(hres))]
    rel_protein_region_matched  <- hres %>% dplyr::select(SNP,all_of(rel_exp_col))
    rel_protein_region_matched  <- rel_protein_region_matched %>% 
      dplyr::select(
        SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,se.exposure,beta.exposure,pval.exposure,samplesize.exposure
      ) %>% 
      dplyr::rename(
        A1 = effect_allele.exposure,
        A2 = other_allele.exposure,
        A1_freq = eaf.exposure,
        beta = beta.exposure,
        se = se.exposure,
        n= samplesize.exposure
      )  %>%
      dplyr::mutate(
        SNP = toupper(SNP),
        p  = get_p_from_bse(beta,se)
      ) %>% dplyr::select(
        SNP,A1,A2,A1_freq,beta,se,p,n
      )


    rel_out_col  <- colnames(hres)[grepl("outcome",colnames(hres))]
    rel_disease_region_matched  <- hres %>% dplyr::select(SNP,all_of(rel_out_col))
    rel_disease_region_matched  <- rel_disease_region_matched %>% dplyr::select(
        SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,se.outcome,beta.outcome,pval.outcome,samplesize.outcome
      ) %>% 
      dplyr::rename(
        A1 = effect_allele.outcome,
        A2 = other_allele.outcome,
        A1_freq = eaf.outcome,
        beta = beta.outcome,
        se = se.outcome,
        n= samplesize.outcome
      ) %>%
      dplyr::mutate(
        SNP = toupper(SNP),
        p  = get_p_from_bse(beta,se)
      ) %>% dplyr::select(
        SNP,A1,A2,A1_freq,beta,se,p,n
      )

    protein_sumstat_file  <- file.path(tdir,glue("coco_{snp}.protein.tsv"))
    fwrite(
      rel_protein_region_matched,
      protein_sumstat_file,
      sep=",",quote=F
    )
    disease_sumstat_file  <- file.path(tdir,glue("coco_{snp}.disease.tsv"))
    fwrite(
      rel_disease_region_matched,
      disease_sumstat_file,
      sep=",",quote=F
    )
    
    bfile  <- file.path(
      refpan_dir,chr
    )
    pwcoco_cmd  <- paste0(
      glue("{pwcoco} "),
      glue("--bfile {bfile} "),
      glue("--sum_stats1 {protein_sumstat_file} "),
      glue("--sum_stats2 {disease_sumstat_file} "),
      glue("--top_snp 10 "),
      glue("--maf 0.01 "),
      glue("--out_cond "),
      glue("--verbose "),
      glue("--out {pwcoco_out_file}")
    )
    system(command = pwcoco_cmd,wait=T)
  }
}
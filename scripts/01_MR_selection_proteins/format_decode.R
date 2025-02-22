#' Obtain the list of proteins to run MR for and process as outcome and exposures
#' 
#' Protein GWAS from DeCode
#' eaf: eaf = effectAlleleFreq (file name)
#' 
library(tidyverse)
library(data.table)
library(glue)
library(yaml)

library(TwoSampleMR)


proj_yml <- "configs/proj_config.yml"
args  <- read_yaml(proj_yml)[['MR_disease_protein']]

source("snp_id_format.R")
source("custom_modifyBuild.R")

outdir  <- file.path(args$rt_dir,'results','deCODE_protein_MR_inputs')



cargs  <- commandArgs(trailingOnly = T)


gene_ids  <- unlist(strsplit(cargs,"\t"))[1]
soma_gene_name_seqID  <- unlist(strsplit(unlist(strsplit(gene_ids,"_SL_"))[2],"_UP_"))[1]
soma_gene_name  <- unlist(strsplit(soma_gene_name_seqID,"@"))[1]
soma_seqID  <- unlist(strsplit(soma_gene_name_seqID,"@"))[2]


if (file.exists(file.path(outdir,glue("{gene_ids}_MR_outcome.tsv.gz")))){
  quit()
}
if (file.exists(file.path(outdir,'results_with_proteins_wo_pQTL',glue("{gene_ids}_MR_outcome.tsv.gz")))){
  quit()
}




decode_sumstat_file  <- list.files(
  args$DeCode_dir,
  pattern=glue("^{soma_seqID}_.*gz$")
)
assertthat::assert_that(length(decode_sumstat_file)==1,msg=glue("{gene_ids}"))


decode_sumstat  <- fread(file.path(args$DeCode_dir,decode_sumstat_file),sep="\t")

updated_allele_info  <- fread('assocvariants.annotated.txt.gz',sep="\t")
# update the allele frequency and effect/other notation as described in readme
decode_sumstat_updated  <- inner_join(
  updated_allele_info,
  decode_sumstat %>%
  dplyr::select(
    Name,Beta,Pval,SE,N
  ),by=c("Name")
)
decode_sumstat_updated  <- decode_sumstat_updated %>% 
  dplyr::filter(rsids != '.' & effectAllele %in% c("A","T","C","G") & otherAllele %in% c("A","T","C","G"))

# update the SNP name to be same as outcome.
# First liftOver because the positions of outcome are in Build 37 (hg19). This is build 38

liftOver  <- "programs/liftOver"
chain  <- "liftover_chain_files/hg38Tohg19.over.chain.gz"
rev_chain  <- "liftover_chain_files/hg19Tohg38.over.chain.gz"
decode_sumstat_updated  <- decode_sumstat_updated %>% 
  rename(chr = Chrom,pos = Pos)
# format the chrom name so it matches up
decode_sumstat_updated$chr  <- gsub("chr","",decode_sumstat_updated$chr)
decode_sumstat_updated  <- own_modifyBuild(
  info_snp = decode_sumstat_updated,
  liftOver = liftOver,
  chain = chain,
  rev_chain = rev_chain,
  from='hg38',
  to='hg19',
  check_reverse=T
)
decode_sumstat_updated  <- decode_sumstat_updated %>% dplyr::filter(!is.na(pos))

decode_sumstat_updated$formatted_name  <- apply(
  decode_sumstat_updated,
  MARGIN = 1,
  FUN = function(x){
    format_snp_id(c(x['chr'],x['pos'],x['effectAllele'],x['otherAllele']),
    args$allele_orders)
  }
)
decode_sumstat_updated  <- decode_sumstat_updated %>% 
  dplyr::filter(!grepl("ERROR:",formatted_name))


decode_sumstat_updated_outcome  <- format_data(
  decode_sumstat_updated,
  type='outcome',
  snp_col = 'formatted_name',
  beta_col = 'Beta',
  se_col = 'SE',
  eaf_col = 'effectAlleleFreq',
  effect_allele_col = 'effectAllele',
  other_allele_col = 'otherAllele',
  pval_col = 'Pval',
  samplesize_col = 'N'
)
decode_sumstat_updated_outcome$outcome  <- gene_ids

fwrite(
  decode_sumstat_updated_outcome,
  file.path(outdir,glue("{gene_ids}_MR_outcome.tsv.gz")),
  quote=F,sep="\t"
)
fwrite(
  decode_sumstat_updated %>% dplyr::select(
    rsids,Name,formatted_name
  )%>%dplyr::rename(
    deCode_original_name = Name,
    formatted_name_b37 = formatted_name
  ),
  file.path(outdir,glue("{gene_ids}_MR_outcome_SNPNAMES.tsv.gz")),sep="\t",quote=F
)
















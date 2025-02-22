#' Make locus zoom plots
library(tidyverse)
library(ggplot2)
library(data.table)
library(patchwork)
library(yaml)
library(glue)
source("get_p_from_bse.R")


make_locus_zoom  <- function(file_title,file_type,pop,width=10,height=7){
  protein_input  <- read.table(
    file.path(coco_inputs,list.files(coco_inputs,pattern=glue("\\.{file_type}\\.tsv"))[1]),header=T,sep=",",quote=""
    )  %>% 
    dplyr::rename(
      effect_allele = A1,
      other_allele = A2,
      eaf = A1_freq
    ) %>% 
    filter(
      eaf <= 0.99 | eaf >= 0.01 # filter for 1%
    )

  snp_info  <- strsplit(protein_input$SNP,':')
  protein_input$chrom  <- lapply(
    snp_info,FUN = function(x){x[[1]]}
  ) %>% unlist()
  protein_input$pos  <- lapply(
    snp_info,FUN = function(x){x[[2]]}
  ) %>% unlist() %>% as.numeric()


  chr  <- unique(protein_input$chrom)
  chr_map_file  <- readRDS(glue("ukb_imp_chr{chr}_v3.list.formatted_name.Rds"))
  name_map  <- setNames(
    chr_map_file$rsid,
    nm = chr_map_file$each_formatted_ids
  )
  protein_input$rsid  <- name_map[toupper(protein_input$SNP)]
  protein_input  <- protein_input %>% dplyr::filter(!is.na(rsid))
  ref_snp  <- gsub("COCO_inputs","",basename(coco_inputs))
  ref_snp_rsid  <- name_map[toupper(ref_snp)]

  # compute own P-values as they are capped by some GWAS
  protein_input  <- protein_input %>% dplyr::select(-p)
  protein_input$p  <- get_neglog10_p_from_bse(
    BETA = protein_input$beta,
    SE = protein_input$se
  )

  locus_zoom_inputs  <- protein_input %>% dplyr::select(rsid,p)




  fwrite(
    protein_input %>% 
    dplyr::select(rsid,p) %>% 
    dplyr::rename(
      MarkerName = rsid,
      'P-value' = p
    ),
    file.path(
      outdir,glue('locuszoom_inputs_{file_title}.tsv')
    ),sep="\t",quote=F
  )


  locuszoom_cmd  <- paste0(
    glue("source ~/.bashrc && conda activate ldsc && cd {outdir} && programs/locuszoom/bin/locuszoom "),
    glue("--meta {file.path(
      outdir,glue('locuszoom_inputs_{file_title}.tsv')
    )} "),
    glue("--build hg19 --pop {pop} --source 1000G_Nov2014 "),
    glue("--no-cleanup --no-transform "),
    glue("width={width} "),
    glue("height={height} "),
    glue("--prefix {file_title} "),
    glue("--flank 500kb "),
    glue("--refsnp {ref_snp_rsid} ")
  )
  system(command = locuszoom_cmd,wait=T)
}


plink  <- "programs/plink"

coco_inputs  <- "COCO_inputs12:57838150:C:T"
outdir_base  <- "INHBC_eGFR"
snp  <- gsub(":","_",gsub("COCO_inputs","",basename(coco_inputs)))
disease  <- "eGFR"
protein  <- "INHBC_eGFR"



outdir  <- file.path(outdir_base,disease)
system(command = glue("mkdir -p {outdir}"))
make_locus_zoom(disease,"disease",pop='EUR')
plot_dir  <- list.dirs(
  file.path(outdir)
)
plot_dir  <- plot_dir[grepl(disease,plot_dir)]
disease_file  <- list.files(
  plot_dir,pattern='^chr.*\\.pdf',full.name=T
)[1]
target_name  <- file.path(
  glue("{disease}_{basename(disease_file)}")
)
system(command = glue("mv {disease_file} {file.path(outdir_base,target_name)}"))
unlink(outdir,recursive=T)

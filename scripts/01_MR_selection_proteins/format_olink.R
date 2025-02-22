#' Common script to format UKB GWAS output
source("snp_id_format.R")
format_single_olink_chr  <- function(chromosome_file,olink_pos_map_dir,chr,mr_type='outcome'){
  chromosome_pos_map_files  <- list.files(
      olink_pos_map_dir,
      pattern=glue("_{chr}_patched_v2\\.tsv\\.gz"),full.names=T
    )
  assertthat::assert_that(
    length(chromosome_pos_map_files) == 1
  )
  chromosome_pos_map_files  <- read.table(
    chromosome_pos_map_files[1],sep="\t",header=T,quote=""
  )
  chromosome_res   <- read.table(
    chromosome_file,sep=" ",header=T,colClasses=c("ALLELE0"="character","ALLELE1"="character")
  )
  chromosome_res_w_grch37_pos  <- inner_join(
    chromosome_res,
    chromosome_pos_map_files %>% dplyr::select(ID,rsid,POS19,POS38),
    by = c("ID")
  )
  assertthat::assert_that(
    identical(chromosome_res_w_grch37_pos$GENPOS,chromosome_res_w_grch37_pos$POS38)
  )
  chromosome_res_w_grch37_pos  <- chromosome_res_w_grch37_pos %>% 
    dplyr::rename(
      beta = BETA,
      se = SE,
      effect_allele = ALLELE1,
      other_allele = ALLELE0,
      eaf=A1FREQ,
      chr=CHROM,
      position = POS19,
      samplesize=N
    ) %>% 
    mutate(
      pval = 10^(-LOG10P)
    )

  doMC::registerDoMC(5)
  snp_df_batch  <- split(seq(1,nrow(chromosome_res_w_grch37_pos)),ceiling(seq_along(seq(1,nrow(chromosome_res_w_grch37_pos))) / 5000))
  chromosome_res_outcome_formatted  <- foreach(batch_idx = seq(1,length(snp_df_batch)),.combine='rbind') %dopar% {
    # print(glue("{batch_idx} / {length(snp_df_batch)}"))
    batch_df  <- chromosome_res_w_grch37_pos[snp_df_batch[[batch_idx]],]
    batch_df$SNP  <- apply(
      X = batch_df,
      MARGIN = 1,
      FUN = function(x){
        ID  <- format_snp_id(
          c(x['chr'],x['position'],x['effect_allele'],x['other_allele']),
          c("A","C","G","T")
        )
        return(ID)
      }
    ) %>% unlist()
    batch_df
  }
  chromosome_res_outcome_formatted  <- chromosome_res_outcome_formatted %>% filter(!grepl("^ERROR",SNP))
  twosample_mr_outcome  <- format_data(chromosome_res_outcome_formatted %>%
    dplyr::select(
        SNP,chr,position,eaf,beta,se,pval,effect_allele,other_allele,samplesize
      ),type=mr_type,header=T
    )
  twosample_mr_outcome$outcome  <- cargs[1]
  return(
    list(
      'twosample_mr_outcome' = twosample_mr_outcome,
      'chromosome_name_map' = chromosome_res_outcome_formatted %>% dplyr::select(ID,SNP,rsid,position,GENPOS,POS38) %>% rename(
        original_name = ID,
        formatted_name = SNP,
        POS19=position
      )
    )
  )
}
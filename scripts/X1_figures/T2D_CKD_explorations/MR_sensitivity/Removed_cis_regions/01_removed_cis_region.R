#' Remove the T2D instruments near INHBC (within 1mb on either side of INHBC) and re-run MR
#' to avoid reverse causation
#' 
#' Results:
#'  No T2D IV are in this region. Nothing to remove.
library(tidyverse)
library(glue)
library(data.table)
library(yaml)


args  <- read_yaml("configs/proj_config.yml")[['proj_resource']]
outdir <- "removed_cis_region_t2d_on_inhbc"

inhbc_chr <- 12
inhbc_transcript <- c(57828543,57844609)

t2d_exp <- readRDS("T2D_processed.exposure_gwas.Rds")
t2d_exp$pos.exposure <- unlist(lapply(
    strsplit(t2d_exp$SNP,":"),
    FUN = function(x){as.numeric(x[[2]])}
))

t2d_IV_inhbc_region_to_remove <- t2d_exp %>% filter(
    # retain variants beyond 1mb region on either side of the transcript
    chr.exposure == inhbc_chr &
    pos.exposure >= inhbc_transcript[1] - 1000000 &
    pos.exposure <= inhbc_transcript[2] + 1000000
)
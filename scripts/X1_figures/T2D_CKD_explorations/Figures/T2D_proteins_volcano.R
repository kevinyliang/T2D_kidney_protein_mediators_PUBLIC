#' Make a volcano plot for all the T2D on proteins results
#' 
#' 
library(tidyverse)
library(glue)
library(data.table)
library(yaml)
library(ggrepel)
yml_file  <-  "configs/proj_config.yml"
args  <- read_yaml(yml_file)[['T2D_CKD_results']]


tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')

final_dir <- "T2D_CKD_results/figures"
final_file <- file.path(final_dir,"Fig1_t2d_proteins_volcano.pdf")
final_file_legend <- file.path(final_dir,"Fig1_volcano_legend.pdf")

sig_t2d_prot_complication_proteins  <- read.table(
  file.path(args$rt_dir,"T2D_CKD_results","Tab_T2D_protein_complications_trios.tsv"),sep="\t",header=T,colClass='character'
) %>% filter(passed_all_QC=='TRUE') %>% pull(Somalogic.seqID) %>% unique()

all_t2d_on_prot_mr_res  <- read.table(
  file.path(args$rt_dir,"T2D_CKD_results","Tab_all_T2D_protein_MR.tsv"),sep="\t",header=T,colClass='character'
)
all_t2d_on_prot_mr_res$bonferroni_p  <- as.numeric(all_t2d_on_prot_mr_res$bonferroni_p)
all_t2d_on_prot_mr_res$b  <- as.numeric(all_t2d_on_prot_mr_res$b)
all_t2d_on_prot_mr_res$pval  <- as.numeric(all_t2d_on_prot_mr_res$pval)
plot_data  <- all_t2d_on_prot_mr_res %>% 
  dplyr::select(
    b,bonferroni_p,passed_all_QC,pval,Somalogic.seqID,Entrez.gene.name
  )
plot_data  <- plot_data %>% 
  dplyr::mutate(
    color = ifelse(
      passed_all_QC == "TRUE",'red',ifelse(
        bonferroni_p < 0.05,'pink', ifelse(
          pval < 0.05,'blue', ifelse(
            pval >= 0.05, 'black','green'
          )
        )
      )
    ),
    shape = ifelse(
      b < 0 , "25", ifelse(
        b > 0, "24","21"
      )
    )
  )
  
most_sig_red  <- plot_data %>% filter(passed_all_QC == "TRUE") %>% slice_min(bonferroni_p) %>% pull(Somalogic.seqID) %>% unique()

assertthat::assert_that(
  plot_data %>% filter(passed_all_QC=='TRUE') %>% pull(bonferroni_p) %>% max() < 0.05 &
  !"green"%in% plot_data$color
)

p  <- ggplot(plot_data) +
  geom_point(aes(
    x = b,y=-log10(pval),color=color,fill=color
  ))+
  scale_color_manual(
    values = c("black"="black",'pink'='pink','red'='red','blue'='blue')
  )+
  scale_fill_manual(
    values = c("black"="black",'pink'='pink','red'='red','blue'='blue')
  )+
  scale_shape_manual(values=c("24"=24,"25"=25,'21'=21))+
  geom_hline(yintercept=-log10(0.05),linetype=2,color='blue')+
  geom_hline(yintercept=-log10(0.05/nrow(plot_data)),linetype=2,color='red')+
  theme_classic() +
  labs(
    x='MR effect size estimate\n(Protein levels are in standardized units)',
    y=bquote(-log[10](p-value))
  )+geom_label_repel(
    data = plot_data %>% filter(Somalogic.seqID %in% c(sig_t2d_prot_complication_proteins,most_sig_red)),
    aes(x = b,y=-log10(pval),label = Entrez.gene.name), box.padding = 0.5, point.padding = 0.5,
    min.segment.length = unit(0, 'lines')
  )+
  theme(
    text=element_text(size=20),
    legend.position = 'none'
  ) 
ggsave(
  filename = tfile,
  plot = p ,
  device = 'png',units='in',height=7,width=6
)
ggsave(
  filename = final_file,
  plot = p ,
  device = 'pdf',units='in',height=7,width=6
)

p  <- ggplot(plot_data) +
  geom_point(aes(
    x = b,y=-log10(pval),color=color,fill=color
  ))+
  scale_color_manual(
    values = c("black"="black",'pink'='pink','red'='red','blue'='blue')
  )+
  scale_fill_manual(
    values = c("black"="black",'pink'='pink','red'='red','blue'='blue')
  )+
  scale_shape_manual(values=c("24"=24,"25"=25,'21'=21))+
  geom_hline(yintercept=-log10(0.05),linetype=2,color='blue')+
  geom_hline(yintercept=-log10(0.05/nrow(plot_data)),linetype=2,color='red')+
  theme_classic() +
  labs(
    x='MR effect size estimate\n(Protein levels are in standardized units)',
    y=bquote(-log[10](p-value))
  )+geom_label_repel(
    data = plot_data %>% filter(Somalogic.seqID %in% c(sig_t2d_prot_complication_proteins,most_sig_red)),
    aes(x = b,y=-log10(pval),label = Entrez.gene.name), box.padding = 0.5, point.padding = 0.5
  )+
  theme(
    text=element_text(size=20)
  ) 
ggsave(
  filename = final_file_legend,
  plot = p ,
  device = 'pdf',units='in',height=7,width=6
)

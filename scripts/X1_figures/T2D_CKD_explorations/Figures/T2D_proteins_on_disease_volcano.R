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

tdir <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,'Fig1_t2d_proteins_on_disease_volcano.pdf')

all_t2d_prot_complication_trios  <- read.table(
  file.path(args$rt_dir,"T2D_CKD_results","Tab_T2D_protein_complications_trios.tsv"),sep="\t",header=T,colClass='character'
)
all_t2d_prot_complication_trios$bonferroni_p  <- as.numeric(all_t2d_prot_complication_trios$bonferroni_p)
all_t2d_prot_complication_trios$b  <- as.numeric(all_t2d_prot_complication_trios$b)
all_t2d_prot_complication_trios$pval  <- as.numeric(all_t2d_prot_complication_trios$pval)
plot_data  <- all_t2d_prot_complication_trios %>% 
  dplyr::select(
    Entrez.gene.name,b,bonferroni_p,passed_all_QC,pval
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
    )
  )

assertthat::assert_that(
  plot_data %>% filter(passed_all_QC=='TRUE') %>% pull(bonferroni_p) %>% max() < 0.05 &
  !"green"%in% plot_data$color
)

p  <- ggplot(plot_data) +
  geom_point(aes(
    x = b,y=-log10(pval),color=color
  ))+
  scale_color_manual(
    values = c("black"="black",'pink'='pink','red'='red','blue'='blue')
  )+
  geom_hline(yintercept=-log10(0.05),linetype=2,color='blue')+
  geom_hline(yintercept=-log10(0.05/nrow(plot_data)),linetype=2,color='red')+
  geom_label_repel(
    data = plot_data %>% filter(passed_all_QC=='TRUE'),
    aes(x = b,y=-log10(pval),label = Entrez.gene.name), box.padding = 0.5, point.padding = 0.5
  ) +
  theme_classic() +
  labs(
    x='MR effect size estimate\n(Protein levels are in standardized units)',
    y=bquote(-log[10](p-value))
  )+
  theme(
    text=element_text(size=20),
    legend.position = 'none'
  )  +
  coord_cartesian(clip = "off")
ggsave(
  filename = tfile,
  plot = p,
  device = 'png',units='in',height=7,width=6
)
ggsave(
  filename = final_file,
  plot = p,
  device = 'pdf',units='in',height=7,width=6
)

#' Make forest plot of significant protein on otucome and T2D on protein
library(ggforestplot)
library(yaml)
library(data.table)
library(glue)
library(tidyverse)
library(patchwork)

tdir  <- tempdir()
tfile  <- file.path(tdir,'tmp.png')
final_dir <- "figures"
final_file <- file.path(final_dir,'FigS2_summary_of_t2d_protein_DKD_results.pdf')


args  <- read_yaml("configs/proj_config.yml")[['T2D_CKD_results']]
 
outdir_base  <- "figures"

rel_t2d_prot_on_outcome  <- read.table(
  "Tab_T2D_protein_complications_trios.tsv",sep="\t",quote="",header=T,colClasses=c("passed_all_QC"="character")
) %>% filter(passed_all_QC=='TRUE')

t2d_on_prot_res  <- read.table(
  "Tab_all_T2D_protein_MR.tsv",sep="\t",header=T,quote="",colClasses=c("passed_all_QC"="character")
) %>% filter(passed_all_QC=='TRUE' & Entrez.gene.name %in% rel_t2d_prot_on_outcome$Entrez.gene.name)

rel_t2d_prot_on_outcome  <- rel_t2d_prot_on_outcome %>% 
    dplyr::select(Entrez.gene.name,complication,b,se,pval) %>% 
    mutate(
      direction = 'Protein on isease'
    )
t2d_on_prot_res <- t2d_on_prot_res  %>% 
    dplyr::select(Entrez.gene.name,b,se,pval) %>% 
    mutate(
      complication = 'T2D',Entrez.gene.name,
      direction = 'Disease on protein'
    )
rel_t2d_prot_on_outcome <- rel_t2d_prot_on_outcome[
  order(
    match(
      rel_t2d_prot_on_outcome$Entrez.gene.name,t2d_on_prot_res$Entrez.gene.name
    )
  ),
]
rel_t2d_prot_on_outcome$Entrez.gene.name <- as_factor(rel_t2d_prot_on_outcome$Entrez.gene.name)
t2d_on_prot_res$Entrez.gene.name <- as_factor(t2d_on_prot_res$Entrez.gene.name)



p1  <- ggforestplot::forestplot(
  df = t2d_on_prot_res,
  name = Entrez.gene.name,
  estimate = b,
  pvalue = pval,
  colour= Entrez.gene.name
) + scale_color_manual(
    values = c(
      "GNPTG"="#e7b323",
      "AGRN"="#7f70ff",
      "LPO"="#17b5b3",
      "INHBC"="#fb303f",
      "CTSD"="#000000"
    )
  )+
  theme_classic() +
  labs(
    x='Estimated casual effects of T2D on proteins',
    y='Outcome',
    colour = 'Outcome'
  ) +
  theme(text=element_text(size=20),legend.position='none')
p2  <- ggforestplot::forestplot(
  df = rel_t2d_prot_on_outcome,
  name = Entrez.gene.name,
  estimate = b,
  pvalue = pval,
  colour= Entrez.gene.name
) + scale_color_manual(
    values = c(
      "GNPTG"="#e7b323",
      "AGRN"="#7f70ff",
      "LPO"="#17b5b3",
      "INHBC"="#fb303f",
      "CTSD"="#000000"
    )
  )+
facet_wrap(~complication,scale='free',ncol=2,nrow=2)+ 
  theme_classic() +
  labs(
    x='Estimated causal effects of protein on traits',
    y='Exposure'
  ) +
  theme(
    text=element_text(size=20),
    legend.position='none',
    plot.margin = unit(c(1,1,1,1),'cm')
  )


ggsave(
  filename = tfile,
  plot = p1/p2,
  device='png',width=11,height=10,units='in'
)
ggsave(
  filename = final_file,
  plot = p1/p2,
  device='pdf',width=11,height=10,units='in'
)

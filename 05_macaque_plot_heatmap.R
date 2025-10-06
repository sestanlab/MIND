###################################################################
# Time: 2025-08-15 15:33:55 EDT                                   #
# Author: Xinyun LI                                               #
# Email: xinyun.li@yale.edu                                       #
# Description: This script includes functions to plot             #
#              normal/distorted heatmaps of                       # 
#              gene modules using macaque micro-array data        #
###################################################################

library(tidyverse)
library(readxl)
library(cowplot)
source("./macaque.plot.func.R")

outdir.heatmap <- paste0("outs/heatmap")
if(!dir.exists(outdir.heatmap)){
  dir.create(outdir.heatmap,recursive = T)
}

# prepare data
rs.order <- rev(c('OFC','MFC','DFC','VFC','M1C','S1C','IPC','V1C','A1C','STC', 'ITC','AMY', "HIP"))
# load data
## zhu18_monkey
fpkm_zhu18 <- read.table('PATH_TO_Zhu2018_RNA_seq_data/nhp_development_RPKM_rmTechRep.txt',header = T) %>%
  mutate(Geneid = rownames(.))
meta_zhu18 <- readRDS('data/meta.zhu2018.rds',header = T) %>% 
  dplyr::filter(Regioncode %in% rs.order)


# load GM
gm.dir <- "outs/Final_GM"
fname <- "shared"
print(fname)
df <- openxlsx::read.xlsx(paste0(gm.dir, "/",fname,".prenatal.GM.xlsx"))%>%
  dplyr::filter(Type%in%c("Af","At","Af n At","S"))%>%
  mutate(Type=factor(Type, levels=c("Af","At","Af n At","S")))%>%arrange(Type)
module.gene.list <- split(df$GeneSymbol, df$Type)

# plotting  #################  
## normal heatmap p4 to p9
n.hp.gm <- plot.GM.normal.heatmap.p4to9.zhu18(
                                    module.gene.list=module.gene.list,
                                    rs.order=rs.order,
                                    zhu18.expr=fpkm_zhu18, 
                                    zhu18.meta=meta_zhu18)
n.n.gm=length(n.hp.gm)
pdf(paste0(outdir.heatmap, "/EF2b_",fname,".p4to9.normal.heatmap.pdf"), width = 6.8*n.n.gm, height = 6)
n.hp.gm.merge = plot_grid(plotlist = n.hp.gm,  ncol = 4)
print(n.hp.gm.merge)
dev.off()

## normal heatmap p4 to p9 without smoothening
n.hp.gm <- plot.GM.normal.heatmap.p4to9.zhu18.nosmooth(
                                    module.gene.list=module.gene.list,
                                    rs.order=rs.order,
                                    zhu18.expr=fpkm_zhu18, 
                                    zhu18.meta=meta_zhu18)
n.n.gm=length(n.hp.gm)
pdf(paste0(outdir.heatmap, "/EF2b_",fname,".p4to9.normal.heatmap.nosmooth.pdf"), width = 6.8*n.n.gm, height = 6)
n.hp.gm.merge = plot_grid(plotlist = n.hp.gm,  ncol = 4)
print(n.hp.gm.merge)
dev.off()

## normal heatmap p4 to p9 without smoothening
n.hp.gm <- plot.GM.normal.heatmap.p4to9.zhu18.nosmooth.noimpute(
                                    module.gene.list=module.gene.list,
                                    rs.order=rs.order,
                                    zhu18.expr=fpkm_zhu18, 
                                    zhu18.meta=meta_zhu18)
n.n.gm=length(n.hp.gm)
pdf(paste0(outdir.heatmap, "/EF2b_",fname,".p4to9.normal.heatmap.nosmooth.noimpute.pdf"), width = 6.8*n.n.gm, height = 6)
n.hp.gm.merge = plot_grid(plotlist = n.hp.gm,  ncol = 4)
print(n.hp.gm.merge)
dev.off()


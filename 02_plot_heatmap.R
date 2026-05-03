###################################################################
# Time: 2025-08-13 14:18:00 EDT                                   #
# Author: Menglei ZHANG, Xinyun LI                                #
# Email: morei.menglei.zhang@yale.edu                             #
#        xinyun.li@yale.edu                                       #
# Description: This script is to plot normal smoothened and       #
#              unsmoothened heatmaps of gene modules,             #
#              and genes using micro-array data                   #
###################################################################

library(tidyverse)
library(openxlsx)
source("./plot.func.R")


outdir.heatmap <- paste0("outs/heatmap")
if(!dir.exists(outdir.heatmap)){
  dir.create(outdir.heatmap,recursive = T)
}

# prepare data
rs.order <- rev(c('OFC','MFC','DFC','VFC','M1C','S1C','IPC','V1C','A1C','STC', 'ITC','AMY', "HIP"))
# load data
# kang11
exon11 <- read.table('PATH_TO_Kang2011_exon_microarray_data/nature105213.GENE1346_nosnp_Log2-expression value.txt',header = T)
meta11 <- readRDS("data/meta.kang2011.rds")

# load GM 
gm.dir <- "outs/Final_GM"
fname <- "shared"
print(fname)
df <- openxlsx::read.xlsx(paste0(gm.dir, "/",fname,".prenatal.GM.xlsx"))%>%
  dplyr::filter(Type%in%c("Af","At","Af n At","S"))%>%
  mutate(Type=factor(Type, levels=c("Af","At","Af n At","S")))%>%arrange(Type)
module.gene.list <- split(df$GeneSymbol, df$Type)

# plotting  #################  
## normal heatmap p3 to p9
### GM 
### smoothened
n.hp.gm <- plot.GM.normal.heatmap.p3to9.kang11(module.gene.list=module.gene.list,
                                               rs.order=rs.order,
                                               kang11.expr=exon11,
                                               kang11.meta=meta11)
n.n.gm <- length(n.hp.gm)
pdf(paste0(outdir.heatmap, "/MF1b_",fname,".p3to9.normal.heatmap.pdf"), width = n.n.gm*6+9, height = 6)
library(cowplot) # Arrange the plots in a custom grid layout
n.hp.gm.merge <- plot_grid(plotlist = n.hp.gm,ncol = 4)  
print(n.hp.gm.merge)
dev.off() # Close the PDF device

### unsmoothened
n.hp.gm <- plot.GM.normal.heatmap.p3to9.kang11.nosmooth(module.gene.list=module.gene.list,
                                               rs.order=rs.order,
                                               kang11.expr=exon11,
                                               kang11.meta=meta11)
n.n.gm <- length(n.hp.gm)
pdf(paste0(outdir.heatmap, "/MF1b_",fname,".p3to9.normal.heatmap.nosmooth.pdf"), width = n.n.gm*6+9, height = 6)
library(cowplot) # Arrange the plots in a custom grid layout
n.hp.gm.merge <- plot_grid(plotlist = n.hp.gm,ncol = 4)  
print(n.hp.gm.merge)
dev.off() # Close the PDF device

### gene 
### smoothened
gs <- c('SEMA7A', 'PLXNC1')
n.hp.g <- plot.gene.normal.heatmap.p3to9.kang11(gene.list=gs,
                                                    rs.order=rs.order,
                                                    kang11.expr=exon11,
                                                    kang11.meta=meta11)
n.n.g <- length(n.hp.g)
pdf(paste0(outdir.heatmap, "/MF2f_SEMA7A_PLXNC1.p3to9.normal.heatmap.pdf"), width = n.n.g*6+3.5, height = 6)
library(cowplot)
n.hp.g.merge <- plot_grid(plotlist = n.hp.g, ncol = 2)  
print(n.hp.g.merge)
dev.off()

### unsmoothened
gs <- c('SEMA7A', 'PLXNC1')
n.hp.g <- plot.gene.normal.heatmap.p3to9.kang11.nosmooth(gene.list=gs,
                                                    rs.order=rs.order,
                                                    kang11.expr=exon11,
                                                    kang11.meta=meta11)
n.n.g <- length(n.hp.g)
pdf(paste0(outdir.heatmap, "/MF2f_SEMA7A_PLXNC1.p3to9.normal.heatmap.nosmooth.pdf"), width = n.n.g*6+3.5, height = 6)
library(cowplot)
n.hp.g.merge <- plot_grid(plotlist = n.hp.g, ncol = 2)  
print(n.hp.g.merge)
dev.off()






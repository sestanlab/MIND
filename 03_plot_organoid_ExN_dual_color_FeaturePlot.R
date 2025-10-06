############################################################################
# Time: 2025-08-14 15:14:04 EDT                                            #
# Author: Menglei ZHANG                                                    #
# Email: morei.menglei.zhang@yale.edu                                      #
# Description: This script is to plot dual color pair in one FeaturePlot.  # 
#              There are two examples:                                #
#              1) GM module pair                                           # 
#              2) gene pair                                                # 
############################################################################

library(tidyverse)
library(Matrix)
library(data.table)
library(Seurat)
source('./plot.func.R')

outDir <- 'outs/Uzquiano2022_Arlotta_6m'
if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

# load single cell data of 6month
data.dir <- 'PATH_TO_Uzquiano2022_Arlotta_data'
exp.flile <- list.files(data.dir,pattern = "^expression")    
# for(x in flist){
expr.id <- str_split_fixed(exp.flile,".txt",2)[1]
id <- str_split_fixed(expr.id,"_",3)[2]

print(id)
seufile <- paste0(data.dir,'/Uzquiano2022_Arlotta_',expr.id,'.rds')
if(!file.exists(seufile)){
  fname=paste0('Uzquiano2022_Arlotta_', id)
  mtx <- readMM(paste0(data.dir,'/',exp.flile))
  umap <-fread(paste0(data.dir,'/umap_',id,".txt"))%>%
    dplyr::filter(NAME!="TYPE")
  meta <- umap%>%dplyr::select(-c(X,Y))
  umap.mtx <- umap%>%dplyr::select(c(X,Y))%>%dplyr::rename("UMAP_1"="X","UMAP_2"="Y")%>%
    sapply(.,as.numeric)%>%as.matrix()
  rownames(umap.mtx) <- umap$NAME
  gs <- fread(paste0(data.dir,'/NormExpression_',id,"_genes.txt"),header = F)$V1
  seu <- CreateSeuratObject(counts = mtx, project = "Uzquiano2022",meta.data = meta)%>%NormalizeData()
  colnames(seu) <- seu$NAME
  rownames(seu) <- gs
  seu[["umap"]] <- CreateDimReducObject(embeddings = umap.mtx, key = "UMAP_", assay = "RNA")
  Idents(seu) <- "CellType"
  saveRDS(seu, seufile)
}else{
  seu <- readRDS(seufile)
}


addmodule.GM.file <- paste0(outDir,'/Uzquiano2022_Arlotta_',id,'_PN_AddModule.GM.seu.rds')
if(!file.exists(addmodule.GM.file)){
  # load GM
  gm.dir <- "outs/Final_GM"
  fname <- "shared"
  print(fname)
  df <- openxlsx::read.xlsx(paste0(gm.dir, "/",fname,".prenatal.GM.xlsx"))%>%
    dplyr::filter(Type%in%c("Af","At","S"))%>%
    mutate(Type=factor(Type, levels=c("Af","At","S")))%>%arrange(Type)
  module.gene.list <- split(df$GeneSymbol, df$Type)
  
  # get ExN cells
  cls.PN <- grep("PN$",unique(seu$CellType),value=T)
  seu.PN <- subset(seu, subset=(CellType%in%cls.PN)&GAD1 ==0)%>%NormalizeData()%>%
    FindVariableFeatures()%>%
    ScaleData()%>%
    RunPCA(., verbose= FALSE)%>%
    RunUMAP(., reduction = "pca", dims = 1:30)
  for(i in 1:length(module.gene.list)){
    if(i==1){
      seu.PN.gm <- seu.PN
    }else{
      seu.PN.gm <- seu.PN.gm
    }
    nm <- names(module.gene.list)[i]
    gm.gs <- module.gene.list[i]
    seu.PN.gm <- AddModuleScore(seu.PN.gm, features = gm.gs, name = paste0("module_",nm))
    colnames(seu.PN.gm@meta.data) <-
      gsub(x = colnames(seu.PN.gm@meta.data)
           , pattern = paste0("module_",nm,1)
           , replacement =paste0("module_",nm))
  }
  
  saveRDS(seu.PN.gm, addmodule.GM.file)
}else{seu.PN.gm <- readRDS(addmodule.GM.file)}

# Dual color FeaturePlot
## GM
fplots.gm <- Dual.Color.FeaturePlot.GM(seu=seu.PN.gm,
                                       outDir=outDir
                                       )

n.fp.gm <- length(fplots.gm)
pdf(paste0(outDir,'/Uzquiano2022_Arlotta_',id, '.Dual.Color.FeaturePlot.GM.pdf'), width = 4*n.fp.gm, height = 4)
library(cowplot)
plot.m <- plot_grid(plotlist = fplots.gm, nrow = 1) 
print(plot.m)
dev.off()


## gene pair
gene.pair=c("SEMA7A","PLXNC1")
fplots.g.pair <-Dual.Color.FeaturePlot.gene.pair(seu=seu.PN.gm, 
                                 gene.pair=gene.pair, 
                                 outDir='outs'
                                 )

pdf(paste0(outDir,'/Uzquiano2022_Arlotta_',id, '.Dual.Color.FeaturePlot.gene.pair',paste(gene.pair,collapse = "_"),'.pdf'), width = 4, height = 4)
print(fplots.g.pair)
dev.off()
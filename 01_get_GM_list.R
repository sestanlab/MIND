###################################################################
# Time: 2025-08-08 17:18:23 EDT                                   #
# Author: Menglei ZHANG                                           #
# Email: morei.menglei.zhang@yale.edu                             #
# Description: This script is to calculate module genes for bulk  # 
#              RNA-seq and micro-array data.                      #
###################################################################

library(tidyverse)
source('./GM.cal.func.R')

# get primary gene module for each data set follows the criteria:

##(1) there is at least one region in that subgroup where the gene is significantly upregulated; 
##(2) the gene is not upregulated in any region of the opposing (S/A) groups; 
##(3) the gene is under-expressed in at least 30% of the areas in the opposing (S/A) group
##(4) the gene is not in a module gene list of any of the opposing (S/A) subgroups

## li18
count_li18 <- read.table('PATH_TO_Li2018_RNA_seq_data/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt',header = T)
meta_li18 <- readRDS("data/meta.li2018.rds")

gm_li18 <- get.GM.RNentropy(count=count_li18, 
                            meta=meta_li18, 
                            periods=c(3:7), 
                            dataset.name="Li2018", 
                            outdir="outs")

## Zhu18
count_zhu18 <- read.table('PATH_TO_Zhu2018_RNA_seq_data/nhp_development_count_rmTechRep.txt',header = T)%>%
  tibble::rownames_to_column(.,var="Geneid")
meta_zhu18 <- readRDS("data/meta.zhu2018.rds")

gm_zhu18 <- get.GM.RNentropy(count=count_zhu18, 
                             meta=meta_zhu18, 
                             periods=c(5:7), 
                             dataset.name="Zhu2018", 
                             outdir="outs")

## kang11
exon_kang11 <- read.table('PATH_TO_Kang2011_exon_microarray_data/nature105213.GENE1346_nosnp_Log2-expression value.txt',header = T)%>%
  mutate(Geneid=paste(ID, Gene.Symbol, sep="|"))  
meta_kang11 <- readRDS("data/meta.kang2011.rds")

gm_kang11 <- get.GM.limma(count=exon_kang11,
                          meta=meta_kang11,
                          periods=c(3:7),
                          dataset.name="Kang2011", 
                          outdir="outs")

# get final gene modules for each group with following steps:

##(1) 2nd round of criteria (4) the gene is not in a module gene list of any of the opposing (S/A) subgroups
##(2) get the union and shared genes within A/S modules 
##(3) Rename modules

## shared
ind.GM.list.shared <- list("kang11"=gm_kang11,
                           "li18"=gm_li18,
                           "zhu18"=gm_zhu18)
gm.final.shared <- get.final.GM(ind.GM.list=ind.GM.list.shared, 
                                fname="shared", 
                                outdir="outs/Final_GM")

## kang11
ind.GM.list.kang11 <- list("kang11"=gm_kang11)
gm.final.kang11 <- get.final.GM(ind.GM.list=ind.GM.list.kang11, 
                                fname="kang11", 
                                outdir="outs/Final_GM")

## li18
ind.GM.list.li18 <- list("li18"=gm_li18)
gm.final.li18 <- get.final.GM(ind.GM.list=ind.GM.list.li18, 
                              fname="li18", 
                              outdir="outs/Final_GM")

## zhu18
ind.GM.list.zhu18 <- list("zhu18"=gm_zhu18)
gm.final.zhu18 <- get.final.GM(ind.GM.list=ind.GM.list.zhu18, 
                               fname="zhu18", 
                               outdir="outs/Final_GM")


###################################################################
# Time: 2025-08-08 17:18:23 EDT                                   #
# Author: Menglei ZHANG                                           #
# Email: morei.menglei.zhang@yale.edu                             #
# Description: This script includes functions to calculate        # 
#              module genes for bulk RNA-seq and micro-array data #
###################################################################

#############################  parameters   ############################
# 1. count should contain a column: "Geneid", and other column names   # 
#       should be same as Sample in meta data                          #
# 2. meta should contain columns: "Sample", "Regioncode", "Period"     # 
#############################  parameters   ############################

###########  Function  ###################
library(tidyverse)
library(RNentropy)
library(edgeR)
library(tweeDEseq)

get.GM.RNentropy <- function(count=NULL, meta=NULL, periods=c(3:7), dataset.name="Li2018", outdir="outs/Individual_dataset_primary_GM"){
  print(paste0("Processing DEGs of dataset: ", dataset.name))
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  
  # prepare data
  module_region <- data.frame(Module=c(rep("Peri-F",4),rep("Peri-T",2),"Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C"),
                                         Regioncode=c(c('DFC','OFC','VFC','MFC'),c('ITC','STC'),'M1C','S1C','V1C','A1C'))
  region.num <- length(unique(module_region$Regioncode))
  print(paste0("There are ",region.num, " Regions from period ", min(periods), " to period ", max(periods)))
  
  # calculate module gene
  GM.file <- paste0(outdir,'/',dataset.name,'.period.',min(periods),'to',max(periods),".",region.num,'Regioncode.GM.rds')
  if(!file.exists(GM.file)){
    res <- lapply(periods, function(period){
      print(period)
      meta.period <- meta%>%dplyr::filter(Period%in%period,Regioncode%in%module_region$Regioncode)
      design.period <- meta.period%>%dplyr::select(Sample,Regioncode)%>%distinct()%>%
        mutate(Regioncode=factor(Regioncode, levels=unique(module_region$Regioncode)))%>%arrange(Regioncode)%>%
        mutate(IN=1)%>%pivot_wider(names_from = Regioncode, values_from = IN)%>%mutate_all(~replace(., is.na(.), 0))%>%
        tibble::column_to_rownames(.,var = 'Sample')%>%as.matrix()
      count.period <- count[rownames(design.period)]
      rownames(count.period) <- count$Geneid
      ## get DEGs among regions using RNentropy
      keep <- edgeR::filterByExpr(count.period)
      count_subset <- count.period[keep, ]
      count_norm.period <- (tweeDEseq::normalizeCounts(count_subset))%>%
        as.data.frame()%>%mutate(Geneid=rownames(.))
      Results.pmi <- RN_pmi(RN_calc(count_norm.period, design = design.period))
      selected.ori <- Results.pmi$selected
      selected.mtx <- selected.ori[,intersect(colnames(selected.ori),module_region$Regioncode)]
      selected <- selected.ori[!apply(selected.mtx, 1, function(x) all(x == 0) | all(is.na(x))), ]
      
      ## get module gene 
      module.gs.df <- lapply(unique(module_region$Module), function(module) {
        print(module)
        rs <- module_region%>%dplyr::filter(Module==module, Regioncode%in%colnames(selected))%>%dplyr::pull(Regioncode)
        rs.Asso <- module_region%>%dplyr::filter(Module%in%c("Peri-F","Peri-T"), Regioncode%in%colnames(selected))%>%dplyr::pull(Regioncode)
        rs.Pri <- module_region%>%dplyr::filter(Module%in%c("Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C"), Regioncode%in%colnames(selected))%>%dplyr::pull(Regioncode)
        # the definition requires two group data
        if((length(rs.Asso)>0)&(length(rs.Pri)>0)){
          if(module%in%c("Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C")){
            gs.rs.Up <- selected[,c("Geneid",rs)][selected[,rs]==1,]%>%drop_na()
            gs.rs.others.Up <- selected[,rs.Asso]%>%mutate(across(everything(), ~replace(., is.na(.)|(.==-1), 0)))
            gs.rs.others.notUp <- gs.rs.others.Up[apply(gs.rs.others.Up, 1, function(row) sum(row)==0),]
            ttt <- selected[rownames(gs.rs.others.notUp),]
            gs.rs.others.Down <- selected[,rs.Asso]%>%mutate(across(everything(), ~replace(., is.na(.)|(.==1), 0)))
            gs.rs.others.30Down <- gs.rs.others.Down[apply(gs.rs.others.Down, 1, function(row) sum(row))<=round((-1*0.3*length(rs.Asso))),]
          }else{
            gs.rs.Up <- selected[,c("Geneid",rs)][apply(selected[,rs], 1, function(row) 1%in%row),]
            gs.rs.others.Up <- selected[,rs.Pri]%>%mutate(across(everything(), ~replace(., is.na(.)|(.==-1), 0)))
            gs.rs.others.notUp <- gs.rs.others.Up[apply(gs.rs.others.Up, 1, function(row) sum(row)==0),]
            gs.rs.others.Down <- selected[,rs.Pri]%>%mutate(across(everything(), ~replace(., is.na(.)|(.==1), 0)))
            gs.rs.others.30Down <- gs.rs.others.Down[apply(gs.rs.others.Down, 1, function(row) sum(row))<=round((-1*0.3*length(rs.Pri))),]
          }
          
          gs.module.df <- gs.rs.Up%>%dplyr::filter(Geneid%in%intersect(rownames(gs.rs.others.notUp), rownames(gs.rs.others.30Down)))%>%
            mutate(Type=module,Period=period)%>%dplyr::select(Type, Period, Geneid)
        }else{
          gs.module.df <- NULL
        }
        return(gs.module.df)
      })%>%do.call(rbind,.)
      module.gs.asso <- module.gs.df%>%dplyr::filter(Type%in%c("Peri-F","Peri-T"))%>%dplyr::pull(Geneid)%>%unique()
      module.gs.pri <- module.gs.df%>%dplyr::filter(Type%in%c("Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C"))%>%dplyr::pull(Geneid)%>%unique()
      # remove shared genes
      module.gs.share <- intersect(module.gs.asso, module.gs.pri)
      module.gs.df.f <- module.gs.df%>%dplyr::filter(!Geneid%in%module.gs.share)
      return(module.gs.df.f)
    })%>%do.call(rbind,.)%>%mutate(GeneSymbol=str_split_fixed(Geneid, "[|]",2)[,2])

    # save results to file
    saveRDS(res, GM.file)
  }else{
    res <- readRDS(GM.file)
  }
  return(res)
}

get.GM.limma <- function(count=NULL, meta=NULL, periods=c(3:7), dataset.name="Kang2011", outdir="outs/Individual_dataset_primary_GM"){
  print(paste0("Processing DEGs of dataset: ", dataset.name))
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  
  # prepare data
  module_region <- data.frame(Module=c(rep("Peri-F",4),rep("Peri-T",2),"Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C"),
                              Regioncode=c(c('DFC','OFC','VFC','MFC'),c('ITC','STC'),'M1C','S1C','V1C','A1C'))
  region.num <- length(unique(module_region$Regioncode))
  print(paste0("There are ",region.num, " Regions from period ", min(periods), " to period ", max(periods)))
  # calculate module gene
  GM.file <- paste0(outdir,'/',dataset.name,'.period.',min(periods),'to',max(periods),".",region.num,'Regioncode.GM.rds')
  if(!file.exists(GM.file)){
    res <- lapply(periods, function(period){
      print(period)
      module.gs.df <- lapply(unique(module_region$Module), function(module){
        print(module)
        rs <- module_region%>%dplyr::filter(Module==module)%>%dplyr::pull(Regioncode)
        rs.Asso <- module_region%>%dplyr::filter(Module%in%c("Peri-F","Peri-T"))%>%dplyr::pull(Regioncode)
        rs.Pri <- module_region%>%dplyr::filter(Module%in%c("Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C"))%>%dplyr::pull(Regioncode)
        if(module%in%c("Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C")){
          rs.ctr <- rs.Asso    # set control group for rs
          rs.ctr.ctr <- rs.Pri # set control group for rs.ctr
        }else{
          rs.ctr <- rs.Pri
          rs.ctr.ctr <- rs.Asso
        }
        print("----get up-regulated genes for module regions-----")
        deg.rs.Up <- lapply(rs, function(rgn){
          print(rgn)
          group.dic <- c(rep("treatment", length(rgn)), rep("control",length(rs.ctr)))%>%setNames(c(rgn, rs.ctr))
          meta.p.rs <- meta%>%dplyr::filter(Period%in%period,Regioncode%in%c(rgn, rs.ctr))%>%
            mutate(group=group.dic[Regioncode])
          design.p.rs <- meta.p.rs%>%dplyr::select(Sample,group)%>%distinct()%>%
            mutate(group=factor(group, levels=c("treatment","control")))%>%arrange(group)%>%
            mutate(IN=1)%>%pivot_wider(names_from = group, values_from = IN)%>%mutate_all(~replace(., is.na(.), 0))%>%
            tibble::column_to_rownames(.,var = 'Sample')%>%as.matrix()
          count.p.rs <- count[rownames(design.p.rs)]
          rownames(count.p.rs) <- count$Geneid
          
          # get up-regulated DEGs using limma
          contrasts_fit.rs.df <- get.DEGs.limma(count=count.p.rs, design=design.p.rs)
          selected.rs <- contrasts_fit.rs.df %>%
            dplyr::filter(adj.P.Val<0.05, logFC>log2(1))%>%
            tibble::rownames_to_column(.,var="Geneid")
          return(selected.rs$Geneid)
        })%>%unlist()%>%unique()
        
        print("----get non-UP genes for module countpart group-----")
        deg.ctr.UP <- list()
        deg.ctr.down <- list()
        for(rgn in rs.ctr){
          print(rgn)
          group.dic <- c(rep("treatment", length(rgn)), rep("control",length(rs.ctr.ctr)))%>%setNames(c(rgn, rs.ctr.ctr))
          meta.p.ctr <- meta%>%dplyr::filter(Period%in%period,Regioncode%in%c(rgn, rs.ctr.ctr))%>%
            mutate(group=group.dic[Regioncode])
          design.p.ctr <- meta.p.ctr%>%dplyr::select(Sample,group)%>%distinct()%>%
            mutate(group=factor(group, levels=c("treatment","control")))%>%arrange(group)%>%
            mutate(IN=1)%>%pivot_wider(names_from = group, values_from = IN)%>%mutate_all(~replace(., is.na(.), 0))%>%
            tibble::column_to_rownames(.,var = 'Sample')%>%as.matrix()
          count.p.ctr <- count[rownames(design.p.ctr)]
          rownames(count.p.ctr) <- count$Geneid
          # get DEGs using limma
          contrasts_fit.df <- get.DEGs.limma(count=count.p.ctr, design=design.p.ctr)
          # get up-regulated DEGs
          selected.up <- contrasts_fit.df%>%
            dplyr::filter(adj.P.Val<0.05, logFC>log2(1))%>%
            tibble::rownames_to_column(.,var="Geneid")
          deg.ctr.UP[[length(deg.ctr.UP)+1]] <- selected.up$Geneid
          
          # get down-regulated DEGs
          selected.down <- contrasts_fit.df%>%
            dplyr::filter(adj.P.Val<0.05, logFC< (-log2(1)))%>%
            mutate(Region=rgn, down=-1)%>%
            tibble::rownames_to_column(.,var="Geneid")%>%
            dplyr::select(Region,Geneid,down)
          deg.ctr.down[[length(deg.ctr.down)+1]] <- selected.down
        }
        deg.ctr.UP <- deg.ctr.UP%>%unlist()%>%unique()
        deg.ctr.noUP <- setdiff(count$Geneid, deg.ctr.UP)
  
        print("----get down-regulated genes-----")
        deg.ctr.down <- deg.ctr.down%>%do.call(rbind,.)
        deg.ctr.down.df <- deg.ctr.down%>%pivot_wider(names_from = "Region", values_from = "down")
        if(length(setdiff(rs.ctr, unique(deg.ctr.down$Region)))>0){
          deg.ctr.down.df[,setdiff(rs.ctr, unique(deg.ctr.down$Region))] <- NA
        }
        deg.ctr.down.df <- deg.ctr.down.df%>%mutate(across(everything(), ~replace(., is.na(.), 0)))
        deg.ctr.30Down <- deg.ctr.down.df[apply(deg.ctr.down.df[rs.ctr], 1, function(row) sum(row))<=round((-1*0.3*length(rs.ctr))),]
        gs.ctr.filter <- intersect(deg.ctr.noUP, deg.ctr.30Down$Geneid)
        gs.module <- intersect(deg.rs.Up,gs.ctr.filter)%>%unique()
        
        if(length(gs.module)>0){
          gs.module.df <- data.frame(Type=module,
                                   Period=period,
                                   Geneid=gs.module)
        }else{
          gs.module.df <- NULL
        }
        return(gs.module.df)
      })%>%do.call(rbind,.)
      # store results
      module.gs.asso <- module.gs.df%>%dplyr::filter(Type%in%c("Peri-F","Peri-T"))%>%dplyr::pull(Geneid)%>%unique()
      module.gs.pri <- module.gs.df%>%dplyr::filter(Type%in%c("Cent-M1C","Cent-S1C","Cent-V1C","Cent-A1C"))%>%dplyr::pull(Geneid)%>%unique()
      module.gs.share <- intersect(module.gs.asso, module.gs.pri)
      module.gs.df.f <- module.gs.df%>%dplyr::filter(!Geneid%in%module.gs.share)
      return(module.gs.df.f)
    })%>%do.call(rbind,.)%>%mutate(GeneSymbol=str_split_fixed(Geneid, "[|]",2)[,2])
    # save results to file
    saveRDS(res, GM.file)
  }else{
    res <- readRDS(GM.file)
  }
  return(res)
}

get.DEGs.limma <- function(count=NULL, design=NULL){
  # filter by Expr
  # Step 1: Filter based on intensity
  threshold <- log2(5)  # Define your expression threshold
  keep_intensity <- rowMeans(count > threshold) >= 0.3  # Keep genes expressed in at least 30% of samples
  
  # Step 2: Filter based on variance (over 10% variable genes)
  gene_variances <- apply(count, 1, var)
  var_threshold <- quantile(gene_variances, 0.1)
  keep_variance <- gene_variances > var_threshold
  
  # Step 3: Combine both filters
  count_subset <- count[keep_intensity & keep_variance, ]
  # DEGs using limma
  # Apply linear model to data
  fit <- limma::lmFit(count_subset, design = design)
  # Apply empirical Bayes to smooth standard errors
  fit <- limma::eBayes(fit)
  # Define the contrast for comparing treatment vs control
  contrast_matrix <- makeContrasts(treatment-control,levels = design)
  
  # Fit the model according to the contrasts matrix
  contrasts_fit <- contrasts.fit(fit, contrast_matrix)
  
  # Re-smooth the Bayes
  contrasts_fit <- eBayes(contrasts_fit)
  
  # View the results for the contrast
  # Apply multiple testing correction and obtain stats
  contrasts_fit.df <- topTable(contrasts_fit, adjust="fdr", number = Inf)
  return(contrasts_fit.df)
}

get.final.GM <- function(ind.GM.list=NULL, fname=NULL, outdir="outs/Final_GM"){
  library(tidyverse)
  library(openxlsx)
  
  if(!dir.exists(outdir)){
    dir.create(outdir,recursive = T)
  }
  
  # get (merged) GM follow Criteria (4)
  merged.gm.df <- NULL
  i <- 1
  for (nm in names(ind.GM.list)) {
    ind.gm <- ind.GM.list[[nm]]%>%
      mutate(dataset=nm)%>%dplyr::select(c(Type,GeneSymbol,dataset))%>%distinct()
    gs.asso <- ind.gm%>%dplyr::filter(grepl("Peri",Type))%>%dplyr::pull(GeneSymbol)%>%unique()
    gs.pri <- ind.gm%>%dplyr::filter(grepl("Cent",Type))%>%dplyr::pull(GeneSymbol)%>%unique()
    gs.share <- intersect(gs.asso, gs.pri)
    ind.gm.f <- ind.gm%>%dplyr::filter(!GeneSymbol%in%gs.share)
    if(i==1){
      merged.gm.df <- ind.gm.f%>%arrange(Type)%>%distinct()
    }else{
      names(ind.gm.f)[which(names(ind.gm.f)=="dataset")] <- "dataset.add"
      merged.gm.df <- merged.gm.df%>%left_join(ind.gm.f)%>%
        dplyr::filter((!is.na(dataset))&(!is.na(dataset.add)))%>%
        mutate(dataset=paste(dataset,dataset.add, sep="."))%>%
        dplyr::select(Type,GeneSymbol,dataset)%>%arrange(Type)%>%distinct()
    }
    i <- i+1
  }
  
  # get the union and shared genes within A/S modules
  m.gm.list <- split(merged.gm.df$GeneSymbol,merged.gm.df$Type)
  # Peri-union FT
  FT.union.gs <- unique(c(m.gm.list$`Peri-F`,m.gm.list$`Peri-T`))
  if(length(FT.union.gs)>0){
    FT.union.gs.df <- data.frame(Type="Peri-union FT",
                                 GeneSymbol=FT.union.gs,
                                 dataset=unique(merged.gm.df$dataset))
  }else(
    FT.union.gs.df <- NULL
  )
  # Peri-shared FT
  FT.shared.gs <- base::Reduce(intersect,m.gm.list[c("Peri-F","Peri-T")])
  if(length(FT.shared.gs)>0){
    FT.shared.gs.df <- data.frame(Type="Peri-shared FT",
                                  GeneSymbol=FT.shared.gs,
                                  dataset=unique(merged.gm.df$dataset))
  }else(
    FT.shared.gs.df <- NULL
  )
  
  # Cent-union SM
  SM.union.gs <- unique(c(m.gm.list$`Cent-A1C`,m.gm.list$`Cent-M1C`,m.gm.list$`Cent-S1C`,m.gm.list$`Cent-V1C`))
  if(length(SM.union.gs)>0){
    SM.union.gs.df <- data.frame(Type="Cent-union SM",
                                 GeneSymbol=SM.union.gs,
                                 dataset=unique(merged.gm.df$dataset))
  }else(
    SM.union.gs.df <- NULL
  )
  
  # Rename modules
  module.nm.dic <- c("Af","Af n At", "Af u At","At","S","Sa1","Sm1","Ss1","Sv1")%>%
    setNames(c("Peri-F","Peri-shared FT", "Peri-union FT","Peri-T","Cent-union SM","Cent-A1C","Cent-M1C","Cent-S1C","Cent-V1C")) 
  merged.gm.final.df <- do.call(rbind,list(merged.gm.df, FT.union.gs.df,FT.shared.gs.df,
                            SM.union.gs.df))%>%
    mutate(Type=factor(module.nm.dic[as.character(Type)], levels=unname(module.nm.dic)))%>%
    arrange(Type)%>%distinct()
  
  # save to file
  if(is.null(fname)){
    fname <-  unique(merged.gm.final.df$dataset)
  }
  openxlsx::write.xlsx(merged.gm.final.df, paste0(outdir, "/",fname,".prenatal.GM.xlsx"), rowNames = F, quote = F,stringsAsFactors = FALSE)
  return(merged.gm.final.df)
}
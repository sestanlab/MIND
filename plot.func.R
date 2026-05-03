###################################################################
# Time: 2025-08-13 15:33:55 EDT                                   #
# Author: Menglei ZHANG, Xinyun LI                                #
# Email: morei.menglei.zhang@yale.edu                             #
#        xinyun.li@yale.edu                                       #
# Description: This script includes functions to plot             #
#              normal 2D smoothened/unsmoothened heatmaps of      # 
#              gene modules, and genes using micro-array data     #
###################################################################


#############################  parameters   ######################################
# 1. kang11.expr should contain a column: "Gene.Symbol", and other column names  # 
#       should be same as Sample in meta data                                    #
# 2. kang11.meta should contain columns: "Sample", "Regioncode", "Period"        # 
#############################  parameters   ######################################


###########  Function  ###################
library(tidyverse)

plot.GM.normal.heatmap.p3to9.kang11 <- function(module.gene.list=NULL, rs.order=NULL, kang11.expr=NULL, kang11.meta=NULL){
  hp.plots <- list()
  for(module in names(module.gene.list)){
    module.gene <-  module.gene.list[[module]]%>%unique()
    if(length(module.gene)>0){
      # set plotting parameters
      ## set color gradient
      if(startsWith(module, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }else{
        hp.colors <- c(colorRampPalette(c("white","#8e0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }
      
      ## set axis labels
      intercepts <- log2(c(70, 91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
      labels <- 3:9
      intercepts.labels <- intercepts[-length(intercepts)]%>%setNames(labels)
      width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>%setNames(labels)
      age_label <- c("10PCW",13, 16, 19, "24PCW", "Birth", "0.5yr")
      midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
      midpoint_labels <- labels[1:length(midpoints)]
      
      rs.order.y <- 1:13 %>%setNames(rs.order)
      gs <- intersect(module.gene, exon11$Gene.Symbol)
      
      
      #kang11
      exon11.module <- exon11%>%dplyr::filter(grepl(paste0(gs,collapse = "|"), Gene.Symbol))%>%
        dplyr::select(grep("X[0-9]*",colnames(.),value = T))%>%
        colMeans(.)%>%
        data.frame()%>%setNames(module)%>%t()%>%as.data.frame()%>%
        pivot_longer(cols = all_of((colnames(.))), names_to = "Sample", values_to = "colMeans_log2expr")
      exon11.df <- kang11.meta%>%left_join(exon11.module,by="Sample",relationship = "many-to-many")%>%
        mutate(Period=as.numeric(Period))%>%
        dplyr::filter(Regioncode%in%rs.order, Period>2, Period<10)%>%
        group_by(Period, Regioncode)%>%
        mutate(log2expr_mean=median(colMeans_log2expr))%>%
        dplyr::select(Period, Regioncode,log2expr_mean)%>%
        distinct()%>%ungroup()%>%
        mutate(Regioncode=factor(Regioncode, levels=rs.order))%>%
        dplyr::select(Period, Regioncode,log2expr_mean)%>%
        mutate(
          Period_x=midpoints[as.character(Period)]
        ) 

      # Smooth
      df_loess <- exon11.df %>%
        group_by(Regioncode) %>%
        do({
          model <- loess(log2expr_mean ~ Period_x, data = ., span = 0.5)
          grid  <- seq(min(.$Period_x), max(.$Period_x), length.out = 1000)
          tibble(Period_x = grid,
                log2expr_smooth = predict(model, newdata = data.frame(Period_x = grid)),
                Regioncode = unique(.$Regioncode))
        }) %>%
        ungroup()
      df_loess.x <- unique(df_loess$Period_x)
      width.x <- c(df_loess.x[-1]-df_loess.x[-length(df_loess.x)], (df_loess.x[length(df_loess.x)]-df_loess.x[length(df_loess.x)-1]))
      
      WH.df <- data.frame(Period_x=df_loess.x,
                            width_x=width.x)
      
      df_loess <- df_loess%>%left_join(WH.df, by="Period_x")

      ## add a blank space between ITC and AMY
      # desired order
      lev <- rs.order
      # make continuous y positions; add a 0.5 gap AFTER "ITC"
      gap_after <- "ITC"
      gap_h <- 0.5
      pos <- seq_along(lev)
      i_gap <- match(gap_after, lev)
      pos[(i_gap):length(lev)] <- pos[(i_gap):length(lev)] + gap_h

      pos_map <- tibble(Regioncode = lev, y_pos = pos)

      df_pos <- df_loess %>%
        mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
        left_join(pos_map, by = "Regioncode")

      # y-axis breaks/labels and limits so nothing gets clipped
      y_breaks <- pos
      y_labels <- lev
      y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)
      
      # plotting
      kang11.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2expr_smooth)) +
        geom_tile(aes(width = width_x, height = 1)) +                  # normal row height
        scale_fill_gradientn(colors = hp.colors, na.value = "white") +
        scale_x_continuous(limits = c(intercepts[1], tail(intercepts,1)),
                          expand = c(0,0.1),
                          labels = midpoint_labels, breaks = midpoints) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                          expand = c(0,0)) +
        annotate(geom = "text", y = Inf, x = intercepts[1:(length(intercepts)-1)],
                label = age_label, size = 2.5, vjust = 1.5, hjust = .25) +
        geom_vline(xintercept = intercepts[-c(1, length(intercepts))],
                  linetype = "dashed", color = "gray25", linewidth = 0.25) +
        geom_vline(xintercept = 8.055282, linetype = "solid",
                  color = "gray25", linewidth = 0.25) +
        theme_bw() +
        labs(title = paste0("kang11: ", module, " (", length(module.gene), " genes)"),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")

      hp.plots[[length(hp.plots)+1]] <- kang11.hp
    }
  }
  all.plots <- hp.plots
  return(all.plots)
}

plot.GM.normal.heatmap.p3to9.kang11.nosmooth <- function(module.gene.list=NULL, rs.order=NULL, kang11.expr=NULL, kang11.meta=NULL){
  hp.plots <- list()
  for(module in names(module.gene.list)){
    module.gene <-  module.gene.list[[module]]%>%unique()
    if(length(module.gene)>0){
      # set plotting parameters
      if(startsWith(module, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }else{
        hp.colors <- c(colorRampPalette(c("white","#8e0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }
      intercepts <- log2(c(70, 91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
      labels <- 3:9
      intercepts.labels <- intercepts[-length(intercepts)]%>%setNames(labels)
      width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>%setNames(labels)
      age_label <- c("10PCW",13, 16, 19, "24PCW", "Birth", "0.5yr")
      midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
      midpoint_labels <- labels[1:length(midpoints)]
      
      rs.order.y <- 1:13 %>%setNames(rs.order)
      gs <- intersect(module.gene, exon11$Gene.Symbol)
      
      
      #kang11
      exon11.module <- exon11%>%dplyr::filter(grepl(paste0(gs,collapse = "|"), Gene.Symbol))%>%
        dplyr::select(grep("X[0-9]*",colnames(.),value = T))%>%
        colMeans(.)%>%
        data.frame()%>%setNames(module)%>%t()%>%as.data.frame()%>%
        pivot_longer(cols = all_of((colnames(.))), names_to = "Sample", values_to = "colMeans_log2expr")
      exon11.df <- kang11.meta%>%left_join(exon11.module,by="Sample",relationship = "many-to-many")%>%
        mutate(Period=as.numeric(Period))%>%
        dplyr::filter(Regioncode%in%rs.order, Period>2, Period<10)%>%
        group_by(Period, Regioncode)%>%
        mutate(log2expr_mean=median(colMeans_log2expr))%>%
        dplyr::select(Period, Regioncode,log2expr_mean)%>%
        distinct()%>%ungroup()%>%
        mutate(Regioncode=factor(Regioncode, levels=rs.order))%>%
        dplyr::select(Period, Regioncode,log2expr_mean)%>%
        mutate(
          Period_x=midpoints[as.character(Period)]
        ) 

      width_tbl <- tibble(
        Period = as.numeric(names(width.labels)),
        width_x = as.numeric(width.labels)
      )

      # join into exon11.df
      exon11.df <- exon11.df %>%
        left_join(width_tbl, by = "Period")

      ## add a blank space between ITC and AMY
      # desired order
      lev <- rs.order
      # make continuous y positions; add a 0.5 gap AFTER "ITC"
      gap_after <- "ITC"
      gap_h <- 0.5
      pos <- seq_along(lev)
      i_gap <- match(gap_after, lev)
      pos[(i_gap):length(lev)] <- pos[(i_gap):length(lev)] + gap_h

      pos_map <- tibble(Regioncode = lev, y_pos = pos)

      df_pos <- exon11.df %>%
        mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
        left_join(pos_map, by = "Regioncode")

      # y-axis breaks/labels and limits so nothing gets clipped
      y_breaks <- pos
      y_labels <- lev
      y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

      kang11.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2expr_mean)) +
        geom_tile(aes(width = width_x, height = 1)) +                  # normal row height
        scale_fill_gradientn(colors = hp.colors, na.value = "white") +
        scale_x_continuous(limits = c(intercepts[1], tail(intercepts,1)),
                          expand = c(0,0.1),
                          labels = midpoint_labels, breaks = midpoints) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                          expand = c(0,0)) +
        annotate(geom = "text", y = Inf, x = intercepts[1:(length(intercepts)-1)],
                label = age_label, size = 2.5, vjust = 1.5, hjust = .25) +
        geom_vline(xintercept = intercepts[-c(1, length(intercepts))],
                  linetype = "dashed", color = "gray25", linewidth = 0.25) +
        geom_vline(xintercept = 8.055282, linetype = "solid",
                  color = "gray25", linewidth = 0.25) +
        theme_bw() +
        labs(title = paste0("kang11: ", module, " (", length(module.gene), " genes)"),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")

    }else{
      kang11.hp <- NULL
    }
    hp.plots[[length(hp.plots)+1]] <- kang11.hp
  }
  return(hp.plots)
}

plot.GM.normal.heatmap.p3to9.li18 <- function(module.gene.list=NULL, rs.order=NULL, li18.expr=NULL, li18.meta=NULL){
  hp.plots <- list()
  for(module in names(module.gene.list)){
    print(module)
    module.gene <-  module.gene.list[[module]]%>%unique()
    if(length(module.gene)>0){
      # set plotting parameters
      if(startsWith(module, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }else{
        hp.colors <- c(colorRampPalette(c("white","#8e0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }
      intercepts <- log2(c(70, 91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
      labels <- 3:9
      intercepts.labels <- intercepts[-length(intercepts)]%>%setNames(labels)
      width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>%setNames(labels)
      age_label <- c("10PCW",13, 16, 19, "24PCW", "Birth", "0.5yr")
      midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
      midpoint_labels <- labels[1:length(midpoints)]
      
      rs.order.y <- 1:13 %>%setNames(rs.order)
      
      gs <- intersect(module.gene, unique(str_split_fixed(li18.expr$Geneid, "[|]",2)[,2]))
      
      
      li18.expr.module <- li18.expr%>%mutate(Gene.Symbol=str_split_fixed(Geneid, "[|]",2)[,2])%>%
        dplyr::filter(Gene.Symbol%in%gs)%>%
        dplyr::select(-c(Geneid,Gene.Symbol))%>%
        colMeans(.)%>%data.frame()%>%setNames(module)%>%t()%>%as.data.frame()%>%
        pivot_longer(cols=colnames(.), names_to = "Sample", values_to = "colMeans_fpkm")
      
      li18.expr.df <- li18.meta%>%left_join(li18.expr.module,by="Sample",relationship = "many-to-many")%>%
        dplyr::filter(Regioncode%in%rs.order, Period>2, Period<10)%>%
        group_by(Period, Regioncode)%>%
        mutate(fpkm_mean=median(colMeans_fpkm))%>%
        dplyr::select(Period, Regioncode,fpkm_mean)%>%
        distinct()%>%ungroup()%>%
        mutate(Regioncode=factor(Regioncode, levels=rs.order),
               log2fpkmPlus=log2(fpkm_mean+1))%>%
        dplyr::select(Period, Regioncode,log2fpkmPlus)
      mfc_p14_exp <- li18.expr.df%>%dplyr::filter(Period==13, Regioncode=="MFC")%>%
        mutate(Period=14)
      li18.expr.df <- rbind(li18.expr.df, mfc_p14_exp)%>%
        mutate(Period_x=midpoints[as.character(Period)]) 
      
      df_loess <- li18.expr.df %>%
        group_by(Regioncode) %>%
        do({
          model <- loess(log2fpkmPlus ~ Period_x, data = ., span = 0.5)
          grid  <- seq(min(.$Period_x), max(.$Period_x), length.out = 1000)
          tibble(Period_x = grid,
                log2expr_smooth = predict(model, newdata = data.frame(Period_x = grid)),
                Regioncode = unique(.$Regioncode))
        }) %>%
        ungroup()
      df_loess.x <- unique(df_loess$Period_x)
      width.x <- c(df_loess.x[-1]-df_loess.x[-length(df_loess.x)], (df_loess.x[length(df_loess.x)]-df_loess.x[length(df_loess.x)-1]))
      
      WH.df <- data.frame(Period_x=df_loess.x,
                            width_x=width.x)
      
      df_loess <- df_loess%>%left_join(WH.df, by="Period_x")

      ## add a blank space between ITC and AMY
      # desired order
      lev <- rs.order
      # make continuous y positions; add a 0.5 gap AFTER "ITC"
      gap_after <- "ITC"
      gap_h <- 0.5
      pos <- seq_along(lev)
      i_gap <- match(gap_after, lev)
      pos[(i_gap):length(lev)] <- pos[(i_gap):length(lev)] + gap_h

      pos_map <- tibble(Regioncode = lev, y_pos = pos)

      df_pos <- df_loess %>%
        mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
        left_join(pos_map, by = "Regioncode")

      # y-axis breaks/labels and limits so nothing gets clipped
      y_breaks <- pos
      y_labels <- lev
      y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

      li18.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2expr_smooth)) +
        geom_tile(aes(width = width_x, height = 1)) +                  # normal row height
        scale_fill_gradientn(colors = hp.colors, na.value = "white") +
        scale_x_continuous(limits = c(intercepts[1], tail(intercepts,1)),
                          expand = c(0,0.1),
                          labels = midpoint_labels, breaks = midpoints) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                          expand = c(0,0)) +
        annotate(geom = "text", y = Inf, x = intercepts[1:(length(intercepts)-1)],
                label = age_label, size = 2.5, vjust = 1.5, hjust = .25) +
        geom_vline(xintercept = intercepts[-c(1, length(intercepts))],
                  linetype = "dashed", color = "gray25", linewidth = 0.25) +
        geom_vline(xintercept = 8.055282, linetype = "solid",
                  color = "gray25", linewidth = 0.25) +
        theme_bw() +
        labs(title = paste0("li18: ", module, " (", length(module.gene), " genes)"),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
    }else{
      li18.hp <- NULL
    }
    hp.plots[[length(hp.plots)+1]] <- li18.hp
  }
  return(hp.plots)
}

plot.GM.normal.heatmap.p3to9.li18.nosmooth <- function(module.gene.list=NULL, rs.order=NULL, li18.expr=NULL, li18.meta=NULL){
  hp.plots <- list()
  for(module in names(module.gene.list)){
    print(module)
    module.gene <-  module.gene.list[[module]]%>%unique()
    if(length(module.gene)>0){
      # set plotting parameters
      if(startsWith(module, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }else{
        hp.colors <- c(colorRampPalette(c("white","#8e0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
      }
      intercepts <- log2(c(70, 91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
      labels <- 3:9
      intercepts.labels <- intercepts[-length(intercepts)]%>%setNames(labels)
      width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>%setNames(labels)
      age_label <- c("10PCW",13, 16, 19, "24PCW", "Birth", "0.5yr")
      midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
      midpoint_labels <- labels[1:length(midpoints)]
      
      rs.order.y <- 1:13 %>%setNames(rs.order)
      
      gs <- intersect(module.gene, unique(str_split_fixed(li18.expr$Geneid, "[|]",2)[,2]))
      
      
      li18.expr.module <- li18.expr%>%mutate(Gene.Symbol=str_split_fixed(Geneid, "[|]",2)[,2])%>%
        dplyr::filter(Gene.Symbol%in%gs)%>%
        dplyr::select(-c(Geneid,Gene.Symbol))%>%
        colMeans(.)%>%data.frame()%>%setNames(module)%>%t()%>%as.data.frame()%>%
        pivot_longer(cols=colnames(.), names_to = "Sample", values_to = "colMeans_fpkm")
      
      li18.expr.df <- li18.meta%>%left_join(li18.expr.module,by="Sample",relationship = "many-to-many")%>%
        dplyr::filter(Regioncode%in%rs.order, Period>2, Period<10)%>%
        group_by(Period, Regioncode)%>%
        mutate(fpkm_mean=median(colMeans_fpkm))%>%
        dplyr::select(Period, Regioncode,fpkm_mean)%>%
        distinct()%>%ungroup()%>%
        mutate(Regioncode=factor(Regioncode, levels=rs.order),
               log2fpkmPlus=log2(fpkm_mean+1))%>%
        dplyr::select(Period, Regioncode,log2fpkmPlus)
      mfc_p14_exp <- li18.expr.df%>%dplyr::filter(Period==13, Regioncode=="MFC")%>%
        mutate(Period=14)
      li18.expr.df <- rbind(li18.expr.df, mfc_p14_exp)%>%
        mutate(Period_x=midpoints[as.character(Period)]) 
      
      width_tbl <- tibble(
        Period = as.numeric(names(width.labels)),
        width_x = as.numeric(width.labels)
      )

      # join into exon11.df
      li18.expr.df <- li18.expr.df %>%
        left_join(width_tbl, by = "Period")

      ## add a blank space between ITC and AMY
      # desired order
      lev <- rs.order
      # make continuous y positions; add a 0.5 gap AFTER "ITC"
      gap_after <- "ITC"
      gap_h <- 0.5
      pos <- seq_along(lev)
      i_gap <- match(gap_after, lev)
      pos[(i_gap):length(lev)] <- pos[(i_gap):length(lev)] + gap_h

      pos_map <- tibble(Regioncode = lev, y_pos = pos)

      df_pos <- li18.expr.df %>%
        mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
        left_join(pos_map, by = "Regioncode")

      # y-axis breaks/labels and limits so nothing gets clipped
      y_breaks <- pos
      y_labels <- lev
      y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

      li18.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2fpkmPlus)) +
        geom_tile(aes(width = width_x, height = 1)) +                  # normal row height
        scale_fill_gradientn(colors = hp.colors, na.value = "white") +
        scale_x_continuous(limits = c(intercepts[1], tail(intercepts,1)),
                          expand = c(0,0.1),
                          labels = midpoint_labels, breaks = midpoints) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                          expand = c(0,0)) +
        annotate(geom = "text", y = Inf, x = intercepts[1:(length(intercepts)-1)],
                label = age_label, size = 2.5, vjust = 1.5, hjust = .25) +
        geom_vline(xintercept = intercepts[-c(1, length(intercepts))],
                  linetype = "dashed", color = "gray25", linewidth = 0.25) +
        geom_vline(xintercept = 8.055282, linetype = "solid",
                  color = "gray25", linewidth = 0.25) +
        theme_bw() +
        labs(title = paste0("li18: ", module, " (", length(module.gene), " genes)"),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
    }else{
      li18.hp <- NULL
    }
    hp.plots[[length(hp.plots)+1]] <- li18.hp
  }
  return(hp.plots)
}

plot.gene.normal.heatmap.p3to9.kang11 <- function(gene.list=NULL, rs.order=NULL, kang11.expr=NULL, kang11.meta=NULL){
  hp.plots <- list()
  for(gene in gene.list){
    print(gene)
    # set plotting parameters
    ## set color gradient
    if(gene%in%c("SEMA7A")){
      hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }else{
      hp.colors <- c(colorRampPalette(c("white","#8e0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }
    
    ## set axis labels
    intercepts <- log2(c(70, 91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
    labels <- 3:9
    intercepts.labels <- intercepts[-length(intercepts)]%>%setNames(labels)
    width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>%setNames(labels)
    age_label <- c("10PCW",13, 16, 19, "24PCW", "Birth", "0.5yr")
    midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
    midpoint_labels <- labels[1:length(midpoints)]
    
    rs.order.y <- 1:13 %>%setNames(rs.order)
    
    #kang11
    gene.ori <- gene
    if(gene=="ZBTB18"){
      gene="ZNF238"
    }
    exon11.gene <- exon11%>%dplyr::filter(grepl(gene, Gene.Symbol))%>%
      pivot_longer(cols=grep("X[0-9]*",colnames(.),value = T), names_to = "Sample", values_to = "log2expr")
    exon11.df <- kang11.meta%>%left_join(exon11.gene,by="Sample",relationship = "many-to-many")%>%
      mutate(Period=as.numeric(Period))%>%
      dplyr::filter(Regioncode%in%rs.order, Period>2, Period<10)%>%
      group_by(Period, Regioncode)%>%
      mutate(log2expr_mean=median(log2expr))%>%
      dplyr::select(Period, Regioncode,log2expr_mean)%>%
      distinct()%>%ungroup()%>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order))%>%
      dplyr::select(Period, Regioncode,log2expr_mean)%>%
      mutate(
        Period_x=midpoints[as.character(Period)]
      )
    df_loess <- exon11.df %>%
        group_by(Regioncode) %>%
        do({
          model <- loess(log2expr_mean ~ Period_x, data = ., span = 0.5)
          grid  <- seq(min(.$Period_x), max(.$Period_x), length.out = 1000)
          tibble(Period_x = grid,
                log2expr_smooth = predict(model, newdata = data.frame(Period_x = grid)),
                Regioncode = unique(.$Regioncode))
        }) %>%
        ungroup()
      df_loess.x <- unique(df_loess$Period_x)
      width.x <- c(df_loess.x[-1]-df_loess.x[-length(df_loess.x)], (df_loess.x[length(df_loess.x)]-df_loess.x[length(df_loess.x)-1]))
      
      WH.df <- data.frame(Period_x=df_loess.x,
                            width_x=width.x)
      
      df_loess <- df_loess%>%left_join(WH.df, by="Period_x")

     ## add a blank space between ITC and AMY
      # desired order
      lev <- rs.order
      # make continuous y positions; add a 0.5 gap AFTER "ITC"
      gap_after <- "ITC"
      gap_h <- 0.5
      pos <- seq_along(lev)
      i_gap <- match(gap_after, lev)
      pos[(i_gap):length(lev)] <- pos[(i_gap):length(lev)] + gap_h

      pos_map <- tibble(Regioncode = lev, y_pos = pos)

      df_pos <- df_loess %>%
        mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
        left_join(pos_map, by = "Regioncode")

      # y-axis breaks/labels and limits so nothing gets clipped
      y_breaks <- pos
      y_labels <- lev
      y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

      kang11.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2expr_smooth)) +
        geom_tile(aes(width = width_x, height = 1)) +                  # normal row height
        scale_fill_gradientn(colors = hp.colors, na.value = "white") +
        scale_x_continuous(limits = c(intercepts[1], tail(intercepts,1)),
                          expand = c(0,0.1),
                          labels = midpoint_labels, breaks = midpoints) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                          expand = c(0,0)) +
        annotate(geom = "text", y = Inf, x = intercepts[1:(length(intercepts)-1)],
                label = age_label, size = 2.5, vjust = 1.5, hjust = .25) +
        geom_vline(xintercept = intercepts[-c(1, length(intercepts))],
                  linetype = "dashed", color = "gray25", linewidth = 0.25) +
        geom_vline(xintercept = 8.055282, linetype = "solid",
                  color = "gray25", linewidth = 0.25) +
        theme_bw() +
        labs(title = paste0("kang11: ", paste0(unique(c(gene.ori,gene)), collapse = "/")),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
    hp.plots[[length(hp.plots)+1]] <- kang11.hp
    
  }
  return(hp.plots)
}

plot.gene.normal.heatmap.p3to9.kang11.nosmooth <- function(gene.list=NULL, rs.order=NULL, kang11.expr=NULL, kang11.meta=NULL){
  hp.plots <- list()
  for(gene in gene.list){
    print(gene)
    # set plotting parameters
    ## set color gradient
    if(gene%in%c("SEMA7A")){
      hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }else{
      hp.colors <- c(colorRampPalette(c("white","#8e0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }
    
    ## set axis labels
    intercepts <- log2(c(70, 91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
    labels <- 3:9
    intercepts.labels <- intercepts[-length(intercepts)]%>%setNames(labels)
    width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>%setNames(labels)
    age_label <- c("10PCW",13, 16, 19, "24PCW", "Birth", "0.5yr")
    midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
    midpoint_labels <- labels[1:length(midpoints)]
    
    rs.order.y <- 1:13 %>%setNames(rs.order)
    
    #kang11
    gene.ori <- gene
    if(gene=="ZBTB18"){
      gene="ZNF238"
    }
    exon11.gene <- exon11%>%dplyr::filter(grepl(gene, Gene.Symbol))%>%
      pivot_longer(cols=grep("X[0-9]*",colnames(.),value = T), names_to = "Sample", values_to = "log2expr")
    exon11.df <- kang11.meta%>%left_join(exon11.gene,by="Sample",relationship = "many-to-many")%>%
      mutate(Period=as.numeric(Period))%>%
      dplyr::filter(Regioncode%in%rs.order, Period>2, Period<10)%>%
      group_by(Period, Regioncode)%>%
      mutate(log2expr_mean=median(log2expr))%>%
      dplyr::select(Period, Regioncode,log2expr_mean)%>%
      distinct()%>%ungroup()%>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order))%>%
      dplyr::select(Period, Regioncode,log2expr_mean)%>%
      mutate(
        Period_x=midpoints[as.character(Period)]
      )
    width_tbl <- tibble(
        Period = as.numeric(names(width.labels)),
        width_x = as.numeric(width.labels)
      )

      # join into exon11.df
      exon11.df <- exon11.df %>%
        left_join(width_tbl, by = "Period")

     ## add a blank space between ITC and AMY
      # desired order
      lev <- rs.order
      # make continuous y positions; add a 0.5 gap AFTER "ITC"
      gap_after <- "ITC"
      gap_h <- 0.5
      pos <- seq_along(lev)
      i_gap <- match(gap_after, lev)
      pos[(i_gap):length(lev)] <- pos[(i_gap):length(lev)] + gap_h

      pos_map <- tibble(Regioncode = lev, y_pos = pos)

      df_pos <- exon11.df %>%
        mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
        left_join(pos_map, by = "Regioncode")

      # y-axis breaks/labels and limits so nothing gets clipped
      y_breaks <- pos
      y_labels <- lev
      y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

      kang11.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2expr_mean)) +
        geom_tile(aes(width = width_x, height = 1)) +                  # normal row height
        scale_fill_gradientn(colors = hp.colors, na.value = "white") +
        scale_x_continuous(limits = c(intercepts[1], tail(intercepts,1)),
                          expand = c(0,0.1),
                          labels = midpoint_labels, breaks = midpoints) +
        scale_y_continuous(breaks = y_breaks, labels = y_labels, limits = y_limits,
                          expand = c(0,0)) +
        annotate(geom = "text", y = Inf, x = intercepts[1:(length(intercepts)-1)],
                label = age_label, size = 2.5, vjust = 1.5, hjust = .25) +
        geom_vline(xintercept = intercepts[-c(1, length(intercepts))],
                  linetype = "dashed", color = "gray25", linewidth = 0.25) +
        geom_vline(xintercept = 8.055282, linetype = "solid",
                  color = "gray25", linewidth = 0.25) +
        theme_bw() +
        labs(title = paste0("kang11: ", paste0(unique(c(gene.ori,gene)), collapse = "/")),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
    hp.plots[[length(hp.plots)+1]] <- kang11.hp
    
  }
  return(hp.plots)
}

Dual.Color.FeaturePlot.GM <- function(seu, outDir='outs'){
  # default gm to plot
  gm.list <- list("S/Af"=c("S","Af"),
                  "S/At"=c("S","At"))
  
  gp.ls <- names(gm.list) %>% map(~{
    print(.x)
    gm.pair.nm <- .x
    gm.pair <- gm.list[[gm.pair.nm]]
    single.plot.tile <- paste0(gm.pair[1],"(green)/",gm.pair[2], "(magenta)")
    # Extract umap coordinates and module score
    umap_data <- as.data.frame(Embeddings(seu, reduction = "umap"))
    umap_data$GM1 <- as.numeric(unlist(FetchData(seu, vars = paste0("module_",gm.pair[1]))))
    umap_data$GM2 <- as.numeric(unlist(FetchData(seu, vars = paste0("module_",gm.pair[2]))))
    umap_data <- umap_data%>%
      mutate(GM1=ifelse(GM1<0,0,GM1),      # only focus on the module score > 0
             GM2=ifelse(GM2<0,0,GM2))%>%
      mutate(Colorcode=ifelse(GM1==0&GM2==0, "no", 
                              ifelse(GM1==0&GM2>0,"gm2", 
                                     ifelse(GM1>0&GM2==0,"gm1", "gm1gm2" ))))
    umap_data_expr <- umap_data%>%dplyr::filter(Colorcode!="no")%>%
      arrange(desc(GM1), GM2,Colorcode)  # sort GM score in two directions
    # set color gradient
    if("gm1"%in%unique(umap_data_expr$Colorcode)){  # set gm1 color gradient darker as expr increase
      color.gm1 <- colorRampPalette(c("#2D9B07", colorRampPalette(c("#2D9B07", "white"))(10)[7]))(table(umap_data$Colorcode)[["gm1"]])
    }else{
      color.gm1 <- NULL
      single.plot.tile <- paste0(single.plot.tile," No ", gm.pair[1])
    }
    if("gm2"%in%unique(umap_data_expr$Colorcode)){  # set gm2 color gradient darker as expr increase
      color.gm2 <- colorRampPalette(c(colorRampPalette(c("#8e0152", "white"))(10)[7],"#8e0152"))(table(umap_data$Colorcode)[["gm2"]])
    }else{
      color.gm2 <- NULL
      single.plot.tile <- paste0(single.plot.tile," No ", gm.pair[2])
    }
    if("gm1gm2"%in%unique(umap_data_expr$Colorcode)){ # set gm1 gm2 shared cell color gradient 
      color.gm1gm2 <- colorRampPalette(colorRampPalette(c(color.gm1[length(color.gm1)],"yellow", color.gm2[1]))(5)[2:4])(table(umap_data$Colorcode)[["gm1gm2"]])
    }else{
      color.gm1gm2 <- NULL
    }
    
    # Create dummy data for the color tile
    colorbar.file <- paste0(outDir,"/colorbar_GM.pdf")
    if(!file.exists(colorbar.file)){
      color_data <- data.frame(
        x = seq(0, 1, length.out = 100),  # Values along the gradient
        y = 1  # Single row for the color tile
      )
      # Plot the customized color tile
      pdf(colorbar.file,width = 1, height=5)
      color_data <- data.frame(x=1, y=1:nrow(umap_data_expr), Colororder=rev(1:nrow(umap_data_expr)))
      label_data <- data.frame(x=1, y=c(0,nrow(umap_data_expr)+1), label=rev(gm.pair))
      gb <- ggplot(color_data, aes(x = x, y = y, fill = Colororder)) +
        geom_tile() +
        geom_text(data = label_data, aes(x = x, y = y, label = label),
                  inherit.aes = FALSE, size = 3, color = "black") +
        scale_fill_gradientn(
          colors = c("#2D9B07", color.gm1[length(color.gm1)],"yellow", color.gm2[1], "#8e0152"),  # Custom gradient colors
          breaks = c(0.1, 1),                     # Breaks at min and max
          labels = gm.pair
        ) +
        theme_void() +  # Remove all unnecessary elements
        theme(
          legend.position = "none"
        ) +
        labs(title = "ColorBar")
      print(gb)
      dev.off()
    }
    
    # umap for each cell as a cot in the plot (similar as FeaturePlot)
    umap_data_expr$Color <- c(color.gm1,color.gm1gm2,color.gm2)
    umap_data2 <- umap_data%>%left_join(umap_data_expr)%>%
      mutate(Color=ifelse(is.na(Color), 'grey85', Color))
    # Plot cell umap with module score levels
    tmp <- ggplot(umap_data2, aes(x = umap_1, y = umap_2)) +
      geom_point(aes(color = Color), shape=16,size = 1, alpha = 1) +
      scale_color_identity() +
      theme_void() + 
      labs(title = single.plot.tile)
    return(tmp)
  })
  
  return(gp.ls)
}

Dual.Color.FeaturePlot.gene.pair <- function(seu, gene.pair=c("SEMA7A","PLXNC1"), outDir='outs'){
  gene.pair.nm <- paste(gene.pair, collapse = "/")
  # Extract umap coordinates and gene expression
  umap_data <- as.data.frame(Embeddings(seu, reduction = "umap"))%>%setNames(toupper(colnames(.)))
  umap_data$Gene1 <- as.numeric(unlist(FetchData(seu, vars = gene.pair[1])))
  umap_data$Gene2 <- as.numeric(unlist(FetchData(seu, vars = gene.pair[2])))
  umap_data <- umap_data%>%mutate(
    Gene1 = scales::rescale(Gene1, to = c(0, 1)),   # scaling the expression
    Gene2 = scales::rescale(Gene2, to = c(0, 1))
  )%>%
    mutate(Colorcode=ifelse(Gene1==0&Gene2==0, "no", 
                            ifelse(Gene1==0&Gene2>0,"g2", 
                                   ifelse(Gene1>0&Gene2==0,"g1", "g1g2" ))))
  umap_data_expr <- umap_data%>%dplyr::filter(Colorcode!="no")%>%
    arrange(desc(Gene1), Gene2,Colorcode)  # sort gene expression in two directions
  color.g1 <- colorRampPalette(c("#2D9B07", colorRampPalette(c("#2D9B07", "white"))(10)[7]))(table(umap_data$Colorcode)[["g1"]])
  color.g2 <- colorRampPalette(c(colorRampPalette(c("#8e0152", "white"))(10)[7],"#8e0152"))(table(umap_data$Colorcode)[["g2"]])
  color.g1g2 <- colorRampPalette(colorRampPalette(c(color.g1[length(color.g1)],"yellow", color.g2[1]))(5)[2:4])(table(umap_data$Colorcode)[["g1g2"]])
  
  
  # umap for each cell as a cot in the plot (similar as FeaturePlot)
  umap_data_expr$Color <- c(color.g1,color.g1g2,color.g2)
  umap_data2 <- umap_data%>%left_join(umap_data_expr)%>%
    mutate(Color=ifelse(is.na(Color), 'grey85', Color))
  # Plot cell umap with expression levels
  tmp <- ggplot(umap_data2, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = Color), shape=16,size = 1, alpha = 1) +
    scale_color_identity() +
    theme_void() + 
    labs(title = paste0(gene.pair[1],"(green)/",gene.pair[2], "(magenta)"))
  return(tmp)
}

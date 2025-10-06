###################################################################
# Time: 2025-08-15 15:33:55 EDT                                   #
# Author: Xinyun LI                                               #
# Email: xinyun.li@yale.edu                                       #
# Description: This script includes functions to plot normal      #
#              smoothened/unsmoothened/unimputed heatmaps of      # 
#              gene modules using macaque micro-array data        #
###################################################################


#############################  parameters   #######################################
# 1. zhu18.expr should contain a column: "Geneid", and other column names         # 
#       should be same as Sample in meta data                                     #
# 2. zhu18.meta should contain columns: "Sample", "Regioncode", "Period"# 
#############################  parameters   #######################################


###########  Function  ###################
library(tidyr)

plot.GM.normal.heatmap.p4to9.zhu18 <- function(module.gene.list=NULL, rs.order=NULL, zhu18.expr=NULL, zhu18.meta=NULL){
  hp.plots <- list()
  for(mod in names(module.gene.list)){
    module.gene <-  module.gene.list[[mod]] %>% unique()
    if (length(module.gene) ==0){
      next
    }
    # set plotting parameters
    ## set color gradient
    if(startsWith(mod, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }else{
        hp.colors <- c(colorRampPalette(c("white","#8E0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }

    ## set axis labels
    intercepts <- log2(c(91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
    labels <- 4:9
    intercepts.labels <- intercepts[-length(intercepts)] %>% setNames(labels)
    width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>% setNames(labels)
    age_label <- c(13, 16, 19, "24PCW", "Birth", "0.5yr")
    midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
    midpoint_labels <- labels[1:length(midpoints)]
      
    rs.order.y <- 1:13 %>% setNames(rs.order)
    
    #zhu18
    zhu18.meta$Sample <- paste(zhu18.meta$Brain, zhu18.meta$Regioncode, sep = ".")
    zhu18.expr.df <- zhu18.expr %>%
      dplyr::filter(grepl(paste0(module.gene,collapse = "|"), Geneid)) %>%
      dplyr::select(-Geneid) %>%
      colMeans(.) %>%
      data.frame() %>% setNames(mod) %>% t() %>% as.data.frame() %>%
      pivot_longer(cols=colnames(.), names_to = "Sample", values_to = "colMeans_fpkm") %>%
      dplyr::filter(Sample %in% zhu18.meta$Sample) %>% left_join(zhu18.meta,by="Sample") %>%
      group_by(Period, Regioncode) %>%
      mutate(fpkm_mean=median(colMeans_fpkm)) %>%
      dplyr::select(Period, Regioncode,fpkm_mean) %>%
      distinct() %>% ungroup() %>%
      dplyr::filter(Regioncode %in% rs.order, Period<10) %>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order),
             log2fpkmPlus=log2(fpkm_mean+1),
             Period=factor(as.character(Period), levels=as.character(c(4:9)))) %>%
      dplyr::select(Period, Regioncode,log2fpkmPlus)
    
    ## Fill in the missing values 
    ## using the mean of the expression in the n-1 and n+1 time point of the same region
    ## Missing values in the earliest time point will be filled with value at the second time point
    ## Missing values in the latest time point will be filled with value at the second to the latest time point
    zhu18_matrix <- zhu18.expr.df %>%
      pivot_wider(names_from = Period, values_from = log2fpkmPlus) %>%
      column_to_rownames("Regioncode")%>%
      dplyr::select(`4`, `5`, `6`, `7`, `8`, `9`)

    zhu18_matrix_filled <- t(apply(zhu18_matrix, 1, function(x) {
      for (i in 2:(length(x)-1)) {
        if (is.na(x[i])) {
          # Replace NA with the mean of the previous and next columns
          x[i] <- mean(c(x[i-1], x[i+1]), na.rm = TRUE)
        }
      }
      # Handle the first and last columns
      if (is.na(x[1])) {
        x[1] <- x[2] + rnorm(1, mean = 0, sd = 0.01)  # Replace the first column with the value of the second column
      }
      if (is.na(x[length(x)])) {
        x[length(x)] <- x[length(x)-1] + rnorm(1, mean = 0, sd = 0.01)  # Replace the last column with the value of the second-last column
      }
      return(x)
    }))

    # Convert to data frame
    df_imputed <- zhu18_matrix_filled %>% 
      as.data.frame() %>%
      rownames_to_column("Regioncode") %>%
      pivot_longer(cols = -Regioncode, names_to = "Period", values_to = "log2fpkmPlus") %>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order),
             Period=factor(as.character(Period), levels=as.character(c(4:9))),
             Period_x = midpoints[as.character(Period)])
      
    #-------------- Smoothing ----------
    df_loess <- df_imputed %>%
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
    
    # plotting
    zhu18.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2expr_smooth)) +
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
        labs(title = paste0("zhu18: ", mod, " (", length(module.gene), " genes)"),
            x = "Period", y = NULL) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.ticks = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              legend.position = "none")
  
    hp.plots <- append(hp.plots, list(zhu18.hp))
  }
  return(hp.plots)
}

plot.GM.normal.heatmap.p4to9.zhu18.nosmooth <- function(module.gene.list=NULL, rs.order=NULL, zhu18.expr=NULL, zhu18.meta=NULL){
  hp.plots <- list()
  for(mod in names(module.gene.list)){
    module.gene <-  module.gene.list[[mod]] %>% unique()
    if (length(module.gene) ==0){
      next
    }
    # set plotting parameters
    ## set color gradient
    if(startsWith(mod, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }else{
        hp.colors <- c(colorRampPalette(c("white","#8E0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }

    ## set axis labels
    intercepts <- log2(c(91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
    labels <- 4:9
    intercepts.labels <- intercepts[-length(intercepts)] %>% setNames(labels)
    width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>% setNames(labels)
    age_label <- c(13, 16, 19, "24PCW", "Birth", "0.5yr")
    midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
    midpoint_labels <- labels[1:length(midpoints)]
      
    rs.order.y <- 1:13 %>% setNames(rs.order)
    
    #zhu18
    zhu18.meta$Sample <- paste(zhu18.meta$Brain, zhu18.meta$Regioncode, sep = ".")
    zhu18.expr.df <- zhu18.expr %>%
      dplyr::filter(grepl(paste0(module.gene,collapse = "|"), Geneid)) %>%
      dplyr::select(-Geneid) %>%
      colMeans(.) %>%
      data.frame() %>% setNames(mod) %>% t() %>% as.data.frame() %>%
      pivot_longer(cols=colnames(.), names_to = "Sample", values_to = "colMeans_fpkm") %>%
      dplyr::filter(Sample %in% zhu18.meta$Sample) %>% left_join(zhu18.meta,by="Sample") %>%
      group_by(Period, Regioncode) %>%
      mutate(fpkm_mean=median(colMeans_fpkm)) %>%
      dplyr::select(Period, Regioncode,fpkm_mean) %>%
      distinct() %>% ungroup() %>%
      dplyr::filter(Regioncode %in% rs.order, Period<10) %>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order),
             log2fpkmPlus=log2(fpkm_mean+1),
             Period=factor(as.character(Period), levels=as.character(c(4:9)))) %>%
      dplyr::select(Period, Regioncode,log2fpkmPlus)
    
    ## Fill in the missing values 
    ## using the mean of the expression in the n-1 and n+1 time point of the same region
    ## Missing values in the earliest time point will be filled with value at the second time point
    ## Missing values in the latest time point will be filled with value at the second to the latest time point
    zhu18_matrix <- zhu18.expr.df %>%
      pivot_wider(names_from = Period, values_from = log2fpkmPlus) %>%
      column_to_rownames("Regioncode") %>%
      dplyr::select(`4`, `5`, `6`, `7`, `8`, `9`)

    zhu18_matrix_filled <- t(apply(zhu18_matrix, 1, function(x) {
      for (i in 2:(length(x)-1)) {
        if (is.na(x[i])) {
          # Replace NA with the mean of the previous and next columns
          x[i] <- mean(c(x[i-1], x[i+1]), na.rm = TRUE)
        }
      }
      # Handle the first and last columns
      if (is.na(x[1])) {
        x[1] <- x[2] #+ rnorm(1, mean = 0, sd = 0.01)  # Replace the first column with the value of the second column
      }
      if (is.na(x[length(x)])) {
        x[length(x)] <- x[length(x)-1] #+ rnorm(1, mean = 0, sd = 0.01)  # Replace the last column with the value of the second-last column
      }
      return(x)
    }))

    # Convert to data frame
    df_imputed <- zhu18_matrix_filled %>% 
      as.data.frame() %>%
      rownames_to_column("Regioncode") %>%
      pivot_longer(cols = -Regioncode, names_to = "Period", values_to = "log2fpkmPlus") %>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order),
             #Period=factor(as.character(Period), levels=as.character(c(4:9))),
             Period = as.numeric(Period),
             Period_x = midpoints[as.character(Period)])
      
    width_tbl <- tibble(
        Period = as.numeric(names(width.labels)),
        width_x = as.numeric(width.labels)
      )

    # join into exon11.df
    df_imputed <- df_imputed %>%
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

    df_pos <- df_imputed %>%
      mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
      left_join(pos_map, by = "Regioncode")

    # y-axis breaks/labels and limits so nothing gets clipped
    y_breaks <- pos
    y_labels <- lev
    y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

    zhu18.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2fpkmPlus)) +
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
      labs(title = paste0("zhu18: ", mod, " (", length(module.gene), " genes)"),
          x = "Period", y = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = "none")
  
    hp.plots <- append(hp.plots, list(zhu18.hp))
  }
  return(hp.plots)
}

plot.GM.normal.heatmap.p4to9.zhu18.nosmooth.noimpute <- function(module.gene.list=NULL, rs.order=NULL, zhu18.expr=NULL, zhu18.meta=NULL){
  hp.plots <- list()
  for(mod in names(module.gene.list)){
    module.gene <-  module.gene.list[[mod]] %>% unique()
    if (length(module.gene) ==0){
      next
    }
    # set plotting parameters
    ## set color gradient
    if(startsWith(mod, "S")){
        hp.colors <- c(colorRampPalette(c("white","#2D9B07"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }else{
        hp.colors <- c(colorRampPalette(c("white","#8E0152"))(50)[c(1:6 ,10,15,20,30,40,50)])
    }

    ## set axis labels
    intercepts <- log2(c(91, 112, 19*7, 24*7, 38*7,(38*7+ c(6*30, 365))))
    labels <- 4:9
    intercepts.labels <- intercepts[-length(intercepts)] %>% setNames(labels)
    width.labels <- (intercepts[-1] - intercepts[-length(intercepts)]) %>% setNames(labels)
    age_label <- c(13, 16, 19, "24PCW", "Birth", "0.5yr")
    midpoints <- ((intercepts[-1] + intercepts[-length(intercepts)]) / 2 )%>% setNames(labels)
    midpoint_labels <- labels[1:length(midpoints)]
      
    rs.order.y <- 1:13 %>% setNames(rs.order)
    
    #zhu18
    zhu18.meta$Sample <- paste(zhu18.meta$Brain, zhu18.meta$Regioncode, sep = ".")
    zhu18.expr.df <- zhu18.expr %>%
      dplyr::filter(grepl(paste0(module.gene,collapse = "|"), Geneid)) %>%
      dplyr::select(-Geneid) %>%
      colMeans(.) %>%
      data.frame() %>% setNames(mod) %>% t() %>% as.data.frame() %>%
      pivot_longer(cols=colnames(.), names_to = "Sample", values_to = "colMeans_fpkm") %>%
      dplyr::filter(Sample %in% zhu18.meta$Sample) %>% left_join(zhu18.meta,by="Sample") %>%
      group_by(Period, Regioncode) %>%
      mutate(fpkm_mean=median(colMeans_fpkm)) %>%
      dplyr::select(Period, Regioncode,fpkm_mean) %>%
      distinct() %>% ungroup() %>%
      dplyr::filter(Regioncode %in% rs.order, Period<10) %>%
      mutate(Regioncode=factor(Regioncode, levels=rs.order),
             log2fpkmPlus=log2(fpkm_mean+1),
             Period=factor(as.character(Period), levels=as.character(c(4:9)))) %>%
      dplyr::select(Period, Regioncode,log2fpkmPlus) %>%
      mutate(Period_x = midpoints[as.character(Period)])
    
      
    width_tbl <- tibble(
        Period = as.factor(names(width.labels)),
        width_x = as.numeric(width.labels)
      )

    # join into exon11.df
    zhu18.expr.df <- zhu18.expr.df %>%
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

    df_pos <- zhu18.expr.df %>%
      mutate(Regioncode = factor(Regioncode, levels = lev)) %>%
      left_join(pos_map, by = "Regioncode")

    # y-axis breaks/labels and limits so nothing gets clipped
    y_breaks <- pos
    y_labels <- lev
    y_limits <- c(min(pos) - 0.5, max(pos) + 0.5)

    zhu18.hp <- ggplot(df_pos, aes(x = Period_x, y = y_pos, fill = log2fpkmPlus)) +
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
      labs(title = paste0("zhu18: ", mod, " (", length(module.gene), " genes)"),
          x = "Period", y = NULL) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            legend.position = "none")
  
    hp.plots <- append(hp.plots, list(zhu18.hp))
  }
  return(hp.plots)
}
















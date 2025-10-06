###################################################################
# Time: 2025-08-17 15:43:25 EDT                                   #
# Author: Xinyun LI                                               #
# Email: xinyun.li@yale.edu                                       #
# Description: This script is used to plot the Main figure 4c.    #
###################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(ggrepel)
library(ggpattern)
library(dplyr)
library(patchwork)
library(scales)
library(openxlsx)

setwd("./outs")

# Read in the conserved GM genes; can be obtained by running 01_get_GM_list.R
GM <- read_excel("shared.prenatal.GM.xlsx") %>%
    as.data.frame()
# Initiate the GM table and its parameters for each GM
Af_gene <- GM[GM$Type == "Af","GeneSymbol"]
At_gene <- GM[GM$Type == "At","GeneSymbol"]
AfAt_shared <- intersect(Af_gene, At_gene)

GM_list <- list(
    "Af GM" = GM %>% subset(Type == "Af"),
    "At GM" = GM %>% subset(Type == "At"),
    "S GM" = GM %>% subset(Type == "S")
)

regions <- c("mPFC", "OFC", "MOs")
colors <- c("#8E0152", "#C75A98", "#2D9B07")

#######################################
# Read in RA dKO DEG result tables [RNA-seq data obtained from Shibata, M. et al., Nature 2021]
dKO_mc_table <- read.csv("DEG_tables/MOs_RAdKO-DEGs-raw.csv") %>%
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    subset(padj < 0.05)
dKO_ofc_table <- read.csv("DEG_tables/OFC_RAdKO-DEGs-raw.csv") %>%
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    subset(padj < 0.05)
dKO_pfc_table <- read.csv("DEG_tables/mPFC_RAdKO-DEGs-raw.csv") %>%
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    subset(padj < 0.05)
dKO_table_list <- list(dKO_pfc_table, dKO_ofc_table, dKO_mc_table)
names(dKO_table_list) <- regions

# Merge the GM table with the dKO tables
final_merged_up_table <- NULL
final_merged_down_table <- NULL
for (gm_name in names(GM_list)){
    # for each of the gene modules..
    gm_table <- GM_list[[gm_name]]
    print("gm_table nrow:")
    print(nrow(gm_table))
    color_name <- colors[which(names(GM_list)==gm_name)]
    merged_up_perGM_table <- NULL
    merged_down_perGM_table <- NULL
    
    for (region in regions){
        print(paste("region: ", region))
        # for each of the dKO tables of different brain region...
        dKO_table <- dKO_table_list[[region]]
        dKO_up_table <- dKO_table[dKO_table$log2FoldChange>0,] %>% na.omit()
        dKO_down_table <- dKO_table[dKO_table$log2FoldChange<0,] %>% na.omit()
        
        # merge GM and dKO table, keeping all genes
        cols_to_replace <- c("log2FoldChange")

        gm_dKO_up <- full_join(dKO_up_table, gm_table, by = c("SYMBOL" = "GeneSymbol")) %>%
            mutate(across(all_of(cols_to_replace), ~ replace_na(.x, 0)), # assign 0 in log2FoldChange to GM genes not differentially expressed
            dot_color = if_else(is.na(Type)|log2FoldChange==0,"grey", color_name)) %>% # grey dots representing non-DE GM genes and non-GM DEGs
            mutate(log2FoldChange = ifelse(log2FoldChange > 3, 3, log2FoldChange)) # cap log2FC to 3

        gm_dKO_down <- full_join(dKO_down_table, gm_table, by = c("SYMBOL" = "GeneSymbol")) %>%
            mutate(across(all_of(cols_to_replace), ~ replace_na(.x, 0)), # assign 0 in log2FoldChange to GM genes not differentially expressed
            dot_color = if_else(is.na(Type)|log2FoldChange==0,"grey", color_name)) %>% # grey dots representing non-DE GM genes and non-GM DEGs
            mutate(log2FoldChange = ifelse(log2FoldChange < -3, -3, log2FoldChange)) # cap log2FC to -3

        if (gm_name == "S GM"){
            # Merged dKO DEGs and S GM genes all have solid dot fill
            gm_dKO_up <- gm_dKO_up  %>% 
                mutate(dot_fill = "solid")

            gm_dKO_down <- gm_dKO_down %>%
                mutate(dot_fill = "solid")

            # calculate overlapping num of GM genes and significant DEGs in dKO
            up_overlap <- nrow(gm_dKO_up[!is.na(gm_dKO_up$Type) & gm_dKO_up$log2FoldChange>0,])
            down_overlap <- nrow(gm_dKO_down[!is.na(gm_dKO_down$Type) & gm_dKO_down$log2FoldChange<0,])

        }else{
            # Merged dKO DEGs and Af & At GM genes have solid dot fill when they are unique to either Af or At GM
            # If the genes are present both in Af and At, their dot fill will be hollow
            gm_dKO_up <- gm_dKO_up  %>%
                mutate(dot_fill = if_else(dot_color == color_name & SYMBOL %in% AfAt_shared, "empty", "solid"))

            gm_dKO_down <- gm_dKO_down %>%
                mutate(dot_fill = if_else(dot_color == color_name & SYMBOL %in% AfAt_shared, "empty", "solid"))

            # calculate overlapping num of GM genes and significant DEGs in dKO
            up_overlap <- nrow(gm_dKO_up[gm_dKO_up$dot_fill=="solid" & !is.na(gm_dKO_up$Type) & gm_dKO_up$log2FoldChange>0,])
            down_overlap <- nrow(gm_dKO_down[gm_dKO_down$dot_fill=="solid" & !is.na(gm_dKO_down$Type) & gm_dKO_down$log2FoldChange<0,])
        }
        # calculate the percentage of overlapping DEG and GM genes over the total number of GMs / DEGs
        up_perc_gm <- up_overlap/nrow(gm_dKO_up[!is.na(gm_dKO_up$Type),])
        down_perc_gm <- down_overlap/nrow(gm_dKO_down[!is.na(gm_dKO_down$Type),])

        up_perc_dKO <- up_overlap/nrow(gm_dKO_up[gm_dKO_up$log2FoldChange>0,])
        down_perc_dKO <- down_overlap/nrow(gm_dKO_down[gm_dKO_down$log2FoldChange<0,])

        # highlight Plxnc1 in blue
        if ("PLXNC1" %in% gm_dKO_up$SYMBOL){
            if (!is.na(gm_dKO_up[gm_dKO_up$SYMBOL=="PLXNC1","Type"]) & gm_dKO_up[gm_dKO_up$SYMBOL=="PLXNC1","log2FoldChange"]>0){
                gm_dKO_up[gm_dKO_up$SYMBOL=="PLXNC1", "dot_color"] <- "blue"
            }
        }
        if ("PLXNC1" %in% gm_dKO_down$SYMBOL){
            if (!is.na(gm_dKO_down[gm_dKO_down$SYMBOL=="PLXNC1","Type"]) & gm_dKO_down[gm_dKO_down$SYMBOL=="PLXNC1","log2FoldChange"]<0){
                gm_dKO_down[gm_dKO_down$SYMBOL=="PLXNC1", "dot_color"] <- "blue"
            }
        }
        
        
        # add columns containing (overlapping) numbers of genes, perc of genes, region, GM name, and total number of GM gene/DEG num
        gm_dKO_up$up_overlap <- up_overlap
        gm_dKO_up$up_perc_gm <- up_perc_gm
        gm_dKO_up$up_perc_dKO <- up_perc_dKO

        gm_dKO_down$down_overlap <- down_overlap
        gm_dKO_down$down_perc_gm <- down_perc_gm
        gm_dKO_down$down_perc_dKO <- down_perc_dKO

        gm_dKO_up$region <- region
        gm_dKO_down$region <- region
        gm_dKO_up$GM <- gm_name
        gm_dKO_down$GM <- gm_name

        gm_dKO_up$GM_num <- nrow(gm_dKO_up[!is.na(gm_dKO_up$Type),])
        gm_dKO_down$GM_num <- nrow(gm_dKO_down[!is.na(gm_dKO_down$Type),])
        gm_dKO_up$dKO_num <- nrow(gm_dKO_up[gm_dKO_up$log2FoldChange>0,])
        gm_dKO_down$dKO_num <- nrow(gm_dKO_down[gm_dKO_down$log2FoldChange<0,])

        # bind the merged table for each regional dKO by row
        if(is.null(merged_up_perGM_table)){
            merged_up_perGM_table <- gm_dKO_up
            merged_down_perGM_table <- gm_dKO_down
        }else{
            merged_up_perGM_table <- rbind(merged_up_perGM_table, gm_dKO_up)
            merged_down_perGM_table <- rbind(merged_down_perGM_table, gm_dKO_down)
        }
        print("merged up table nrow:")
        print(nrow(merged_up_perGM_table))
        print("merged down table nrow:")
        print(nrow(merged_down_perGM_table))
    }

    final_up_perGM_table <- merged_up_perGM_table %>% mutate(color = color_name)
    final_down_perGM_table <- merged_down_perGM_table %>% mutate(color = color_name)

    # bind the merged table for each GM by row
    if(is.null(final_merged_up_table)){
        final_merged_up_table <- final_up_perGM_table
        final_merged_down_table <- final_down_perGM_table
    }else{
        final_merged_up_table <- rbind(final_merged_up_table, final_up_perGM_table)
        final_merged_down_table <- rbind(final_merged_down_table, final_down_perGM_table)
    }

}

########### Additional plotting parameters #############
# Define the three groups manually
region_order <- c("mPFC-Af GM", "OFC-Af GM", "MOs-Af GM", 
                  "mPFC-At GM", "OFC-At GM", "MOs-At GM", 
                  "mPFC-S GM", "OFC-S GM", "MOs-S GM")

# Assign x positions to make groups visually clustered
x_positions <- c(1, 2, 3,    # First group (Af)
                 5, 6, 7,    # Second group (At)
                 9, 10, 11)  # Third group (S)
x_map <- data.frame(region_GM = factor(region_order, levels = region_order), x = x_positions)

# Preparing table columns for plotting
table_modifier <- function(table, direction){
    table$region_GM <- factor(paste0(table$region,"-",table$GM), 
        levels = unique(paste0(table$region,"-",table$GM)))
    table$color <- factor(table$color, levels = colors)

    if (direction == "up"){
        table <- table %>% subset(log2FoldChange>0)
    }else{
        table <- table %>% subset(log2FoldChange<0)
    }
    table <- table %>%
    mutate(dot_alpha = if_else(dot_color=="grey",0.2,0.7)) %>%
    left_join(x_map, by = "region_GM")

    return(table)
}

final_merged_up_table <- table_modifier(final_merged_up_table, "up")

final_merged_down_table <- table_modifier(final_merged_down_table, "down")

















########## Significance test ##############
# whether the proportions of overlapping genes in two groups (e.g., mPFC-Af GM vs. mPFC-S GM) are significantly different.
# functions getting the fishers test values and performing fishers test
test_nums <- function(table, direction, region_GM, total_num){
    overlap <- unique(table[table$region_GM == region_GM, paste0(direction,"_overlap")])
    total <- unique(table[table$region_GM == region_GM, total_num])
    return(list(overlap = overlap, total = total))
}
fishers_test <- function(overlap1, total1, overlap2, total2){
    contingency_table <- matrix(c(overlap1, total1 - overlap1, 
                                    overlap2, total2 - overlap2), nrow = 2, byrow = TRUE)
    res <- fisher.test(contingency_table)
    return(res$p.value)
}
######### GM_num - total number of GM genes in each GM as the denominator ############
test_num_list_up <- list()
test_num_list_down <- list()
# loop through the region_GMs and calculate the entries to the fishers test
for(region in region_order){
    overlap <- test_nums(final_merged_up_table, "up", region, "GM_num")$overlap
    total <- test_nums(final_merged_up_table, "up", region, "GM_num")$total
    test_num_list_up[[region]] <- c(overlap, total)
    
    overlap <- test_nums(final_merged_down_table, "down", region, "GM_num")$overlap
    total <- test_nums(final_merged_down_table, "down", region, "GM_num")$total
    test_num_list_down[[region]] <- c(overlap, total)
}
# perform fishers test
sig_pvalues_gm_num <- list()
# comparison within GM
for (gm in c("Af", "At", "S")){
    print(gm)
    region_gms <- names(test_num_list_up)[grep(gm, names(test_num_list_up))]
    sublist_up <- test_num_list_up[names(test_num_list_up) %in% region_gms]
    sublist_down <- test_num_list_down[names(test_num_list_up) %in% region_gms]
    for(i in 2:length(region_gms)){
        print("Up:")
        res <- fishers_test(sublist_up[[1]][1], sublist_up[[1]][2],
                            sublist_up[[i]][1], sublist_up[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_gm_num[[paste0(region_gms[1], "_", region_gms[i],"_up")]] <- res
        }
        print("Down:")
        res <- fishers_test(sublist_down[[1]][1], sublist_down[[1]][2],
                            sublist_down[[i]][1], sublist_down[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_gm_num[[paste0(region_gms[1], "_", region_gms[i],"_down")]] <- res
        }
    }
}

# comparison between GM
for (region in c("mPFC", "OFC", "MOs")){
    print(region)
    region_gms <- names(test_num_list_up)[grep(region, names(test_num_list_up))]
    sublist_up <- test_num_list_up[names(test_num_list_up) %in% region_gms]
    sublist_down <- test_num_list_down[names(test_num_list_up) %in% region_gms]
    for(i in 2:length(region_gms)){
        print("Up:")
        res <- fishers_test(sublist_up[[1]][1], sublist_up[[1]][2],
                            sublist_up[[i]][1], sublist_up[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_gm_num[[paste0(region_gms[1], "_", region_gms[i],"_up")]] <- res
        }
        print("Down:")
        res <- fishers_test(sublist_down[[1]][1], sublist_down[[1]][2],
                            sublist_down[[i]][1], sublist_down[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_gm_num[[paste0(region_gms[1], "_", region_gms[i],"_down")]] <- res
        }
    }
}

######### dKO_num - total number of DEGs in each dKO as the denominator ############
test_num_list_up <- list()
test_num_list_down <- list()
# loop through the region_GMs and calculate the entries to the fishers test
for(region in region_order){
    overlap <- test_nums(final_merged_up_table, "up", region, "dKO_num")$overlap
    total <- test_nums(final_merged_up_table, "up", region, "dKO_num")$total
    test_num_list_up[[region]] <- c(overlap, total)
    
    overlap <- test_nums(final_merged_down_table, "down", region, "dKO_num")$overlap
    total <- test_nums(final_merged_down_table, "down", region, "dKO_num")$total
    test_num_list_down[[region]] <- c(overlap, total)
}
# perform fishers test
sig_pvalues_dKO_num <- list()
# comparison within GM
for (gm in c("Af", "At", "S")){
    print(gm)
    region_gms <- names(test_num_list_up)[grep(gm, names(test_num_list_up))]
    sublist_up <- test_num_list_up[names(test_num_list_up) %in% region_gms]
    sublist_down <- test_num_list_down[names(test_num_list_up) %in% region_gms]
    for(i in 2:length(region_gms)){
        print("Up:")
        res <- fishers_test(sublist_up[[1]][1], sublist_up[[1]][2],
                            sublist_up[[i]][1], sublist_up[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_dKO_num[[paste0(region_gms[1], "_", region_gms[i],"_up")]] <- res
        }
        print("Down:")
        res <- fishers_test(sublist_down[[1]][1], sublist_down[[1]][2],
                            sublist_down[[i]][1], sublist_down[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_dKO_num[[paste0(region_gms[1], "_", region_gms[i],"_down")]] <- res
        }
    }
}

# comparison between GM
for (region in c("mPFC", "OFC", "MOs")){
    print(region)
    region_gms <- names(test_num_list_up)[grep(region, names(test_num_list_up))]
    sublist_up <- test_num_list_up[names(test_num_list_up) %in% region_gms]
    sublist_down <- test_num_list_down[names(test_num_list_up) %in% region_gms]
    for(i in 2:length(region_gms)){
        print("Up:")
        res <- fishers_test(sublist_up[[1]][1], sublist_up[[1]][2],
                            sublist_up[[i]][1], sublist_up[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_dKO_num[[paste0(region_gms[1], "_", region_gms[i],"_up")]] <- res
        }
        print("Down:")
        res <- fishers_test(sublist_down[[1]][1], sublist_down[[1]][2],
                            sublist_down[[i]][1], sublist_down[[i]][2])
        print(paste(region_gms[1], "vs", region_gms[i], ":", res))
        if(res < 0.05){
            sig_pvalues_dKO_num[[paste0(region_gms[1], "_", region_gms[i],"_down")]] <- res
        }
    }
}

# format the p-values for plotting
for (i in 1:length(sig_pvalues_dKO_num)){
    sig_pvalues_dKO_num[[i]] <- formatC(sig_pvalues_dKO_num[[i]], format = "e", digits = 2)
}
for (i in 1:length(sig_pvalues_gm_num)){
    sig_pvalues_gm_num[[i]] <- formatC(sig_pvalues_gm_num[[i]], format = "e", digits = 2)
}






############### Visualization ##############
# plotting parameters
axis_margin <- 2
x_label <- rep(regions,3)
max_log2FC <- max(final_merged_up_table$log2FoldChange, abs(final_merged_down_table$log2FoldChange), na.rm = TRUE)

# Plotting the up-regulated overlapping genes
# df1 for the bar plot showing the percentage of overlapping genes over GM genes/DEGs
df1 <- final_merged_up_table %>% 
    dplyr::select(region_GM, up_perc_gm, up_perc_dKO, color, x) %>% 
    distinct() %>%
    pivot_longer(
        cols = c(up_perc_gm, up_perc_dKO),
        names_to = "source",       # Indicates whether it's from GM or dKO
        values_to = "up_perc"      # Combines both into one column
    ) %>%
    mutate(x_adjusted = if_else(row_number() %% 2 == 1, x - 0.25, x + 0.25)) %>% # bar width
    mutate(bar_fill = if_else(source == "up_perc_gm", "solid", "hollow")) %>%
    mutate(fill_final = if_else(bar_fill == "solid", color, NA_character_))
# df2 for the dot plot showing the log2FC of the overlapping GM genes
df2 <- final_merged_up_table %>% 
    dplyr::select(region_GM, SYMBOL, log2FoldChange, dot_color, dot_alpha, dot_fill, x) %>% 
    distinct() %>%
    mutate(scaled_log2FC = log2FoldChange / (max_log2FC / 0.4))  # Scale log2FC to match percentage range
df2$dot_fill <- factor(df2$dot_fill, levels = c("solid", "empty"))

# combine the bar plot and dot plot together
p1 <- ggplot() +
    geom_col(data = df1, aes(y = up_perc, x = x_adjusted, fill = fill_final, color = color), alpha=0.5, width = 0.5) + # bar plot
    geom_point(data=df2, aes(x = x, y = scaled_log2FC, color = dot_color, alpha = dot_alpha, shape = dot_fill), 
                size = 2, position = position_jitter(width = 0.2, height = 0)) + # dot plot
    scale_shape_manual(values = c("solid" = 16, "empty" = 1), name = "Unique in Af or At", labels = c("Yes", "No"), drop = FALSE) + # dot shape
    scale_x_continuous(
        breaks = x_positions,  # Define exact positions for x-ticks
        labels = x_label  # Assign correct labels
    ) +
    scale_y_continuous(
        name = "Percentage of GM genes upregulated in dKO",
        limits = c(0,0.4), # Primary y-axis range
        sec.axis = sec_axis(~ . * (max_log2FC / 0.4), name = "Log2 Fold Change")  # Secondary y-axis for log2FC
    ) +
    theme_minimal() +
    theme(
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.margin = margin(t = axis_margin, r= axis_margin, b= 0, l=axis_margin),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.5, "cm")) + 
    # legend parameters
    guides(alpha = "none") +
    scale_fill_identity(
        name = "Gene Module",
        labels = c("Af only", "At only", "S union"),
        guide = "legend",
        breaks = c("#8E0152", "#C75A98", "#2D9B07")
        ) + 
    scale_color_identity() 

# Plotting the down-regulated overlapping genes
# df1 for the bar plot showing the percentage of overlapping genes over GM genes/DEGs
df1 <- final_merged_down_table %>% 
    dplyr::select(region_GM, down_perc_gm, down_perc_dKO, color, x) %>% 
    distinct() %>%
    pivot_longer(
        cols = c(down_perc_gm, down_perc_dKO),
        names_to = "source",       # Indicates whether it's from GM or dKO
        values_to = "down_perc"      # Combines both into one column
    ) %>%
    mutate(x_adjusted = if_else(row_number() %% 2 == 1, x - 0.25, x + 0.25)) %>% # bar width
    mutate(bar_fill = if_else(source == "down_perc_gm", "solid", "hollow")) %>%
    mutate(fill_final = if_else(bar_fill == "solid", color, NA_character_))
# df2 for the dot plot showing the log2FC of the overlapping GM genes
df2 <- final_merged_down_table %>% 
    dplyr::select(region_GM, SYMBOL, log2FoldChange, dot_color, dot_alpha, dot_fill, x) %>% 
    distinct() %>%
    mutate(scaled_log2FC = log2FoldChange / (max_log2FC / 0.4))  # Scale log2FC to match percentage range

p2 <- ggplot() +
    geom_col(data = df1, aes(y = down_perc, x = x_adjusted, fill = fill_final, color = color), alpha=0.5, width = 0.5) + # bar plot
    geom_point(data = df2, aes(x = x, y = -scaled_log2FC, color = dot_color, alpha = dot_alpha, shape = dot_fill), 
                size = 2, show.legend = FALSE, position = position_jitter(width = 0.2, height = 0)) + # dot plot
    scale_shape_manual(values = c("solid" = 16, "empty" = 1)) +
    scale_x_continuous(
        breaks = x_positions  # Define exact positions for x-ticks
    ) +
    scale_y_continuous(
        name = "Percentage of GM genes downregulated in dKO",
        limits = c(0.4,0), # Primary y-axis range
        trans = "reverse", # reverse the y-axis to flip the plot upside down
        sec.axis = sec_axis(~ -(. * (max_log2FC / 0.4)), name = "Log2 Fold Change")  # Secondary y-axis for log2FC
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t=0, r=axis_margin, b=axis_margin, l=axis_margin))+
    labs(x = "Regions", fill = "Gene module") +
    scale_fill_identity() +
    scale_color_identity() +
    # adding text for significant p-values in the fishers test
    # Bracket for mPFC-Af GM vs. mPFC-At GM
    geom_segment(aes(x = 1.25, xend = 5.25, y = 0.3, yend = 0.3), size = 0.5, color = "grey40") +  # Horizontal line for bracket
    geom_segment(aes(x = 1.25, xend = 1.25, y = 0.28, yend = 0.3), size = 0.5, color = "grey40") +  # Left vertical line
    geom_segment(aes(x = 5.25, xend = 5.25, y = 0.28, yend = 0.3), size = 0.5, color = "grey40") +  # Right vertical line
    annotate("text", x = 3.25, y = 0.31, label = paste("p =", sig_pvalues_dKO_num[[1]]), size = 2, color = "grey40") +  # P-value text
    # Bracket for mPFC-Af GM vs. mPFC-S GM
    geom_segment(aes(x = 1.25, xend = 9.25, y = 0.34, yend = 0.34), size = 0.5, color = "grey40") +
    geom_segment(aes(x = 1.25, xend = 1.25, y = 0.32, yend = 0.34), size = 0.5, color = "grey40") +
    geom_segment(aes(x = 9.25, xend = 9.25, y = 0.32, yend = 0.34), size = 0.5, color = "grey40") +
    annotate("text", x = 5.25, y = 0.35, label = paste("p =", sig_pvalues_dKO_num[[2]]), size = 2, color = "grey40") +
    # Bracket for mPFC-Af GM vs. OFC-Af GM
    geom_segment(aes(x = 0.75, xend = 1.75, y = 0.22, yend = 0.22), size = 0.5) +
    geom_segment(aes(x = 0.75, xend = 0.75, y = 0.2, yend = 0.22), size = 0.5) +
    geom_segment(aes(x = 1.75, xend = 1.75, y = 0.2, yend = 0.22), size = 0.5) +
    annotate("text", x = 1.25, y = 0.23, label = paste("p =", sig_pvalues_gm_num[[1]]), size = 2) +
    # Bracket for mPFC-Af GM vs. MOs-Af GM
    geom_segment(aes(x = 0.75, xend = 2.75, y = 0.26, yend = 0.26), size = 0.5) +
    geom_segment(aes(x = 0.75, xend = 0.75, y = 0.24, yend = 0.26), size = 0.5) +
    geom_segment(aes(x = 2.75, xend = 2.75, y = 0.24, yend = 0.26), size = 0.5) +
    annotate("text", x = 1.75, y = 0.27, label = paste("p =", sig_pvalues_gm_num[[2]]), size = 2)

pdf("RA-dKO-bar-dot-plot.pdf", width = 5.5, height = 6)
print(p1/p2)
dev.off()






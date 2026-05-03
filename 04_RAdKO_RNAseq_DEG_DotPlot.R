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
Af_gene <- GM[GM$Type == "Peri-F","GeneSymbol"]
At_gene <- GM[GM$Type == "Peri-T","GeneSymbol"]
AfAt_shared <- intersect(Af_gene, At_gene)

GM_list <- list(
    "Af GM" = GM %>% subset(Type == "Peri-F"),
    "At GM" = GM %>% subset(Type == "Peri-T"),
    "S GM" = GM %>% subset(Type == "Cent-union S")
)

regions <- c("mPFC", "OFC", "MOs")
colors <- c("#8E0152", "#C75A98", "#2D9B07")

#######################################
# Read in RA dKO DEG result tables [RNA-seq data obtained from Shibata, M. et al., Nature 2021]
dKO_mc_table <- read.csv("DEG_tables/MOs_RAdKO-DEGs-raw.csv") 
total_gene_num_mc <- nrow(dKO_mc_table)
dKO_mc_table %>% dKO_mc_table
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    subset(padj < 0.05, abs(log2FoldChange) >= 0.6)

dKO_ofc_table <- read.csv("DEG_tables/OFC_RAdKO-DEGs-raw.csv") 
total_gene_num_ofc <- nrow(dKO_ofc_table)
dKO_ofc_table %>% dKO_ofc_table
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    subset(padj < 0.05, abs(log2FoldChange) >= 0.6)

dKO_pfc_table <- read.csv("DEG_tables/mPFC_RAdKO-DEGs-raw.csv") 
total_gene_num_pfc <- nrow(dKO_pfc_table)
dKO_pfc_table %>% dKO_pfc_table
    mutate(SYMBOL = toupper(SYMBOL)) %>%
    subset(padj < 0.05, abs(log2FoldChange) >= 0.6)
dKO_table_list <- list(dKO_pfc_table, dKO_ofc_table, dKO_mc_table)
names(dKO_table_list) <- regions

total_gene_num_list <- list(total_gene_num_pfc, total_gene_num_ofc, total_gene_num_mc)
names(total_gene_num_list) <- regions

# Merge the GM table with the dKO tables
final_merged_up_table <- NULL
final_merged_down_table <- NULL
for (gm_name in names(GM_list)){
    # for each of the gene modules..
    gm_table <- GM_list[[gm_name]]
    print("gm_table nrow:")
    print(nrow(gm_table))
    color_name <- colors[which(names(GM_list)==gm_name)]
    merged_up_table <- NULL
    merged_down_table <- NULL
    
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
            mutate(log2FoldChange = ifelse(log2FoldChange > 3, 3, log2FoldChange)) %>% # cap log2FC to 3
            mutate(dot_fill = "solid")

        gm_dKO_down <- full_join(dKO_down_table, gm_table, by = c("SYMBOL" = "GeneSymbol")) %>%
            mutate(across(all_of(cols_to_replace), ~ replace_na(.x, 0)), # assign 0 in log2FoldChange to GM genes not differentially expressed
            dot_color = if_else(is.na(Type)|log2FoldChange==0,"grey", color_name)) %>% # grey dots representing non-DE GM genes and non-GM DEGs
            mutate(log2FoldChange = ifelse(log2FoldChange < -3, -3, log2FoldChange)) %>% # cap log2FC to -3
            mutate(dot_fill = "solid")

        # calculate overlapping num of GM genes and significant DEGs in dKO
        up_overlap <- nrow(gm_dKO_up[!is.na(gm_dKO_up$Type) & gm_dKO_up$log2FoldChange>0,])
        down_overlap <- nrow(gm_dKO_down[!is.na(gm_dKO_down$Type) & gm_dKO_down$log2FoldChange<0,])
        
        
        # add columns containing (overlapping) numbers of genes, perc of genes, region, GM name, and total number of GM gene/DEG num
        gm_dKO_up$up_overlap <- up_overlap
        gm_dKO_down$down_overlap <- down_overlap

        gm_dKO_up$region <- region
        gm_dKO_down$region <- region
        gm_dKO_up$GM <- gm_name
        gm_dKO_down$GM <- gm_name

        gm_dKO_up$GM_num <- nrow(gm_dKO_up[!is.na(gm_dKO_up$Type),])
        gm_dKO_down$GM_num <- nrow(gm_dKO_down[!is.na(gm_dKO_down$Type),])
        gm_dKO_up$dKO_num <- nrow(gm_dKO_up[gm_dKO_up$log2FoldChange>0,])
        gm_dKO_down$dKO_num <- nrow(gm_dKO_down[gm_dKO_down$log2FoldChange<0,])

        gm_dKO_up$total_gene_num <- total_gene_num_list[[region]]
        gm_dKO_down$total_gene_num <- total_gene_num_list[[region]]

        # bind the merged table for each regional dKO by row
        if(is.null(merged_up_table)){
            merged_up_table <- gm_dKO_up
            merged_down_table <- gm_dKO_down
        }else{
            merged_up_table <- rbind(merged_up_table, gm_dKO_up)
            merged_down_table <- rbind(merged_down_table, gm_dKO_down)
        }
        print("merged up table nrow:")
        print(nrow(merged_up_table))
        print("merged down table nrow:")
        print(nrow(merged_down_table))
    }

    final_up_table <- merged_up_table %>% mutate(color = color_name)
    final_down_table <- merged_down_table %>% mutate(color = color_name)

    # bind the merged table for each GM by row
    if(is.null(final_merged_up_table)){
        final_merged_up_table <- final_up_table
        final_merged_down_table <- final_down_table
    }else{
        final_merged_up_table <- rbind(final_merged_up_table, final_up_table)
        final_merged_down_table <- rbind(final_merged_down_table, final_down_table)
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

final_up_table_DEG <- final_merged_up_table %>%
                        dplyr::select(ensembl_gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, SYMBOL, region, dot_color) %>%
                        arrange(padj) %>%
                        rename(in_which_GM = dot_color) %>%
                        distinct() %>%
                        mutate(DEG_direction = "up")
final_up_table_DEG$in_which_GM[final_up_table_DEG$in_which_GM == "grey"] <- "None"
final_up_table_DEG$in_which_GM[final_up_table_DEG$in_which_GM == "#8E0152"] <- "Af"
final_up_table_DEG$in_which_GM[final_up_table_DEG$in_which_GM == "#C75A98"] <- "At"
final_up_table_DEG$in_which_GM[final_up_table_DEG$in_which_GM == "#2D9B07"] <- "S"

final_down_table_DEG <- final_merged_down_table %>%
                        dplyr::select(ensembl_gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, SYMBOL, region, dot_color) %>%
                        arrange(padj)%>%
                        rename(in_which_GM = dot_color)%>%
                        distinct() %>%
                        mutate(DEG_direction = "down")
final_down_table_DEG$in_which_GM[final_down_table_DEG$in_which_GM == "grey"] <- "None"
final_down_table_DEG$in_which_GM[final_down_table_DEG$in_which_GM == "#8E0152"] <- "Af"
final_down_table_DEG$in_which_GM[final_down_table_DEG$in_which_GM == "#C75A98"] <- "At"
final_down_table_DEG$in_which_GM[final_down_table_DEG$in_which_GM == "#2D9B07"] <- "S"

final_merged_table <- rbind(final_up_table_DEG, final_down_table_DEG)
write.csv(final_merged_table, "MF3c_dKO_DEG_table.csv")


###### hypergeometric test ########
up_results <- final_merged_up_table %>%
  select(region, GM, up_overlap, GM_num, dKO_num, total_gene_num) %>%
  distinct() %>%
  mutate(
    p_val = phyper(
      q = up_overlap - 1, 
      m = GM_num,
      n = total_gene_num - GM_num,
      k = dKO_num,
      lower.tail = FALSE
    ),
    padj = p.adjust(p_val, method = "BH"),
    a = up_overlap,
    c = GM_num - up_overlap,
    b = dKO_num - up_overlap,
    d = (total_gene_num - GM_num) - b,
    odds_ratio = (a * d) / (b * c)
  ) %>%
  dplyr::select(-c(a,b,c,d))

down_results <- final_merged_down_table %>%
  select(region, GM, down_overlap, GM_num, dKO_num, total_gene_num) %>%
  distinct() %>%
  mutate(
    p_val = phyper(
      q = down_overlap - 1, 
      m = GM_num,
      n = total_gene_num - GM_num,
      k = dKO_num,
      lower.tail = FALSE
    ),
    padj = p.adjust(p_val, method = "BH"),
    a = down_overlap,
    c = GM_num - down_overlap,
    b = dKO_num - down_overlap,
    d = (total_gene_num - GM_num) - b,
    odds_ratio = (a * d) / (b * c)
  )%>%
  dplyr::select(-c(a,b,c,d))

up_res_modified <- up_results %>% 
            rename(overlap_num = up_overlap, dKO_DEG_num = dKO_num, total_bg_gene_num = total_gene_num) %>%
            mutate(DEG_direction = "up")
down_res_modified <- down_results %>% 
            rename(overlap_num = down_overlap, dKO_DEG_num = dKO_num, total_bg_gene_num = total_gene_num) %>%
            mutate(DEG_direction = "down")
full_results <- rbind(up_res_modified, down_res_modified)
write.csv(full_results,"MF3c_GM_enrichment_test_statistics.csv")



############### Visualization ##############
# plotting parameters
axis_margin <- 2
x_label <- rep(regions,3)
format_padj <- formatC(min(down_results$padj), format = "e", digits = 2)


df <- final_merged_up_table %>% 
    dplyr::select(region_GM, SYMBOL, log2FoldChange, dot_color, dot_alpha, dot_fill, x) %>% 
    distinct()


p1 <- ggplot() +
    geom_point(data=df, aes(x = x, y = log2FoldChange, color = dot_color, alpha = dot_alpha), 
                size = 2, position = position_jitter(width = 0.2, height = 0)) + 
    scale_x_continuous(
        breaks = x_positions,  # Define exact positions for x-ticks
        labels = x_label  # Assign correct labels
    ) +
    scale_y_continuous(
        name = "Log2FC of upregulated genes in dKO",
        limits = c(0,3)
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
    scale_color_identity(
        name = "Gene Module",
        labels = c("Af only", "At only", "S union"),
        guide = "none",
        breaks = c("#8E0152", "#C75A98", "#2D9B07")
        )



df <- final_merged_down_table %>% 
    dplyr::select(region_GM, SYMBOL, log2FoldChange, dot_color, dot_alpha, dot_fill, x) %>% 
    distinct() 

p2 <- ggplot() +
    geom_point(data = df, aes(x = x, y = -log2FoldChange, color = dot_color, alpha = dot_alpha), 
                size = 2, show.legend = TRUE, position = position_jitter(width = 0.2, height = 0)) +
    scale_x_continuous(
        breaks = x_positions  # Define exact positions for x-ticks
    ) +
    scale_y_continuous(
        name = "Log2FC of downregulated genes in dKO",
        limits = c(3,0), # Primary y-axis range
        trans = "reverse", # reverse the y-axis to flip the plot upside down
        breaks = seq(0, 3, by = 1),       # Sets the tick positions
        labels = function(x) paste0("-", x)
    ) +
    theme_minimal() +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        plot.margin = margin(t=0, r=axis_margin, b=axis_margin, l=axis_margin))+
    labs(x = "Regions") +
    scale_color_identity(
        name = "Gene Module",
        labels = c("Af", "At", "S union"),
        guide = "legend",
        breaks = c("#8E0152", "#C75A98", "#2D9B07")
        )+
    guides(alpha = "none") +
    annotate("text", x = 1.5, y = 2.5, label = paste("padj =", format_padj), size = 2, color = "grey40")

pdf("MF3c-RA-dKO-dot-plot-hypergeometric-padj.pdf", width = 5.5, height = 6)
print(p1/p2)
dev.off()
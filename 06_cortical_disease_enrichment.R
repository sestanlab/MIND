###################################################################
# Time: 2025-08-17 15:36:50 EDT                                   #
# Author: Xinyun LI                                               #
# Email: xinyun.li@yale.edu                                       #
# Description: This script is used to plot the Extended figure 5. #
###################################################################

library(tidyverse)
library(readxl)
library(purrr)
library(ggplot2)
library(ggtext) 
library(qvalue)
library(scales)

# function to modify the dot size in the grid plot
power_trans <- function(power = 2) {
    trans_new(
        name = paste0("power", power),
        transform = function(x) x^power,
        inverse = function(x) x^(1/power)
    )
}
# function to extract background gene list for hypergeometric test
kang11 <- function(){
    exon11 <- read.table('PATH_TO_Kang2011_exon_microarray_data/nature105213.GENE1346_nosnp_Log2-expression value.txt',header = T)%>%
        mutate(Geneid=paste(ID, Gene.Symbol, sep="|"))  
    meta11 <- readRDS("data/meta.kang2011.rds")%>%
        mutate(Sample=rownames(.))
    prenatal_samples <- unique(meta11[meta11$Period %in% 3:7, "Sample"])
    indi <- prenatal_samples[which(prenatal_samples %in% colnames(exon11))]
    exon11_prenatal <- exon11[, indi]
    all_genes <- sub("^[^|]*\\|", "", exon11$Geneid)
    print(paste("# of all_genes:", length(all_genes)))

    return(all_genes)
}
li18 <- function(){
    fpkm_li18 <- read.table('PATH_TO_Li2018_RNA_seq_data/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.count.txt',header = T)
    meta18 <- readRDS("data/meta.li2018.rds")
    prenatal_samples <- unique(meta18[meta18$Period %in% 3:7, "Braincode"])
    fpkm_li18_prenatal <- fpkm_li18[, grep(paste(prenatal_samples, collapse = "|"), colnames(fpkm_li18))]
    all_genes <- sub("^[^|]*\\|", "", fpkm_li18$Geneid)
    print(paste("# of all_genes:", length(all_genes)))

    return(all_genes)
}
zhu18 <- function(){
    fpkm_zhu18 <- read.table('PATH_TO_Zhu2018_RNA_seq_data/nhp_development_count_rmTechRep.txt',header = T)
    fpkm_zhu18$Geneid <- rownames(fpkm_zhu18)
    regions <- c("HIP",'OFC','MFC','DFC','VFC','M1C','S1C','IPC','V1C','A1C','STC', 'ITC','AMY')
    meta_zhu18 <- readRDS("data/meta.zhu2018.rds") %>%
        dplyr::filter(Regioncode%in%regions)
    prenatal_samples <- unique(meta_zhu18[meta_zhu18$Period %in% 5:7, "Brain"])
    fpkm_zhu18_prenatal <- fpkm_zhu18[, grep(paste(prenatal_samples, collapse = "|"), colnames(fpkm_zhu18))]
    non_zero_rows <- which(apply(fpkm_zhu18_prenatal, 1, function(x) any(x != 0)))
    nonzero_genes <- rownames(fpkm_zhu18_prenatal[non_zero_rows,])
    all_genes <- sub("^[^|]*\\|", "", nonzero_genes)
    print(paste("# of all_genes:", length(all_genes)))

    return(all_genes)
}

# ------------ collect genes associated with diseases ------------
outdir <- "./outs"
## Read in supplementary table 1 from Mato-Blanco, X., Kim, SK., Jourdon, A. et al. Nat Common 2025
file_path <- paste0(outdir,"/supp_table.1.Gene_disease_associations.xlsx")
df <- read_excel(file_path, col_names = FALSE)
colnames(df) <- df[5, ]
colnames(df)[1] <- "Gene"
df <- df[-c(1:5), ]
disease_gene_list <- list()
# Iterate through each disease (column) and extract genes (row indices) associated with the disease (whose entry value is 1)
for (col in colnames(df)[2:ncol(df)]) {
  disease_gene_list[[col]] <- df[which(df[[col]] == 1),]$Gene
}

# Extract/merge cortical diseases of interest
merged_disease_gene_list <- list(
  "AD" = disease_gene_list[["AD_2019"]],
  "ADHD" = disease_gene_list[["ADHD_2019"]],
  "AN" = disease_gene_list[["AN_2019"]],
  "ASD_HC(S1)" = disease_gene_list[["SFAR_S1"]],
  "ASD_SFAR123" = unique(c(disease_gene_list[["SFAR_S1"]], disease_gene_list[["SFAR_S2"]], disease_gene_list[["SFAR_S3"]])),
  "BD" = disease_gene_list[["BD_2019"]],
  "DD" = disease_gene_list[["DD"]],
  "IQ" = disease_gene_list[["IQ_2018"]],
  "MDD" = disease_gene_list[["MDD_2018"]],
  "NEUROT" = disease_gene_list[["NEUROT_2018"]],
  "PD" = disease_gene_list[["PD_2014"]],
  "SCZ" = disease_gene_list[["SCZ_2020"]]
)
# Read in non-brain (control) diseases
non_brain <- read_excel(paste0(outdir, "/non-brain-disease-gene-lists.xlsx"))
merged_disease_gene_list[['CAD']] <- unique(na.omit(non_brain$`CAD-st22-3methods`))
merged_disease_gene_list[['Cohen']] <- unique(na.omit(non_brain$Cohen))
merged_disease_gene_list[['Lupus']] <- unique(na.omit(non_brain$lupus))
merged_disease_gene_list[['MetS']] <- unique(gsub("`","",na.omit(non_brain$MetS)))
merged_disease_gene_list[['T2D']] <- unique(na.omit(non_brain$T2D))

# print out the number of genes associated with each disease
for (i in 1:length(merged_disease_gene_list)) {
  element_length <- length(merged_disease_gene_list[[i]])
  print(paste("Length of", names(merged_disease_gene_list)[i], "is", element_length))
}

# merge all disease genes into a vector
all_disease_genes <- unique(unlist(merged_disease_gene_list))
length(unique(unlist(merged_disease_gene_list[1:12]))) # 2502 cortical disease genes
length(unique(unlist(merged_disease_gene_list[13:17]))) # 697 non-brain disease genes

dataset <- c("Human.Kang2011", "Human.Li2018", "Macaque.Zhu2018", "shared")
for (ds in dataset){
    # ------ Read in prenatal GM genes identified in each study; can be obtained via running 01_get_GM_list.R ---------
    file_path <- paste0(outdir, "/", ds, ".prenatal.GM.xlsx")
    GM_df <- read_excel(file_path)
    GM_list <- GM_df %>%
        group_by(Type) %>%
        summarise(values = list(GeneSymbol), .groups = "drop") %>%
        deframe()
    GM_list <- GM_list[c("Af", "At", "Af u At", "Af n At", "Ss1", "Sm1", "Sv1", "Sa1", "S")]
    
    # ----- Read in the background genes (genes expressed in prenatal periods; p3-p7 for human and p4-p7 for macaque) for the enrichment test -------
    if(ds == "Human.Kang2011"){
        all_genes <- kang11()
    }
    if(ds == "Human.Li2018"){
        all_genes <- li18()
    }
    if(ds == "Macaque.Zhu2018"){
        all_genes <- zhu18()
    }
    if(ds == "shared"){
        kang_genes <- kang11()
        li_genes <- li18()
        zhu_genes <- zhu18()
        all_genes <- unique(intersect(intersect(li_genes, kang_genes),zhu_genes))
    }

    # ------- select disease genes that are present in the background gene set (all_genes) ------
    sum(all_disease_genes %in% all_genes) < length(all_disease_genes) # True: not all disease genes are present in the background gene set

    get_common_genes <- function(merged_disease_gene_list) {
        for (i in names(merged_disease_gene_list)) {
            current_genes <- merged_disease_gene_list[[i]]
            shared_genes <- intersect(current_genes, all_genes)
            merged_disease_gene_list[[i]] <- shared_genes
        }
        return(merged_disease_gene_list)
    }
    shared_disease_list <- get_common_genes(merged_disease_gene_list)
    print("# of cortical disease genes in BG gene set: ")
    print(length(unique(unlist(shared_disease_list[1:12]))))
    print("# of non-brain disease genes in BG gene set: ")
    print(length(unique(unlist(shared_disease_list[13:17]))))

    # ----- select genes present in the background gene set but not in any GM ----
    GM_list[["Remaining Genes"]] <- unique(setdiff(all_genes, unique(unlist(GM_list))))

    # ---- perform hypergeometric test to calculate the enrichment of disease genes in each GM gene set ----
    p_values <- numeric(length(GM_list))
    gene_list <- shared_disease_list
    # enrichment test function
    enrichment_results <- map_df(names(gene_list), function(dg_name) {
        disease_genes <- gene_list[[dg_name]]
        
        map_df(names(GM_list), function(gm_name) {
            gm_genes <- GM_list[[gm_name]]
            
            # Define sizes
            m <- length(disease_genes)  # Size of the disease gene set
            n <- length(all_genes) - m  # Size of the background gene set
            k <- length(gm_genes)  # Size of the GM gene set of interest

            # Create contingency table
            overlap <- length(intersect(disease_genes, gm_genes))
            only_disease <- length(setdiff(disease_genes, gm_genes))
            only_gm <- length(setdiff(gm_genes, disease_genes))
            neither <- length(setdiff(all_genes, union(disease_genes, gm_genes)))
            OR <- (overlap*neither)/(only_disease*only_gm)
        
            # Hypergeometric test
            p_values[i] <- phyper(overlap - 1, m, n, k, lower.tail = FALSE)

            # overlapped GM and disease gene names
            overlapping_genes <- intersect(gm_genes, disease_genes)
            
            tibble(
            Disease = dg_name,
            GM = gm_name,
            pvalue = p_values[i],
            OddsRatio = OR,
            Count = overlap,
            GeneRatio = paste0(overlap,"/", length(gm_genes)),
            BgRatio = paste0(length(disease_genes), "/", length(all_genes)), 
            OverlappingGenes = paste(overlapping_genes, collapse = ", ")
            )
        })
    })

    # adjust hypergeometric test p-value using BH
    enrichment_results <- enrichment_results %>%
    mutate(p.adj = p.adjust(pvalue, method = "BH"),
            qvalue = qvalue(pvalue)$qvalues)
    enrichment_results <- enrichment_results[, c("Disease","GM","GeneRatio","BgRatio","OddsRatio","pvalue","p.adj","qvalue","Count","OverlappingGenes")]
    # save the enrichment test results
    write.csv(enrichment_results, paste0(outdir, "/cortical-disease-enrichment-table-",ds,".csv"), row.names = FALSE)


    # ------------- visualization ------------
    # Define the desired order of GMs
    custom_order <- rev(c(
    "Af", "At", "Af u At", "Af n At",
    "Ss1", "Sm1", "Sv1", "Sa1", "S", "Remaining Genes"
    ))
    module_rank <- match(enrichment_results$GM, custom_order)

    # format the GM label for plotting
    enrichment_results$GM_label <- ifelse(
    grepl("^A", enrichment_results$GM),
    paste0("<span style='color:#8e0152;'>", enrichment_results$GM, "</span>"),
    ifelse(
        grepl("^S", enrichment_results$GM),
        paste0("<span style='color:#2D9B07;'>", enrichment_results$GM, "</span>"),
        enrichment_results$GM
    )
    )
    enrichment_results_sorted <- enrichment_results[order(module_rank), ]
    enrichment_results_sorted$GM_label <- factor(enrichment_results_sorted$GM_label, levels=unique(enrichment_results_sorted$GM_label))
    enrichment_results_sorted$Disease <- factor(enrichment_results_sorted$Disease, levels = unique(names(merged_disease_gene_list)))

    # make non-significant OR 0 and p.adj 1
    enrichment_results_sorted_sig_adj <- enrichment_results_sorted %>%
    mutate(
        OddsRatio = ifelse(p.adj > 0.05, 0, OddsRatio),
        p.adj = ifelse(p.adj > 0.05, 1, p.adj)
    )

    # Plotting
    pdf(paste0(outdir, "/cortical-disease-dotplot-",ds,".pdf"), width = 17, height = 5)
    p <- ggplot(enrichment_results_sorted_sig_adj, aes(x = Disease, y = GM_label)) +
        geom_point(
            aes(
            size = OddsRatio,
            color = -log10(p.adj),
            alpha = OddsRatio > 0
            )
        ) +
        scale_size_continuous(
            range = c(4, 12),
            trans = power_trans(2),  # power transformation here to adjust the dot size
            breaks = c(0, 2, 4, 6, 8, 10, 12), # legend ticks
            labels = c(0, 2, 4, 6, 8, 10, 12), # legend labels
            name = "Odds Ratio"
        ) +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0)) +
        scale_color_gradient(low = "white", high = "red") +
        theme_minimal() +
        theme(
            axis.text.y = element_markdown()
        ) +
        labs(
            title = "Gene Set Enrichment",
            x = "Disease",
            y = "Gene module",
            color = "-log10(p-adj)",
            alpha = NULL
        )
    print(p)
    dev.off()

}

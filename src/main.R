suppressPackageStartupMessages(library("DESeq2"));
suppressPackageStartupMessages(library("tidyverse"));
suppressPackageStartupMessages(library("EnhancedVolcano"));
suppressPackageStartupMessages(library("pheatmap"));
suppressPackageStartupMessages(library("clusterProfiler"));
suppressPackageStartupMessages(library("org.Hs.eg.db"));
suppressPackageStartupMessages(library("enrichplot"));
suppressPackageStartupMessages(library("plotly"));
suppressPackageStartupMessages(library("ReactomePA"));
suppressPackageStartupMessages(library("rWikiPathways"));
suppressPackageStartupMessages(library("DT"));
suppressPackageStartupMessages(library("pathview"));

# set reactome database for pathview


# Set the warning behavior to "suppress"
options(warn = -1)

display_table <- function(table, caption) {
    suppressMessages(datatable(table,
                    caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:200% ;',
                    caption)))
}

#' Reads the count dataset and metadata, and preprocesses them.
#'
#' This function reads the count dataset from the specified file path
#' and stores it as a matrix. It removes the length column from the count matrix.
#' It reads the metadata from the specified file path and converts the "condition"
#' column to a factor. The metadata is sorted by sample, and the row names of
#' the metadata are set to be the order of the samples. It checks if all the samples
#' in the metadata are present in the count matrix. The function returns a list
#' containing the count matrix and the processed metadata.
#'
#' @param count_path File path to the count dataset.
#' @param metadata_path File path to the metadata.
#'
#' @return A list containing the count matrix and processed metadata.
#'
#' @examples
#' count_path <- "data/count.csv"
#' metadata_path <- "data/metadata.txt"
#' result <- read_and_preprocess(count_path, metadata_path)
#' count_matrix <- result$count_matrix
#' coldata <- result$coldata
#'
read_and_preprocess <- function(count_path, metadata_path) {

    # get count dataset
    count_matrix <- as.matrix(read.csv(count_path, row.names = "gene_id"))

    coldata <- read.table(metadata_path, header=TRUE)
    coldata$condition <- as.factor(coldata$condition)

    # sort by sample
    coldata <- coldata[order(coldata$sample), ]

    # set rowname to be the order samples
    rownames(coldata) <- coldata[,1]
    coldata[,1] <- NULL

    # check to be sure all samples in metadata are in the count matrix
    stopifnot(all(rownames(coldata) %in% colnames(count_matrix)))
    stopifnot(all(rownames(coldata) == colnames(count_matrix)))

    display_table(count_matrix, caption = "Subset Selection")
    display_table(coldata, caption = "test")

    list(count_matrix = count_matrix, coldata = coldata)
}

#' Performs differential expression analysis using DESeq2.
#'
#' This function takes a count matrix and metadata as input, performs
#' differential expression analysis using DESeq2, and returns the DESeq2 object (dds).
#' It first creates a DESeqDataSet from the count matrix and metadata, filters genes
#' with low counts, and then performs differential expression analysis using DESeq.
#'
#' @param count_matrix Count matrix of gene expression data.
#' @param coldata Metadata containing information about samples.
#'
#' @return DESeq2 object containing differential expression analysis results.
#'
#' @examples
#' count_path <- "data/count.csv"
#' metadata_path <- "data/metadata.txt"
#' result <- read_and_preprocess(count_path, metadata_path)
#' count_matrix <- result$count_matrix
#' coldata <- result$coldata
#' dds <- perform_differential_expression(count_matrix, coldata)
#'
differential_expression <- function(count_matrix, coldata) {

    # Create a DESeqDataSet from the count matrix
    suppressMessages(dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition))

    # Pre-filter the genes which have low counts. Low count genes may not have sufficient evidence for differential gene expression.
    # Furthermore, removing low count genes reduce the load of multiple hypothesis testing corrections.
    # Here, I will remove the genes which have < 10 reads (this can vary based on research goal) in total across all the samples.
    # Pre-filtering helps to remove genes that have very few mapped reads, reduces memory, and increases the speed of the DESeq2 analysis
    dds <- dds[rowSums(counts(dds)) >= 10,]

    # Perform differential expression analysis
    suppressMessages(dds <- DESeq(dds))

    dds
}

#' Extracts and writes the results of differential gene expression analysis.
#'
#' This function takes a DESeq2 object (dds) and a file path (result_file_path),
#' and extracts the differential gene expression analysis results from the dds object.
#' The results are ordered by adjusted p value (Benjamini-Hochberg FDR method),
#' and then written to a CSV file specified by the result_file_path.
#' The function returns the extracted results.
#'
#' @param dds DESeq2 object containing differential expression analysis results.
#' @param result_file_path File path to write the results CSV file.
#'
#' @return The extracted differential gene expression analysis results.
#'
#' @examples
#' result_file_path <- "results.csv"
#' extracted_results <- extract_and_write_results(dds, result_file_path)
#'
extract_and_write_results <- function(dds, p_value, out_file="results/condition_infected_vs_control_dge.csv") {
    # Get gene expression table
    res <- results(dds)

    # see all comparisons (here there is only one)
    resultsNames(dds)

    # keep gene names as rownames
    rownames(res) <- sub(".*\\|", "", rownames(res))

    # Order gene expression table by adjusted p value (Benjamini-Hochberg FDR method)
    res = res[order(res$padj),]

    # create results directory
    dir.create("results", recursive = TRUE)

    # Export differential gene expression analysis table to CSV file
    write.csv(as.data.frame(res), file = out_file)

    # Get summary of differential gene expression with adjusted p value cut-off at 0.05,
    summary(results(dds, alpha=p_value))

    res
}

# Normalizes counts
get_normalized_counts <- function(dds){
    normalized_counts <- counts(dds, normalized=TRUE)
    rownames(normalized_counts) <- sub(".*\\|", "", rownames(normalized_counts))

    normalized_counts
}

# The shrinkage of effect size (LFC) helps to remove the low count genes (by shrinking towards zero).
# The low or highly variable read count genes can give large estimates of LFCs which may not represent
# true difference in changes in gene expression between two conditions.
lfc_shrink <- function(dds){
    resLFC <- lfcShrink(dds, coef="condition_infected_vs_control", type="apeglm")
    rownames(resLFC) <- sub(".*\\|", "", rownames(resLFC))

    resLFC
}

MA_plots <- function(res, resLFC) {

    # Visualize the shrinkage estimation of LFCs with MA plot and compare it without shrinkage of LFCs,
    par(mfrow = c(1, 2))
    plotMA(resLFC, main="Shrinkage of LFCs", ylim=c(-4,4))
    plotMA(res, main="No shrinkage of LFCs", ylim=c(-4,4))

}

PCA_plot <- function(dds){

    # Calculate the PCA
    vsd <- vst(dds, blind=FALSE)
    pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    # Create the PCA plot
    ggplot(pcaData, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      coord_fixed() +
      theme_bw() +
      scale_color_discrete(name="Condition") +
      ylim(-40,40)
      # xlim(-20, 20) +
}

volcano_plot <- function(res, p_value_cutoff, logFC_cutoff) {
    EnhancedVolcano(toptable = res,              # We use the shrunken log2 fold change as noise associated with low count genes is removed
                x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "padj",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(res),
                pCutoff = p_value_cutoff,
                title = paste("(fold change cutoff = ", logFC_cutoff, "p-value cutoff = ", p_value_cutoff, ")"),
                FCcutoff = logFC_cutoff,
                )
}

heatmap_plot <- function(normalized_counts, logFC_cutoff) {

    counts_scaled =
      normalized_counts %>%
      t(.) %>%                              # transpose to have the genes in columns
      scale() %>%                           # scale(x, center = TRUE, scale = TRUE)
      t(.)                                  # back in original shape

    genes_differential_fold =
     res %>%
     as.data.frame() %>%
     rownames_to_column("gene") %>%
     filter(log2FoldChange > logFC_cutoff) %>%
     dplyr::select(gene)

    counts_scaled_filtered_high_fold_change =
     counts_scaled[row.names(counts_scaled) %in% genes_differential_fold$gene, ]

    # cluster genes and samples
    pheatmap(counts_scaled_filtered_high_fold_change,
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            show_rownames = FALSE,
            show_colnames = TRUE,
            clustering_method = "average",
            main = "Clustering of genes and samples (high fold only)")

    counts_scaled_filtered_high_fold_change
}

run_gseGO <- function(res, organism, p_value) {
    # we want the log2 fold change
    original_gene_list <- res$log2FoldChange

    # name the vector
    names(original_gene_list) <- rownames(res)

    # omit any NA values
    gene_list<-na.omit(original_gene_list)

    # sort the list in decreasing order (required for clusterProfiler)
    gene_list = sort(gene_list, decreasing = TRUE)

    gseBP <- suppressMessages(gseGO(geneList=gene_list,
                 ont ="BP",
                 keyType = "SYMBOL",
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = p_value,
                 verbose = TRUE,
                 OrgDb = organism,
                 pAdjustMethod = "none"))

    gseMF <- suppressMessages(gseGO(geneList=gene_list,
                 ont ="MF",
                 keyType = "SYMBOL",
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = p_value,
                 verbose = TRUE,
                 OrgDb = organism,
                 pAdjustMethod = "none"))

    gseCC <- suppressMessages(gseGO(geneList=gene_list,
                 ont ="CC",
                 keyType = "SYMBOL",
                 nPerm = 10000,
                 minGSSize = 3,
                 maxGSSize = 800,
                 pvalueCutoff = p_value,
                 verbose = TRUE,
                 OrgDb = organism,
                 pAdjustMethod = "none"))

    list(original_gene_list=original_gene_list, gseBP=gseBP, gseMF=gseMF, gseCC=gseCC)
}

plot_gseGO <- function(gseGO_results) {

    print(dotplot(gseGO_results$gseBP, showCategory=10, title = "Biological Processes", split=".sign") + facet_grid(.~.sign))
    print(dotplot(gseGO_results$gseMF, showCategory=10, title = "Molecular Function", split=".sign") + facet_grid(.~.sign))
    print(dotplot(gseGO_results$gseCC, showCategory=10, title = "Cell Components", split=".sign") + facet_grid(.~.sign))

    bp = pairwise_termsim(gseGO_results$gseBP)
    bp_plot = emapplot(bp, showCategory = 20) + ggtitle("Biological Processes")
    ggsave(bp_plot, file = "results/biological_function.pdf", width=10, height=10)

    mf = pairwise_termsim(gseGO_results$gseMF)
    mf_plot = emapplot(mf, showCategory = 20) + ggtitle("Molecular Function")
    ggsave(mf_plot, file = "results/molecular_function.pdf", width=10, height=10)

    cc = pairwise_termsim(gseGO_results$gseCC)
    cc_plot = emapplot(cc, showCategory = 20) + ggtitle("Cell Components")
    ggsave(cc_plot, file = "results/cellular_components.pdf", width=10, height=10)
}

run_gseKEGG <- function(original_gene_list, organism, kegg_organism, reactome_organism, wikiPathways_organism, p_value) {
    # Convert gene IDs for gseKEGG function
    # We will lose some genes here because not all IDs will be converted
    ids<-bitr(names(original_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

    # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
    dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

    # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
    df2 = res[rownames(res) %in% dedup_ids$SYMBOL,]
    df2 = df2[!duplicated(rownames(df2)),]

    # Create a new column in df2 with the corresponding ENTREZ IDs
    df2$Y = dedup_ids$ENTREZID

    # Create a vector of the gene unuiverse
    kegg_gene_list <- df2$log2FoldChange

    # Name vector with ENTREZ ids
    names(kegg_gene_list) <- df2$Y

    # omit any NA values
    kegg_gene_list<-na.omit(kegg_gene_list)

    # sort the list in decreasing order (required for clusterProfiler)
    kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

    kegg <- suppressMessages(gseKEGG(geneList = kegg_gene_list,
                   organism = kegg_organism,
                   nPerm = 10000,
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff = p_value,
                   pAdjustMethod = "none",
                   verbose = TRUE))

    reactome <- suppressMessages(gsePathway(geneList = kegg_gene_list,
                   organism     = reactome_organism,
                   nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = p_value,
                   pAdjustMethod = "none",
                   verbose = TRUE))

    wikiPathways <- suppressMessages(gseWP(geneList = kegg_gene_list,
                   organism     = wikiPathways_organism,
                   nPerm        = 10000,
                   minGSSize    = 3,
                   maxGSSize    = 800,
                   pvalueCutoff = p_value,
                   pAdjustMethod = "none",
                   verbose = TRUE))

    list(kegg_gene_list=kegg_gene_list, kegg=kegg, reactome=reactome, wikiPathways=wikiPathways)
}


plot_gseKEGG <- function(gseKEGG_results) {

    print(dotplot(gseKEGG_results$kegg, showCategory = 10, title = "KEGG" , split=".sign") + facet_grid(.~.sign))
    print(dotplot(gseKEGG_results$reactome, showCategory = 10, title = "Reactome" , split=".sign") + facet_grid(.~.sign))
    print(dotplot(gseKEGG_results$wikiPathways, showCategory = 10, title = "WikiPathways" , split=".sign") + facet_grid(.~.sign))

    kegg = pairwise_termsim(gseKEGG_results$kegg)
    kegg_plot = emapplot(kegg, showCategory = 10, cex_label_category = 1, group_category=TRUE, group_legend=TRUE)
    ggsave(kegg_plot, file = "results/kegg_emap.pdf", width=10, height=10)

    react = pairwise_termsim(gseKEGG_results$reactome)
    react_plot = emapplot(react, showCategory = 10, cex_label_category = 1, group_category=TRUE, group_legend=TRUE)
    ggsave(react_plot, file = "results/reactome_emap.pdf", width=10, height=10)

    wp = pairwise_termsim(gseKEGG_results$wikiPathways)
    wp_plot = emapplot(wp, showCategory = 10, cex_label_category = 1, group_category=TRUE, group_legend=TRUE)
    ggsave(react_plot, file = "results/wikiPathways_emap.pdf", width=10, height=10)
}

pathview_plots <- function(kegg_gene_list, kegg_organism, kegg, num_hits) {

    kegg_top_hits  <- head(kegg, num_hits)

    # combine ID and Description of pathway for renaming files.
    combined_columns <- unite(kegg_top_hits, name, ID, Description, sep = "_")
    hsa_names = head(kegg_top_hits$ID, num_hits)
    kegg_names <- head(combined_columns$name, num_hits)

    # visualize KEGG pathways using pathview
    for (count in seq_along(hsa_names)) {
        hsa = hsa_names[count]
        name = kegg_names[count]
        suppressMessages(pathview(gene.data=kegg_gene_list, pathway.id=hsa, species = kegg_organism, file.format = "png"))
        move_rename_files(hsa, name, source="kegg")
    }

}


move_rename_files <- function(hsa, name, source) {
    # Specify the source directory where the files are located
    source_directory <- getwd()

    # Specify the destination directory where the files will be moved
    destination_directory <- paste0("results/", source, "/", name)
    dir.create(destination_directory, recursive = TRUE)

    # Specify the wildcard pattern to match the files to be moved
    wildcard_pattern <- paste0(hsa, "*")  # Example: Move all text files

    # Get a list of files matching the wildcard pattern in the source directory
    file_list <- list.files(path = source_directory, pattern = wildcard_pattern, full.names = TRUE)

    # Create an empty string array
    file_renamed_list <- character()

    # Iterate over each file and rename it
    for (file_path in file_list) {
      new_file_path <- gsub(hsa, name, file_path)
      file.rename(file_path, new_file_path)
      file_renamed_list <- c(file_renamed_list, new_file_path)
    }

    for (file in file_renamed_list) {
        new_file_path <- file.path(destination_directory, basename(file))
        if (file.exists(file)) {
          file.rename(from = file, to = new_file_path)
        } else {
          warning(paste("File not found:", file))
        }
    }

}

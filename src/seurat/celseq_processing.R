library("Seurat", quietly = TRUE)
library("SeuratDisk", quietly = TRUE)
library("dplyr", quietly = TRUE)


load_celseq_data <- function(expression_mtx, metadata_mtx) {
  #' Loads CelSeq data 
  #' 
  #' Loads CelSeq data into a Seurat object for analysis. 
  #' Assumes that in the expression matrix, the feature (gene) names are in their own column named "gene"
  #' Assumes that in the metadata matrix, the cell names are in their own column named "cell_name"
  #' 
  #' @param expression_mtx the single-cell expression matrix to load into a Seurat object
  #' @param metadata_mtx the single-cell metadata matrix to load into a Seurat object
  
  gene_df <- read.table(expression_mtx, sep="\t", header=TRUE)
  gene_df[is.na(gene_df)] <- 0
  
  # Seurat object needs features as rows and not features as their own column
  rownames(gene_df) <- gene_df$gene
  gene_df <- gene_df %>% dplyr::select(-gene)
  
  metadata_df <- read.table(metadata_mtx, sep="\t", header=TRUE)
  rownames(metadata_df) <- metadata_df$cell_name
  metadata_df <- metadata_df %>% dplyr::select(-cell_name)
  
  celseq_object = CreateSeuratObject(counts=gene_df, meta.data=metadata_df)
  
  return(celseq_object)
}


add_cell_scores <- function(celseq_object, cell_score_dir) {
  #' Add Cell Scores
  #' 
  #' Load the scores calculated by scDRS into the Seurat object's metadata
  #' Scores are generated on a per-cell basis based on the representation of trait-associated genes in a cell's expression matrix
  #' 
  #' @param celseq_object the Seurat object to add the cell scores to
  #' @param cell_score_dir the directory that the cell score files are located in. these scores are assumed to be in the same 
  #'    format as output from scDRS, with the cell IDs located in a column called "cell_name"

  object_metadata <- celseq_object@meta.data
  
  # create a new column using the cell names since this will be merged on to 
  # associate cell-scores with the appropriate cell
  object_metadata$cell_name <- rownames(object_metadata)
  
  # add each set of cell scores found in the cell score directory to the Seurat
  # object metadata
  for (cell_score_file in list.files(cell_score_dir, full.names=TRUE)) {
    cell_score_df <- read.table(cell_score_file, sep="\t", header=TRUE)
    
    suffix <- sub(pattern=".tsv", "", basename(cell_score_file))
    disease_score_colname <- sprintf("norm_score_%s", suffix)
    pval_colname <- sprintf("disease_nlog10_pval_%s", suffix)
    colnames(cell_score_df)[colnames(cell_score_df) == "norm_score"] <- disease_score_colname
    colnames(cell_score_df)[colnames(cell_score_df) == "nlog10_pval"] <- pval_colname
    
    object_metadata <- merge(object_metadata, cell_score_df[c("cell_name", disease_score_colname, pval_colname)], by="cell_name")
  }
  
  cell_scores <- dplyr::select(object_metadata, matches("cell_name|nlog|norm_score"))
  cell_scores <- distinct(cell_scores)
  rownames(cell_scores) <- cell_scores$cell_name
  cell_scores <- cell_scores %>% dplyr::select(-cell_name)
  
  celseq_object <- AddMetaData(celseq_object, cell_scores)
  
  return(celseq_object)
}

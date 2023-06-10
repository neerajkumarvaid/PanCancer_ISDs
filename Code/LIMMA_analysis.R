library(Biobase)
library(BiocManager)
library(data.table)
library(dplyr)
library(FactoMineR)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(limma)
library(MTLR)
library(openxlsx)
library(org.Hs.eg.db)
library(parallel)
library(pheatmap)
library(randomForestSRC)
library(RColorBrewer)
library(survival)
library(svglite)
library(umap)


projectDir <- paste0('cancer-topics')
setwd(projectDir)

types <- c(
  BLCA = "Bladder urothelial carcinoma",
  BRCA = "Breast invasive carcinoma",
  CESC = "Cervical and endocervical cancers",
  COAD = "Colon adenocarcinoma",
  GBM  = "Glioblastoma multiforme",
  HNSC = "Head and neck squamous cell carcinoma",
  KIRC = "Kidney renal clear cell carcinoma",
  KIRP = "Kidney renal papillary cell carcinoma",
  LAML = "Acute myeloid leukemia",
  LGG  = "Brain lower grade glioma",
  LIHC = "Liver hepatocellular carcinoma",
  LUAD = "Lung adenocarcinoma",
  LUSC = "Lung squamous cell carcinoma",
  OV   = "Ovarian serous cystadenocarcinoma",
  PAAD = "Pancreatic adenocarcinoma",
  PRAD = "Prostate adenocarcinoma",
  SKCM = "Skin cutaneous melanoma",
  STAD = "Stomach adenocarcinoma",
  THCA = "Thyroid carcinoma",
  UCEC = "Uterine corpus endometrial carcinoma"
)

# Data prep ====================================================================
# Gene expression data ---------------------------------------------------------
# Gene expression data file retrieved from cbioportal
exprs_mat <- fread(
  file = 'EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv',
  sep  = '\t',
  header = T,
  na.strings = NULL,
  stringsAsFactors = F,
  logical01 = F,
  data.table = F
)
gene_names <- exprs_mat[ ,1]
exprs_mat[[1]] <- NULL
exprs_mat <- as.data.frame(lapply(exprs_mat, as.numeric))
row.names(exprs_mat) <- gene_names
colnames(exprs_mat) <- gsub("\\.", "-", colnames(exprs_mat))
exprs_mat <- log2(1 + exprs_mat)

# Clinical data ----------------------------------------------------------------
# Data file retrieved from cbioportal
clinical_data <- fread(
  file = "clinical_PANCAN_patient_with_followup.tsv",
  sep = "\t",
  header = T,
  na.strings = NULL,
  stringsAsFactors = F,
  logical01 = T,
  data.table = F
)
sort(table(
  unlist(clinical_data)[grepl("\\[[A-Za-z0-9 ]*\\]", unlist(clinical_data))]
))
sel <- clinical_data == "[Not Available]"  |
       clinical_data == "[Not Applicable]" |
       clinical_data == "[Not Available]"  |
       clinical_data == "[Not Evaluated]"  |
       clinical_data == "[Not Reported]"   |
       clinical_data == "[Unknown]"        |
       clinical_data == "[Discrepancy]"
clinical_data[sel] <- NA
# There are a bunch of columns whose data are delimited on pipes ("|") but
# without some documentation the significance of the delimited data is unclear.
sel <- sapply(clinical_data, function(x) any(grepl("\\|", x)))
clinical_data <- clinical_data[ ,!sel]

# Survival data ----------------------------------------------------------------
surv_data <- fread(
  file = "survivaldata_7922cases.csv",
  sep = ",",
  header = T,
  na.strings = NULL,
  stringsAsFactors = F,
  logical01 = T,
  data.table = F
)
row.names(surv_data) <- surv_data$Sample_ID
surv_data$Sample_ID <- NULL

# NNMF gene-factor matrix ------------------------------------------------------
gf_mat <- fread(
  file = "nmf_basis_rank100_PanCancer.csv",
  sep = ",",
  header = T,
  na.strings = NULL,
  stringsAsFactors = F,
  logical01 = T,
  data.table = F
)
row.names(gf_mat) <- gf_mat$Gene_Id
gf_mat$Gene_Id <- NULL
gf_mat$V1 <- NULL  # Not sure what this column is, maybe indexing from 0
anyNA(match(row.names(gf_mat), row.names(exprs_mat)))

# NNMF patient-factor matrix ---------------------------------------------------
pf_mat <- fread(
  file = "nmf_coefficient_rank100_PanCancer_10274cases.csv",
  sep = ",",
  header = T,
  na.strings = NULL,
  stringsAsFactors = F,
  logical01 = T,
  data.table = F
)
row.names(pf_mat) <- pf_mat$Sample_ID
pf_mat$Sample_ID <- NULL
m <- match(row.names(pf_mat), colnames(exprs_mat))
anyNA(m)  # FALSE: all pf_mat sample IDs map to exprs_mat sample IDs
exprs_mat <- exprs_mat[ ,m]
# # There are a bunch of patient IDs that do not map to
# # clinical_data$bcr_patient_barcode. Save these in a separate data frame.
# sel <- is.na(match(pf_mat$Patient_ID, clinical_data$bcr_patient_barcode))
# pf_mat_no_cd <- pf_mat[sel, ]   # i.e. subset w/ no clinical data available
# pf_mat_cd    <- pf_mat[!sel, ]  # i.e. subset w/ clinical data available
# # Reorder the clincal data to match the order of pf_mat_cd
# m <- match(pf_mat_cd$Patient_ID, clinical_data$bcr_patient_barcode)
# clinical_data <- clinical_data[m, ]
# pf_mat_no_cd$Patient_ID <- NULL
# pf_mat_cd$Patient_ID <- NULL
pf_mat$Patient_ID <- NULL

# Combine patient-factor and survival data into a working set ------------------
m <- match(row.names(surv_data), row.names(pf_mat))
surv_data$Type2 <- types[surv_data$Type]
set <- cbind(surv_data, pf_mat[m, ])
table(set$Type2)
table(set$vital_status)
set <- set[set$vital_status %in% c("Alive", "Dead"), ]

# Map features to names, IDs ---------------------------------------------------
re <- "(?<=\\|)[0-9]+"
m <- regexpr(re, row.names(gf_mat), perl = T)
entrez <- regmatches(row.names(gf_mat), m)
feat_mapping <- data.frame(
  feature     = row.names(gf_mat),
  ENTREZID    = entrez,
  ENSEMBL     = mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
                       column = "ENSEMBL", multiVals = "first"),
  ENSEMBLPROT = mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
                       column = "ENSEMBLPROT", multiVals = "first"),
  SYMBOL      = mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
                       column = "SYMBOL", multiVals = "first"),
  GENENAME    = mapIds(org.Hs.eg.db, keys = entrez, keytype = "ENTREZID",
                       column = "GENENAME", multiVals = "first")
)

# ==============================================================================

# Analysis =====================================================================

# Top genes deferentially expressed in relation to patient scores --------------
# The interpretation of the model is gene expression per unit increase in score
# Useful because the NNMF basis values do not tell you the nature of a factor's
# relationship with a gene, it only tells you how important it the gene is. With
# differential expression analysis we can see the direction of the relationship
# between factors and gene expression.

# Note: To replicate the top_genes data in the github repository for GSEA 
# analysis the number should be changed to 16335
if (any(dir("data") %in% "nnmf_pancancer_f100_limma_top500genes.xlsx")) {
  message("What, you expected something to happen?")
} else {
  top_genes <- mclapply(colnames(pf_mat), function(x) {
    mm <- model.matrix(as.formula(paste("~", x)), data = pf_mat)
    fit <- lmFit(setNames(exprs_mat, row.names(pf_mat)), design = mm)
    topTable(
      eBayes(fit),
      coef          = 2,
      sort.by       = "B",
      adjust.method = "fdr",
      number        = 500  #Change to 16335 for GSEA purposes
    )
  }, mc.cores = 2)
  names(top_genes) <- paste0("F", 1:length(top_genes))
  
  wb <- createWorkbook()
  invisible(mapply(
    function(data, name) {
      addWorksheet(wb, sheetName = name)
      m <- match(row.names(data), feat_mapping$feature)
      data$GENENAME <- feat_mapping$GENENAME[m]
      writeData(wb, sheet = name, data, rowNames = T)
    },
    top_genes, names(top_genes)
  ))
  saveWorkbook(wb, file = "nnmf_pancancer_f100_limma_top500genes.xlsx",
               overwrite = T)
}

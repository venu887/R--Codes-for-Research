## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(progress = FALSE)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)

## ----results = 'hide', message=FALSE, warning=FALSE, eval = F-----------------
  query <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Gene expression",
    data.type = "Gene expression quantification",
    platform = "Illumina HiSeq",
    file.type  = "normalized_results",
    experimental.strategy = "RNA-Seq",
#    barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
    access = "open",
    legacy = TRUE
  )
setwd("/Users/user/Documents/R/TCGAbiolinks")
  GDCdownload(query, method = "api", files.per.chunk = 10)
  data <- GDCprepare(query)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
data <- gbm.exp.legacy

## ----message=FALSE, warning=FALSE---------------------------------------------
# Gene expression aligned against hg19.
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
# Only first 100 to make render faster
datatable(assay(data)[1:4,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

rowRanges(data)

## ----results = 'hide', message=FALSE, warning=FALSE, eval=FALSE---------------
#  # Gene expression aligned against hg38
  query <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
#    barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
     access = "open",
#     sample.type= c("Primary Tumor","Solid Tissue Normal","Metastatic")
    sample.type= c("Primary Tumor","Solid Tissue Normal")
#    sample.type= c("Metastatic")
 )
  GDCdownload(query)
  data <- GDCprepare(query)

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
data <- gbm.exp.harmonized

## ----message=FALSE, warning=FALSE---------------------------------------------
datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

datatable(assay(data)[1:4,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)

rowRanges(data)

## ----message=FALSE, warning=FALSE, eval =FALSE--------------------------------
#  #-------------------------------------------------------
#  # Example to idat files from TCGA projects
#  #-------------------------------------------------------
#  projects <- TCGAbiolinks:::getGDCprojects()$project_id
#  projects <- projects[grepl('^TCGA',projects,perl=T)]
#  match.file.cases.all <- NULL
#  for(proj in projects){
#    print(proj)
#    query <- GDCquery(
#      project = proj,
#      data.category = "Raw microarray data",
#      data.type = "Raw intensities",
#      experimental.strategy = "Methylation array",
#      legacy = TRUE,
#      file.type = ".idat",
#      platform = "Illumina Human Methylation 450"
#    )
#    match.file.cases <- getResults(query,cols=c("cases","file_name"))
#    match.file.cases$project <- proj
#    match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
#    tryCatch(
#      GDCdownload(query, method = "api", files.per.chunk = 20),
#      error = function(e) GDCdownload(query, method = "client")
#    )
#  }
#  # This will create a map between idat file name, cases (barcode) and project
#  readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
#  # code to move all files to local folder
#  for(file in dir(".",pattern = ".idat", recursive = T)){
#    TCGAbiolinks::move(file,basename(file))
#  }

## ---- eval = FALSE------------------------------------------------------------
#  query_meth.hg19 <- GDCquery(
#    project= "TCGA-LGG",
#    data.category = "DNA methylation",
#    platform = "Illumina Human Methylation 450",
#    barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05"),
#    legacy = TRUE
#  )
#  GDCdownload(query_meth.hg19)
#  data.hg19 <- GDCprepare(query_meth.hg19)

## ---- eval = FALSE------------------------------------------------------------
#  query <- GDCquery(
#    project = "TCGA-GBM",
#    data.category = "Protein expression",
#    legacy = TRUE,
#    barcode = c("TCGA-OX-A56R-01A-21-A44T-20","TCGA-08-0357-01A-21-1898-20")
#  )
#  GDCdownload(query)
#  data <- GDCprepare(
#    query, save = TRUE,
#    save.filename = "gbmProteinExpression.rda",
#    remove.files.prepared = TRUE
#  )

## ---- eval = FALSE------------------------------------------------------------
#  # Aligned against Hg19
#  query.exp.hg19 <- GDCquery(
#    project = "TCGA-GBM",
#    data.category = "Gene expression",
#    data.type = "Gene expression quantification",
#    platform = "Illumina HiSeq",
#    file.type  = "normalized_results",
#    experimental.strategy = "RNA-Seq",
#    barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
#    legacy = TRUE
#  )
#  GDCdownload(query.exp.hg19)
#  data <- GDCprepare(query.exp.hg19)

## ---- eval = FALSE------------------------------------------------------------
#  query <- GDCquery(
#    project = "TCGA-ACC",
#    data.category = "Copy Number Variation",
#    data.type = "Copy Number Segment",
#    barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01")
#  )
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ---- eval = FALSE------------------------------------------------------------
#  query <- GDCquery(
#    project = "TCGA-ACC",
#    data.category = "Copy Number Variation",
#    data.type = "Gene Level Copy Number Scores",
#    access = "open"
#  )
#  GDCdownload(query)
#  data <- GDCprepare(query)

## ---- eval = FALSE------------------------------------------------------------
#  # mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
#  query.exp.hg38 <- GDCquery(
#    project = "TCGA-GBM",
#    data.category = "Transcriptome Profiling",
#    data.type = "Gene Expression Quantification",
#    workflow.type = "HTSeq - FPKM-UQ",
#    barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01")
#  )
#  GDCdownload(query.exp.hg38)
#  expdat <- GDCprepare(
#    query = query.exp.hg38,
#    save = TRUE,
#    save.filename = "exp.rda"
#  )

## ---- eval = FALSE------------------------------------------------------------
#  library(TCGAbiolinks)
#  query.mirna <- GDCquery(
#    project = "TARGET-AML",
#    experimental.strategy = "miRNA-Seq",
#    data.category = "Transcriptome Profiling",
#    barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
#    data.type = "miRNA Expression Quantification"
#  )
#  GDCdownload(query.mirna)
#  mirna <- GDCprepare(
#    query = query.mirna,
#    save = TRUE,
#    save.filename = "mirna.rda"
#  )
#  
#  
#  query.isoform <- GDCquery(
#    project = "TARGET-AML",
#    experimental.strategy = "miRNA-Seq",
#    data.category = "Transcriptome Profiling",
#    barcode = c("TARGET-20-PATDNN","TARGET-20-PAPUNR"),
#    data.type = "Isoform Expression Quantification"
#  )
#  GDCdownload(query.isoform)
#  
#  isoform <- GDCprepare(
#    query = query.isoform,
#    save = TRUE,
#    save.filename = "mirna-isoform.rda"
#  )

## ---- eval = FALSE------------------------------------------------------------
#  #--------------------------------------
#  # DNA methylation data
#  #--------------------------------------
#  # DNA methylation aligned to hg38
#  query_met.hg38 <- GDCquery(
#    project= "TCGA-LGG",
#    data.category = "DNA Methylation",
#    platform = "Illumina Human Methylation 450",
#    barcode = c("TCGA-HT-8111-01A-11D-2399-05","TCGA-HT-A5R5-01A-11D-A28N-05")
#  )
#  GDCdownload(query_met.hg38)
#  data.hg38 <- GDCprepare(query_met.hg38)

## ---- eval = FALSE------------------------------------------------------------
#  # Using sesame 	http://bioconductor.org/packages/sesame/
#  # Please cite 10.1093/nar/gky691 and doi: 10.1093/nar/gkt090.
#  library(TCGAbiolinks)
#  proj <- "TCGA-ACC"
#  query <- GDCquery(
#    project = proj,
#    data.category = "Raw microarray data",
#    data.type = "Raw intensities",
#    experimental.strategy = "Methylation array",
#    legacy = TRUE,
#    barcode = c("TCGA-OR-A5JT","CGA-OR-A5LG","TCGA-OR-A5JX"),
#    file.type = ".idat",
#    platform = "Illumina Human Methylation 450"
#  )
#  tryCatch(
#    GDCdownload(query, method = "api", files.per.chunk = 20),
#    error = function(e) GDCdownload(query, method = "client")
#  )
#  betas <- GDCprepare(query)

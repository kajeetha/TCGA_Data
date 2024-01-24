library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

#getting clinical data for TCGA-BRCA cohort
clinical_brca <- GDCquery_clinic("TCGA-BRCA")

any(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death")
    
#index those colnames:
which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death")
clinical_brca[ ,c(9,39,44)]

#from vital_status checking to see how many are dead and alive: 
table(clinical_brca$vital_status)

#need time(overall_survival), status(death, or alive) and event(or strate low high TP53) info to perform survival analysis. Need to censor patients who have not experienced death:
clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE) #if it is alive then set it to false, else set it to true.

#create an "overall survival" variable that is equal to days_to_death for dead patients, and to days_to_last_follow_up for patients who are still alive:
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive", clinical_brca$days_to_last_follow_up, clinical_brca$days_to_death)

#gather gene expression data:
query_brca_all = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq", 
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification", 
  sample.type = "Primary Tumor",
  access = "open")

#we want to store the output of the query into a dataframe object using getResults from the TCGA biolinks package:

output_brca <- getResults(query_brca_all)

#get only the 20 primary tissue sample barcodes for small sample size. This will give a vector of those specific patient samples
tumor <- output_brca[output_brca$sample_type == "Primary Tumor", "cases"][1:20]

query_brca = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq", 
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification", 
  sample.type = "Primary Tumor",
  access = "open",
  barcode = tumor)

#download the data
GDCdownload(query_brca)

#get the read counts
tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = TRUE)
brca_matrix <- assay(tcga_brca_data, "unstranded")

gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(colData(tcga_brca_data))

#create a countdata object for differential gene expression analysis: 1 group implying all samples belong to the same group
dds <- DESeqDataSetFromMatrix(countData = brca_matrix, colData = coldata, design = ~ 1) 

#lets remove genes with less then 10 reads across all samples: 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

#A variance-stabilizing transformation is applied to count data, such as RNA-seq read counts, to stabilize the variance across the range of expression levels. In genomics, the count data is often characterized by higher variability at lower expression levels and a trend of stabilization as expression levels increase.
#this will return another summarized experiment object and save those counts in another matrix:

vsd <- vst(dds, blind=FALSE) #this accounts for the replicates and metadata
brca_matrix_vst <- assay(vsd)

#add gene symbol information from metadata into the matrix: 
brca_tp53 <- brca_matrix_vst %>%
  as.data.frame() %>%
  rownames_to_column(var = 'gene_id') %>%
  gather(key = 'case_id', value = 'counts', -gene_id) %>%
  left_join(brca_matrix_vst, gene_metadata, by = "gene_id") %>%
  filter(gene_name == "TP53")

#divide the cohort into high expressing and low expressing using median value vst count: 
#you can also alternatively divide with quartile ranges, can also use a t-test or mannwhitney u test to bifurcate the cohort. 

median_value <- median(brca_tp53$counts)
brca_tp53$strata <- ifelse(brca_tp53$counts >= median_value, "HIGH", "LOW")

#combine the table with the survival info and status info post adjustment of identifier: 
brca_tp53$case_id <- gsub('-01[AB]. *$', "", brca_tp53$case_id)
brca_tp53 <- merge(brca_tp53, clinical_brca, by.x = "case_id", by.y = "submitter_id")

#compute survival curve
fit <- survfit(Surv(overall_survival, deceased) ~strata, data = brca_tp53)
#plot survival curve
ggsurvplot(fit, data = brca_tp53, pval = T, risk.table = T)
#survdiff( for categorical variables)
  

 
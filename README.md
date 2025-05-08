# scRNA-seq
This is a pipeline for single-cell RNA sequencing (scRNA-seq) data quality control, applicable to 3' scRNA-seq, 5' scRNA-seq, and snoRNA-seq.

## Part I Introduction
### i. Workflow
As illustrated in the figure,
(i) yellow circles represent the steps where commands need to be entered;
(ii) pink dashed rectangular boxes represent the output results after processing at each step.
![image](https://github.com/user-attachments/assets/ccef0d00-c140-4422-ba21-ebdd44a7b063)
### ii.Conda Environment
```
conda install hcc::cellranger
```
### iii.Index Prepared
```
cellranger mkref --genome=GRCh38 \
   --fasta=/mnt/share/FileTransfer/GRCh38.primary_assembly.genome.fa \
   --genes=/mnt/share/FileTransfer/gencode.v47.annotation.gtf \
   --nthreads=30 >GRCh38.mkref.log &
```

## Part II Generation of Data for Analysis
### i.Raw Data Quality Check(QC) & Mapping Data to Genome
```
cellranger count --id=Control1 --localcores=30 \
   --fastqs=/mnt/project/scASI/10X/E-MATB-10378/1.rawdata/ \
   --sample=Control1 \
   --transcriptome=/mnt/project/scASI/GRCh38 \
   --nosecondary \
   --chemistry SC5P-PE > Control1.log &
```
#### Filter the data based on the QC report from CellRanger, web_summary.html
Qualified scRNA-seq data should adhere to the following standards. Any data that does not meet these criteria should be excluded.  
• Summary part  
![image](https://github.com/user-attachments/assets/99200b58-8963-4d65-b4b8-01646423e525)  
(i) Estimated Number of Cells: >3000  
#the discrepancy between Estimated Number of Cells and the expected number of cells in the experiment: >50% & <200%  
(ii) Mean Reads per Cell: >20,000  
(iii) Median Genes per Cell: >800  
• Sequencing part  
![image](https://github.com/user-attachments/assets/dcc69310-811e-45c5-af02-67f7def20bd0)  
(i) Q30 Bases in Barcode: >80%  
(ii) Q30 score of RNA Read: >60%  
• Mapping part  
![image](https://github.com/user-attachments/assets/9e280989-4a27-4c14-9be7-8f9e9db66c48)  
(i) Reads Mapped Confidently to Exonic Regions: >60%  

### ii.Remove low-quality cells and genes
```
dir<-" "
pbmc.data<-Read10X(data.dir=dir)
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
expr_mat <- GetAssayData(pbmc, slot = "counts")
genes_use <- rowSums(expr_mat > 0) > 0.05 * ncol(pbmc)
pic<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
![image](https://github.com/user-attachments/assets/fe68f2e3-bcde-42cd-89fd-0a4910966ce4)  

• Number of counts: >1500  
• Number of genes: >700  
• Percent of mitochondrial transcripts: <15%  
Filtering criteria for some exception cell types  
Cardiac muscle cells, Skeletal muscle cells, Aged cells: <50%  
Neurons, Cancer cells: No filtering  
• Percent of ribosomal transcripts: No filtering  
• Number of cells in which gene is present: >5%  

```
pbmc <- subset(pbmc, subset = nCount_RNA > 1500 &
                 nFeature_RNA > 700 & 
                 percent.mt < 15)
pbmc <- subset(pbmc, features = rownames(expr_mat)[genes_use])
```
• Remove doublets
| Total Cells | Expected Doublet Rate |
|-------------|------------------------|
| ~500        | 0.004 (0.4%)           |
| ~1000       | 0.008 (0.8%)           |
| ~2000       | 0.016 (1.6%)           |
| ~3000       | 0.024 (2.4%)           |
| ~4000       | 0.032 (3.2%)           |
| ~5000       | 0.040 (4.0%)           |
| ~6000       | 0.048 (4.8%)           |
| ~7000       | 0.056 (5.6%)           |
| ~8000       | 0.064 (6.4%)           |
| ~9000       | 0.072 (7.2%)           |
| ~10000      | 0.080 (8.0%)           |
```
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

sweep.res.list <- paramSweep(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)  #choose the best pK
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.056 * nrow(pbmc@meta.data))  # 0.056 is the expected doublet rate for ~7000 cells
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

pbmc <- doubletFinder(pbmc, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
pANN_col <- grep("pANN", colnames(pbmc@meta.data), value = TRUE)
pbmc <- doubletFinder(pbmc, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = FALSE)
doublet_col <- grep("DF.classifications", colnames(pbmc@meta.data), value = TRUE)
colnames(pbmc@meta.data)[which(colnames(pbmc@meta.data) == doublet_col)] <- "Doublet"

pbmc <- subset(pbmc, subset = Doublet == "Singlet")
```









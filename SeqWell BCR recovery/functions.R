#setwd("L:/Duncan/eoepaper")

reticulate::use_python('C:/Users/dmorgan/AppData/Local/Continuum/anaconda3/python.exe', required = TRUE)

library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(feather)
library(dplyr)
library(reshape2)
library(viridis)
library(tidyr)

# update: 
# seurat is renormalizing the data (since python data is already normalized)
# i'm going to implement a brute-forcey work around, for now
# analysis-wise, this requires getting rid of poisson scaling
# read in files transferred from python
pyImport = function(name, rawfile = 'cellsAll.feather') {
  normdata = read_feather(paste0(name, '.feather')) %>% as.data.frame()
  #rownames(data) = data[,1]
  #data = data[,-1]
  
  print('reading in raw data')
  rawdata = read_feather(rawfile) %>% as.data.frame()
  rownames(rawdata) = rawdata[,1]
  rawdata = rawdata[,-1]
  rawdata = rawdata[,colnames(rawdata) %in% colnames(normdata)]
  
  
  metadata = read.csv(paste0(name, '_meta.txt'), row.names = 1, stringsAsFactors = FALSE)
  seurat = CreateSeuratObject(rawdata, min.cells = 5)
  seurat@meta.data = metadata
  seurat
}

# standard Seurat processing
pyProcess = function(seurat) {
  seurat = NormalizeData(seurat) 
  seurat = FindVariableGenes(seurat, do.plot= FALSE)
  seurat = ScaleData(seurat, genes.use =seurat@var.genes, vars.to.regress = c('n_genes'), model.use = 'poisson')
  seurat = RunPCA(seurat, dims.use = seurat@var.genes, do.print = FALSE)
  seurat = RunUMAP(seurat, dims.use = 1:20)
  seurat@meta.data$UMAP1 = seurat@dr$umap@cell.embeddings[,1]
  seurat@meta.data$UMAP2 = seurat@dr$umap@cell.embeddings[,2]
  seurat
}

# take the files from the exportSeurat function, assemble them into a seurat object, and process
pyToSeurat = function(name, rawfile = 'cellsAll.feather') {
  print('reading in data')
  seurat = pyImport(name, rawfile)
  print('processing data')
  seurat = pyProcess(seurat)
  #print('converting to sparse format')
  #seurat = MakeSparse(seurat)
  seurat
}


addUMAP = function(seurat) {
    seurat@meta.data$UMAP1 = seurat@dr$umap@cell.embeddings[,1]
    seurat@meta.data$UMAP2 = seurat@dr$umap@cell.embeddings[,2]
    seurat
} 

seuratToPython = function(seurat, name) {
    write_feather(as.data.frame(as.matrix(seurat@raw.data[,seurat@cell.names])), paste0(name, '_raw.feather'))
    write_feather(seurat@meta.data, paste0(name, '_meta.feather'))
}
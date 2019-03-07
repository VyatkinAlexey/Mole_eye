normalisation_qc_f1000  <- function(file, runQuickCluster=T, output="sce_normalised_object.rds"){
    
    # Аннотацию к этой функции смотри в описании
    
    library(scater)
    library(SingleCellExperiment)
    library(scran)
    
    count_matrix = read.csv(file, row.names=1)
    sce = SingleCellExperiment(assays=list(counts=as.matrix(count_matrix)))
    
    is.spike <- grepl("^ERCC", rownames(sce))
    is.mito <- grepl("^mt-", rownames(sce))
    
    sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))
    head(colnames(colData(sce)))

    pdf("Genes and library sizes statistics.pdf")
    par(mfrow=c(1,2))
    hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",
         breaks=20, col="grey80", ylab="Number of cells")
    hist(sce$total_features, xlab="Number of expressed genes", main="",
         breaks=20, col="grey80", ylab="Number of cells")
    dev.off()
    
    libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
    feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

    mito.drop <- isOutlier(sce$pct_counts_feature_controls_Mt, nmads=3, type="higher")
    spike.drop <- isOutlier(sce$pct_counts_feature_controls_ERCC, nmads=3, type="higher")
    
    if(length(mito.drop) == 0) {
        mito.drop = rep(FALSE, length(libsize.drop))
    }

    if(length(spike.drop) == 0) {
        spike.drop = rep(FALSE, length(libsize.drop))
    }
    
    sce <- sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
    
    ave.counts <- rowMeans(counts(sce))
    keep <- ave.counts >= 1
    
    pdf("Average count number histogram.pdf")
    hist(log10(ave.counts), breaks=100, main="", col="grey80",
         xlab=expression(Log[10]~"average count"))
    abline(v=log10(1), col="blue", lwd=2, lty=2)
    dev.off()
    
    numcells <- nexprs(sce, byrow=TRUE)
    alt.keep <- numcells >= 10

    pdf("Smooth scatter distribution of expression by cells.pdf")
    smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"),
         ylab="Number of expressing cells")
    dev.off()
    
    # normalization
    message("Running normalization. It can take a while, depends on the 
    number of cells.")
    if(runQuickCluster){
        cl <- tryCatch(scran::quickCluster(sce), error=function(e) NULL)
    } else {
        cl <- NULL
    }

    print("compute sizeFactors which will be used for normalization")
                       
    # compute sizeFactors which will be used for normalization
    sceNorm <- scran::computeSumFactors(sce, sizes = sizes, clusters = cl)

    message("summary(sizeFactors(sceObject)):")
    print(summary(sizeFactors(sceNorm)))
    if(length(sizeFactors(sceNorm)[sizeFactors(sceNorm) <= 0]) > 0){
        message("Cells with negative sizeFactors will be deleted before the 
    downstream analysis.")
    }
    sceNorm <- sceNorm[, sizeFactors(sceNorm) > 0]
    sceNorm <- scater::normalize(sceNorm)
    rm(sce)
                       
    saveRDS(sceNorm,file=output)
                       
    return("Succesfull run!")
}
                       
normalisation_qc_seurat <- function(dir, output="seurat_normalised_object.rds"){
    
    library(Seurat)
    library(dplyr)
    
    message("Reading raw data")
    
    raw_data <- Read10X(data.dir = dir)
    
    message("Raw data has read")
    
    
    # creating seurat object 
    seurat_object <- CreateSeuratObject(raw.data = raw_data, min.cells = 3, min.genes = 200, 
    project = "single_cell")
    
    # QC and selecting cells for further analysis
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = seurat_object@data), value = TRUE)
    percent.mito <- Matrix::colSums(seurat_object@raw.data[mito.genes, ])/Matrix::colSums(seurat_object@raw.data)

    # AddMetaData adds columns to object@meta.data, and is a great place to
    # stash QC stats
    seurat_object <- AddMetaData(object = seurat_object, metadata = percent.mito, col.name = "percent.mito")

    # plotting qc metrics
    pdf("Violin plot of Genes, UMI, mito percentage.pdf")
    VlnPlot(object = seurat_object, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, do.return = T)
    dev.off()
    
    # Plotting correlation between charasteristics
    pdf("Feature correlation: nUMI vs percent.mito or nGene.pdf")
    par(mfrow = c(1, 2))
    GenePlot(object = seurat_object, gene1 = "nUMI", gene2 = "percent.mito")
    GenePlot(object = seurat_object, gene1 = "nUMI", gene2 = "nGene")
    dev.off()
    
    
    # filtering and normalisation
    seurat_object <- FilterCells(object = seurat_object, subset.names = c("nGene", "percent.mito"), 
        low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
    seurat_object <- NormalizeData(object = seurat_object, normalization.method = "LogNormalize", 
        scale.factor = 10000)
    
    
    # plotting variable genes
    pdf("Variable genes.pdf")
    seurat_object <- FindVariableGenes(object = seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, 
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    dev.off()
    
    
    # scaling the data
    seurat_object <- ScaleData(object = seurat_object, vars.to.regress = c("nUMI", "percent.mito"))
    
    # converting to SingleCellExperiment
    sce <- Convert(from = seurat_object, to = "sce")
    
    # saving the object
    saveRDS(sce,file=paste(output, "sce.rds", sep="_"))
    saveRDS(seurat_object,file=paste(output, "seurat.rds", sep="_"))
    
    return("Succesfull run!")
}
                       
args = commandArgs(trailingOnly=TRUE)
                       
if(args[1] == "1"){
    normalisation_qc_seurat(dir=args[2], output=args[4])
} else {
    runQuickCluster = (args[3] == "1")
    normalisation_qc_f1000(file=args[2], runQuickCluster=runQuickCluster, output=args[4]) 
}

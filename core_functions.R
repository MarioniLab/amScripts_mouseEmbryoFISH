# load embryo 8.5 data
load_embryo_8.5 = function(filterNullArea = TRUE, threshTotalRNA, filterBigClumps = TRUE){
  
  require(SingleCellExperiment)
  require(scater)
  data.dir = "/Users/alsu/Develop/spatial/mouse_embryo/data/8_5/source/sce_all.Rds"
  
  sce = readRDS( data.dir )
  meta = colData(sce)
  # add libsize
  counts = assay(sce, "counts")
  libsize = colSums( counts )
  meta = cbind(meta, libsize)
  # to sce: add log2 normalization by cell area and libsize
  sizeFactors.libsize <- meta$libsize/mean(meta$libsize)
  sizeFactors.area <- meta$Area/mean(meta$Area)
  assay(sce, "logcounts.libsize") <- log2(t(t(counts)/sizeFactors.libsize) + 1)
  assay(sce, "logcounts.area") <- log2(t(t(counts)/sizeFactors.area) + 1)

  # filter cells with NULL area
  if (filterNullArea) {
    idx = meta$Area > 0
    print( paste0( "Discarded ", sum(meta$Area == 0), " null area cells (out of ", 
                   dim(meta)[1], "): ", round( mean(meta$Area == 0)*100, 2), "%" ))
    meta = meta[idx, ]
    sce = sce[, idx]
  }
  # filter out cells with low number of mRNA molecules
  idx = meta$libsize > threshTotalRNA
  print( paste0( "Discarded ", sum(meta$libsize <= threshTotalRNA), " empty cells (out of ", 
                   dim(meta)[1], "): ", round( mean(meta$libsize <= threshTotalRNA)*100, 2), "%" ))
  meta = meta[idx,]
  sce = sce[, idx]
  # filter big clumps (too big of the area)
  if (filterBigClumps){
    meta$sqrt.area = sqrt(meta$Area)
    pnorm.fit = pnorm(meta$sqrt.area, mean = median(meta$sqrt.area), 
                      sd = mad(meta$sqrt.area) , lower.tail = F)
    thresh.sqrt = min(meta$sqrt.area[which(p.adjust(pnorm.fit, method = "fdr") < 0.05)])
    print( paste0( "Discarded ", sum(meta$sqrt.area > thresh.sqrt), " big clumps (out of ", 
                   dim(meta)[1], "): ", round( mean(meta$sqrt.area > thresh.sqrt)*100, 2), "%" ))
    idx = meta$sqrt.area <= thresh.sqrt
    meta = meta[idx,]
    sce = sce[,idx]
  }
  
  assign("sce", sce, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  
  invisible(0)
}

load_data_atlas = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, load_corrected = FALSE, rownames_ensembl = TRUE){
  
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  
  counts = readRDS("/nfs/research1/marioni/jonny/embryos/data/raw_counts.rds")
  genes = read.table("/nfs/research1/marioni/jonny/embryos/data/genes.tsv", stringsAsFactors = F)
  meta = read.table("/nfs/research1/marioni/jonny/embryos/data/meta.tab", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  if (!rownames_ensembl){
    rownames(counts) = genes[,2] # systematic name
  } else {
    rownames(counts) = genes[,1] #ensembl
  }
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table("/nfs/research1/marioni/jonny/embryos/data/sizefactors.tab", stringsAsFactors = F)[,1]
    sizeFactors(sce) = sfs
    sce = scater::normalize(sce)
  }
  
  if(remove_doublets){
    sce = scater::normalize(sce[,!meta$doublet])
    meta = meta[!meta$doublet,]
  }
  
  if(remove_stripped){
    sce = scater::normalize(sce[,!meta$stripped])
    meta = meta[!meta$stripped, ]
  }
  
  if(load_corrected){
    corrected = readRDS("/nfs/research1/marioni/jonny/embryos/data/corrected_pcas.rds")
    assign("corrected", corrected, envir = .GlobalEnv)
  }
  
  assign("genes.atlas", genes, envir = .GlobalEnv)
  assign("meta.atlas", meta, envir = .GlobalEnv)
  assign("sce.atlas", sce, envir = .GlobalEnv)
  
  invisible(0)
}



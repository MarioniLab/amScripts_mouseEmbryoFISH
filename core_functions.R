# load embryo 8.5 data
load_embryo_8.5 = function(filterNullArea = TRUE, threshTotalRNA, filterBigClumps = TRUE, 
                           dir = "local"){
  
  require(SingleCellExperiment)
  require(scater)
  
  if (dir == "cluster") {
    data.dir = "/nfs/research1/marioni/alsu/spatial/mouse_embryo/data/8_5/source/sce_all.Rds"
  } else {
    data.dir = "/Users/alsu/Develop/spatial/mouse_embryo/data/8_5/source/sce_all.Rds"
  }
  
  sce = readRDS( data.dir )
  meta = colData(sce)
  
  # rename Cavin3 --> Prkcdbp
  rownames.sce = rownames(sce)
  rownames.sce[rownames.sce == "Cavin3"] = "Prkcdbp"
  rownames(sce) = rownames.sce
  
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

# fikter empty cells (for the reduced gene subset)
filterEmptyCells = function(sce, meta){
  counts = counts(sce)
  reduced.libsize = colSums( counts )
  idx = reduced.libsize > 0
  meta = meta[idx,]
  sce = sce[,idx]
  
  assign("sce", sce, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  
  invisible(0)
}

adjust_atlas_to_seq = function(sce.atlas, meta.atlas, sce.seq){
  
  # get rid of doublets and stripped - we don't want to map to them
  idx = meta.atlas$doublet == F & meta.atlas$stripped == F
  meta.atlas = meta.atlas[idx,]
  sce.atlas = sce.atlas[,idx]
  # keep only those genes that are in seqFISH 
  idx = rownames(sce.atlas) %in% rownames(sce.seq)
  sce.atlas = sce.atlas[idx,]
  
  assign("sce.atlas", sce.atlas, envir = .GlobalEnv)
  assign("meta.atlas", meta.atlas, envir = .GlobalEnv)
  invisible(0)
}

reduceAtlas = function(sce.atlas, meta.atlas, thresh){
  celltypes = unique(meta.atlas$celltype)
  celltypes = celltypes[!is.na(celltypes)]
  idx = sapply(celltypes, function(x){
    idx.CT = which(meta.atlas$celltype == x)
    if (length(idx.CT) > thresh) {
      return(sample(idx.CT, thresh))
    } else {
      return(idx.CT)
    }
  })
  idx = unlist(idx)
  meta.atlas = meta.atlas[idx,]
  sce.atlas = sce.atlas[,idx]
  
  assign("sce.atlas", sce.atlas, envir = .GlobalEnv)
  assign("meta.atlas", meta.atlas, envir = .GlobalEnv)
  invisible(0)
}

seqBatchCorrect = function(counts.seq, counts.atlas, meta.atlas){
  require(batchelor)
  set.seed(42)
  unq.samples = unique(meta.atlas$sample)
  counts.seq = t(cosineNorm(counts.seq))
  counts.atlas.bySamples = lapply(unq.samples, function(x){
    out = counts.atlas[,meta.atlas$sample == x]
    return(t(cosineNorm(out)))
  })
  # batch correct atlas
  atlas.corrected = do.call(reducedMNN, c(counts.atlas.bySamples))
  atlas.corrected = atlas.corrected$corrected
  # correction for corrected atlas + seq
  joint.corrected = reducedMNN(atlas.corrected, counts.seq)
  joint.corrected = joint.corrected$corrected
  
  atlas = 1:nrow(atlas.corrected)
  out = list("atlas" = joint.corrected[atlas,], "seq" = joint.corrected[-atlas,])
  return(out)
}

seqPCAandBatchCorrect = function(counts.seq, counts.atlas, meta.atlas, nPC){
  require(irlba)
  require(batchelor)
  set.seed(42)
  unq.samples = unique(meta.atlas$sample)
  counts.seq = cosineNorm(counts.seq)
  counts.atlas = cosineNorm(counts.atlas)
  counts.joint = cbind(counts.atlas, counts.seq)
  counts.pca = prcomp_irlba(t( counts.joint ), n = nPC)$x
  rownames(counts.pca) = colnames(counts.joint) 
  atlas.pca = counts.pca[1:dim(counts.atlas)[2],]
  seq.pca = counts.pca[-(1:dim(counts.atlas)[2]),]
  
  atlas.pca.bySamples = lapply(unq.samples, function(x){
    out = atlas.pca[meta.atlas$sample == x,]
    return(out)
  })
  # batch correct atlas
  atlas.corrected = do.call(reducedMNN, c(atlas.pca.bySamples))
  atlas.corrected = atlas.corrected$corrected
  # correction for corrected atlas + seq
  joint.corrected = reducedMNN(atlas.corrected, seq.pca)
  joint.corrected = joint.corrected$corrected
  
  atlas = 1:nrow(atlas.corrected)
  out = list("atlas" = joint.corrected[atlas,], "seq" = joint.corrected[-atlas,])
  return(out)
}

map_KNN = function(atlas, seq, meta.atlas, meta, k.neigh, mcparam){
  knns = queryKNN(atlas, seq, k = k.neigh, get.index = TRUE, get.distance = TRUE, BPPARAM = mcparam)
  k.mapped = t(apply(knns$index, 1, function(x) meta.atlas$cell[x]))
  celltypes = t(apply(k.mapped, 1, function(x) meta.atlas$celltype[match(x, meta.atlas$cell)]))
  stages = t(apply(k.mapped, 1, function(x) meta.atlas$stage[match(x, meta.atlas$cell)]))
  theiler = t(apply(k.mapped, 1, function(x) meta.atlas$theiler[match(x, meta.atlas$cell)]))
  samples = t(apply(k.mapped, 1, function(x) meta.atlas$sample[match(x, meta.atlas$cell)]))
  mapping = lapply(1:dim(knns$index)[1], function(x){
    out = list(cells.mapped = k.mapped[x,],
               celltypes.mapped = celltypes[x,],
               distances.mapped = knns$distance[x,], 
               stages.mapped = stages[x,],
               theiler.mapped = theiler[x,],
               samples.mapped = samples[x,])
    return(out)
  })
  names(mapping) = meta$uniqueID
  return(mapping)
}

# add regressed by CT assay
addCTregressedLogcounts = function(sce, assay.type, meta){
  # make sure meta and sce are aligned
  meta = meta[order(meta$uniqueID),]
  sce = sce[, order(colnames(sce))]
  
  counts = assay(sce, assay.type)
  counts.perGene = apply(counts, 1, function(x){
    df = data.frame(CT = meta$CT, counts = as.numeric(x), uniqueID = meta$uniqueID, embryo = meta$embryo)
    stat = df%>% group_by(CT, embryo) %>% 
      dplyr::summarize(mean.CT = mean(counts, na.rm=T) , sd.CT = sd(counts, na.rm=T))
    # change sd 0/NA to sd 1 (so we don't divide to 0)
    stat$sd.CT[stat$sd.CT == 0 | is.na(stat$sd.CT)] <- 1
    df = merge(df, stat, by = c("CT", "embryo"), all.x=T)
    df = df[order(df$uniqueID),]
    return( as.numeric( (df$counts - df$mean.CT)/df$sd.CT ) )
  })
  regressed.counts = as.data.frame( t(counts.perGene) )
  rownames(regressed.counts) = rownames(counts)
  colnames(regressed.counts) = meta$uniqueID
  assay(sce, paste0("CT_regressed.", assay.type)) <- regressed.counts
  
  return(sce)
}


load_data_atlas = function(normalise = TRUE, remove_doublets = FALSE, remove_stripped = FALSE, 
                           load_corrected = FALSE, rownames_ensembl = TRUE){
  if(load_corrected & (!remove_doublets | !remove_stripped)){
    message("Using corrected PCs, also removing doublets + stripped now.")
    remove_doublets = TRUE
    remove_stripped = TRUE
  }
  require(scran)
  require(scater)
  require(SingleCellExperiment)
  require(Matrix)
  
  root.dir = "/nfs/research1/marioni/jonny/embryos/data/"
  #root.dir = "/Users/alsu/Develop/FetalAlcoholSyndrome/atlas_Jonny/data/"
  counts = readRDS(paste0( root.dir, "raw_counts.rds"))
  genes = read.table(paste0( root.dir, "genes.tsv"), stringsAsFactors = F)
  meta = read.table(paste0( root.dir, "meta.tab"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "$")
  
  if (!rownames_ensembl){
    rownames(counts) = genes[,2] # systematic name
  } else {
    rownames(counts) = genes[,1] #ensembl
  }
  colnames(counts) = meta$cell
  
  sce = SingleCellExperiment(assays = list("counts" = counts))
  
  if(normalise){
    sfs = read.table(paste0(root.dir, "sizefactors.tab"), stringsAsFactors = F)[,1]
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
    corrected = readRDS(paste0(root.dir, "corrected_pcas.rds"))
    assign("corrected", corrected, envir = .GlobalEnv)
  }
  assign("genes.atlas", genes, envir = .GlobalEnv)
  assign("meta.atlas", meta, envir = .GlobalEnv)
  assign("sce.atlas", sce, envir = .GlobalEnv)
  assign("sfs.atlas", sfs, envir = .GlobalEnv)
  invisible(0)
}

getCorrForOneCell = function(uniqueID, counts, meta, cell.locations){
  require(fields)
  
  transcriptional.dist = rdist( t(counts[, uniqueID]) , t(counts) )
  physical.dist = rdist( cell.locations[meta$uniqueID == uniqueID, ], cell.locations )
  
  out = as.numeric(cor(t(physical.dist), t(transcriptional.dist), method="pearson"))
  return(out)
}


getmode <- function(v, dist) {
  tab = table(v)
  #if tie, break to shortest distance
  if(sum(tab == max(tab)) > 1){
    tied = names(tab)[tab == max(tab)]
    sub = dist[v %in% tied]
    names(sub) = v[v %in% tied]
    return(names(sub)[which.min(sub)])
  } else {
    return(names(tab)[which.max(tab)])
  }
}

getHVGs = function(sce, min.mean = 1e-3, gene_df = genes){
  require(biomaRt)
  trend = scran::trendVar(sce, loess.args = list(span = 0.05))
  decomp = scran::decomposeVar(sce, fit = trend)
  decomp = decomp[decomp$mean > min.mean,]
  
  #exclude sex genes
  xist = "ENSMUSG00000086503"
  xist = "Xist"
  # mouse_ensembl = useMart("ensembl")
  # mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
  # gene_map = getBM(attributes=c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
  # ychr = gene_map[gene_map[,2] == "Y", 1]
  ychr = read.table("/nfs/research1/marioni/jonny/embryos/data/ygenes.tab", stringsAsFactors = FALSE)[,2]
  #ychr = read.table("/Users/alsu/Develop/FetalAlcoholSyndrome/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
  
  other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, ychr, other),]
  
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}




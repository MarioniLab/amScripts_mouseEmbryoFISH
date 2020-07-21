# load embryo 8.5 data
load_embryo_8.5 = function(dir = "local"){
  
  require(SingleCellExperiment)
  require(scater)
  
  if (dir == "cluster") {
    data.dir = "/nfs/research1/marioni/alsu/spatial/mouse_embryo/data/8_5/source/"
  } else {
    data.dir = "/Users/alsu/Develop/spatial/mouse_embryo/data/8_5/source/"
  }
  
  sce = readRDS( paste0( data.dir , "E8.5_sce_filt_unlabelled.Rds"))
  # add normalization by libsize 
  assay(sce, "cpm") <- logcounts(scater::logNormCounts(sce, size_factors = sce$total))
  
  meta = colData(sce)
  meta = data.frame(meta)
  
  # rename Cavin3 --> Prkcdbp
  rownames.sce = rownames(sce)
  rownames.sce[rownames.sce == "Cavin3"] = "Prkcdbp"
  rownames(sce) = rownames.sce
  
  assign("sce", sce, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
  invisible(0)
}

getNumGenes = function(sce){
  counts = counts(sce)
  n.genes = sapply(1:ncol(counts), function(x){
    current.counts = counts[,x]
    return(sum(current.counts > 0))
  })
  return(n.genes)
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

seqPCAandBatchCorrect = function(counts.seq, counts.atlas, meta.atlas, nPC){
  require(irlba)
  require(batchelor)
  set.seed(2020)
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
                           load_corrected = FALSE, rownames_ensembl = TRUE, stages = c("E8.25", "E8.5") ){
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
  # only keep relevant stages
  idx = meta$stage %in% stages
  meta = meta[idx,]
  sce = sce[, idx]
  sfs = sfs[idx]
  
  assign("genes.atlas", genes, envir = .GlobalEnv)
  assign("meta.atlas", meta, envir = .GlobalEnv)
  assign("sce.atlas", sce, envir = .GlobalEnv)
  assign("sfs.atlas", sfs, envir = .GlobalEnv)
  invisible(0)
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
  trend = scran::trendVar(sce, loess.args = list(span = 0.05), use.spikes = F)
  decomp = scran::decomposeVar(sce, fit = trend)
  decomp = decomp[decomp$mean > min.mean,]
  
  #exclude sex genes
  xist = "ENSMUSG00000086503"
  xist = "Xist"
  # mouse_ensembl = useMart("ensembl")
  # mouse_ensembl = useDataset("mmusculus_gene_ensembl", mart = mouse_ensembl)
  # gene_map = getBM(attributes=c("ensembl_gene_id", "chromosome_name"), filters = "ensembl_gene_id", values = rownames(decomp), mart = mouse_ensembl)
  # ychr = gene_map[gene_map[,2] == "Y", 1]
  #ychr = read.table("/nfs/research1/marioni/jonny/embryos/data/ygenes.tab", stringsAsFactors = FALSE)[,2]
  #ychr = read.table("/Users/alsu/Develop/FetalAlcoholSyndrome/data/ygenes.tab", stringsAsFactors = FALSE)[,1]
  
  other = c("tomato-td") #for the chimera
  decomp = decomp[!rownames(decomp) %in% c(xist, other),]
  
  decomp$FDR = p.adjust(decomp$p.value, method = "fdr")
  return(rownames(decomp)[decomp$p.value < 0.05])
}



denoisingCounts = function(counts, meta, neigh.net, n.neigh, avg_metric, mcparam){
  # subgraph neigh.net to the vertices we have
  vertices.neigh.net = as_ids( V(neigh.net) )
  vertices.neigh.net = intersect(vertices.neigh.net, colnames(counts))
  neigh.net = subgraph(neigh.net, vertices.neigh.net)
  adjacency.list = adjacent_vertices(neigh.net, V(neigh.net), mode = "all")
  
  denoised.counts.perCell = bplapply(colnames(counts), function(cell){
    current.CT = meta$celltype_mapped_denoised[meta$uniqueID == cell]
    current.adjacency.list = adjacency.list[names(adjacency.list) == cell]
    if (!is_empty(current.adjacency.list)){
      
      if (n.neigh == 1){
        neighbor.cells = as_ids(current.adjacency.list[[1]])
      } else if (n.neigh == 2){
        neighbor.cells = as_ids(current.adjacency.list[[1]])
        second.neighbor.cells = sapply(neighbor.cells, function(x){
          current.second.adjacency.list = adjacency.list[names(adjacency.list) == x]
          return(as_ids(current.second.adjacency.list[[1]]))
        })
        neighbor.cells = unique(c(neighbor.cells, unlist(second.neighbor.cells)))
      }
      neighbor.CT = sapply(neighbor.cells, function(x) meta$celltype_mapped_denoised[meta$uniqueID == x])
      cells.sameCT = c( cell , names(neighbor.CT)[neighbor.CT == current.CT])
      if (length(cells.sameCT) > 1){
        current.counts = counts[,cells.sameCT]
        if (avg_metric == "mean"){
          current.counts = apply(current.counts, 1, mean)
        } else if (avg_metric == "median"){
          current.counts = apply(current.counts, 1, median)
        } 
        return(current.counts)
      } else {
        return(counts[,cell])
      }
    } else {
      return(counts[,cell])
    }
  }, BPPARAM = mcparam)
  denoised.counts.perCell = do.call(cbind, denoised.counts.perCell)
  colnames(denoised.counts.perCell) = colnames(counts)
  return(denoised.counts.perCell)
}

binaryDenoisingCell = function(x, original.counts){
  if (sum(x) > length(x)/2){
    return(1)
  } else if (sum(x) < length(x)/2){
    return(0)
  } else {
    return(original.counts)
  }
}

denoisingBinaryCounts = function(counts, meta, neigh.net, n.neigh, mcparam){
  # subgraph neigh.net to the vertices we have
  vertices.neigh.net = as_ids( V(neigh.net) )
  vertices.neigh.net = intersect(vertices.neigh.net, colnames(counts))
  neigh.net = subgraph(neigh.net, vertices.neigh.net)
  adjacency.list = adjacent_vertices(neigh.net, V(neigh.net), mode = "all")
  
  denoised.counts.perCell = bplapply(colnames(counts), function(cell){
    current.CT = meta$celltype_mapped_denoised[meta$uniqueID == cell]
    current.adjacency.list = adjacency.list[names(adjacency.list) == cell]
    if (!is_empty(current.adjacency.list)){
      
      if (n.neigh == 1){
        neighbor.cells = as_ids(current.adjacency.list[[1]])
      } else if (n.neigh == 2){
        neighbor.cells = as_ids(current.adjacency.list[[1]])
        second.neighbor.cells = sapply(neighbor.cells, function(x){
          current.second.adjacency.list = adjacency.list[names(adjacency.list) == x]
          return(as_ids(current.second.adjacency.list[[1]]))
        })
        neighbor.cells = unique(c(neighbor.cells, unlist(second.neighbor.cells)))
      }
      neighbor.CT = sapply(neighbor.cells, function(x) meta$celltype_mapped_denoised[meta$uniqueID == x])
      cells.sameCT = c( cell , names(neighbor.CT)[neighbor.CT == current.CT])
      if (length(cells.sameCT) > 1){
        current.counts = counts[,cells.sameCT]
        original.counts = counts[,cell]
        current.denoised.counts = sapply(1:nrow(current.counts), function(current.gene){
          return(binaryDenoisingCell(current.counts[current.gene,], original.counts[current.gene]))
        })
        return(current.denoised.counts)
      } else {
        return(counts[,cell])
      }
    } else {
      return(counts[,cell])
    }
  }, BPPARAM = mcparam)
  denoised.counts.perCell = do.call(cbind, denoised.counts.perCell)
  colnames(denoised.counts.perCell) = colnames(counts)
  rownames(denoised.counts.perCell) = rownames(counts)
  return(denoised.counts.perCell)
}




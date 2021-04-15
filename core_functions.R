# load embryo 8.5 data
load_embryo_8.5 = function(dir = "local"){
  require(SingleCellExperiment)
  require(scater)
  
  # load sce for seqFISH
  if (dir == "cluster") {
    data.dir = "/nfs/research1/marioni/alsu/spatial/mouse_embryo/data/8_5/source/"
  } else {
    data.dir = "/Users/alsu/Develop/spatial/mouse_embryo/data/8_5/source/"
  }
  
  sce = readRDS( paste0( data.dir , "E8.5_sce_filt_unlabelled.Rds"))
  # add normalization by libsize 
  assay(sce, "cpm") <- logcounts(scater::logNormCounts(sce, size_factors = sce$total))
  assay(sce, "cpm_wo_xist") <- logcounts(scater::logNormCounts(sce, size_factors = as.numeric( sce$total - counts( sce["Xist"] )) ))
  
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


getSegmentationVerticesDF = function(DF,
                                     xname = "segmentation_vertices_x_global",
                                     yname = "segmentation_vertices_y_global",
                                     othercols = c("uniqueID","z")) {
  # DF is a DataFrame object
  # othercols is the others to keep
  
  long_x = unlist(DF[,xname])
  long_y = unlist(DF[,yname])
  
  if (length(long_x) != length(long_y)) stop("x and y need to be same length")
  
  long_xy = data.frame(
    long_x,
    long_y
  )
  colnames(long_xy) <- c(xname, yname)
  
  long_DF = cbind(
    rep(DF[,othercols], times = unlist(lapply(DF[,xname], length))),
    long_xy
  )
  
  return(as.data.frame(long_DF))
}

add_scalebar = function(dist_um = 250, x = 2.75, y = -3.20, ...) {
  # this is a quantity to add to an existing ggplot
  
  # dist_um is the distance for the scalebar in um, default 250um
  # x and y are the coordinates to place the scalebar
  # usage: to add to a ggplot object like "g + add_scalebar()"
  
  # useful optional argument is box.col = "white"
  
  # need to make sure that ggplot of interest doesn't have group aesthetic
  # defined in the ggplot, but in the geom_polygon itself
  
  require(ggsn)
  add = ggsn::scalebar(location = "bottomright",
                       dist = dist_um/227.74, dist_unit = "units", 
                       transform = FALSE,
                       x.min = x, x.max = x,
                       y.min = y, y.max = y,
                       height = 0.2, 
                       box.fill = "black",
                       st.size = 0,
                       inherit.aes = FALSE,
                       ...)
  return(add)
}

# function to rotate the embryos
rotateDF = function(DF, 
                    xname = "segmentation_vertices_x_global_affine", 
                    yname = "segmentation_vertices_y_global_affine", 
                    ang = 0) {
  # ang is a numeric vector named three values corresponding to embryo1, embryo2, and embryo3
  
  ang_long = ang[as.character(DF$embryo)]
  ang_rad = ang_long/180
  
  x = DF[,xname]
  y = DF[,yname]
  
  x_turn = x*cos(ang_rad) - y*sin(ang_rad)
  y_turn = x*sin(ang_rad) + y*cos(ang_rad)
  
  # reset the columns and then return the DF
  
  DF[,xname] <- x_turn
  DF[,yname] <- y_turn
  
  return(DF)
  
}
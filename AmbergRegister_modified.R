#' @importFrom Rvcg vcgClean vcgClost vcgUpdateNormals
#' @importFrom Morpho meshcube applyTransform computeTransform pcAlign
AmbergRegister_modified <- function(x, mesh2, lm1=NULL, lm2=NULL, k=1, lambda=1, iterations=15, rho=pi/2, dist=2, border=FALSE, smooth=TRUE, smoothit=1, smoothtype="t", tol=1e-10, useiter=TRUE, minclost=50, distinc=1,affine=NULL, rigid=NULL, similarity=NULL, tps=FALSE, pcAlign=FALSE,nn=20, silent=FALSE, useConstrained=TRUE, forceLM=FALSE,visualize=FALSE, folder=NULL,noinc=FALSE,
                                    bboxCrop=NULL,threads=0,iterations4pcAlign=150,  mc.cores4pcAlign=2)
{
  if (inherits(x, "mesh3d")) {
    mesh1 <- x
    Bayes <- NULL
  } else if (inherits(x, "BayesDeform"))
    Bayes <- x
  else
    stop("x must be an object of class mesh3d or BayesDeform")
  
  if (!is.null(Bayes)) {
    if (!requireNamespace("RvtkStatismo"))
      stop("for using the option Bayes, please install RvtkStatismo from https://github.com/zarquon42b/RvtkStatismo")
    mesh1 <- RvtkStatismo::DrawMean(Bayes$model)
  }
  #mesh1 <- rmUnrefVertex(mesh1, silent=TRUE)
  mesh1 <- vcgUpdateNormals(mesh1)
  mesh2 <- vcgUpdateNormals(mesh2)
  meshbord <- vcgBorder(mesh2)
  count <- 0
  if (iterations < 1)
    iterations <- 1e10
  if (length(lambda) == 1)
    lambda <- rep(lambda,iterations)
  else if (length(lambda) != iterations)
    stop("lambda must be vector of length 'iterations'")
  k <- round(k)# make sure k is integer - otherwise RAM overkill
  if (length(k) == 1)
    k <- rep(k,iterations)
  else if (length(k) != iterations)
    stop("k must be vector of length 'iterations'")
  
  affinemat <- NULL
  meshorig <- mesh1
  stopit <- FALSE
  hasLM <- FALSE
  if (!is.null(lm1) && !is.null(lm2)) {
    hasLM <- TRUE
    bary <- vcgClost(lm1,mesh1,barycentric = T)
  }
  if (!is.null(Bayes$initparams)) {
    mesh1 <- RvtkStatismo::DrawSample(Bayes$model,Bayes$initparams)
    if (hasLM)
      lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
    #wire3d(mesh1);spheres3d(lm1)
    #return(1)
  }
  if (hasLM) {## case: landmarks are provided
    
    if (!is.null(Bayes) && hasLM) {
      ##register landmarks on model and constrain reference
      lm2tmp <- rotonto(lm1,lm2,scale=Bayes$model@scale,reflection=FALSE)$yrot
      constMod <- RvtkStatismo::statismoConstrainModel(Bayes$model,lm2tmp,lm1,Bayes$ptValueNoise)
      if (useConstrained) {
        Bayes$model <- constMod
        mesh1 <- vcgUpdateNormals(RvtkStatismo::DrawMean(Bayes$model))
        lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
      }
    }                
    if (tps) {
      mesh1 <- tps3d(mesh1,lm1,lm2,threads=threads)
      lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
    } else {
      if (is.null(affine) && is.null(rigid) && is.null(similarity) && !pcAlign) {
        if (is.null(Bayes)) {                
          rigid <- list(iterations=0)
          if (!silent)
            cat("\n landmarks but no transform specified, performing rigid transform\n")
        } else if (Bayes$align) {
          rigid <- list(iterations=0)
          if (!silent)
            cat("\n landmarks but no transform specified, performing rigid transform\n")
        }
      }
      if (pcAlign) {
      mesh1 <- pcAlign(mesh1,mesh2,optim = TRUE,iterations = iterations4pcAlign,mc.cores = mc.cores4pcAlign)
        lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
      }
      if (!is.null(affine)) {##similarity matching
        if (!is.null(rigid) && is.null(similarity)) {
          affine$lm1 <- lm1
          affine$lm2 <- lm2
        }
        mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent,threads=threads)
        lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
      }
      if (!is.null(rigid)) { ##perform rigid icp-matching
        if (!pcAlign) {
          rigid$lm1 <- lm1
          rigid$lm2 <- lm2
        }
        mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent,threads=threads)
        if (hasLM)
          lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
      }
      if (!is.null(similarity)) {##similarity matching
        if (is.null(rigid)) {
          similarity$lm1 <- lm1
          similarity$lm2 <- lm2
        }
        mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent,threads=threads)
        lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
      }
    }
    
    affinemat <- computeTransform(vert2points(mesh1),vert2points(meshorig))
    tmp <- list()
    tmp$mesh <- mesh1
    if (!useiter && !forceLM)
      tmp$S <- createS(mesh1)
    
    if (forceLM && hasLM) {
      lm1 <- bary2point(bary$barycoords,bary$faceptr,mesh1)
      tmp <- AmbergDeformSpam(mesh1,lm1,lm2,k0=k[1],lambda=lambda[1])
      count <- count+1
      if (iterations == 1)
        stopit <- TRUE
    }
    verts0 <- vert2points(mesh1)
    
    
  } else if (pcAlign || !is.null(affine) || !is.null(similarity) || !is.null(rigid)) {
    if (pcAlign) {
      mesh1 <- pcAlign(mesh1,mesh2,optim = TRUE,iterations = iterations4pcAlign,mc.cores = mc.cores4pcAlign)
    }
    if (!is.null(affine)) {##similarity matching
      mesh1 <- rigSimAff(mesh1,mesh2,affine,type="a",silent = silent)
    }
    if (!is.null(similarity)) {##similarity matching
      mesh1 <- rigSimAff(mesh1,mesh2,similarity,type="s",silent = silent)
    }
    if (!is.null(rigid)){ ##perform rigid icp-matching
      mesh1 <- rigSimAff(mesh1,mesh2,rigid,type="r",silent = silent)
    }
    affinemat <- computeTransform(vert2points(mesh1),vert2points(meshorig))
    tmp <- list(mesh=mesh1)
    if (!useiter)
      tmp$S <- createS(mesh1)
    verts0 <- vert2points(mesh1)
    
  } else {   ## case: meshes are already aligned
    affinemat <- diag(4)
    tmp <- list()
    tmp$mesh <- mesh1
    if (!useiter)
      tmp$S <- createS(mesh1)
    verts0 <- vert2points(mesh1)
  }
  if (!is.null(bboxCrop)) {
    mesh2 <- cropOutsideBBox(mesh1,mesh2,extend=bboxCrop)
    if (!silent)
      cat("cropping target mesh\n")
  } 
  if (visualize) {
    rglid <- NULL
    if (!length(rgl.ids()$id)) 
      open3d(windowRect=c(80,80,800,800))
    else {
      rgl.bringtotop()
      rgl.clear()
    }
    points3d(meshcube(tmp$mesh),col="white",alpha=0)
    shade3d(mesh2,col=2,specular=1)
    if (!is.null(rglid))
      rgl.pop(id=rglid)
    rglid <- shade3d(tmp$mesh,col="white",front="lines", back="lines")
    
    if (!is.null(folder)) {
      if (substr(folder,start=nchar(folder),stop=nchar(folder)) != "/") 
        folder <- paste(folder,"/",sep="")
      dir.create(folder,showWarnings=F)
      movie <- paste(folder,"deformation",sep="")
      
      npics <- nchar(iterations+1)
      ndec <- paste0("%s%0",npics,"d.png")
    }
    #if (interactive())
    #readline("please select viewpoint\n")
    
    
    if (!is.null(folder)) {
      filename <- sprintf("%s%04d.png", movie, 1)
      rgl.snapshot(filename,fmt="png")
      movcount <- 2
    }
  }
  
  if (stopit) {
    ## set error and counter appropriately
    distance <- 1e12
    error <- 1e12
    count <- count+1
    while (count <= iterations && error > tol) {
      time0 <- Sys.time()
      if (useiter) {
        verts0 <- vert2points(tmp$mesh)
        mesh1 <- tmp$mesh
      }
      vert_old <- vert2points(tmp$mesh)
      clost <- vcgClostKD(tmp$mesh,mesh2,k=nn,threads=threads)
      verts1 <- vert2points(clost)
      nc <- normcheck(clost,tmp$mesh,threads = threads)                        
      
      ## find valid hits
      normgood <- as.logical(nc < rho)
      distgood <- as.logical(abs(clost$quality) <= dist)
      bordergood <- 1
      if (!border) 
        bordgood <- as.logical(!meshbord$borderit[clost$faceptr])
      #dupes <- !(as.logical(vcgClean(clost)$remvert))
      dupes <- TRUE
      good <- sort(which(as.logical(normgood*distgood*bordergood*dupes)))
      
      
      
      
      ### in case no good hit is found within the given distance we increase the distance by 1mm until valid references are found:
      increase <- distinc
      while (length(good) < minclost) {
        distgood <- as.logical(abs(clost$quality) <= (dist+increase))
        good <- sort(which(as.logical(normgood*distgood*bordergood)))
        increase <- increase+distinc
        cat(paste("distance increased to",dist+increase,"\n"))
      }
      
      
      ## update reference points
      lmtmp1 <- verts0[good,]
      lmtmp2 <- verts1[good,]
      ## map it according to new reference points
      #points3d(lmtmp2,col=count)
      if (useiter)
        tmp$S <- NULL
      
      tmpold <- tmp
      chk <- try(tmp <- AmbergDeformSpam(mesh1,lmtmp1,lmtmp2,k0=k[count],lambda=lambda[count],S=tmp$S),silent = TRUE)
      if (inherits(chk,"try-error")) {
        tmp <- tmpold
        cat(paste("iteration failed with:",chk,"previous iteration used\n"))
      }
      gc()
      if (smooth)
        tmp$mesh <- vcgSmooth(tmp$mesh,iteration = smoothit,type=smoothtype)
      ## calculate error
      if (!is.null(Bayes) && length(Bayes$sdmax) >= count) {
        if (!is.null(Bayes$wt)) {
          wt <- Bayes$wt[count]
          wts <- c(1,wt)
          wts <- wts/sum(wts)
          tmpmesh <- RvtkStatismo::PredictSample(Bayes$model,tmp$mesh,TRUE, sdmax=Bayes$sdmax[count],align=Bayes$align,mahaprob=Bayes$mahaprob)
          tmp$mesh$vb[1:3,] <- wts[1]*tmp$mesh$vb[1:3,]+wts[2]*tmpmesh$vb[1:3,]
        } else {
          tmp$mesh <- RvtkStatismo::PredictSample(Bayes$model,tmp$mesh,TRUE, sdmax=Bayes$sdmax[count],align=Bayes$align,mahaprob=Bayes$mahaprob)
        }
        
      }
      distance_old <- distance
      distance <- mean(vcgClostKD(mesh2,tmp$mesh,k0=10,sign=F,threads=threads)$quality)
      if (distance > distance_old && !is.null(Bayes) && noinc) {
        cat("\n=========================================\n")
        message(paste(" Info: Distance is increasing, matching stopped after ",count,"iterations\n"))
        count <- 1e10
        tmp <- tmpold
      }
      tmp$mesh <- vcgUpdateNormals(tmp$mesh)
      if (visualize) {
        
        if (!is.null(rglid))
          rgl.pop(id=rglid)
        rglid <- shade3d(tmp$mesh,col="white",front="lines",back="lines")
        if (!is.null(folder)) {
          filename <- sprintf("%s%04d.png", movie, movcount)
          movcount <- movcount+1
          rgl.snapshot(filename,fmt="png")
        }
      }
      error <- sum((vert2points(tmp$mesh)-vert_old)^2)/nrow(vert_old)
      
      time1 <- Sys.time()
      if (!silent && count < 1e10) {
        cat(paste("-> finished iteration",count,"in",round(as.numeric(time1-time0,unit="secs"),2), "seconds\n"))
        cat(paste(" Info: MSE between iterations:",error,"\n"))
        cat(paste(" Info: Average distance to target:",distance,"\n"))
        if (error < tol)
          cat(paste("***\n==> Convergence threshold reached after",count,"iterations\n"))
      }
      count <- count+1
    }
  }
  lm1map <- NULL
  if (!is.null(lm1))
    lm1map <- lm1 <- bary2point(bary$barycoords,bary$faceptr,tmp$mesh)
  return(list(mesh=tmp$mesh,affine=affinemat,lm1=lm1map))
}


rigSimAff <- function(mesh1,mesh2,args,type="r",silent=TRUE,threads=1) {
  iterations <- args$iterations; if (is.null(iterations)) iterations <- 3
  lm1=args$lm1
  lm2 <- args$lm2
  uprange <- args$uprange; if (is.null(uprange)) uprange <- 0.9
  maxdist <- args$maxdist
  minclost <- args$minclost; if (is.null(minclost)) minclost <- 50
  distinc <- args$distinc;
  rhotol <- args$rhotol; if (is.null(rhotol)) rhotol <- pi
  k <- args$k; if (is.null(k)) k <- 50
  reflection <- args$reflection;  if (is.null(reflection)) reflection <- FALSE
  subsample <- args$subsample
  out <- icp(mesh1, mesh2, iterations=iterations, lm1=lm1, lm2=lm2, uprange=uprange ,maxdist=maxdist, minclost=minclost, distinc=distinc, rhotol=rhotol, k=k, reflection=reflection, silent = silent,subsample=subsample, type=type,threads=threads)
  return(out)
}

#============================================================================
# you need to have a folder with .nii files of CT/µCT scans for your sample!
#============================================================================

# steps for mesh registration
#===============================
# load libraries
# note that many libraries need compilation and may take some time to install
library(rgl)
library(Rvcg)
library(shapes)
library(abind)
library(Morpho)
library(geomorph)
library(parallel)
library(doParallel)
library(mesheR)
library(Arothron)
library(ctrlR)
library(deformetrics)
library(devtools)
library(RNifti)
library(ANTsR)
library(RANTs)
mc.cores<-parallel::detectCores() # detect how many core the computer has
#============================================================================
# set your working directory with setwd()
# provide a path to .nii files (eg mouse crania)
meshes_list=list.files("crania_nii",full.names = T)
cbind(meshes_list)
# have a look at one file
antsImageRead(meshes_list[1]) # ANTsR ants uses mm scale
antsGetSpacing(antsImageRead(meshes_list[1])) # voxel size in mm
readNifti(meshes_list[1]) # to check voxel size, do antsGetSpacing(refmeshim)*1000 to get back to µm.
antsGetSpacing(antsImageRead(meshes_list[1]))*1000
# costumized function to create a 3d mesh surface, super useful! From RANTs package.
# function to create a 3d mesh surface from RANTs package
isosurfANTsR<-function(x,threshold,as.int=TRUE,scale2microns=1000,...) {
    segsurf <- vcgIsosurface(as.array(x),threshold=threshold,spacing = antsGetSpacing(x)*scale2microns,
    origin = antsGetOrigin(x),direction = antsGetDirection(x),
    as.int = as.int,IJK2RAS = diag(c(1, 1, 1, 1)))
    return(segsurf)
} # set scale2microns=1 for mm, µm otherwise
# have a look at one specimen
# 0.051665 is vozel size in my case. resample by 3 and dilate by 1, threshold between 20-60 for bone, depending on CT.
# check
?iMath
?resampleImage
shade3d(vcgSmooth(mesh = isosurfANTsR(iMath(resampleImage(antsImageRead(meshes_list[1]),0.051665*3,0,1), "GD",1),threshold = 40),lambda = .6),col="ivory")
#
rgl.clear()
##============================================================================
### make a list of all 3d meshes
system.time({ target_meshes<-mclapply(seq(length(meshes_list)), function(i) 
	vcgSmooth(mesh = isosurfANTsR(iMath(resampleImage(antsImageRead(meshes_list[1]),0.051665*3,0,1), "GD",1),threshold = 40),lambda = .6), mc.cores = mc.cores-1) }) # also tells you the time it took..
# check
length(target_meshes)
shade3d(target_meshes[[1]],col="white") # specimen #1
rgl.clear();rgl.close()
###============================================================================
# decimate meshes to low resolution
###============================================================================
# how many vertices?
dim(target_meshes[[1]]$vb)[2]
wire3d(target_meshes[[30]],col="white")
shade3d(target_meshes[[30]],col="red")
replicate(length(rgl.dev.list()),rgl.close())
# decimate the targets to get lowres. meshes regulation the percent param.
?vcgQEdecim
system.time({ target_meshes_lowres<-mclapply(seq(length(meshes_list)), function(i) rmUnrefVertex(vcgIsolated(vcgQEdecim(target_meshes[[i]],percent = .4))),mc.cores = mc.cores-1) })
# check
# how many vertices? ideally between 15-20.000
dim(target_meshes_lowres[[1]]$vb)[2]
wire3d(target_meshes_lowres[[1]],col="white")
shade3d(target_meshes_lowres[[1]],col="red",alpha=.7)
replicate(length(rgl.dev.list()),rgl.close())
###============================================================================
# low resolution meshes uniformely remeshed
# regulate the voxelSize param
?vcgUniformRemesh
target_meshes_lowres_remeshed<-mclapply(seq(length(meshes_list)), function(i) rmUnrefVertex(vcgIsolated(vcgUniformRemesh(x = rmUnrefVertex(vcgIsolated(target_meshes_lowres[[i]])),voxelSize = 180,multiSample = T))),mc.cores = mc.cores)
dim(target_meshes_lowres_remeshed[[1]]$vb)[2] # number of vertices
# load the landmarks that belong to one arbitrary specimen, eg specimen #1
# there are different ways to import your landmarks, it depends in what format they are.
# eg. in .csv format
landmarks<-as.matrix.data.frame(read.csv("your_landmarks.csv")) # just for 1 specimen
# check them
spheres3d(landmarks,radius = .3) # play with radius, it depends on the scale of your landmarks, eg for µm the radius is larger then 2
# project these lanmarks onto the surface of the remeshed specimen #1 using KD-Tree Search
landmarks<-vert2points(vcgClostKD(landmarks,vcgSmooth(target_meshes_lowres_remeshed[[1]])))
# prepare refmesh which is our specimen #1
dim(rmUnrefVertex(vcgSmooth(target_meshes_lowres_remeshed[[1]]))$vb) # how many vertices initially
refmesh1<-target_meshes_lowres_remeshed[[1]] # assing a mesh
refmesh1$vb<-cbind(refmesh1$vb,t(cbind(landmarks,1.000))) # add or hide your landmarks at the end of the mesh vertices list
dim(refmesh1$vb)
# have a look
shade3d(vcgSmooth(refmesh1),col="ivory")
# plot lms
spheres3d(landmarks,radius = .3)
# close all open rgl windonws
replicate(length(rgl.dev.list()),rgl.close())

# prepare your onw indices of right, left and saggital landmarks
# such that they would start from dim(refmesh1$vb)+1
dim(refmesh1$vb)+1 # this is the first lm from your landmark matrix
dim(refmesh1$vb) # is the last vertex that belongs to the mesh
# now define what belongs to what
sagittal_lms_inds_refmesh1
surlmsL_refmesh1
surlmsR_refmesh1


#============================================================================
start_time_icpAmbergregistGaussmatch_refmesh1<-Sys.time()
# load the modified gaussMatch function (no hidden landmarks removal and larger rgl window, and added landmarks)
source("AmbergRegister_modified.R")
source("gaussMatch_modified.R")
# check original functions
?AmbergRegister
?gaussMatch
affine4Ambergregist<-list(iterations=50,pcAlign = T,type="affine") # params for affine registration, play with pcAlign=F and more iterations
# align by principal axes and then fine-tine with affine registration
icpAmbergRegist_refmesh1<-lapply(seq(length(meshes_list)), function(i) 
	AmbergRegister_modified(vcgSmooth(refmesh1),
	vcgSmooth(target_meshes_lowres[[i]]),iterations = 20,pcAlign = T,affine = affine4Ambergregist,
	mc.cores4pcAlign = mc.cores-1))
# check how it looks before alignement
shade3d(vcgSmooth(refmesh1),col="white",front="lines",back="lines")
shade3d(vcgSmooth(target_meshes_lowres[[1]]),col="red")
spheres3d(landmarks,radius = .3)
###
# check after principal axes aligment and affine registration
open3d()
shade3d(icpAmbergRegist_refmesh1[[1]]$mesh,col="white",front="lines",back="lines")
shade3d(vcgSmooth(target_meshes_lowres[[1]]),col="red")
# close all open rgl windonws
replicate(length(rgl.dev.list()),rgl.close())
# now use the output from affine registration for further # elastic matching using smooth displacement fields
gaussMatch_refmesh1<-lapply(seq(length(meshes_list)), function(i) 
	gaussMatch_modified(icpAmbergRegist_refmesh1[[i]]$mesh,vcgSmooth(target_meshes_lowres[[i]]), iterations = 4, 
	smoothit = 20, smooth = 10, sigma = 20, pro = "kd", gamma = 2, nh = 5, threads = 0, displacementsmooth = "Gauss", 
	visualize=T, alpha = 1,radius4spheres3d = .3, mesh1lms_indicesR = surlmsR_refmesh1, mesh1lms_indicesL = 
	surlmsL_refmesh1,mesh1lms_indicesSagittal = sagittal_lms_inds_refmesh1))
# close all windows
replicate(length(rgl.dev.list()),rgl.close())
end_time_icpAmbergregistGaussmatch_refmesh1<-Sys.time()
time_taken4icpAmbergregistGaussmatch_refmesh1<-end_time_icpAmbergregistGaussmatch_refmesh1 - start_time_icpAmbergregistGaussmatch_refmesh1
time_taken4icpAmbergregistGaussmatch_refmesh1
#============================================================================
# check the mesh distances
meshDist(vcgSmooth(target_meshes_lowres[[1]]),gaussMatch_refmesh1[[1]],rampcolors = c("blue","ivory","red"),tolcol = "ivory",alpha=1)

# project landmarks on the targets and returm them as a list
registered_lms<-mclapply(seq(length(complete_specs_names)), function(i) 
	vert2points(vcgClostKD(vert2points(gaussMatch_refmesh1_mean[[i]])[c(c(surlmsR_refmesh1),c(surlmsL_refmesh1),sagittal_lms_inds_refmesh1),], 
	vcgSmooth(target_meshes_lowres[[i]]))),mc.cores = mc.cores-1)
# convert list to array
registered_lms_arr<-list2array(registered_lms)
# now use your registered array of landmarks for morphometric analysis!
#eg
?procSym



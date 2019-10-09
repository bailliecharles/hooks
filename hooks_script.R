########################################
####          HOOKS SCRIPT          ####
####        Charles Baillie         ####
####              2019              ####
####   bailliecharles@gmail.com     ####
########################################

# This script was written to quantify and analyse shape variation of cymothoid dactyli 
# between species exhibiting different ecological traits (attachment location).
# The script uses the package geomorph, landmarks were plotted using tpsDig2.


###-------------------------------
###        1. Preliminaries
###-------------------------------

# 1. set up environment 
  rm(list=ls())

  setwd("~/Dropbox/projects/paper_geomorph")
  ### you don't need all of these packages to run this script!!!
  library(geomorph)
  library(Rphylopars)
  library(plyr)
  library(abind)
  library(ggplot2)
  library(grid)
  library(reshape2)
  library(Cairo)
  library(jtools)
  library(MASS)
  library(RRPP)
  library(ggtree)
  library(geiger)
  library(caper)
  library(nlme)
  library(phangorn)
  library(phytools)
  library(scales)
  library(psych)
  library(StereoMorph)
  #library(phylobase)

  source("scripts/hooks_funcs.R")
  source("scripts/btShapes.R")


# 2.  Read in data and check it out
  land <-  readland.tps("hooks.tps", specID="ID", readcurves=T, warnmsg=T)

  
# 3. Sort names out a bit
  nam.a <- dimnames(land)[[3]]
  nam.b <- sub(".*\\\\", "", nam.a) # regex to select everything up to end of last double slash (if whole path included in name)
  #y <- gsub(" ", "_", x, fixed=T) # replace all spaces with an underscore 
  nam.c <- gsub(" (slightly broken)", "", nam.b,fixed=T)
  nam.d <- gsub('.JPG', '', nam.c)
  nam.e <- gsub('.jpg', '', nam.d)
  unique(nam.e) # check 
  dimnames(land)[[3]] <- nam.e

# 4. define semilandmarks and their bounds for the procrustes
  #define.sliders(a_data_1, 39)  #plot the point: remember the slider needs to be point in the middle!!! Will write .csv
  curves <-  as.matrix(read.csv("hooks_curves.csv", header=T))


# 5. read tree
  phy = read.tree("trees/hooks.tre") 
  
  phy=drop.tip(phy,c(16)) # remove outgroup
  phy = root(phy, "Elthusa_raynaudii")
  
  phy$tip.label[5]<- "Nerocila_depressa"

  phy_chron = chronos(phy, model="correlated",lambda=0)
  plot(phy_chron)
  # make phy_chron into class phylo
  phy_chron = list(edge=phy_chron$edge,edge.length=phy_chron$edge.length, Nnode=phy_chron$Nnode, node.label=phy_chron$node.label, tip.label=phy$tip.label)
  mostattributes(phy_chron)<-list(order= "cladewise",  names = list("edge","edge.length","Nnode","node.label","tip.label"))
  class(phy_chron)<- "phylo"
  str(phy_chron)
  phy_chron$node.label=c(1:phy_chron$Nnode)
  class(phy_chron)

# To find optimal value for smoothing parameter in chronos():
# l= c(0:500)
# LL.out=NULL
# for (i in 1:length(l)) {
#   LL.out[i]=attributes(chronos(phy, lambda=l[i]))$ploglik
# }
# l[LL.out==max(LL.out)]
# plot(log(l+1), LL.out, type="l", lwd=3, col="grey", xlab=expression(paste("Log(",lambda,")")), ylab="log likelihood")





###----------------------------
###   2. GPA, Sample Sizes etc. 
###----------------------------

#1. split into p1 and p7 
  a_P7 <-  land[,,-grep("P1", dimnames(land)[[3]])]
  a_P1 <-  land[,,-grep("P7", dimnames(land)[[3]])]

#2.Make a factor as a grouping variable (by parasitic mode)
  a_P1_mode <- as.factor(substr(dimnames(a_P1)[[3]],1,4))
  levels(a_P1_mode) <-  c("E","M", "M","M", "G","G","E", "E","G")   

  a_P7_mode <- as.factor(substr(dimnames(a_P7)[[3]],1,4))
  levels(a_P7_mode) <-  c("E","M", "M","M", "G","G", "E", "E", "G")
  
#3. perform superimposition of data using generalized procrustes analysis
  a_P1_gpa <-  gpagen(a_P1, curves = curves)
  a_P7_gpa <-  gpagen(a_P7, curves = curves)

#4. check for outliers  
  plotOutliers(a_P1_gpa$coords, groups=a_P1_mode) 
  plotOutliers(a_P7_gpa$coords, groups=a_P7_mode) 

  plot(a_P1_gpa) #check consensus shape
  plot(a_P7_gpa)
  
#5. and make geomorph style DFs for laters
a_P1_gdf <-  geomorph.data.frame(shape = a_P1_gpa$coords, mode = a_P1_mode, CS=log(a_P1_gpa$Csize) ) # geomorph way of making a DF
a_P7_gdf <-  geomorph.data.frame(shape = a_P7_gpa$coords, mode = a_P7_mode, CS=log(a_P7_gpa$Csize) )


#6. Get sample sizes
sum(colSums(table(dimnames(a_P1_gdf$shape)[[3]], a_P1_gdf$mode))) # sample sizes
sum(colSums(table(dimnames(a_P7_gdf$shape)[[3]], a_P7_gdf$mode))) # sample sizes

samples.P1 = plyr::count(as.factor(gsub('.{7}$', '', dimnames(a_P1_gdf$shape)[[3]])))
samples.P7 =plyr::count(as.factor(gsub('.{7}$', '', dimnames(a_P7_gdf$shape)[[3]])))



###----------------------------
###   3. Mean Shapes
###----------------------------

### P1
P1.MEAN = a_P1

dimnames(P1.MEAN)[[3]] = gsub("(\\ .*?)\\ (.*)", "\\1",dimnames(P1.MEAN)[[3]])
dimnames(P1.MEAN)[[3]] = gsub(" ", "_", dimnames(P1.MEAN)[[3]])

P1.MEAN.names = as.factor(dimnames(P1.MEAN)[[3]])

P1.MEAN.RAW_listBYspp = list() 
for(i in 1:nlevels(P1.MEAN.names)){
  P1.MEAN.RAW_listBYspp[[i]] = P1.MEAN[,,which(dimnames(P1.MEAN)[[3]]==levels(P1.MEAN.names)[i])]
}

names(P1.MEAN.RAW_listBYspp) = levels(P1.MEAN.names)
P1.MEAN.FIT_listBYspp = lapply(P1.MEAN.RAW_listBYspp, function(x) gpagen(x, curves=curves))


# mshape coords
P1.MEAN.MSHAPE_listBYspp = list()
for (i in 1:length(P1.MEAN.FIT_listBYspp)){
  P1.MEAN.MSHAPE_listBYspp[[i]] = P1.MEAN.FIT_listBYspp[[i]]$consensus
}
names(P1.MEAN.MSHAPE_listBYspp) = names(P1.MEAN.RAW_listBYspp)
P1.MEAN.MSHAPE = abind(P1.MEAN.MSHAPE_listBYspp, along=3)
P1.MEAN.MSHAPE = P1.MEAN.MSHAPE[,,order(dimnames(P1.MEAN.MSHAPE)[[3]])]


#### If you need to drop spp not in tree do it here
P1.MEAN.MSHAPE_phy = P1.MEAN.MSHAPE[,,order(match(dimnames(P1.MEAN.MSHAPE)[[3]],phy$tip.label))]
P1.MEAN.MSHAPE_phy = P1.MEAN.MSHAPE_phy[,,-c(16:18)] 
dimnames(P1.MEAN.MSHAPE_phy)[[3]]

# centroid sizes
P1.MEAN.CS=list()
for (i in 1:length(P1.MEAN.FIT_listBYspp)){
  P1.MEAN.CS[[i]] = mean(P1.MEAN.FIT_listBYspp[[i]]$Csize)
}
names(P1.MEAN.CS) = names(P1.MEAN.RAW_listBYspp)
P1.MEAN.CS = unlist(P1.MEAN.CS)
P1.MEAN.CS = P1.MEAN.CS[order(names(P1.MEAN.CS))]


#### If you need to drop spp not in tree do it here
P1.MEAN.CS_phy = P1.MEAN.CS[order(match(names(P1.MEAN.CS),phy$tip.label))]
P1.MEAN.CS_phy = P1.MEAN.CS_phy[-c(16:18)]

#second procrustes OLS
P1.MEAN.FIT2=gpagen(P1.MEAN.MSHAPE, curves=curves) # procrustes of the mean shapes for each spp.
P1.MEAN.coords= P1.MEAN.FIT2$coords
P1.MEAN.coords= P1.MEAN.coords[,,order(dimnames(P1.MEAN.coords)[[3]])]
P1.MEAN.MODE <- as.factor(substr(dimnames(P1.MEAN.coords)[[3]],1,4))
levels(P1.MEAN.MODE) <-  c("E","M", "M","M", "G","G", "E","E","G")
names(P1.MEAN.MODE) = dimnames(P1.MEAN.coords)[[3]]

#second procrustes PGLS
P1.MEAN.FIT2_phy=gpagen(P1.MEAN.MSHAPE_phy, curves=curves) # procrustes of the mean shapes for each spp.
P1.MEAN.coords_phy= P1.MEAN.FIT2_phy$coords
#P1.MEAN.coords_phy= P1.MEAN.coords[,,order(match(dimnames(P1.MEAN.coords_phy)[[3]], phy_chron$tip.label))]
P1.MEAN.MODE_phy <- as.factor(substr(dimnames(P1.MEAN.coords_phy)[[3]],1,4))
levels(P1.MEAN.MODE_phy) <-  c("E","M", "M","M", "G","G", "E","G")
names(P1.MEAN.MODE_phy) = dimnames(P1.MEAN.coords_phy)[[3]]


P1.MEAN_gdf = geomorph.data.frame(shape=P1.MEAN.coords, CS=as.vector(log(P1.MEAN.CS)), mode=P1.MEAN.MODE)
P1.MEAN_gdf_phy = geomorph.data.frame(shape=P1.MEAN.coords_phy, CS=as.vector(log(P1.MEAN.CS_phy)), mode=P1.MEAN.MODE_phy)



### P7
P7.MEAN = a_P7

dimnames(P7.MEAN)[[3]] = gsub("(\\ .*?)\\ (.*)", "\\1",dimnames(P7.MEAN)[[3]])
dimnames(P7.MEAN)[[3]] = gsub(" ", "_", dimnames(P7.MEAN)[[3]])

P7.MEAN.names = as.factor(dimnames(P7.MEAN)[[3]])

P7.MEAN.RAW_listBYspp = list() 
for(i in 1:nlevels(P7.MEAN.names)){
  P7.MEAN.RAW_listBYspp[[i]] = P7.MEAN[,,which(dimnames(P7.MEAN)[[3]]==levels(P7.MEAN.names)[i])]
}

names(P7.MEAN.RAW_listBYspp) = levels(P7.MEAN.names)
P7.MEAN.FIT_listBYspp = lapply(P7.MEAN.RAW_listBYspp, function(x) gpagen(x, curves=curves))


# mshape coords
P7.MEAN.MSHAPE_listBYspp = list()
for (i in 1:length(P7.MEAN.FIT_listBYspp)){
  P7.MEAN.MSHAPE_listBYspp[[i]] = P7.MEAN.FIT_listBYspp[[i]]$consensus
}
names(P7.MEAN.MSHAPE_listBYspp) = names(P7.MEAN.RAW_listBYspp)
P7.MEAN.MSHAPE = abind(P7.MEAN.MSHAPE_listBYspp, along=3)
P7.MEAN.MSHAPE = P7.MEAN.MSHAPE[,,order(dimnames(P7.MEAN.MSHAPE)[[3]])]


#### If you need to drop spp not in tree do it here
P7.MEAN.MSHAPE_phy = P7.MEAN.MSHAPE[,,order(match(dimnames(P7.MEAN.MSHAPE)[[3]],phy$tip.label))]
P7.MEAN.MSHAPE_phy = P7.MEAN.MSHAPE_phy[,,-c(16:19)] 
dimnames(P7.MEAN.MSHAPE_phy)[[3]]

# centroid sizes
P7.MEAN.CS=list()
for (i in 1:length(P7.MEAN.FIT_listBYspp)){
  P7.MEAN.CS[[i]] = mean(P7.MEAN.FIT_listBYspp[[i]]$Csize)
}
names(P7.MEAN.CS) = names(P7.MEAN.FIT_listBYspp)
P7.MEAN.CS = unlist(P7.MEAN.CS)
P7.MEAN.CS = P7.MEAN.CS[order(names(P7.MEAN.CS))]

#### If you need to drop spp not in tree do it here
P7.MEAN.CS_phy = P7.MEAN.CS[order(match(names(P7.MEAN.CS),phy$tip.label))]
P7.MEAN.CS_phy = P7.MEAN.CS_phy[-c(16:19)]

#second procrustes (OLS)
P7.MEAN.FIT2=gpagen(P7.MEAN.MSHAPE, curves=curves) # procrustes of the mean shapes for each spp.
P7.MEAN.coords= P7.MEAN.FIT2$coords
P7.MEAN.coords= P7.MEAN.coords[,,order(dimnames(P7.MEAN.coords)[[3]])]
P7.MEAN.MODE <- as.factor(substr(dimnames(P7.MEAN.coords)[[3]],1,4))
levels(P7.MEAN.MODE) <-  c("E","M", "M","M", "G","G", "E","E","G")
names(P7.MEAN.MODE) = dimnames(P7.MEAN.coords)[[3]]

#second procrustes (PGLS)
P7.MEAN.FIT2_phy=gpagen(P7.MEAN.MSHAPE_phy, curves=curves) # procrustes of the mean shapes for each spp.
P7.MEAN.coords_phy= P7.MEAN.FIT2_phy$coords
#P7.MEAN.coords_phy= P7.MEAN.coords[,,order(match(dimnames(P7.MEAN.coords_phy)[[3]], phy$tip.label))]
P7.MEAN.MODE_phy <- as.factor(substr(dimnames(P7.MEAN.coords_phy)[[3]],1,4))
levels(P7.MEAN.MODE_phy) <-  c("E","M", "M","M", "G","G", "E","G")
names(P7.MEAN.MODE_phy) = dimnames(P7.MEAN.coords_phy)[[3]]


P7.MEAN_gdf = geomorph.data.frame(shape=P7.MEAN.coords, CS=as.vector(log(P7.MEAN.CS)), mode=P7.MEAN.MODE)
P7.MEAN_gdf_phy = geomorph.data.frame(shape=P7.MEAN.coords_phy, CS=as.vector(log(P7.MEAN.CS_phy)), mode=P7.MEAN.MODE_phy)





###----------------------------------
###           4. Modelling
###----------------------------------


# 1.  Make some better DFs from PCA ###

### P1 OLS DF
P1_2d = two.d.array(P1.MEAN.coords) #shapes to 2d array

# Project into tangent space
P1_pca = prcomp(P1_2d) #PCA
summary(P1_pca) #48.7 28.6 (77.3) 17.7
P1_pca_multi = P1_pca$x*10000   # multiply to avoid IEEE error

colnames(P1_pca_multi ) <- sub("PC", "PCm", colnames(P1_pca_multi)) # rename mulitplied cols

P1_DF = data.frame(P1_2d, P1_pca$x[,1:12], P1_pca_multi[,1:12], log(P1.MEAN.CS), P1.MEAN.MODE) # make into df and change column names
colnames(P1_DF)[seq(from=2,84,2)] <- paste0("Y", 1:42)
colnames(P1_DF)[seq(from=1,83,2)] <- paste0("X", 1:42)
colnames(P1_DF)[ncol(P1_DF)-1] = "CS"
colnames(P1_DF)[ncol(P1_DF)] = "mode"

# check variance explained by the retained PCAs
sum(diag(P1.MEAN.FIT2$points.VCV)) #trace of covariance matrix 
sum(P1_pca$sdev[1:12]^2) # eigenvalues 
(sum(P1_pca$sdev[1:12]^2)/sum(diag(P1.MEAN.FIT2$points.VCV)))*100 # 99.9% of data is described by these PCs

####
#From Bolzan et al. CHECK pairwise diffs of allometric trajectories
p1Epca = P1_pca$x[c(1:3,16,17),]
p1Gpca = P1_pca$x[c(11:15,18),]
p1Mpca = P1_pca$x[c(4:10),]
sqrt(crossprod(c(p1Epca[,1],p1Gpca[,1]))) #implies no diff between modes
sqrt(crossprod(c(p1Epca[,1],p1Mpca[,1]))) #implies no diff between modes
sqrt(crossprod(c(p1Gpca[,1],p1Mpca[,1]))) #implies no diff between modes
####

### P7 OLS DF
P7_2d = two.d.array(P7.MEAN.coords) #shapes to 2d array

# Project into tangent space
P7_pca = prcomp(P7_2d) 
summary(P7_pca) #
P7_pca_multi = P7_pca$x*10000   # multiply to avoid IEEE error

colnames(P7_pca_multi) <- sub("PC", "PCm", colnames(P7_pca_multi)) # rename mulitplied cols

P7_DF = data.frame(P7_2d, P7_pca$x[,1:12], P7_pca_multi[,1:12], log(P7.MEAN.CS), P7.MEAN.MODE) # make into df and change column names
colnames(P7_DF)[seq(from=2,84,2)] <- paste0("Y", 1:42)
colnames(P7_DF)[seq(from=1,83,2)] <- paste0("X", 1:42)
colnames(P7_DF)[ncol(P7_DF)-1] = "CS"
colnames(P7_DF)[ncol(P7_DF)] = "mode"

# check variance explained by the retained PCAs
sum(diag(P7.MEAN.FIT2$points.VCV)) #trace of covariance matrix 
sum(P7_pca$sdev[1:8]^2) # eigenvalues 
(sum(P7_pca$sdev[1:12]^2)/sum(diag(P7.MEAN.FIT2$points.VCV)))*100 # 99.9% of data is described by these PCs

####
#From Bolzan et al. CHECK pairwise diffs of allometric trajectories
P7Epca = P7_pca$x[c(1:3,16,17,18),]
P7Gpca = P7_pca$x[c(11:15,19),]
P7Mpca = P7_pca$x[c(4:10),]
sqrt(crossprod(c(P7Epca[,1],P7Gpca[,1]))) #implies no diff between modes 
sqrt(crossprod(c(P7Epca[,1],P7Mpca[,1]))) #implies no diff between modes
sqrt(crossprod(c(P7Gpca[,1],P7Mpca[,1]))) #implies no diff between modes


### P1 PGLS DF
P1_pca_phy = phyl.pca(two.d.array(P1.MEAN.coords_phy), tree=phy_chron) #PCA
summary(P1_pca_phy) #71.9 16.4 (88.3) 9.0
P1_pca_multi_phy = P1_pca_phy$S*10000   # multiply to avoid IEEE error

colnames(P1_pca_multi_phy ) <- sub("PC", "PCm", colnames(P1_pca_multi_phy)) # rename mulitplied cols

P1_DF_phy = data.frame(row.names(P1_pca_phy$S),P1_pca_phy$S, P1_pca_multi_phy[,1:12], log(P1.MEAN.CS_phy), P1.MEAN.MODE_phy) # make into df and change column names
colnames(P1_DF_phy)[ncol(P1_DF_phy)-1] = "CS"
colnames(P1_DF_phy)[ncol(P1_DF_phy)] = "mode"
colnames(P1_DF_phy)[1] ="species"
colnames(P1_DF_phy)[16:27]=LETTERS[1:12]

P1_DF_phy_gdf = geomorph.data.frame(shape=as.matrix(P1_DF_phy[,16:27]), CS=log(P1.MEAN.CS_phy), pmode=P1.MEAN.MODE_phy)
row.names(P1_DF_phy_gdf$shape) = row.names(P1_DF_phy)

### P7 PGLS DF
P7_pca_phy = phyl.pca(two.d.array(P7.MEAN.coords_phy), tree=phy_chron) #PCA
summary(P7_pca_phy) #71.9 16.4 (88.3) 9.0
P7_pca_multi_phy = P7_pca_phy$S*10000   # multiply to avoid IEEE error

colnames(P7_pca_multi_phy ) <- sub("PC", "PCm", colnames(P7_pca_multi_phy)) # rename mulitplied cols

P7_DF_phy = data.frame(row.names(P7_pca_phy$S),P7_pca_phy$S, P7_pca_multi_phy[,1:12], log(P7.MEAN.CS_phy), P7.MEAN.MODE_phy) # make into df and change column names
colnames(P7_DF_phy)[ncol(P7_DF_phy)-1] = "CS"
colnames(P7_DF_phy)[ncol(P7_DF_phy)] = "mode"
colnames(P7_DF_phy)[1] ="species"
colnames(P7_DF_phy)[16:27]=LETTERS[1:12]

P7_DF_phy_gdf = geomorph.data.frame(shape=as.matrix(P7_DF_phy[,16:27]), CS=log(P7.MEAN.CS_phy), pmode=P7.MEAN.MODE_phy)
row.names(P7_DF_phy_gdf$shape) = row.names(P7_DF_phy)


# Define OLS lhs formula
fmla.LHS <-paste(colnames(P1_DF[97:108]),collapse="  , ")
fmla.full <- as.formula(paste("cbind(", fmla.LHS, ")" , "~","CS*mode")) 
fmla.size.plus.mode <- as.formula(paste("cbind(", fmla.LHS, ")" , "~","CS+mode")) 
fmla.size <- as.formula(paste("cbind(", fmla.LHS, ")" , "~","CS")) 
fmla.mode <- as.formula(paste("cbind(", fmla.LHS, ")" , "~","mode")) 
fmla.int <- as.formula(paste("cbind(", fmla.LHS, ")" , "~","1")) 

fmla.LHS.phy <-paste(LETTERS[1:12],collapse="  , ")
fmla.full.phy  <- as.formula(paste("cbind(", fmla.LHS.phy, ")" , "~","CS*mode")) 
fmla.size.plus.mode.phy  <- as.formula(paste("cbind(", fmla.LHS.phy, ")" , "~","CS+mode")) 
fmla.size.phy  <- as.formula(paste("cbind(", fmla.LHS.phy, ")" , "~","CS")) 
fmla.mode.phy  <- as.formula(paste("cbind(", fmla.LHS.phy, ")" , "~","mode")) 
fmla.int.phy  <- as.formula(paste("cbind(", fmla.LHS.phy, ")" , "~","1")) 

# 2. 
#################
# Ordinary Least Squares (OLS)
# 1. Full model
# 2. Group Means Type I SS
# 3. Group Slopes Type I SS
# 4. Allo-free Group Means Type I SS

# P1
P1_mod_full = procD.lm(fmla.full, data=P1_DF)
summary(P1_mod_full)

P1_mode_pairwise_mean = advanced.procD.lm(fmla.size, fmla.size.plus.mode, 
                                     groups=~mode,data=P1_DF, angle.type = "deg", iter=1000)
P1_mode_pairwise_slope = advanced.procD.lm(fmla.full, fmla.size.plus.mode,slope=~CS, 
                                          groups=~mode,data=P1_DF, angle.type = "deg", iter=1000)


P1.size = procD.lm(fmla.size, data=P1_DF)
P1_DF_size.free = geomorph.data.frame(shape=P1.size$residuals, CS= P1_DF$CS, mode=P1_DF$mode)
P1.size_free.mode = procD.lm(shape~mode, data=P1_DF_size.free)
P1.size_free.mode_pair=advanced.procD.lm(shape~CS,shape~CS+mode, groups=~mode, data=P1.MEAN_gdf, angle.type = "deg", iter=1000)


# P7
P7_mod_full = procD.lm(fmla.full, data=P7_DF)
summary(P7_mod_full) 

P7_mode_pairwise_slope = advanced.procD.lm(fmla.full, fmla.size.plus.mode, slope=~CS,
                                           groups=~mode,data=P7_DF, angle.type = "deg", iter=1000)
P7_mode_pairwise_mean = advanced.procD.lm(fmla.size, fmla.size.plus.mode,
                                          groups=~mode,data=P7_DF, angle.type = "deg", iter=1000)

P7.size = procD.lm(fmla.size, data=P7_DF) 
P7_DF_size.free = geomorph.data.frame(shape=P7.size$residuals, mode=P7_DF$mode)
P7.size_free.mode = procD.lm(shape~mode, data=P7_DF_size.free)
P7.size_free.mode_pair=advanced.procD.lm(shape~mode,shape~1, groups=~mode, data=P7.MEAN_gdf, angle.type = "deg", iter=1000)


# 3.
###########
# Phylogentic Generalised Least Squares (PGLS)
# 1. Estimate lambda while fitting the model sensu Revell (2010) using nlme::gls
# 2. Use estimate to make var-covar correction matrix
# 3. procD.lm with correction matrix


P1.sig_df = data.frame(species=row.names(P1_pca_multi_phy),P1_pca_multi_phy[,1:12])
P1.sig_df2 =data.frame(P1.sig_df,CS=log(P1.MEAN.CS_phy), mode=P1.MEAN.MODE_phy)
colnames(P1.sig_df2)[2:13]=LETTERS[1:12]

P7.sig_df = data.frame(species=row.names(P7_pca_multi_phy),P7_pca_multi_phy[,1:12])
P7.sig_df2 =data.frame(P7.sig_df,CS=log(P7.MEAN.CS_phy), mode=P7.MEAN.MODE_phy)
colnames(P7.sig_df2)[2:13]=LETTERS[1:12]


### P1
# First check loglik profile of lambda for full model
lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode,
                                                        correlation = corPagel(value = lambda, phy = phy_chron, fixed = TRUE),
                                                        data = P1.sig_df2)))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                       lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
fitPagel <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(value=0, phy=phy_chron, fixed=T), data=P1.sig_df2)
fitPagel2 <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(value=0.5, phy=phy_chron, fixed=T), data=P1.sig_df2)
fitPagel3 <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(value=1, phy=phy_chron, fixed=T),data=P1.sig_df2)
fitPagel4 <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(0, phy=phy_chron, fixed=F),method="REML",data=P1.sig_df2)
abline(v = fitPagel4$modelStruct$corStruct, col = "red")
# lambda = 0.19

### Using scores from phylogenetic PCA

PCov=vcv(corPagel(0.19, phy=phy_chron))

p1_pgls_phy = procD.lm(fmla.full.phy,data=P1_DF_phy,Cov=PCov, SS.type="I")

p1_pgls_pairwise_slope_ppca = advanced.procD.lm(shape~CS*pmode,shape~CS+pmode,
                                                group=~pmode, slope=~CS,Cov=PCov,data=P1_DF_phy_gdf, angle.type = "deg")
p1_pgls_pairwise_mean_ppca = advanced.procD.lm(shape~CS,shape~CS+pmode, group=~pmode, Cov=PCov,data=P1_DF_phy_gdf)


#### P7
# First check loglik profile of lambda for full model
lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda, function(lambda) logLik(nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode,
                                                        correlation = corPagel(value = lambda, phy = phy_chron, fixed = TRUE),
                                                        data = P7.sig_df2)))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
                                                       lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
fitPagel <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corMartins(value=0, phy=phy_chron, fixed=T), data=P7.sig_df2)
fitPagel2 <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(value=0.5, phy=phy_chron, fixed=T), data=P7.sig_df2)
fitPagel3 <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(value=1, phy=phy_chron, fixed=T),data=P7.sig_df2)
fitPagel4 <- nlme::gls(A+B+C+D+E+F+G+H+I+J+K+L~CS*mode, correlation=corPagel(value=0, phy=phy_chron,fixed=F), method="REML",data=P7.sig_df2)
# REML estimate is 0.17

PCov7=vcv(corPagel(0.17,phy_chron))

P7_pgls_phy = procD.lm(fmla.full.phy, Cov=PCov7,data = P7_DF_phy, SS.type="I")

P7_pgls_pairwise_slope = advanced.procD.lm(shape~CS*pmode,shape~CS+pmode, group=~pmode, slope=~CS,Cov=PCov7,
                                           iter = 999, effect.type = "F", data = P7_DF_phy_gdf, angle.type = "deg")
P7_pgls_pairwise_mean = advanced.procD.lm(shape~CS,shape~CS+mode, group=~mode, Cov=PCov7,
                                          iter = 999, effect.type = "F", data = P7.MEAN_gdf_phy, angle.type = "deg")



###-----------------------------
###       5. Plots
###-----------------------------


# 1. Phylomorphospaces
# Cairo(file="2019-17-02_phylogeomorph-P1.png",
#       type="png",
#       units="cm",
#       width=30,
#       height=25,
#       pointsize=1,
#       dpi=500)
plotGMPhyloMorphoSpace(phy, P1.MEAN.coords_phy, tip.labels = T, node.labels = F, plot.param=list(lwd=2.5, t.pch=NA, n.pch=NA))
#dev.off()

# Cairo(file="figures/2019-17-02_phylogeomorph-P7.png",
#        type="png",
#        units="cm",
#        width=30,
#        height=25,
#        pointsize=1,
#        dpi=500)
plotGMPhyloMorphoSpace(phy, P7.MEAN.coords_phy, tip.labels = T, node.labels = F, plot.param=list(lwd=3, t.pch=NA, n.pch=NA))
#dev.off()


# 2. Tree
# Can't directly access tree from chronos object. So can make our own with attributes: 

# plot(phy_chron)
# Cairo(file="phylogeny.png",
#       type="png",
#       units="cm",
#       width=25,
#       height=30,
#       pointsize=1,
#       dpi=600)
plot.phylo(ladderize(phy_chron, right=F), edge.width=2.5, label.offset=0.1, adj= 0,cex=15, align.tip.label = 0)
add.scale.bar(x=0.1,y=12, cex=10, lwd=1.5)
#dev.off()



# 3.  PCA 

########
## P1 ##
########
a_P1_pca_plot <-  data.frame(one=P1_pca$x[,1]*-1,two=P1_pca$x[,2], three=P1_pca$x[,3], mode=P1.MEAN.MODE) # flipped PC1 to match rotation of backtransforms
P1_PC1_2 <- myPCA2(a_P1_pca_plot,col=c(2),ylim_max = c(0.06), ylim_min = c(-0.07), xlim_max = c(0.08), xlim_min = c(-0.07))
#P1_PC1_3 <- myPCA2(a_P1_pca_plot,col=c(3),ylim_max = c(0.05), ylim_min = c(-0.05), xlim_max = c(0.08), xlim_min = c(-0.06))
P1_PC1_2
#P1_PC1_3


########
## P7 ##
########

a_P7_pca_plot <-  data.frame(P7_pca$x[,1],P7_pca$x[,2], P7_pca$x[,3], P7.MEAN.MODE)
P7_PC1_2 <-myPCA2(a_P7_pca_plot,col=c(2),ylim_max = c(0.1), ylim_min = c(-0.09), xlim_max = c(0.1), xlim_min = c(-0.11))
#P7_PC1_3 <- myPCA2(a_P7_pca_plot,col=c(3),ylim_max = c(0.04), ylim_min = c(-0.05), xlim_max = c(0.11), xlim_min = c(-0.11))
P7_PC1_2
#P7_PC1_3


#PRINT

# print.names= c("P7_PC1_2","P7_PC1_3","P1_PC1_2","P1_PC1_3")
# print.list= list(P7_PC1_2,P7_PC1_3, P1_PC1_2,P1_PC1_3)
# for (i in 1:length(print.names)){ 
# Cairo(file=paste(print.names[[i]],".png"), 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=18, 
#       pointsize=1, 
#       dpi=600)
# print(print.list[[i]])
# dev.off()}


# 4. Allometric trajectories OLS and PGLS - not included in paper but coold to see diffs after pgls

# PD1 = procD.lm(shape ~ mode*CS, Cov = PCov, data = P1.MEAN_gdf_phy, iter = 999, RRPP = TRUE)
# rrr1 = plot(PD1, predictor = P1.MEAN_gdf_phy$CS, type="regression",reg.type="PredLine")
# plot(rrr1$PredLine~P1.MEAN_gdf_phy$CS )
# PD1 = procD.lm(shape ~ mode*CS, data = P1.MEAN_gdf_phy, iter = 999, RRPP = TRUE)
# rrr1 = plot(PD1, predictor = P1.MEAN_gdf_phy$CS, type="regression",reg.type="PredLine")
# plot(rrr1$PredLine~P1.MEAN_gdf_phy$CS )
# 
# PD7 = procD.lm(shape ~ mode*CS,Cov = PCov7, data = P7.MEAN_gdf_phy)
# rrr7 = plot(PD7, predictor = P7.MEAN_gdf_phy$CS, type="regression",reg.type="PredLine")
# plot(rrr7$PredLine~P7.MEAN_gdf_phy$CS )
# PD7 = procD.lm(shape ~ mode*CS, data = P7.MEAN_gdf_phy)
# rrr7 = plot(PD7, predictor = P7.MEAN_gdf_phy$CS, type="regression",reg.type="PredLine")
# plot(rrr7$PredLine~P7.MEAN_gdf_phy$CS )

# PD7phy = procD.lm(shape ~ pmode*CS, Cov = PCov7, data = P7_DF_phy_gdf, iter = 999, RRPP = TRUE)
# PD7phy_data = plot(PD7phy, predictor = P7_DF_phy_gdf$CS, type="regression",reg.type="PredLine")
# P7_LINES_phy=data.frame(Pred=(PD7phy_data$PredLine), mode=PD7phy_data$groups, CS= P7_DF_phy_gdf$CS)
# 
# PD7 = procD.lm(fmla.full, data=P7_DF_phy)
# PD7_data= plot(PD7, predictor = P7_DF_phy$CS, type="regression",reg.type="PredLine")
# P7_LINES = data.frame(Pred=(PD7_data$PredLine), mode=PD7_data$groups, CS= P7_DF_phy$CS)
# 
# Cairo(file="sP7_slopes.png",
#       type="png",
#       units="cm",
#       width=12,
#       height=10,
#       pointsize=1,
#       dpi=500)
# ggplot() + 
#   geom_line(data=P7_LINES_phy, aes(CS, Pred, group = pmode, color=pmode),size=1 ) +
#   geom_line(data=P7_LINES, aes(CS, Pred, group = mode, color=mode), alpha=0.4,size=1)+
#   scale_color_manual(values=c("#1b9e90","#E84A5F", "#F3B54A")) +
#   theme_bw()+
#   theme( panel.background = element_rect(fill = "transparent"), # bg of the panel
#          plot.background = element_rect(fill = "transparent", color = NA),
#          panel.grid.major = element_blank(), panel.grid.minor = element_line(size=0.05))+
#   theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
#   theme(legend.position="none")+
#   scale_y_continuous(breaks=pretty_breaks(n=4))+
#   theme(axis.text.x=element_text(colour="gray27", size=4), 
#         axis.text.y=element_text(colour = "gray27", angle = 90,hjust=0.5, size=4),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.line = element_line(colour = 'gray27', size = 0.07),
#         axis.ticks = element_line(colour = "black", size = 0.07))
# dev.off()
# 
# 
# PD1phy = procD.lm(shape ~ pmode*CS, Cov = PCov, data = P1_DF_phy_gdf, iter = 999, RRPP = TRUE)
# PD1phy_data = plot(PD1phy, predictor = P1_DF_phy_gdf$CS, type="regression",reg.type="PredLine")
# P1_LINES_phy=data.frame((Pred=(PD1phy_data$PredLine))[c(3,13,1,2,5:12,14,4,15)], mode=PD1phy_data$groups, CS= P1_DF_phy_gdf$CS)
# 
# 
# PD1 = procD.lm(fmla.full,data=P1_DF_phy)
# PD1_data= plot(PD1, predictor = P1_DF_phy$CS, type="regression",reg.type="PredLine")
# P1_LINES = data.frame(Pred=(PD1_data$PredLine), mode=PD1_data$groups, CS= P1_DF_phy$CS)
# 
# 
# Cairo(file="P1_slopes.png",
#       type="png",
#       units="cm",
#       width=12,
#       height=10,
#       pointsize=1,
#       dpi=500)
# ggplot() + 
#   geom_line(data=P1_LINES_phy, aes(CS, Pred, group = pmode, color=pmode),size=1 ) +
#   geom_line(data=P1_LINES, aes(CS, Pred, group = mode, color=mode), alpha=0.4,size=1)+
#   scale_color_manual(values=c("#1b9e90","#E84A5F", "#F3B54A")) +
#   theme_bw()+
#   theme( panel.background = element_rect(fill = "transparent"), # bg of the panel
#          plot.background = element_rect(fill = "transparent", color = NA),
#          panel.grid.major = element_blank(), panel.grid.minor = element_line(size=0.05))+
#   theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm")) +
#   theme(legend.position="none")+
#   scale_y_continuous(breaks=pretty_breaks(n=4))+
#   theme(axis.text.x=element_text(colour="gray27", size=4), 
#         axis.text.y=element_text(colour = "gray27", angle = 90,hjust=0.5, size=4),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.line = element_line(colour = 'gray27', size = 0.07),
#         axis.ticks = element_line(colour = "black", size = 0.07))
# dev.off()

#5. Backtransforms

##### Backtransforms #####

resEig <-eigen(cov(two.d.array(P1.MEAN.FIT2$coords)))
scores <- two.d.array(P1.MEAN.FIT2$coords) %*% resEig$vectors
pcs <- 1:2
lm_array <- P1.MEAN.FIT2$coords[, 1:2, ]
dimnames(lm_array)<- list(c(1:42), dimnames(lm_array)[[2]], 
                          dimnames(lm_array)[[3]])

# Cairo(file="p1_back.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=18, 
#       pointsize=1, 
#       dpi=600)
plot(scores[, pcs], type='n', ylab = "", xlab = "", xaxt='n', yaxt="n", frame.plot = F)
btShapes(scores=scores, vectors=resEig$vectors, fcn=plot_beak_lateral, 
         pcs=pcs, n=c(5,6), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.08,0.08), phy.mean=NULL, size=0.05, col=gray(0.7),centroid.size=NULL)
#dev.off()


##### P7


resEig7 <- eigen(cov(two.d.array(P7.MEAN.FIT2$coords)))
scores7 <- (two.d.array(P7.MEAN.FIT2$coords) %*% (resEig7$vectors))
pcs <- 1:2
lm_array7 <- P7.MEAN.FIT2$coords[, 1:2, ]
dimnames(lm_array7)<- list(c(1:42), dimnames(lm_array7)[[2]], 
                           dimnames(lm_array7)[[3]])

# Cairo(file="p7_back.png",
#       type="png",
#       units="cm",
#       width=20,
#       height=18,
#       pointsize=1,
#       dpi=600)
plot(scores7[, pcs], type='n', ylab = "", xlab = "", xaxt='n', yaxt="n", frame.plot=F)
btShapes(scores=scores7, vectors=(resEig7$vectors), fcn=plot_beak_lateral, 
         pcs=pcs, n=c(5,6), m=dim(lm_array7)[2], row.names=dimnames(lm_array7)[[1]], 
         pc.margin=c(0.08,0.08), phy.mean=NULL, size=0.065, col=gray(0.7),centroid.size=NULL)
#dev.off()



#6. Warps  - Another way to visualise shape change across morphospace, not as nice as backtrans

#P1
# a_P1_pca = geomorph::plotTangentSpace(a_P1_gpa$coords)
# ref <- as.data.frame(mshape(a_P1_gpa$coords))
# p1_PC1_min <-  as.data.frame(a_P1_pca$pc.shapes$PC1min)
# p1_PC1_max <- as.data.frame(a_P1_pca$pc.shapes$PC1max)
# p1_PC2_min <- as.data.frame(a_P1_pca$pc.shapes$PC2min)
# p1_PC2_max <-as.data.frame(a_P1_pca$pc.shapes$PC2max)
# p1_PC3_min <- as.data.frame(a_P1_pca$pc.shapes$PC3min)
# p1_PC3_max <-as.data.frame(a_P1_pca$pc.shapes$PC3max)
# P1_minMax <-  cbind(ref,p1_PC1_min,p1_PC1_max, p1_PC2_min ,p1_PC2_max, p1_PC3_min,p1_PC3_max)
# 
# colnames(P1_minMax)=c("ref_x","ref_y","p1min_x","p1min_y","p1max_x","p1max_y","p2min_x","p2min_y","p2max_x","p2max_y","p3min_x","p3min_y","p3max_x","p3max_y")
# ord <-  c(1,4:16,2,17:42,3)
# 
# plotRefToTarget(ref,P1_minMax[,11:12], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,11][ord],P1_minMax[,12][ord], col="#E84A5F", lwd=2)
# 
# 
# Cairo(file="pc1min.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=600)
# plotRefToTarget(ref,P1_minMax[,3:4], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,3][ord],P1_minMax[,4][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc1max.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=600)
# plotRefToTarget(ref,P1_minMax[,5:6], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,5][ord],P1_minMax[,6][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc2min.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref,P1_minMax[,7:8], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,7][ord],P1_minMax[,8][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc2max.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref,P1_minMax[,9:10], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,9][ord],P1_minMax[,10][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc3min.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref,P1_minMax[,11:12], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,11][ord],P1_minMax[,12][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc3max.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref,P1_minMax[,13:14], method="vector", gridPars = gridPar(pt.size=0))
# lines(P1_minMax[,1][ord],P1_minMax[,2][ord], col="dark grey", lty=2)
# lines(P1_minMax[,13][ord],P1_minMax[,14][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# 
# 
# ### P7
# a_P7_pca = geomorph::plotTangentSpace(a_P7_gpa$coords)
# ref7 <- as.data.frame(mshape(a_P7_gpa$coords))
# p7_PC1_min <-  as.data.frame(a_P7_pca$pc.shapes$PC1min)
# p7_PC1_max <- as.data.frame(a_P7_pca$pc.shapes$PC1max)
# p7_PC2_min <- as.data.frame(a_P7_pca$pc.shapes$PC2min)
# p7_PC2_max <-as.data.frame(a_P7_pca$pc.shapes$PC2max)
# p7_PC3_min <- as.data.frame(a_P7_pca$pc.shapes$PC3min)
# p7_PC3_max <-as.data.frame(a_P7_pca$pc.shapes$PC3max)
# P7_minMax <-  cbind(ref7,p7_PC1_min,p7_PC1_max, p7_PC2_min ,p7_PC2_max, p7_PC3_min,p7_PC3_max)
# 
# colnames(P7_minMax)=c("ref_x","ref_y","p1min_x","p1min_y","p1max_x","p1max_y","p2min_x","p2min_y","p2max_x","p2max_y","p3min_x","p3min_y","p3max_x","p3max_y")
# ord <-  c(1,4:16,2,17:42,3)
# 
# Cairo(file="pc1min7.png", 
#       type="png",
#       units="cm", 
#       width=30, 
#       height=25, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref7,P7_minMax[,3:4], method="vector", gridPars = gridPar(pt.size=0))
# lines(P7_minMax[,1][ord],P7_minMax[,2][ord], col="dark grey", lty=2)
# lines(P7_minMax[,3][ord],P7_minMax[,4][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc1max7.png", 
#       type="png",
#       units="cm", 
#       width=30, 
#       height=25, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref7,P7_minMax[,5:6], method="vector", gridPars = gridPar(pt.size=0))
# lines(P7_minMax[,1][ord],P7_minMax[,2][ord], col="dark grey", lty=2)
# lines(P7_minMax[,5][ord],P7_minMax[,6][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc2min7.png", 
#       type="png",
#       units="cm", 
#       width=30, 
#       height=25, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref7,P7_minMax[,7:8], method="vector", gridPars = gridPar(pt.size=0))
# lines(P7_minMax[,1][ord],P7_minMax[,2][ord], col="dark grey", lty=2)
# lines(P7_minMax[,7][ord],P7_minMax[,8][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc2max7.png", 
#       type="png",
#       units="cm", 
#       width=30, 
#       height=25, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref7,P7_minMax[,9:10], method="vector", gridPars = gridPar(pt.size=0))
# lines(P7_minMax[,1][ord],P7_minMax[,2][ord], col="dark grey", lty=2)
# lines(P7_minMax[,9][ord],P7_minMax[,10][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# 
# Cairo(file="pc3min7.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref,P7_minMax[,11:12], method="vector", gridPars = gridPar(pt.size=0))
# lines(P7_minMax[,1][ord],P7_minMax[,2][ord], col="dark grey", lty=2)
# lines(P7_minMax[,11][ord],P7_minMax[,12][ord], col="#E84A5F", lwd=2)
# dev.off()
# 
# Cairo(file="pc3max7.png", 
#       type="png",
#       units="cm", 
#       width=20, 
#       height=17, 
#       pointsize=1, 
#       dpi=800)
# plotRefToTarget(ref,P7_minMax[,13:14], method="vector", gridPars = gridPar(pt.size=0))
# lines(P7_minMax[,1][ord],P7_minMax[,2][ord], col="dark grey", lty=2)
# lines(P7_minMax[,13][ord],P7_minMax[,14][ord], col="#E84A5F", lwd=2)
# dev.off()




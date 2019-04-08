#!/usr/bin/env Rscript
# functions for geomorph paper

# distances for each group to the overall mean of all groups combined
myDist_groups2overalMean <- function(x,shape, grp, grpName){
  mn <- matrix(t(mshape(x[[shape]])), nrow = 1)
  A.d = two.d.array(x[[shape]][,,which(x[[grp]]==grpName)])
  d <- NULL
  for (i in 1:nrow(A.d)) {
    d <- c(d, as.vector(dist(rbind(mn, A.d[i, ]))))
  }
  DF= data.frame(d)
  row.names(DF) = dimnames(x[[shape]][, , which(x[[grp]]==grpName)])[[3]]
  DF = DF[order(d,decreasing=T),,drop=F]
  return(DF)
}

# distances for each group to the mean of each group
myDist_groups2groupMean <- function(x,shape, grp, grpName){
  mn <- matrix(t(mshape(x[[shape]][,,which(x[[grp]]==grpName)])), nrow = 1)
  A.d = two.d.array(x[[shape]][,,which(x[[grp]]==grpName)])
  d <- NULL
  for (i in 1:nrow(A.d)) {
    d <- c(d, as.vector(dist(rbind(mn, A.d[i, ]))))
  }
  DF= data.frame(d)
  row.names(DF) = dimnames(x[[shape]][, , which(x[[grp]]==grpName)])[[3]]
  DF = DF[order(d,decreasing=T),,drop=F]
  return(DF)
}

# plot PCAs nicely
myPCA <- function(x,col,ylim_min,ylim_max){
  a1pca <-  ggplot(x, aes(x=x[,1],y=x[,col]))
  a1pca+geom_point(size=3, alpha=0.5, stroke=1,aes(color=x[,4], shape=x[,4], fill=x[,4]))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA))+
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))+
    #geom_vline(xintercept=c(0), linetype="dotted")+
    #geom_hline(yintercept=c(0), linetype="dotted")+
    scale_shape_manual(values=c(21, 22, 24), name="",labels=c("External","Buccal","Branchial"))+ 
    scale_color_manual(values=c("#1b9e90","#E84A5F", "#F3B54A"),name="", labels=c("External","Buccal","Branchial"))+
    scale_fill_manual(values=c("#1b9e90","#E84A5F", "#F3B54A"), name="",labels=c("External","Buccal","Branchial"))+
    scale_x_continuous(breaks=seq(-0.2,0.2, 0.05), labels=seq(-0.2,0.2, 0.05))+
    scale_y_continuous(breaks=seq(-0.2,0.2, 0.05), labels=seq(-0.2,0.2, 0.05))+
    ylim(ylim_min,ylim_max)+
    #xlim(xlim_min,xlim_max)+
    theme(axis.text.x=element_text(colour="gray27", size=10), 
          axis.text.y=element_text(colour = "gray27", angle = 90,hjust=0.5, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line = element_line(colour = 'gray27', size = 0.1))+
    #theme(axis.title.x=element_text(colour="gray27"), axis.title.y=element_text(colour="gray27"))+
    #ggtitle("PCA: Pereopod 1") + 
    #xlab("Principal Component 1") + ylab(ylabel)+
    theme(legend.position="none")+
    #       legend.text = element_text(colour="grey27", size=5),
    #      legend.key = element_blank(),
    #      legend.title = element_blank(),
    #      legend.position = c(0.96, 0.87),
    #      legend.key.width = unit(0.5,"cm"),
    #      legend.key.height = unit(0.5,"cm"),
    #      legend.background = element_rect("transparent"))+
    # guides(color = guide_legend(override.aes = list(size=2,
    #                                                 alpha=1,shape=c(21, 22, 24),  fill=c(c("#1b9e90","#E84A5F", "#F3B54A")))))+
    stat_ellipse(type="t", mapping=aes(color=x[,4]), show.legend = F, size=0.4) }

# mean shape coords for each species present in tree
mshapePerTipSpp <- function(x,y){ 
  myList = list()
  for (i in 1:length(y$tip.label)){
    myList[[i]] <- mshape(x[,,which(dimnames(x)[[3]]==y$tip.label[i])])
  }
  names(myList)=y$tip.label
  myArray = abind(myList, along=3)
  return(myArray)
}



#################################



##### from Olsen (2016)

plot_beak_lateral <- function(xy, coor, size=1, col='black'){

	# If 3D, rotate points about x-axis using 3D rotation matrix
	if(ncol(coor) == 3){
		coor <- coor %*% matrix(c(1,0,0, 0,cos(-pi/2),sin(-pi/2), 
			0,-sin(-pi/2),cos(-pi/2)), nrow=3, ncol=3)
	}

	# Get just x,y coordinates (orthographic projection into xy-plane)
	coor <- coor[, 1:2]

	# Get plot aspect ratio
	w <- par('pin')[1]/diff(par('usr')[1:2])
	h <- par('pin')[2]/diff(par('usr')[3:4])
	asp <- w/h

	#Correct for plot aspect ratio not necessarily being 1:1
	coor[, 1] <- coor[, 1] * (1/asp)

	# Scale points and place back in position
	coor <- coor*size

	# Center about zero based on range of coordinates
	coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
		nrow=nrow(coor), ncol=ncol(coor), byrow=TRUE)

	# Move shape to PC score
	coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow=TRUE)

	# Set order in which to draw points to create polygon
	polygon_order <- c(1,4:16,2,17:42,3)

	# Create filled polygon
	polygon(coor[polygon_order, ], col=col, border=col)
}



### myPCA2

myPCA2 <-function(x,col,ylim_min,ylim_max,xlim_min,xlim_max){
  a1pca <-  ggplot(x, aes(x=x[,1],y=x[,col]))
  a1pca+geom_point(size=3, alpha=1, stroke=1,aes(color=x[,4], shape=x[,4], fill=x[,4]))+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent"), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA))+
    theme(plot.margin = unit(c(1.5,1.5,1.5,1.5), "cm"))+
    scale_shape_manual(values=c(21, 22, 24), name="",labels=c("External","Buccal","Branchial"))+ 
    scale_color_manual(values=c("#1b9e90","#E84A5F", "#F3B54A"),name="", labels=c("External","Buccal","Branchial"))+
    scale_fill_manual(values=c("#1b9e90","#E84A5F", "#F3B54A"), name="",labels=c("External","Buccal","Branchial"))+
    scale_x_continuous(breaks=seq(-0.2,0.2, 0.05), labels=seq(-0.2,0.2, 0.05))+
    scale_y_continuous(breaks=seq(-0.2,0.2, 0.05), labels=seq(-0.2,0.2, 0.05))+
    ylim(ylim_min,ylim_max)+
    xlim(xlim_min,xlim_max)+
    theme(axis.text.x=element_text(colour="gray27", size=10), 
          axis.text.y=element_text(colour = "gray27", angle = 90,hjust=0.5, size=10),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line = element_line(colour = 'gray27', size = 0.1))+
    theme(legend.position="none")
}



	

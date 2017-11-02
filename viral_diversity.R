set.seed(1234)


# Install the required packages:
# install.packages( c("fossil","vegan","FactoClass","scatterplot3d","MASS","cluster","gplots","RcolorBrewer","pvclust","splines", "stats4", "survival", "mvtnorm", "modeltools", "coin","fields") )
# source("https://bioconductor.org/biocLite.R")
# biocLite("metagenomeSeq")

####### FUNCTIONS  #######
library(fossil)
library(vegan)
source("./supplementary_diversity_functions.R")
library(RCurl)
###############################################################################
############## Metadata and sampling ##########################################
###############################################################################

## Read the data table of the relative abundaces of the viral "species". 
file.i <- "https://raw.githubusercontent.com/jorgevazcast/Viromic-diversity/master/simulated_viral_data.tsv"
abundance.table<-read.table(textConnection(getURL(file.i)), header=T, row.names = 1 ,dec=".",sep="\t") ## Function to read the file that is in the table format
# In case that you download the file
# abundance.table<-read.table(file="simulated_viral_data.tsv", header=T, row.names = 1 ,dec=".",sep="\t") ## Function to read the file that is in the table format
head(abundance.table)

## Create the metadata table
metadata<-cbind(colnames(abundance.table),c(rep("H",100) , rep("C",100) )) # The metadata simulates two conditions, the "C" and the "H" condition
rownames(metadata)<-colnames(abundance.table) # Change the name of the 
metadata<-as.data.frame(metadata)
colnames(metadata)<-c("samle_name","HC")

### Exclude rare species that could be taxonomic binning artifacts
zero.vec<-apply(  abundance.table,1,function(x){ length(x[x==0])/length(x)}  )
### Exclude all those virus like particles whose viral like particles is above the 10% of the samples
abundance.table.nr.sp<-abundance.table[zero.vec<0.1,]

### Compare 
dim(abundance.table)
dim(abundance.table.nr.sp)

### Plot the diversity curves using the function rarecurve from the vegan libary
min.number.secs<-min(apply(abundance.table.nr.sp,2,sum)) ### Estimates the minimum number of sequences
#png(filename = "rarefaction_curves.png", width = 600, height = 600)
pdf(file = "Fig1_rarefaction_curves.pdf")
	rarecurve(t(abundance.table.nr.sp), step = 200, sample = min.number.secs, col = "blue", cex = 0.6) ## The 
dev.off()

#### Rarefaction and normalization
##### Rarefaction: performs a rarefaction to the sample with the minimum number of reads
viromic.rar<-rrarefy(t(abundance.table.nr.sp), min(apply(abundance.table.nr.sp,2,sum))) # Performs the rarefaction
apply(viromic.rar,1,sum) # Verify the output

##### Standardization: performs the CSS normalization
library("metagenomeSeq") # Use the metagenomeSeq library, where there is the CSS normalization function 
viromic.CSS<-t(cumNormMat(as.matrix(abundance.table.nr.sp))) # It standardizes each one of the samples by its respective CSS value

###############################################################################
############## Alpha diversity ################################################
###############################################################################
#### Alpha diversity indexes 
Shannon <- diversity(viromic.rar) # Shannon index, the function is coded in the "vegan" library
Specnumber <- specnumber(viromic.rar) # Expected number of species, the function is coded in the "vegan" library
Pielou <- pielou(viromic.rar) # Pielou index, the function is coded in the supplementary functions file "supplementary_diversity_functions.R"
Chao1<-apply(viromic.rar,1, chao1) # Chao1 index, the function is coded in the "fossil" library  
ACE<-apply(viromic.rar,1, ACE) # ACE index, the function is coded in the "fossil" library  

pdf("Fig2_Alpha_diversity_index_comparison.pdf") #  Plots the boxplot of the different diversity index. The significance between condition was given by the Wilcoxon test
par(mfrow=c(3,2)) # Divide the plotting window in trhee rows and two columns
p.val<-round(wilcox.test(Shannon~metadata[,2])$p.value, digits=3) # Calculate the p-value of the Wilcoxon test
 boxplot(Shannon ~ metadata[,2],col=c("red","darkgreen"), main="Shannon index",sub=paste("WT p-val =",p.val), ylab="Shannon index") # Represents the data in a boxplot
legend("bottomleft", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=0.7,box.lty=0, title="Legend") # Add the plot legend
p.val<-round(wilcox.test(Specnumber~metadata[,2])$p.value, digits=3)
 boxplot(Specnumber ~ metadata[,2],col=c("red","darkgreen"), main="Expected number of species",sub=paste("WT p-val =",p.val), ylab="Num species")
legend("bottomleft", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=0.7,box.lty=0, title="Legend") # Add the plot legend
p.val<-round(wilcox.test(Pielou~metadata[,2])$p.value, digits=3)
 boxplot(Pielou ~ metadata[,2],col=c("red","darkgreen"), main="Pielou index",sub=paste("WT p-val =",p.val), ylab="Pielou index")
legend("bottomleft", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=0.7,box.lty=0, title="Legend") # Add the plot legend
p.val<-round(wilcox.test(Chao1~metadata[,2])$p.value, digits=3)
 boxplot(Chao1 ~ metadata[,2],col=c("red","darkgreen"), main="Chao1 index",sub=paste("WT p-val =",p.val), ylab="Chao1 index")
legend("bottomleft", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=0.7,box.lty=0, title="Legend") # Add the plot legend
p.val<-round(wilcox.test(ACE~metadata[,2])$p.value, digits=3)
 boxplot(ACE ~ metadata[,2],col=c("red","darkgreen"), main="ACE index",sub=paste("WT p-val =",p.val), ylab="ACE index")
legend("bottomleft", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=0.7,box.lty=0, title="Legend") # Add the plot legend
dev.off()

###############################################################################
############## Beta diversidad  ###############################################
###############################################################################

### Quantitative data distance matrix
PCoA<-capscale(viromic.CSS ~ 1, metaMDS = TRUE,sqrt.dist= TRUE) # Estimates a PCoA using the capscale function encoded in the vegan library
pdf("Fig3_Default_PCoA.pdf")  
 plot(PCoA) ## Plot the PCoA using default parameters
dev.off() # Close the plotting window

### Dissimilarity index to calculate the beta diversity
# Bray-Curtis
rar.bray.dist<-vegdist(viromic.rar)
CSS.bray.dist<-vegdist(viromic.CSS)
# Hellinger norm and euclidean
rar.hell.dist<-vegdist(decostand(viromic.rar,method="hell"),method="euclidean")
CSS.hell.dist<-vegdist(decostand(viromic.CSS,method="hell"),method="euclidean")
## Jaccard distance
viromic.rar.binary<-viromic.rar
viromic.rar.binary[viromic.rar.binary>0]=1 # Jaccard distance is for binary data, then the abundance matrix must be transformed into a presence/absence matrix
rar.jaccard.dist<-vegdist(viromic.rar.binary, method = "jaccard", binary = TRUE)

#### Compare the group categories using the Adonis test. Note It also can be used the ANOSIM 
CH.cond<-metadata[rownames(as.matrix((rar.bray.dist))),2]
adonis.rar.bray<-adonis(   CSS.bray.dist ~ CH.cond    ) # Perform the ADONIS test
adonis.rar.bray$aov.tab[ 6]$Pr[1] # Extract the ADONIS p-value

#### Compare the variance homogeneity between groups
betadisper.res<-betadisper(CSS.bray.dist, CH.cond, type = c("median")) # Perform the betadisper test
permutest(betadisper.res, pairwise = FALSE, permutations = 9999) # Look for the significance for the test

### Represents the PCoA using different function approximations
PCoA.scores<-scores(PCoA) ### Function that decouples data from a vegan function
eig <- eigenvals(PCoA) #  calculate the eigenvalues
eigVariance<- round(eig / sum(eig) * 100,digits=3) # proportion of variance explained
x.lab<-paste("Comp1",paste(eigVariance[1],"%",sep="")) # proportion of variance explained of the first component
y.lab<-paste("Comp2",paste(eigVariance[2],"%",sep="")) # proportion of variance explained of the second component
z.lab<-paste("Comp3",paste(eigVariance[3],"%",sep="")) # proportion of variance explained of the third component

### Different ways of representing a PCoA
pdf("Fig4_Different_PCoA_representations.pdf",width=12, height=10)
par(mfrow=c(2,2))

## Plot the first two components of the PCoA
plot(PCoA.scores$sites, cex = 1, col= c(rep("red",100),rep("darkgreen",100)), pch = 16,### Plot the PCoA
 cex.main=1,cex.sub=0.7, xlab=x.lab, ylab=y.lab, main="PCoA B-C index 2D")

## Plot the first three components of the PCoA, using the libraries: FactoClass and scatterplot3d
library("FactoClass") 
library("scatterplot3d") 
PCoA.scores<-scores(PCoA,choices=c(1,2,3))
scatterplot3d(PCoA.scores$sites[, 1:3], pch = 16, grid=FALSE, box=FALSE,color=c(rep("red",100),rep("blue",100)),
 xlab=x.lab, ylab=y.lab,zlab=z.lab,main="PCoA B-C index 3D") # Plot the first trhee components
addgrids3d(PCoA.scores$sites[, 1:3], grid = c("xy", "xz", "yz")) # Add the grid

## Plot the first two components of the PCoA and the sample density unsing the MASS library
library("MASS")
x <- PCoA.scores$sites[, 1]
y <- PCoA.scores$sites[, 2]
z <- PCoA.scores$sites[, 3]
dens3d <- kde2d(x, y,n = 50)
res<-persp(dens3d, box=TRUE, theta=30,phi=30, xlab=x.lab, ylab=y.lab, zlab="Sample density", main="PCoA B-C index 3D density plot")
## install.packages("fields")
library(fields) ### Library to plot the sample points into the density plot, library fields
density.z<- interp.surface( dens3d, cbind(x,y))
points(trans3d(x, y, density.z, pmat = res), col = c(rep("red",100),rep("darkgreen",100)), pch = 16, cex=1)
legend("bottomleft", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=1.5,box.lty=0, title="Legend") # Add the plot legend
## Heat density plot
image(dens3d, xlab=x.lab, ylab=y.lab, main="PCoA B−C index 3D heatmap density plot") # kernel density estimation
dev.off()

#### Different PCoA_using_different distance index and matrices
distance_matrix_comparison(viromic.CSS,viromic.rar,metadata)

#### Use a diferet aprpximation a NMDS
################## Create the NMDS ####################
pdf("Fig6_NMDS_comparison.pdf",width=17)
par(mfrow=c(1,3))
nmds.CSS.bray<-metaMDS(CSS.bray.dist, trymax=2000,pc=T) # Creates NMDS. Note It could take some time
plot(nmds.CSS.bray$points, cex = 1.5, col= c(rep("darkgreen",100),rep("red",100)), sub=paste("Stress",round(nmds.CSS.bray$stress,digits=3),sep="="),
 cex.main=1.7,cex.sub=1,pch=16, main=" Non-metric multidimensional scaling B-C index")
legend("bottomright", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=1.5,box.lty=0, title="Legend") # Add the plot legend	
CSS.hell.dist<-vegdist(decostand(viromic.CSS,method="hell"),method="euclidean")
nmds.CSS.hell<-metaMDS(CSS.hell.dist, trymax=2000,pc=T) # Crea NMDS
plot(nmds.CSS.hell$points, cex = 1.5, col= c(rep("darkgreen",100),rep("red",100)), sub=paste("Stress",round(nmds.CSS.hell$stress,digits=3),sep="="),
 cex.main=1.7,cex.sub=1,pch=16, main=" Non-metric multidimensional scaling Hellinger distance")
legend("bottomright", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=1.5,box.lty=0, title="Legend") # Add the plot legend
prtest<-protest(nmds.CSS.bray,nmds.CSS.hell) #Perform the protest to verify if the cluster configurations are congruent
plot(prtest)

dev.off()
###############################################################################
############## Clustering analysis   ##########################################
###############################################################################

##### PAM Clusters
#### Estimate the silhouette index and finds the configuration that maximizes the index
library(cluster)
obs.silhouette.temp=NULL
nclusters.silhouette=NULL
pdf("Fig7_NMDS_PAM_clustering.pdf")
# Estimate the value of the silhouette index for clusters configurations from 2 to 20
num.clusters<-(length(colnames(as.matrix(CSS.bray.dist)))-180)
for (k in 1:num.clusters) { 
	if (k==1) {
		obs.silhouette.temp[k]=NA 
	} else {
		data.cluster.temp<-as.vector(pam(CSS.bray.dist, k, diss=TRUE)$clustering)
		obs.silhouette.temp<-mean(silhouette(data.cluster.temp, CSS.bray.dist)[,3])
		nclusters.silhouette[k]=obs.silhouette.temp
	}
}
silhouette.num.clusters<-grep("TRUE",max(nclusters.silhouette, na.rm = T)==nclusters.silhouette) # Obtain the cluster configuration that maximize the silhouette index
data.cluster<-as.vector(pam(CSS.bray.dist, silhouette.num.clusters, diss=TRUE)$clustering) # Convert into a vector
adonis.clusters<- adonis(   CSS.bray.dist ~ data.cluster    )$aov.tab[ 6]$Pr[1] # Obtain the ADONIS p-value
clust.title<-paste("Adonis pval Clusters",round(adonis.clusters,digits=3),sep=" = ")

# Plot the NMDS
plot(nmds.CSS.bray$points, cex = 1, pch=20, col= c(rep("darkgreen",100),rep("red",100)),  
 main=clust.title, sub=paste("Stress",round(nmds.CSS.bray$stress,digits=3),sep="="),
 cex.main=0.7,cex.sub=0.7)
legend("bottomright", legend=c("C", "H"), pch=16, col=c("red", "darkgreen"), cex=1,box.lty=0, title="Legend") # Add the plot legend
# Plot an ellispes for each of the conditions, H and C
treat<-data.cluster
for(i in unique(treat)) {
	print(i)
	if(i==1){ # set the colors for each condition
		color="darkgreen"
	}
	if(i==2){
		color="red"

	}
	# The ellipes rerepsetn the 80% of viraicne of the group, this can be modified by changing the parameter "conf"
	ordiellipse(nmds.CSS.bray$point[grep(i,treat),], groups=treat[treat==i],col=color,draw="polygon",alpha=40, conf=0.80, ### Its confidence is set to the 
	show.groups=i) 
	# The ordispider function draw a line from the centroid of the ellipse to each of the samples	
	ordispider(nmds.CSS.bray$point[grep(i,treat),], groups=treat[treat==i],col=color,lty=2,lwd=0.7) 
} 

dev.off()

##### Hierarchical clustering and heatmaps
library("gplots")
library (RColorBrewer)

# Estimate the hierarchical clustering
hcluster.col <- hclust(CSS.bray.dist, method ="ward.D")

### Only plot the most abuntand species
prop.viromic.CSS<-prop.table(viromic.CSS,1) * 100
apply(prop.viromic.CSS,1,sum)
vec.above<-apply(prop.viromic.CSS,2,function(x){mean(x)>0.5})
virus.aboce<-names(vec.above[grep(TRUE,vec.above)])

distance.row<-vegdist(t(viromic.CSS[,virus.aboce]))
hcluster.row <- hclust(distance.row, method ="ward.D")

ColSideColors<-c(rep("red",100),rep("blue",100))
heatcolors <- colorRampPalette(brewer.pal(9, "Spectral"))

plot.heatmap.matrix<-t(prop.viromic.CSS)
rownames(plot.heatmap.matrix)<-sapply(rownames(plot.heatmap.matrix),function(x){unlist(strsplit(x, ";"))[5]})


pdf("Fig8_Heatmap_hclust.pdf")
# Calculate and represent the heatmap
heatmap.2( plot.heatmap.matrix,  
	#main = paste( "test"), 
	col=heatcolors,
	trace="none",          
	margins =c(5,7),
	ColSideColors=ColSideColors,   
	colCol =   ColSideColors,    
	dendrogram="both",      
	Rowv = as.dendrogram(hcluster.row),  
	Colv = as.dendrogram(hcluster.col), 
	key.xlab = "% Relative abundance",
	cexRow =1,
	cexCol = 0.3,
	na.rm = TRUE ) 
dev.off()
	
##### Add p-values to the clustering methods
library("pvclust")
pdf("Fig9_Bootstrapping_HC.pdf",width=24)
pv.hclust <- pvclust(t(viromic.CSS), method.dist="cor", method.hclust="ward.D",parallel=T, nboot=1000)
plot(pv.hclust)
pvrect(pv.hclust, alpha=0.95) 
dev.off()

###############################################################################
############## Biomarker dyscovery   ##########################################
###############################################################################

##### Implement the fitFeatureModel to computes  a differential abundance analysis using a zero-inflated log-normal model
# Create a MRexperiment object
Metadata = AnnotatedDataFrame(metadata)
taxVir<-sapply(rownames(t(viromic.CSS)),function(x){unlist(strsplit(x, ";"))})
taxVir<-t(taxVir)
taxVir<-cbind(rownames(taxVir),taxVir)
colnames(taxVir)=c("Taxonomy","superkingdom","order","family","genus","species")
taxVir=data.frame(taxVir)
OTUdata = AnnotatedDataFrame(taxVir)
obj = newMRexperiment(as.matrix(abundance.table.nr.sp),phenoData=Metadata,featureData=OTUdata) # Create the input matrix

### Normilize data
obj <-cumNorm(obj, p = 0.5)
pd <-pData(obj)
mod <-model.matrix(~1 + HC, data = pd)

# Computes de model
vir.mod =fitFeatureModel(obj, mod) # The function fitFeatureModel is the one that computes de model, it is allocated into the 
df.mod.res<-MRcoefs(vir.mod, number = dim(taxVir)[1]) # Extract the results
df.mod.sigVirtax<-subset(df.mod.res,adjPvalues<0.05)  # Select those coefficients which are statistically significative
df.mod.sigVirtax

#### Creates the infiles for the LEfSe 
##### LEfSe infile
LEfSe.infile<-t(viromic.CSS)
rownames(LEfSe.infile)<-gsub(";","|",rownames(LEfSe.infile)) # Change the ";" to "|" So that LEfSe can use them in the hierarchical classification of its analysis
Classes<-as.character(metadata[colnames(LEfSe.infile),2]) # First row, contain the classes
samples<-colnames(LEfSe.infile) # Second row, the samples name
LEfSe.infile<-rbind(Classes,samples,LEfSe.infile) # Bind all rows
write.table(LEfSe.infile, "LEfSe.infile" ,col.names=F,row.names = TRUE,quote=FALSE,sep = "\t") # Create the data table


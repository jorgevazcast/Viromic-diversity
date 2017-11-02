set.seed(1234)

library(fossil) ### Estimate the diversity index
library(vegan)

## supplementary_diversity_functions
######## Functions
# pielou function: Estimate the Pileu index
pielou<-function(matrix){
 matriz<-t(matrix)
 H <- diversity(matrix)
 S <- specnumber(matrix)
 J <- H/log(S)
 return(J)
}


# Dist.statistics function: Performs a PCoA analysis, the ADONIS test, and the betadisper test
# Inputs: 
# int.matrix: It is the abundance matrix in which the columns are the species counts and the rows are the samples
# dist.name: the desired index/distance to estimate the dissimilarity matrix. Distance availabe: bray (Bray-Curtis), jac (Jaccard), and hell (Hellinger)
# stand.method: the method used to standardize the abundance matrix
# metadata: a data.frame dindicating the group of which each sample belongs.
# Outputs:
# PCoA.stats: a vector whith different statistics: the variance explained by the first three components, the ADONIS p-value, and the betadisper p-value
Dist.statistics<-function(int.matrix="", dist.name="bray",  stand.method="",metadata=""){

 ### Estimate the distance matrix and the PCoA
 if(dist.name=="jac"){
  dist.veg<-vegdist(int.matrix,method=dist.name, binary = TRUE)
  PCoA<-capscale(int.matrix ~ 1, dist=dist.name, binary = TRUE,metaMDS = TRUE) 
 }else if(dist.name=="hell"){
  dist.veg<-vegdist(decostand(int.matrix,method="hell"),method="euclidean")
  hell.norm.matrix<-decostand(int.matrix,method="hell")
  PCoA<-capscale( hell.norm.matrix~ 1, dist="euclidean", metaMDS = TRUE)
 } else{
  dist.veg<-vegdist(int.matrix,method=dist.name)
  PCoA<-capscale(int.matrix ~ 1, dist=dist.name,metaMDS = TRUE)
 }

 ### Estimate the eigenvector values
 eig <- eigenvals(PCoA)
 eigVariance<- round(eig / sum(eig) * 100,digits=3)
 x.lab<-paste("Comp1",paste(eigVariance[1],"%",sep=""))
 y.lab<-paste("Comp2",paste(eigVariance[2],"%",sep=""))
 z.lab<-paste("Comp3",paste(eigVariance[3],"%",sep=""))
 sum.three.comp<-eigVariance[1]+eigVariance[2]+eigVariance[3] ## Sum of the explained variance of the first three components

 ### Perform the ADONIS test
 adonis.pvalue<-adonis(   dist.veg ~ metadata[rownames(as.matrix((dist.veg))),2]    )$aov.tab[ 6]$Pr[1] # Obtain the p-value

 ### Perform the betadisper test
 betadisper.res<-betadisper(dist.veg, metadata[rownames(as.matrix((dist.veg))),2], type = c("median")) # Perform the test
 permutest.p<-permutest(betadisper.res, pairwise = FALSE, permutations = 9999)$tab[6]$Pr[1] # Estimate the significance of the test using 999 permutations, then obtain the p-value

 ### Plot the PCoA
 PCoA.scores<-scores(PCoA,choices=c(1,2,3))
 plot(PCoA.scores$sites, cex = 1, col= c(rep("darkgreen",100),rep("red",100)), main=paste("PCoA",stand.method,dist.name,"distance"),
 sub=paste("ADONIS p-val =",adonis.pvalue), cex.main=1, cex.sub=0.7, ylab=y.lab, xlab=x.lab, pch=16)
legend("bottom", legend=c("C", "H"), pch=15, col=c("red", "darkgreen"), cex=0.7,box.lty=0, title="Legend") # Add the plot legend
 PCoA.stats<-c(sum.three.comp, adonis.pvalue, permutest.p) # Return statistics vector
 names(PCoA.stats)<-c("Sum variance three comp", "ADONIS p-value","betadisper p-value")
 return(PCoA.stats)
}

# Perform different PCoA using the Bray-Curtis, Hellinger and Jaccard distances over the rarefied and CSS normalized datasets. The function is a wrapper to plot 6 different PCoA. 
# Input: The CSS and the rarefied abundance tables as well as the metadata data frame
# Output: A plot with the six PCoA and a matrix with the sum of the three first components, the ADONIS p-value and the betadisp p-value
distance_matrix_comparison <- function(viromic.CSS,viromic.rar,metadata){
	pdf("Fig5_PCoA_using_different_distance_matrix.pdf",width=11)
	par(mfrow=c(2,3))
	matrix.dist.stats<-matrix(0,6,3)
	name.index<-c()
	index<-1
	for(i in c("bray","hell")){
		print(paste("CSS",i))
		matrix.dist.stats[index,]<-Dist.statistics(int.matrix=viromic.CSS, stand.method="CSS", dist.name=i, metadata=metadata)
		name.index[index]<-paste("CSS",i)
		index<-index+1

	}
	for(i in c("bray","hell")){
		print(paste("rar",i))
		matrix.dist.stats[index,]<-Dist.statistics(int.matrix=viromic.rar, stand.method="rar", dist.name=i, metadata=metadata)
		name.index[index]<-paste("rar",i)
		index<-index+1
	}
	viromic.rar.binary<-viromic.rar
	viromic.rar.binary[viromic.rar.binary>0]=1
	name.index[index]<-paste("rar","jac")
	matrix.dist.stats[index,]<-Dist.statistics(int.matrix=viromic.rar.binary, stand.method="rar", dist.name="jac", metadata=metadata)
	index<-index+1

	viromic.CSS.binary<-viromic.CSS
	viromic.CSS.binary[viromic.CSS.binary>0]=1
	name.index[index]<-paste("CSS","jac")
	matrix.dist.stats[index,]<-Dist.statistics(int.matrix=viromic.CSS.binary, stand.method="CSS", dist.name="jac", metadata=metadata)

	dev.off()

	colnames(matrix.dist.stats)<-c("Sum variance three comp", "ADONIS p-value","betadisper p-value")
	rownames(matrix.dist.stats)<-name.index
	return(matrix.dist.stats[sort(name.index),])
}


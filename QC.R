#!/bin/env Rscript

rm(list = ls())

args<-commandArgs(TRUE)
input=toString(args[1])
output=toString(args[2])

setwd(input)

# ****************************************************

# Quality control

# ****************************************************

#source('http://www.bioconductor.org/biocLite.R')
#biocLite('cummeRbund', ask=FALSE)
library(cummeRbund)

setwd(input)
cuff<-readCufflinks()

cuff

disp<-dispersionPlot(genes(cuff))
pdf(paste(output,"/dispersionPlot.pdf", sep=""))
disp
dev.off()
rm(disp)

genes.scv<-fpkmSCVPlot(genes(cuff))
pdf(paste(output,"/genes.scv.pdf", sep=""))
genes.scv
dev.off()
rm(genes.scv)

isoforms.scv<-fpkmSCVPlot(isoforms(cuff))
pdf(paste(output,"/isoforms.scv.pdf", sep=""))
isoforms.scv
dev.off()
rm(isoforms.scv)

dens<-csDensity(genes(cuff))
pdf(paste(output,"/csDensity.pdf",sep=""))
dens
dev.off()
rm(dens)

densRep<-csDensity(genes(cuff),replicates=T)
pdf(paste(output,"/csDensRep.pdf",sep=""))
densRep
dev.off()
rm(densRep)

b<-csBoxplot(genes(cuff))
pdf(paste(output,"/csBoxplot.pdf",sep=""))
b
dev.off()
rm(b)

brep<-csBoxplot(genes(cuff),replicates=T)
# brep
brepY<-brep+coord_cartesian(ylim=c(-5.15,5.15))
pdf(paste(output,"/csBoxplotrep.pdf", sep=""))
brepY
dev.off()
rm(brepY)

s<-csScatterMatrix(genes(cuff))
pdf(paste(output,"/csScatterMatrix.pdf",sep=""))
s
dev.off()
rm(s)

# srep<-csScatterMatrix(genes(cuff),replicates=T)
# srep

pdf(paste(output,"/csDendro.pdf",sep=""))
dend<-csDendro(genes(cuff))
dev.off()
rm(dend)

pdf(paste(output,"/csDendrorep.pdf",sep=""))
dend.rep<-csDendro(genes(cuff),replicates=T)
dev.off()
rm(dend.rep)

genes.PCA<-PCAplot(genes(cuff),"PC1","PC2")
pdf(paste(output,"/genes.PCA.pdf",sep=""))
genes.PCA
dev.off()
rm(genes.PCA)

genes.PCA.rep<-PCAplot(genes(cuff),"PC1","PC2",replicates=T)
pdf(paste(output,"/genes.PCA.rep.pdf",sep=""))
genes.PCA.rep
dev.off()
rm(genes.PCA.rep)

cond.PCoA <- MDSplot(genes(cuff))
pdf(paste(output,"/cond.PCoA.pdf",sep=""))
cond.PCoA
dev.off()

cond.PCoA.rep <- MDSplot(genes(cuff), replicates = T)
pdf(paste(output,"/cond.PCoA.rep.pdf",sep=""))
cond.PCoA.rep
dev.off()

mySigGeneIds<-getSig(cuff,alpha=0.05,level='genes')
mySigGenes<-getGenes(cuff,mySigGeneIds)

h<-csHeatmap(mySigGenes,cluster='both')
h$scales$scales[[2]]$labels <- sub("\\|.*", "",h$scales$scales[[2]]$labels)
h$scales$scales[[2]]$labels <- sub(",.*", "",h$scales$scales[[2]]$labels)
pdf(paste(output,"/heatmap.pdf",sep=""))
h
dev.off()
#tiff(filename = paste(output,"/heatmap.tiff",sep=""),width=2400,height=1350)
#h
#dev.off()
rm(h)

h.rep<-csHeatmap(mySigGenes,cluster='both',replicates=T)
h.rep$scales$scales[[2]]$labels <- sub("\\|.*", "",h.rep$scales$scales[[2]]$labels)
h.rep$scales$scales[[2]]$labels <- sub(",.*", "",h.rep$scales$scales[[2]]$labels)
pdf(paste(output,"/heatmaprep.pdf",sep=""))
h.rep
dev.off()
rm(h.rep)
#tiff(filename = paste(output,"/heatmaprep.tiff",sep=""),width=2400,height=1350)
#h.rep
#dev.off()
#rm(h.rep)

sigMatrix(cuff)
pdf(paste(output,"/Significance_matrix.pdf",sep=""))
sigMatrix(cuff)
dev.off()


rm(list = ls())

q(save="no")

# ****************************************************

setwd("/raid1/Liuw/zebrafish/");
library(WGCNA);
enableWGCNAThreads(6)
options(stringsAsFactors = FALSE)
dat0=read.csv("zebrafish.csv",header=TRUE)
datSummary=dat0[,1:1];
dim(dat0)
datExpr = t(dat0[,2: ncol(dat0)]);
no.samples = dim(datExpr)[[1]];
dim(datExpr)
a=t(duplicated(datExpr))
b=t(datExpr[a==FALSE,]) 
write.table(data.frame(datSummary,b),"zebrafish2.csv",row.names=FALSE,sep=',')
write.csv(data.frame(colnames(dat0)[2:1496],a[1,]),"duplicated-samples.csv") 

dat0=read.csv("zebrafish.csv",header=TRUE)
library(preprocessCore)
datExpr=t(normalize.quantiles(as.matrix(dat0[,-1])))
GeneName= dat0$datSummary
ArrayName= names(data.frame(dat0[,-1]))
powers=c(seq(1,10,by=1),seq(12,14,by=2));
sft=pickSoftThreshold(datExpr, powerVector=powers,networkType = "signed")
RpowerTable=sft[[2]]
sizeGrWindow(9, 5);
pdf('choosing power.pdf');
par(mfrow = c(1,2));cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red");
dev.off()
# Mean connectivity as a function of the soft-thresholding power
sizeGrWindow(9, 5);
pdf('mean connectivity.pdf');
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");
dev.off()
softPower =16
Connectivity=softConnectivity(datExpr,corFnc = "cor", corOptions = "use ='p'",
                              power=softPower,type="signed")
pdf("scale-free.pdf");
scaleFreePlot(Connectivity,nBreaks = 10,truncated = FALSE,removeFirst = FALSE, main = "");
dev.off()
adjacency = adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",
                      type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType="signed");dissTOM = 1-TOM
#method="complete"  ?
geneTree = hclust(as.dist(dissTOM), method = "average")#高版本已经用hclust
minModuleSize =30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 4, 
                            pamRespectsDendro = FALSE,minClusterSize = minModuleSize,
                            cutHeight=0.99);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
sizeGrWindow(12, 9)
pdf("DendroAndColors.pdf")
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, 
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
moduleColors = mergedColors
colorOrder = c("grey", standardColors(unique(moduleColors)));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
pdf("METree.pdf")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
nSamples=nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = cbind.data.frame(datSummary,corPvalueStudent(as.matrix(geneModuleMembership), 
                                                        nSamples));
write.table(data.frame(ArrayName,MEs),"MEs.csv",row.name=F)
kMEdat=data.frame(geneModuleMembership,MMPvalue)
write.table(data.frame(datSummary,kMEdat),"kME-MMPvalue.csv",row.names=FALSE)
k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'", 
                                        type = "signed", power = softPower), 
                              moduleColors,scaleByMax = FALSE)
datout=data.frame(datSummary, colorNEW=moduleColors, k.in)
write.table(datout, file="OutputCancerNetwork.csv", sep=",", row.names=F)
hubs    = chooseTopHubInEachModule(datExpr, moduleColors)
write.csv(data.frame(module=names(hubs),moduleColor=labels2colors(names(hubs)),hub=hubs),
          "num2color.csv",row.names=F)

library(GEOquery)
Data <- getGEO("GSE66688")

##module stability test
modules=unique(moduleColors[moduleColors !='grey']) #一个一个模块来,不用grey
n=length(modules)
nSamples=ncol(dat0)-1
pb <- txtProgressBar(min = 0, max = n, style = 3)
for (p in 1:n){
  inModule = is.finite(match(moduleColors, modules[p]));
  dat2=dat0[inModule,]
  resamples=lapply(1:1000,function(i) a=t(sample(dat2[,2:nSamples],nSamples/8,replace=F)))
  K1=sapply(resamples,softConnectivity,power=softPower,type="signed") #,type="signed"?
  K=softConnectivity(t(dat2[,2:nSamples]),power=softPower,type="signed") #,type="signed"
  #outfile=paste(modules[p],"-edit.txt",sep="")
  write.table(data.frame(mean(cor(K,K1)),apply(cor(K,K1),1,sd)), file = "module-stability-8.csv", row.names = modules[p], append = TRUE, col.names = FALSE, sep = ", ")
  setTxtProgressBar(pb, p)}#所有loop结果写到一个文件
close(pb) 
#write.table(data.frame(mean(cor(K,K1)),apply(cor(K,K1),1,sd)),quote=FALSE,sep=", ",outfile)

gene=read.csv("OutputCancerNetwork.csv",header=T)
library(gProfileR)
for (i in unique(gene$colorNEW)){
  genes=subset(gene$datSummary,gene$colorNEW==i)
  go=gprofiler(genes, 
               organism = "drerio",numeric_ns="")
  write.table(data.frame(mod=i,go),"moduel_enrichment.csv",append =T)}


#PGE Analysis
# perl pge.pl -r entrez -q list
# perl SFDRv166.pl -assoc input.txt -SFDR -out test
setwd("/raid1/Liuw/pge")
arg1 <- "-r entrez"
arg2 <- "-q list"
arg3 <- "-out test"
cmd <- paste("perl", "pge.pl", arg1, arg2)
list=read.csv("zebrafish.csv",header=F)
write(list[list$V1=="65",2],"list")
system(cmd)

for (i in unique(list$V1)) {write(list[list$V1==as.character(i),2],"list");sink(paste("out",i,".txt"));system(cmd);sink()}

#Intermodule proportion
for (i in no){da=(data[data$colorNEW==i,]$kOut/data[data$colorNEW==i,]$kTotal);
d=data.frame(mean=mean(da),sd=sd(da),row.names=i);write.table(d,"intermoduleProp.csv",append=T,col.names = F,sep=",")}
#higher order network
dir(".")
library(WGCNA);options(stringsAsFactors = FALSE);
ME=read.csv("MEs.csv",header=T,sep=" ",row.names=1)[,1:50]
datExpr = ME;
no.samples = dim(datExpr)[[1]];
dim(datExpr)
powers=c(seq(1,10,by=1),seq(12,14,by=2));
sft=pickSoftThreshold(datExpr, powerVector=powers,networkType = "signed")
RpowerTable=sft[[2]]
softPower =7
Connectivity=softConnectivity(datExpr,corFnc = "cor", corOptions = "use ='p'",
                              power=softPower,type="signed")
adjacency = adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'",
                      type = "signed", power = softPower)
TOM = TOMsimilarity(adjacency,TOMType="signed");dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")#高版本已经用hclust
minModuleSize =2;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 4, 
                            pamRespectsDendro = FALSE,minClusterSize = minModuleSize,
                            cutHeight=0.99);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicMods)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
MEDissThres = 0.2
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3);
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
moduleColors = mergedColors
colorOrder = c("grey", standardColors(unique(moduleColors)));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");#
dir.create("higherorder");setwd("./higherorder")
pdf("METree.pdf")
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
sizeGrWindow(12, 9)
pdf("DendroAndColors.pdf")
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = names(ME), 
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()
MEList = moduleEigengenes(datExpr, colors = dynamicMods)
nSamples=nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
datSummary=names(ME)
MMPvalue = cbind.data.frame(datSummary,corPvalueStudent(as.matrix(geneModuleMembership), 
                                                        nSamples));
ArrayName=row.names(ME)
write.table(data.frame(ArrayName,MEs),"MEs.csv",row.name=F)
kMEdat=data.frame(geneModuleMembership,MMPvalue)
write.table(data.frame(datSummary,kMEdat),"kME-MMPvalue.csv",row.names=FALSE)
k.in=intramodularConnectivity(adjacency(datExpr,corFnc = "cor", corOptions = "use ='p'", 
                                        type = "signed", power = softPower), 
                              moduleColors,scaleByMax = FALSE)
datout=data.frame(datSummary, colorNEW=moduleColors, k.in)
write.table(datout, file="OutputCancerNetwork.csv", sep=",", row.names=F)
hubs    = chooseTopHubInEachModule(datExpr, moduleColors)
write.csv(data.frame(module=names(hubs),moduleColor=labels2colors(names(hubs)),hub=hubs),
          "num2color.csv",row.names=F)

##############module preservation
setwd("/raid1/Liuw/zebrafish/preservation/")
datFemale = read.csv("zebrafish-3.csv")
datMale = read.csv("GSE83466-3.csv")
setLabels = c("Female", "Male");
datSummaryFemale=datFemale[,1]
datSummaryMale=datMale[,1]
dim(datMale)
dim(datFemale)
datExprFemale= t(datFemale[,2:(ncol(datFemale)-1)])
no.samplesFemale <- dim(datExprFemale)[[1]]
dim(datExprFemale)
datExprMale= t(datMale[,2:ncol(datMale)])
colorsFemale = datFemale$Module
colnames(datExprMale)=toupper(datMale[,1])
colnames(datExprFemale)=toupper(datFemale[,1])
nSets = 2
ref = 1
test = 2
female2male = match(colnames(datExprFemale), colnames(datExprMale));
table(is.finite(female2male))
datExprMale = datExprMale[, female2male];
all.equal(colnames(datExprFemale), colnames(datExprMale))
multiExpr = list(Female = list(data = datExprFemale), Male = list(data = datExprMale));
multiColor = list(Female = colorsFemale);
mp = modulePreservation(multiExpr, multiColor,referenceNetworks = 1,networkType="signed",nPermutations = 100,randomSeed = 1,parallelCalculation=F,quickCor = 0,verbose = 3)
save(mp, file = "modulePreservation.RData");
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];


# ??ͼleave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
#text = modColors[plotMods];
labs = match(modColors[plotMods], standardColors(unique(modColors)-2))
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
pdf(file="FemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5,onefile=TRUE)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labs = match(modColors[plotMods], standardColors(length(unique(modColors))))
  write.table(data.frame(mod))
  #replace text to labs as number labeling: labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08); 
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], labs, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


# Re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
# Exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
# Create numeric labels for each module
labs = match(modColors[plotMods], standardColors(length(unique(modColors))));  #50 should larger than module number
# Start the plot: open a suitably sized graphical window and set sectioning and margins. Alternatively,
# plot into a pdf file.
sizeGrWindow(10, 9);
pdf(file="PreservationZStatistics.pdf", w=10, h=9)
par(mfrow = c(4,4))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
for (s in 1:ncol(statsZ))
{
  min = min(statsZ[plotMods, s], na.rm = TRUE);
  max = max(statsZ[plotMods, s], na.rm = TRUE);
  if (min > -max/12) min = -max/12
  plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
       main = colnames(statsZ)[s],
       cex = 2.2,
       ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
       ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
       xlim = c(30, 1200),
       cex.lab = 1.2, cex.axis = 1.2)
  labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 1, offs = 0.06);
  #text(moduleSizes[-1], statsZ[-c(1:2), s], labels = letter[-c(1:2)], col = "black"); #modColors[-2]);
  abline(h=0)
  abline(h=2, col = "blue", lty = 2)
  abline(h=10, col = "darkgreen", lty = 2)
}
# If plotting into a file, close it, otherwise it is unreadable.
dev.off();

# ??ͼplot the clusterRepro results.,Find the maximum and minimum
max = 0; min = 1;
for (ref in c(1:nSets))
{
  p = 1;
  for (test in 1:nSets)
  {
    stats = cr[[ref]][[test]]$Actual.IGP;
    max = max(max, stats, na.rm = TRUE);
    min = min(min, stats, na.rm = TRUE);
  }
}
max = max + 0.05*(max-min);
min = min - 0.03*(max-min);
refNames = setLabels
# Plot the results on a uniform scale
sizeGrWindow(10,9);
#pdf(file="Plots/BxHLiverFemaleOnly-clusterRepro-IGP.pdf", w=10, h=9)
par(mfrow = c(2,2))
par(mar = c(3.5, 3.5, 3, 0.4))
par(mgp = c(2.0, 0.5, 0));
for (ref in 1:nSets)
{
  p = 1;
  for (test in 1:nSets)
  {
    stats = cr[[ref]][[test]]$Actual.IGP;
    moduleSizes = table(multiColor[[ref]])
    labelsX = names(moduleSizes);
    #modNumbers = match(labelsX, standardColors(20))
    plotMods = !(labelsX %in% c("grey", "gold"));
    labelsX = labelsX[plotMods]
    xmin = min(moduleSizes[plotMods]);
    xmax = max(moduleSizes[plotMods]);
    xlim = c(xmin * (xmin/xmax)^.20, xmax * (xmax/xmin)^0.15);
    plot(moduleSizes[plotMods], stats, bg = labelsX, pch = 21,
         main = spaste(LETTERS[ref], p, ". Modules: ", refNames[ref], "\n Test data: ", setLabels[test]),
         cex = 2, cex.axis = 1.2, cex.lab = 1.2,
         ylab = "Actual IGP", xlab = "Module size", log = "x", ylim = c(min, max), xlim = xlim)
    #abline(h = stats[colors=="orange",s], col = "grey30", lty = 2)
    abline(h=0)
    labelPoints(moduleSizes[plotMods], stats, labels = labelsX, offs = 0.070,
                jiggle = 0, cex = 1)
    p = p+ 1;
  }
}
# If plotting into a file, close it.
dev.off()


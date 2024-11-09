```{r}
options(stringAsFactors=F)
tumor.cells<-row.names(tumorlist$pml2@meta.data)[which(tumorlist$pml2@meta.data$celltype=="Epithellial cell")]
length(tumor.cells)
tu=as.data.frame(GetAssayData(subset(tumorlist$pml2,cells=tumor.cells),slot = "count"))

immune.cells<-row.names(tumorlist$pml2@meta.data)[which(tumorlist$pml2@meta.data$celltype=="Immune cell")]
length(immune.cells)
imm=as.data.frame(GetAssayData(subset(tumorlist$pml2,cells=immune.cells),slot = "count"))

str.cells<-row.names(tumorlist$pml2@meta.data)[which(tumorlist$pml2@meta.data$celltype=="Stromal cell")]
length(str.cells)
str=as.data.frame(GetAssayData(subset(tumorlist$pml2,cells=str.cells),slot = "count"))

dat=cbind(tu,imm,str)
groupinfo=data.frame(v1=colnames(dat),v2=c(rep("Tumor",ncol(tu)),rep("Immune",ncol(imm)),rep("Stromal",ncol(str))))

geneInfor=annoGene(rownames(dat),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

geneInfor$chr<-factor(geneInfor$chr,levels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9", "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19", "chr20", "chr21", "chr22", "chrM", "chrX", "chrY"))
geneInfor[order(geneInfor$chr),]

dat=dat[rownames(dat) %in% geneInfor[,1],]
write.table(dat,file = "../expFile.txt",sep = '\t',quote = F)
write.table(groupinfo,file = "../groupFile.txt",sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneInfor,file = "../geneFile.txt",sep = '\t',quote = F,col.names = F,row.names = F)
rm(tu,str,imm,dat,geneInfor,groupinfo,str.cells,immune.cells,tumor.cells)
```




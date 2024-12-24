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
```{r}
library(Seurat)
library(ggplot2)
library(infercnv)
options(stringsAsFactors = F)
expFile='../expFile.txt'
groupFile='../groupFile.txt'
geneFile='../geneFile.txt'
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file=groupFile,
                                    delim="\t",
                                    gene_order_file= geneFile,
                                    ref_group_names= c("Immune","Stromal"))
infercnv_obj2 = infercnv::run(infercnv_obj,
                              cutoff=0.1,
                              out_dir=  './inferCNV_crc10' ,
                              cluster_by_groups=TRUE,analysis_mode="subclusters",output_format="pdf",
                              denoise=TRUE, HMM = TRUE,num_threads=4,tumor_subcluster_partition_method ="random_trees")

```

```{r,find malingant}
library(infercna)
library(ggplot2)
library(ggsci)

infercnv.observations <- read.csv("../crc5/infercnv.observations.txt", sep="")
infercnv.references <- read.csv("../crc5/infercnv.references.txt", sep="")

infercnv.observations<-infercnv.observations[rownames(infercnv.references),]
crc5cna<-cbind(infercnv.observations,infercnv.references)
colnames(crc5cna)<-gsub("\\.","-",colnames(crc5cna))

refCells=list()

crc5<-subset(metadata,subset=sample=="crc5")
Immune<-row.names(crc5)[which(crc5$celltype=='Immune cell')]
Stromal<-row.names(crc5)[which(crc5$celltype=='Stromal cell')]
refCells$Immune<-Immune
refCells$Stromal<-Stromal


crc5cna<-cbind(infercnv.observations,infercnv.references)
colnames(crc5cna)<-gsub("\\.","-",colnames(crc5cna))

crc5_signal<-cnaSignal(cna=crc5cna,refCells = refCells)
crc5_cor<-cnaCor(cna=crc5cna,refCells = refCells)
crc5_signal<-as.data.frame(crc5_signal)
crc5_cor<-as.data.frame(crc5_cor)
crc5_signal<-crc5_signal[rownames(crc5_cor),,drop=FALSE]
crc5cs<-cbind(crc5_cor,crc5_signal)

group<-as.data.frame(crc5[,c("cells","celltype")])
group<-group[rownames(crc5cs),]
crc5cs<-cbind(crc5cs,group)
crc5cs<-crc5cs[,-3]
colnames(crc5cs)<-gsub("crc5","cna",colnames(crc5cs))
ep<-subset(crc5cs,subset=celltype=="Epithellial cell")
tumorcrc5<-row.names(ep)[which(ep$cna_cor>0.3)]
rm(crc5,crc5_cor,crc5_signal,crc5cna)
```

```{r}

p<-ggplot(crc5cs,aes(x=cna_signal,y=cna_cor,color=celltype))+theme_classic()+ggtitle("crc5")+geom_point()+geom_hline(yintercept = 0.3,linetype="longdash")+scale_fill_nejm()
p
ggsave(p,filename = "../crc5tumor.tiff")
save(crc5cs,file = "../crc5data.RData")
```







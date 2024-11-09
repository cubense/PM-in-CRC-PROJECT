```{r}
library(Seurat)
library(ggsci)
library(scales)
mypal<-c(ggsci::pal_npg("nrc",alpha = 0.5)(6),ggsci::pal_simpsons("springfield",alpha = 0.5)(7),ggsci::pal_nejm("default",alpha = 0.5)(7),ggsci::pal_aaas("default",alpha = 0.5)(9),ggsci::pal_igv("default",alpha = 0.5)(7))
```


```{r}
filnoep$maincellcluster<-gsub("Epithellial cell",'Epithelial cell',filnoep$maincellcluster)
filnoep$maincellcluster<-gsub("TNK cell",'T/NK cell',filnoep$maincellcluster)
filnoep$maincelltype<-filnoep$maincellcluster
filnoep$maincelltype<-gsub("Tumor",'Epithelial cell',filnoep$maincelltype)
filnoep$maincelltype<-factor(filnoep$maincelltype,levels = c("T/NK cell","B cell","Plasma cell","Myeloid cell","Neutrophil cell","Mast cell","Stromal cell","Epithelial cell"))

sceasy::convertFormat(filnoep, from="seurat", to="anndata", outFile='./filnoep.h5ad')


```


```{r}
library(Seurat)
library(SeuratWrappers)

load("D:/PROJECT/pm/2/2-fastmnn/ep/ep.RData")
normalep<-subset(ep,subset = malignant=="normal")
```

```{r}
ep@meta.data$tissuegroup = paste(ep@meta.data$colon,ep@meta.data$tissue, sep = "_")
ep@meta.data$tumorgroup = paste(ep@meta.data$colon,ep@meta.data$malignant, sep = "_")
table(ep$tissuegroup)
table(ep$tumorgroup)
ep@meta.data$tumorsite = paste(ep@meta.data$site,ep@meta.data$malignant, sep = "_")
table(ep$tumorsite)
```

```{r}

Idents(ep)<-"tumorsite"
eptumorsitemark<-FindAllMarkers(ep,only.pos = T,test.use = "MAST")
write.table(eptumorsitemark, file=paste("D:/PROJECT/pm/2/2-fastmnn/ep/eptumorsitemark.csv", sep="/t"))
Idents(ep)<-"malignant"
eptumormark<-FindAllMarkers(ep,only.pos = T,test.use = "MAST")
write.table(eptumormark, file=paste("D:/PROJECT/pm/2/2-fastmnn/ep/eptumormark.csv", sep="/t"))

```

```{r}

normalep <- NormalizeData(normalep)
normalep <- FindVariableFeatures(normalep )
normalep  <- RunFastMNN(object.list = SplitObject(normalep , split.by = "orig.ident"))
normalep  <- RunUMAP(normalep , reduction = "mnn", dims = 1:30)
normalep  <- FindNeighbors(normalep , reduction = "mnn", dims = 1:30)
normalep <- FindClusters(normalep ,resolution =0.5)
DimPlot(normalep,label = T,repel = T)
```

```{r}
DimPlot(normalep,label = T,repel = T)
DimPlot(normalep,label = T,repel = T,split.by = "site")
DimPlot(normalep,label = T,repel = T,split.by = "colon")
DimPlot(normalep,label = T,repel = T,group.by = "patient")
```


```{r}
DimPlot(normalep,label = T,repel = T,split.by = "site")
DimPlot(normalep,label = T,repel = T,split.by = "colon")
DimPlot(normalep,label = T,repel = T,group.by = "patient")
VlnPlot(normalep,features = c("percent.ribo", "percent.mt","nFeature_RNA", "nCount_RNA"),stack = T,flip = T,pt.size = 0)+ggsci::scale_fill_npg()
```

```{r}
library(magrittr)
library(tidyverse)
clusternep<-FindAllMarkers(normalep,test.use = "MAST",only.pos = T)
genes<-clusternep %>% group_by(cluster) %>% top_n(5, avg_log2FC)
genes<-unique(genes$gene)
DotPlot(normalep,features = genes,cols = c("blue","red"))+coord_flip()
VlnPlot(normalep,features = genes,stack = T,flip = T)
```



```{r}
DotPlot(normalep, features = c("PTPRC","CD3G","CD3E","IGKC","MZB1" ,"CD14","CD68","EPCAM","MS4A2","TPSAB1","CSF3R","MNDA","COL1A1","TAGLN","KLRF1","CD79A", "MS4A1","PECAM1","VWF"),cols ="RdBu") + coord_flip()
DotPlot(normalep, features = c("LGR5","ASCL2","OLFM4","MKI67","PCNA" ,"GUCA2B","SLC26A3","TFF3","SPINK1","REG4","AGR2","MEP1A","FGF15","CLEC2H","LCT","CBR1","EPHX2","LRMP","DCLK1","CD24A","LYZ","DEFA7P"),cols ="RdBu") + coord_flip()
DotPlot(normalep, features =grep("^DEFA",rownames(normalep),value = T),cols = c("blue","red"))+RotatedAxis()
DotPlot(normalep, features = c("NUPR1","RARRES2","MLEC","TMSB10","HES1","SNHG5","TFF3","WFDC2","MUC2","ITLN1","SPINK4","MGST1","ADH1C","UGT2B17","CRYBA2","SCGN","PCSK1N","PTMS","TUBA1A","CEACAM1","CEACAM7","SELENBP1","CA1","GUCA2A","SLC26A3","AQP8","BEST4","OTOP2"),cols ="RdBu") + coord_flip()
DotPlot(normalep, features = c("LGR5","ASCL2","SLC12A2","AXIN2","OLFM4","GKN3","ALPI","APOA1","APOA4","FABP1","MUC2","CLCA3","TFF3","AGR2","LYZ","DEFA6","ANG4","CHGA","CHGB","TAC1","TPH1","NEUROG3","DCLK1","TRPM5","GFI1B","IL25"),cols ="RdBu") + coord_flip()


```



```{r,anno_example}
library(ggsci)
#5,9  Goblet:MUC2,ITLN1,SPINK4,TFF3
#0,2,3,6,7  Undifferentiated cell
#10 BEST4/OTOP2 cell
#1  colonocyte
#4,8 Enteroendocrine cell
#11 Tuft cell


#add column 'main_seurat_cluster' to store cluster ids from this step.
normalep@meta.data$ep_cluster <- normalep@meta.data$RNA_snn_res.0.5
# Change the column of the resolution if you ended up using a different one than 0.5 
cluster.ids <- sort(as.numeric(unique(as.character(normalep@meta.data$RNA_snn_res.0.5))))
 
ep_cluster <- c("Undifferentiated cell","Colonocyte","Undifferentiated cell","Undifferentiated cell","Enteroendocrine cell","Goblet cell","Undifferentiated cell","Undifferentiated cell","Enteroendocrine cell","Goblet cell","BEST4/OTOP2 cell","Tuft cell")
# Add annotation to the Seurat object 
normalep@meta.data$ep_cluster  <- plyr::mapvalues(x = normalep@meta.data$ep_cluster , from = cluster.ids, to = ep_cluster)
# Make a table 
table(normalep@meta.data$ep_cluster)
table(normalep@meta.data$ep_cluster, normalep@meta.data$RNA_snn_res.0.5)
DimPlot(normalep,group.by = "ep_cluster",repel = T,label = T,label.size = 5)+scale_color_igv()


```







```{r}
DotPlot(filnoep, features = c("PTPRC","CD3G","CD3E","CD79A", "MS4A1","IGKC","MZB1","CD14","CD68","CSF3R","MNDA" ,"MS4A2","TPSAB1","COL1A1","TAGLN","EPCAM"),group.by = "maincelltype",cols ="RdBu") + theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
```

```{r}
FeaturePlot(filnoep,features = c("program11","program13","program22"),cols = c("black","red"),max.cutoff = "q90",min.cutoff = "q10",ncol = 3)

FeaturePlot(filnoep,features = c("program11"),cols = c("black","red"),max.cutoff = "q90",min.cutoff = "q10")
FeaturePlot(filnoep,features = c("program13"),cols = c("black","red"),max.cutoff = "q90",min.cutoff = "q10")
FeaturePlot(filnoep,features = c("program22"),cols = c("black","red"),max.cutoff = "q90",min.cutoff = "q10")
```




```{r}
#neu
FeaturePlot(nffil,"NLRP3",split.by = "tumorsite",min.cutoff = "q10",max.cutoff = "q90",cols = c("grey","red"))
VlnPlot(nffil,features = "NLRP3",group.by = "tumorsite",pt.size = 0)+stat_compare_means(comparisons = list(c("metastasis_normal","primary_tumor")),label = "p.signif")+ylim(0,6)+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

VlnPlot(nffil,features = "NLRP3",group.by = "tumorsite",pt.size = 0)+stat_compare_means(comparisons = list(c("primary_tumor","metastasis_normal"),c("primary_tumor","metastasis_tumor"),c("primary_tumor","primary_normal")),label = "p.signif")+ylim(0,7)+theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))


```

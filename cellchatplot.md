```{r}
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
options(stringsAsFactors = F)
```
```{r}
data.dir <- '../cellchat2/'
dir.create(data.dir)
setwd(data.dir)

```

```{r}

metanormal<-subset(filnoep,subset=tumorsite=="metastasis_normal")
metatumor<-subset(filnoep,subset=tumorsite=="metastasis_tumor")
primarynormal<-subset(filnoep,subset=tumorsite=="primary_normal")
primarytumor<-subset(filnoep,subset=tumorsite=="primary_tumor")
```

```{r}
mtinput<-metatumor@assays$RNA@data
mt=data.frame(group=metatumor$celltype20,row.names = names(metatumor$celltype20))
unique(mt$group)
mtmd<-metatumor@meta.data
```

```{r}
mtcc<-createCellChat(object = mtinput,meta = mtmd,group.by = "celltype20")
mtcc<-setIdent(mtcc,ident.use = "celltype20")
levels(mtcc@idents)
groupSize<-as.numeric(table(mtcc@idents))

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
mtcc@DB <- CellChatDB

mtcc<-subsetData(mtcc)
future::plan("multiprocess",workers=10)
mtcc<-identifyOverExpressedGenes(mtcc)
mtcc <- identifyOverExpressedInteractions(mtcc)
mtcc <- projectData(mtcc, PPI.human)

library(future)
plan("multisession",workers=6)
options(future.globals.maxSize=6000*1024^2)
mtcc <- computeCommunProb(mtcc)
mtcc <- computeCommunProbPathway(mtcc)
mtcc <- aggregateNet(mtcc)


mtcc<-netAnalysis_computeCentrality(mtcc,slot.name = "netP")

saveRDS(mtcc, file = "mtcc.rds")
rm(mtinput,mtmd)
```
```{r}

mninput<-metanormal@assays$RNA@data
mn=data.frame(group=metanormal$celltype20,row.names = names(metanormal$celltype20))
unique(mn$group)
mnmd<-metanormal@meta.data


mncc<-createCellChat(object = mninput,meta = mnmd,group.by = "celltype20")
mncc<-setIdent(mncc,ident.use = "celltype20")
levels(mncc@idents)
groupSize<-as.numeric(table(mncc@idents))

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
mncc@DB <- CellChatDB

mncc<-subsetData(mncc)
future::plan("multiprocess",workers=10)
mncc<-identifyOverExpressedGenes(mncc)
mncc <- identifyOverExpressedInteractions(mncc)
mncc <- projectData(mncc, PPI.human)

library(future)
plan("multisession",workers=6)
options(future.globals.maxSize=6000*1024^2)
mncc <- computeCommunProb(mncc)
mncc <- computeCommunProbPathway(mncc)
mncc <- aggregateNet(mncc)


mncc<-netAnalysis_computeCentrality(mncc,slot.name = "netP")

saveRDS(mncc, file = "mncc.rds")
```

```{r}

ptinput<-primarytumor@assays$RNA@data
pt=data.frame(group=primarytumor$celltype20,row.names = names(primarytumor$celltype20))
unique(pt$group)
ptmd<-primarytumor@meta.data


ptcc<-createCellChat(object = ptinput,meta = ptmd,group.by = "celltype20")
ptcc<-setIdent(ptcc,ident.use = "celltype20")
levels(ptcc@idents)
groupSize<-as.numeric(table(ptcc@idents))

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
ptcc@DB <- CellChatDB

ptcc<-subsetData(ptcc)
future::plan("multiprocess",workers=10)
ptcc<-identifyOverExpressedGenes(ptcc)
ptcc <- identifyOverExpressedInteractions(ptcc)
ptcc <- projectData(ptcc, PPI.human)

library(future)
plan("multisession",workers=6)
options(future.globals.maxSize=6000*1024^2)
ptcc <- computeCommunProb(ptcc)
ptcc <- computeCommunProbPathway(ptcc)
ptcc <- aggregateNet(ptcc)


ptcc<-netAnalysis_computeCentrality(ptcc,slot.name = "netP")

saveRDS(ptcc, file = "ptcc.rds")
```

```{r}

pninput<-primarynormal@assays$RNA@data
pn=data.frame(group=primarynormal$celltype20,row.names = names(primarynormal$celltype20))
unique(pn$group)
pnmd<-primarynormal@meta.data


pncc<-createCellChat(object = pninput,meta = pnmd,group.by = "celltype20")
pncc<-setIdent(pncc,ident.use = "celltype20")
levels(pncc@idents)
groupSize<-as.numeric(table(pncc@idents))

CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
pncc@DB <- CellChatDB

pncc<-subsetData(pncc)
future::plan("multiprocess",workers=10)
pncc<-identifyOverExpressedGenes(pncc)
pncc <- identifyOverExpressedInteractions(pncc)
pncc <- projectData(pncc, PPI.human)

library(future)
plan("multisession",workers=6)
opnions(future.globals.maxSize=6000*1024^2)
pncc <- computeCommunProb(pncc)
pncc <- computeCommunProbPathway(pncc)
pncc <- aggregateNet(pncc)


pncc<-netAnalysis_computeCentrality(pncc,slot.name = "netP")

saveRDS(pncc, file = "pncc.rds")
```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt

ptvspn<-mergeCellChat(list(ptcc,pncc),add.names = c("primary_tumor","primary_normal"))
rankNet(ptvspn,mode = "comparison",do.stat = T,color.use = c("#ff7f0e","#1f77b4"))
rankNet(ptvspn,mode = "comparison",stacked = T,do.stat = T,color.use = c("#ff7f0e","#1f77b4"))

mtvsmn<-mergeCellChat(list(mtcc,mncc),add.names = c("metastasis_tumor","metastasis_normal"))
rankNet(mtvsmn,mode = "comparison",do.stat = T,color.use = c("#d62728","#2ca02c"))
rankNet(mtvsmn,mode = "comparison",stacked = T,do.stat = T,color.use = c("#d62728","#2ca02c"))
```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt

mtvspt<-mergeCellChat(list(mtcc,ptcc),add.names = c("metastasis_tumor","primary_tumor"))
rankNet(mtvspt,mode = "comparison",do.stat = T,color.use = c("#d62728","#ff7f0e"))
rankNet(mtvspt,mode = "comparison",stacked = T,do.stat = T,color.use = c("#d62728","#ff7f0e"))

mnvspn<-mergeCellChat(list(mncc,pncc),add.names = c("metastasis_normal","primary_normal"))
rankNet(mnvspn,mode = "comparison",do.stat = T,color.use = c("#2ca02c","#1f77b4"))
rankNet(mnvspn,mode = "comparison",stacked = T,do.stat = T,color.use = c("#2ca02c","#1f77b4"))
```

```{r}
save(mnvspn,mtvsmn,mtvspt,ptvspn,file = "comparison.RData")

```





```{r}
list<-list(mtcc,mncc,ptcc,pncc)
weight.max <- getMaxWeight(list, slot.name = c("netP"), attribute = "FN1") # control the edge weights across different datasets
for (i in 1:length(list)) {
  netVisual_aggregate(list[[i]], signaling = "FN1", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste("FN1", names(list)[i]))
}
```

```{r}
list<-list(mtcc,mncc,ptcc,pncc)
weight.max <- getMaxWeight(list, slot.name = c("netP"), attribute = "IFN-II") # control the edge weights across different datasets
for (i in 1:length(list)) {
  netVisual_aggregate(list[[i]], signaling = "IFN-II", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste("IFN-II", names(list)[i]))
}
```

```{r}
list<-list(mtcc,mncc,ptcc,pncc)
weight.max <- getMaxWeight(list, slot.name = c("netP"), attribute = "CCL") # control the edge weights across different datasets
for (i in 1:length(list)) {
  netVisual_aggregate(list[[i]], signaling = "CCL", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste("CCL", names(list)[i]))
}
```


```{r}
list<-list(mtcc,mncc,ptcc,pncc)
weight.max <- getMaxWeight(list, slot.name = c("netP"), attribute = "CXCL") # control the edge weights across different datasets
for (i in 1:length(list)) {
  netVisual_aggregate(list[[i]], signaling = "CXCL", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste("CXCL", names(list)[i]))
}
```


```{R}
netAnalysis_signalingRole_scatter(mtcc)
netAnalysis_signalingRole_scatter(mncc)
netAnalysis_signalingRole_scatter(ptcc)
netAnalysis_signalingRole_scatter(pncc)

```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt
netAnalysis_signalingChanges_scatter(mtvsmn,idents.use = "Mesothelial",color.use = c("grey10","#d62728","#2ca02c"),point.shape = c(15,16,17,18))

```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt
netAnalysis_signalingChanges_scatter(mtvsmn,idents.use = "MSC/ASC",color.use = c("grey10","#d62728","#2ca02c"),point.shape = c(15,16,17,18))

```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt
netAnalysis_signalingChanges_scatter(mtvsmn,idents.use = "Fibroblast",color.use = c("grey10","#d62728","#2ca02c"),point.shape = c(15,16,17,18))

```

```{r}
netVisual_diffInteraction(mtvsmn,measure = "weight",sources.use = 4,targets.use = c(5,6),comparison = c(1,2))

```
```{R}
netAnalysis_diff_signalingRole_scatter(mtvsmn)
netAnalysis_diff_signalingRole_scatter(mtvspt)
netAnalysis_diff_signalingRole_scatter(mnvspn)
netAnalysis_diff_signalingRole_scatter(ptvspn)
```
```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt
table(mncc@idents)
netVisual_bubble(ptvspn, sources.use = 17, targets.use = c(2,3,4),  comparison = c(1, 2), angle.x = 90,color.text = c("#ff7f0e","#1f77b4"),color.heatmap = "viridis")
```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt
table(mncc@idents)
netVisual_bubble(ptvspn, sources.use = 11, targets.use = c(2,3,4,8,9,14,15,16,18),  comparison = c(1, 2), angle.x = 90,color.text = c("#ff7f0e","#1f77b4"),color.heatmap = "viridis")
```

```{r}
#color=c("#1f77b4","#ff7f0e","#2ca02c","#d62728")
#pn pt mn mt
table(mncc@idents)
netVisual_bubble(mtvsmn, sources.use = 13, targets.use = c(8,14),  comparison = c(2, 1), angle.x = 90,color.text = c("#2ca02c","#d62728"),color.heatmap = "viridis")
```

```{r}

immune<-subset(metadata,subset=maincelltype!="Epithellial cell")
immune<-subset(immune,subset=maincelltype!="Stromal cell")
immune2<-immune[,c("patient","tumorsite","maincelltype","xicluster")]
immune3<-immune2 %>% group_by(patient,tumorsite,xicluster) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup() %>% select(patient,tumorsite,xicluster,proportion) 

immune3_2<-immune2 %>% group_by(patient,tumorsite,xicluster) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup()%>% select(patient,tumorsite,xicluster,count)




immune4<-subset(immune3,subset=tumorsite!="metastasis_normal")
immune4<-subset(immune4,subset=tumorsite!="primary_normal")

immune4_2<-subset(immune3_2,subset=tumorsite!="metastasis_normal")
immune4_2<-subset(immune4_2,subset=tumorsite!="primary_normal")
immune5_2<-tidyr::spread(immune4_2,xicluster,count)
immune5_2[is.na(immune5_2)]<-0
immune6_2<-reshape2::melt(immune5_2,id.vars=c("patient","tumorsite"),variable.name="xicluster",value.name ="count")
immune7_2<-immune6_2%>%group_by(tumorsite,xicluster)%>% mutate(Mean=mean(count))%>%ungroup()
immune8_2<-immune7_2[,c(2,3,5)]
immune8_2<-unique(immune8_2)
immune8_3<-tidyr::spread(immune8_2,tumorsite,Mean)
immune8_3<-immune8_3 %>%  mutate(Enrichment=(metastasis_tumor-primary_tumor)/sqrt(primary_tumor)) %>% ungroup()
immune8_4<-immune8_3[,c(1,4)]
immune8_4$primary_tumor<-0
colnames(immune8_4)<-c("xicluster","metastasis_tumor","primary_tumor")
immune9_2<-reshape2::melt(immune8_4,id.vars="xicluster",variable.name="tumorsite",value.name = "Enrichment")

immune5<-tidyr::spread(immune4,xicluster,proportion)
immune5[is.na(immune5)]<-0

immune6<-reshape2::melt(immune5,id.vars=c("patient","tumorsite"),variable.name="xicluster",value.name ="proportion")

immune7<-immune6%>%group_by(tumorsite,xicluster)%>% mutate(Mean=mean(proportion))%>%ungroup()
immune7<-immune7%>%group_by(tumorsite,xicluster)%>% mutate(Median=median(proportion))%>%ungroup()
immune8<-immune7[,c(2,3,5,6)]
immune8<-unique(immune8)


immune9_2$xicluster<-factor(immune9_2$xicluster,levels =c("IEL_GDT","CD8_TEFF",'CD8_Tcycling',"CD8_TEX","CD8_TPEX","CD8_TRM","CD4_TH17","CD4_TREG","CD4_TNAIVE","CD4_TFH","NK","ILC","GC_B","Follicular_B","Memory_B","Breg","IGHA_PC","IGHM_PC","IGHG4_PC","IGHG1G3_PC","IGHG2_PC","Neu-IL1B-1","Neu-IL1B-2","Neu-DUSP6","Neu-MMP9-1","Neu-MMP9-2","Neu-HSPA1A","Neu-ISG15-1","Neu-ISG15-2","Macrophage-FN1","Macrophage-IL1B","Macrophage-SELENOP","Macrophage-TREM2","Mono-CD14","Mono-CD16","cDC1","cDC2","cDC-Mature","monolike-FCN1","pDC","Mast-EGR1","Mast-NFKB1"))

immune8$xicluster<-factor(immune8$xicluster,levels =c("IEL_GDT","CD8_TEFF",'CD8_Tcycling',"CD8_TEX","CD8_TPEX","CD8_TRM","CD4_TH17","CD4_TREG","CD4_TNAIVE","CD4_TFH","NK","ILC","GC_B","Follicular_B","Memory_B","Breg","IGHA_PC","IGHM_PC","IGHG4_PC","IGHG1G3_PC","IGHG2_PC","Neu-IL1B-1","Neu-IL1B-2","Neu-DUSP6","Neu-MMP9-1","Neu-MMP9-2","Neu-HSPA1A","Neu-ISG15-1","Neu-ISG15-2","Macrophage-FN1","Macrophage-IL1B","Macrophage-SELENOP","Macrophage-TREM2","Mono-CD14","Mono-CD16","cDC1","cDC2","cDC-Mature","monolike-FCN1","pDC","Mast-EGR1","Mast-NFKB1"))
```



```{r}
metapro<-immune5[immune5$tumorsite=="metastasis_tumor","IEL_GDT"]
pripro<-immune5[immune5$tumorsite=="primary_tumor","IEL_GDT"]
result<-wilcox.test(metapro$IEL_GDT,pripro$IEL_GDT,paired = T)

```

```{r}
metapro<-immune5[immune5$tumorsite=="metastasis_tumor","IGHA_PC"]
pripro<-immune5[immune5$tumorsite=="primary_tumor","IGHA_PC"]

result<-wilcox.test(metapro$IGHA_PC,pripro$IGHA_PC,paired = T)
result$p.value
```

```{r}
metapro<-immune5[immune5$tumorsite=="metastasis_tumor","NK"]
pripro<-immune5[immune5$tumorsite=="primary_tumor","NK"]
  
result<-wilcox.test(metapro$NK,pripro$NK,paired = T)
result$p.value

p<-immune5%>%wilcox.test(formula=.data[[colnames(immune5)[3]]]~tumorsite)
```


```{r}
d=list()
for (a in 3:44){
  
 var=colnames(immune5)[a]
 p<-wilcox.test(immune5[[var]]~tumorsite,data = immune5,paired=T)
 d[[a-2]]=p$p.value
  
}
d=do.call(cbind,d)  
colnames(d)<-colnames(immune5)[3:44]


value=round(p.adjust(d[1,],"BH"),4 )
d<-rbind(d,value)

```

```{r}
p<-ggplot()+geom_tile(data=immune9_2,aes(xicluster,tumorsite,fill=Enrichment),colour="white",size=1)+scale_fill_gradientn(colours = gplots::bluered(128),limit=c(-50,50),breaks=c(-50,0,50),labels=c("-50","0","50"),name="Enrichment \nMetastasis vs Primary")+geom_point(data=immune8,aes(xicluster,tumorsite,size=Mean),shape=1)+scale_size_area(breaks=c(0.0001,0.001,0.01,0.03,0.05,0.07,0.1),labels=c(0.0001,0.001,0.01,0.03,0.05,0.07,">0.1"),name="among\nall immune cells")+
  labs(x="",y="")+scale_x_discrete(position = "top")+theme_bw()+theme(panel.grid.major = element_blank(),axis.text.x = element_text(angle = 90),legend.position = "bottom",legend.direction = "vertical",legend.title = element_text(angle=90),legend.title.align = 0.5,legend.box.just="left")

```

```{r}

#immune<-subset(metadata,subset=maincelltype!="Epithellial cell")
#immune<-subset(immune,subset=maincelltype!="Stromal cell")
immune2<-immune[,c("patient","sample","tumorsite","ctypes")]
immune3<-immune2 %>% group_by(patient,tumorsite,ctypes) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup() %>% select(patient,tumorsite,ctypes,proportion) 



immune3_2<-immune2 %>% group_by(patient,tumorsite,ctypes) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup()%>% select(patient,tumorsite,ctypes,count)




immune4<-subset(immune3,subset=tumorsite!="metastasis_normal")
immune4<-subset(immune4,subset=tumorsite!="primary_normal")

immune4_2<-subset(immune3_2,subset=tumorsite!="metastasis_normal")
immune4_2<-subset(immune4_2,subset=tumorsite!="primary_normal")

immune5_2<-tidyr::spread(immune4_2,ctypes,count)
immune5_2[is.na(immune5_2)]<-0
immune6_2<-reshape2::melt(immune5_2,id.vars=c("patient","tumorsite"),variable.name="ctypes",value.name ="count")
immune7_2<-immune6_2%>%group_by(tumorsite,ctypes)%>% mutate(Mean=mean(count))%>%ungroup()
immune8_2<-immune7_2[,c(2,3,5)]
immune8_2<-unique(immune8_2)
immune8_3<-tidyr::spread(immune8_2,tumorsite,Mean)
immune8_3<-immune8_3 %>%  mutate(Enrichment=(metastasis_tumor-primary_tumor)/sqrt(primary_tumor)) %>% ungroup()
immune8_4<-immune8_3[,c(1,4)]
immune8_4$primary_tumor<-0
colnames(immune8_4)<-c("ctypes","metastasis_tumor","primary_tumor")
immune9_2<-reshape2::melt(immune8_4,id.vars="ctypes",variable.name="tumorsite",value.name = "Enrichment")

immune5<-tidyr::spread(immune4,ctypes,proportion)
immune5[is.na(immune5)]<-0

immune6<-reshape2::melt(immune5,id.vars=c("patient","tumorsite"),variable.name="ctypes",value.name ="proportion")

immune7<-immune6%>%group_by(tumorsite,ctypes)%>% mutate(Mean=mean(proportion))%>%ungroup()
immune7<-immune7%>%group_by(tumorsite,ctypes)%>% mutate(Median=median(proportion))%>%ungroup()
immune8<-immune7[,c(2,3,5,6)]
immune8<-unique(immune8)


immune9_2$ctypes<-factor(immune9_2$ctypes,levels =c("IEL_GDT","CD8T",'CD4T',"NK","ILC","B cell","Plasma cell","Neutrophil","Mono/Macro","cDC","pDC","Mast cell"))

immune8$ctypes<-factor(immune8$ctypes,levels =c("IEL_GDT","CD8T",'CD4T',"NK","ILC","B cell","Plasma cell","Neutrophil","Mono/Macro","cDC","pDC","Mast cell"))
```

```{r}
d=list()
for (a in 3:14){
  
 var=colnames(immune5)[a]
 p<-wilcox.test(immune5[[var]]~tumorsite,data = immune5,paired=T)
 d[[a-2]]=p$p.value
  
}
d=do.call(cbind,d)  
colnames(d)<-colnames(immune5)[3:14]
value=round(p.adjust(d[1,],"BH"),4 )
d<-rbind(d,value)
```

```{r}
p<-ggplot()+geom_tile(data=immune9_2,aes(ctypes,tumorsite,fill=Enrichment),colour="white",size=1)+scale_fill_gradientn(colours = gplots::bluered(128),limit=c(-50,50),breaks=c(-50,0,50),labels=c("-50","0","50"),name="Enrichment \nMetastasis vs Primary")+geom_point(data=immune8,aes(ctypes,tumorsite,size=Mean),shape=1)+scale_size_area(breaks=c(0.002,0.001,0.01,0.05,0.1,0.2,0.3),labels=c(0.002,0.001,0.01,0.05,0.1,0.2,">0.3"),name="among\nall immune cells")+
  labs(x="",y="")+scale_x_discrete(position = "top")+theme_bw()+theme(panel.grid.major = element_blank(),axis.text.x = element_text(angle = 90),legend.position = "bottom",legend.direction = "vertical",legend.title = element_text(angle=90),legend.title.align = 0.5,legend.box.just="left")

```


```{r}

meta2<-metadata[,c("patient","sample","tumorsite","celltype")]
meta3<-meta2 %>% group_by(patient,tumorsite,celltype) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup() %>% select(patient,tumorsite,celltype,proportion) 
meta3_2<-meta2 %>% group_by(patient,tumorsite,celltype) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup()%>% select(patient,tumorsite,celltype,count)
meta4<-subset(meta3,subset=tumorsite!="metastasis_normal")
meta4<-subset(meta4,subset=tumorsite!="primary_normal")
meta4_2<-subset(meta3_2,subset=tumorsite!="metastasis_normal")
meta4_2<-subset(meta4_2,subset=tumorsite!="primary_normal")
meta5_2<-tidyr::spread(meta4_2,celltype,count)
meta5_2[is.na(meta5_2)]<-0
meta6_2<-reshape2::melt(meta5_2,id.vars=c("patient","tumorsite"),variable.name="celltype",value.name ="count")
meta7_2<-meta6_2%>%group_by(tumorsite,celltype)%>% mutate(Mean=mean(count))%>%ungroup()
meta8_2<-meta7_2[,c(2,3,5)]
meta8_2<-unique(meta8_2)
meta8_3<-tidyr::spread(meta8_2,tumorsite,Mean)
meta8_3<-meta8_3 %>%  mutate(Enrichment=(metastasis_tumor-primary_tumor)/sqrt(primary_tumor)) %>% ungroup()
meta8_4<-meta8_3[,c(1,4)]
meta8_4$primary_tumor<-0
colnames(meta8_4)<-c("celltype","metastasis_tumor","primary_tumor")
meta9_2<-reshape2::melt(meta8_4,id.vars="celltype",variable.name="tumorsite",value.name = "Enrichment")
meta5<-tidyr::spread(meta4,celltype,proportion)
meta5[is.na(meta5)]<-0
meta6<-reshape2::melt(meta5,id.vars=c("patient","tumorsite"),variable.name="celltype",value.name ="proportion")
meta7<-meta6%>%group_by(tumorsite,celltype)%>% mutate(Mean=mean(proportion))%>%ungroup()
meta7<-meta7%>%group_by(tumorsite,celltype)%>% mutate(Median=median(proportion))%>%ungroup()
meta8<-meta7[,c(2,3,5,6)]
meta8<-unique(meta8)


meta9_2$celltype<-factor(meta9_2$celltype,levels =c("IEL_GDT","CD8T",'CD4T',"NK","ILC","B cell","Plasma cell","Neutrophil","Mono/Macro","cDC","pDC","Mast cell"))

meta8$celltype<-factor(meta8$celltype,levels =c("IEL_GDT","CD8T",'CD4T',"NK","ILC","B cell","Plasma cell","Neutrophil","Mono/Macro","cDC","pDC","Mast cell"))
```


```{r}
d=list()
for (a in 3:5){
  
 var=colnames(meta5)[a]
 p<-wilcox.test(meta5[[var]]~tumorsite,data = meta5,paired=T)
 d[[a-2]]=p$p.value
  
}
d=do.call(cbind,d)  
colnames(d)<-colnames(meta5)[3:5]


value=round(p.adjust(d[1,],"BH"),4 )
d<-rbind(d,value)
for (i in 1:6){if (meta8[i,"tumorsite"]=="metastasis_tumor"){meta8[i,"Mean"]<-(-meta8[i,"Mean"])}}
p2<-ggplot(data=meta8, aes(x=maincelltype,y=Mean,fill=tumorsite))+geom_bar(stat="identity",position="identity",color="black",size=0.25)+theme_light()
```

```{r}
p<-ggplot()+geom_tile(data=meta9_2,aes(celltype,tumorsite,fill=Enrichment),colour="white",size=1)+scale_fill_gradientn(colours = gplots::bluered(128),limit=c(-50,50),breaks=c(-50,0,50),labels=c("-50","0","50"),name="Enrichment \nMetastasis vs Primary")+geom_point(data=meta8,aes(celltype,tumorsite,size=Mean),shape=1)+scale_size_area(breaks=c(0.05,0.1,0.2,0.3,0.4,0.5,0.6),labels=c(0.05,0.1,0.2,0.3,0.4,0.5,">0.6"),name="among\nall meta cells")+
  labs(x="",y="")+scale_x_discrete(position = "top")+theme_bw()+theme(panel.grid.major = element_blank(),axis.text.x = element_text(angle = 90),legend.position = "bottom",legend.direction = "vertical",legend.title = element_text(angle=90),legend.title.align = 0.5,legend.box.just="left")

save(metadata,meta9_2,meta8,d,p,file="D:/PROJECT/pm/3/新建文件夹/celltype.RData")
```
```{r}
meta5<-tidyr::spread(meta3,celltype,proportion)
meta5[is.na(meta5)]<-0
ggplot(meta5,aes(tumorsite,`Epithellial cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`Stromal cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`Immune cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
```
```{r}


meta<-subset(metadata,subset=maincelltype!="Epithellial cell")
meta<-subset(meta,subset=maincelltype!="Stromal cell")
meta2<-meta[,c("patient","sample","tumorsite","maincelltype")]
meta3<-meta2 %>% group_by(patient,tumorsite,maincelltype) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup() %>% select(patient,tumorsite,maincelltype,proportion) 



meta3_2<-meta2 %>% group_by(patient,tumorsite,maincelltype) %>% summarize(count=n()) %>% group_by(patient,tumorsite) %>% mutate(total_count=sum(count),proportion=count/total_count) %>% ungroup()%>% select(patient,tumorsite,maincelltype,count)




meta4<-subset(meta3,subset=tumorsite!="metastasis_normal")
meta4<-subset(meta4,subset=tumorsite!="primary_normal")

meta4_2<-subset(meta3_2,subset=tumorsite!="metastasis_normal")
meta4_2<-subset(meta4_2,subset=tumorsite!="primary_normal")

meta5_2<-tidyr::spread(meta4_2,maincelltype,count)
meta5_2[is.na(meta5_2)]<-0
meta6_2<-reshape2::melt(meta5_2,id.vars=c("patient","tumorsite"),variable.name="maincelltype",value.name ="count")
meta7_2<-meta6_2%>%group_by(tumorsite,maincelltype)%>% mutate(Mean=mean(count))%>%ungroup()
meta8_2<-meta7_2[,c(2,3,5)]
meta8_2<-unique(meta8_2)
meta8_3<-tidyr::spread(meta8_2,tumorsite,Mean)
meta8_3<-meta8_3 %>%  mutate(Enrichment=(metastasis_tumor-primary_tumor)/sqrt(primary_tumor)) %>% ungroup()
meta8_4<-meta8_3[,c(1,4)]
meta8_4$primary_tumor<-0
colnames(meta8_4)<-c("maincelltype","metastasis_tumor","primary_tumor")
meta9_2<-reshape2::melt(meta8_4,id.vars="maincelltype",variable.name="tumorsite",value.name = "Enrichment")

meta5<-tidyr::spread(meta4,maincelltype,proportion)
meta5[is.na(meta5)]<-0

meta6<-reshape2::melt(meta5,id.vars=c("patient","tumorsite"),variable.name="maincelltype",value.name ="proportion")

meta7<-meta6%>%group_by(tumorsite,maincelltype)%>% mutate(Mean=mean(proportion))%>%ungroup()
meta7<-meta7%>%group_by(tumorsite,maincelltype)%>% mutate(Median=median(proportion))%>%ungroup()
meta8<-meta7[,c(2,3,5,6)]
meta8<-unique(meta8)

for (i in 1:12){if (meta8[i,"tumorsite"]=="metastasis_tumor"){meta8[i,"Mean"]<-(-meta8[i,"Mean"])}}

```
```{r,calculate fdr}
d=list()
for (a in 3:8){
  
 var=colnames(meta5)[a]
 p<-wilcox.test(meta5[[var]]~tumorsite,data = meta5,paired=T)
 d[[a-2]]=p$p.value
  
}
d=do.call(cbind,d)  
colnames(d)<-colnames(meta5)[3:8]


value=round(p.adjust(d[1,],"BH"),4 )
d<-rbind(d,value)

```


```{r,boxplot}

meta5<-tidyr::spread(meta3,celltype,proportion)
meta5[is.na(meta5)]<-
  
ggplot(meta5,aes(tumorsite,`B cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`Mast cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`Myeloid cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`Neutrophil cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`Plasma cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())
ggplot(meta5,aes(tumorsite,`TNK cell`))+geom_boxplot(aes(fill=tumorsite),notch = FALSE)+geom_jitter(position = "jitter",size=1)+scale_fill_manual(values=c(brewer.pal(7,"Set2")[c(1,2,3,4)]))+theme(panel.grid = element_blank())

```


```{r}
d=list()
for (a in 3:16){
  
 var=colnames(meta5)[a]
 p<-wilcox.test(meta5[[var]]~tumorsite,data = meta5,paired=T)
 d[[a-2]]=p$p.value
  
}
d=do.call(cbind,d)  
colnames(d)<-colnames(meta5)[3:16]


value=round(p.adjust(d[1,],"BH"),4 )
d<-rbind(d,value)

```
```{r,proportion scater plot}



p2<-ggplot()+geom_tile(data=meta9_2,aes(factor(xicluster,levels = c("Pericyte-RGS5","SMC-DES","SMC-RERGL","Fibroblast-INHBA","Fibroblast-ACTA2","Fibroblast-APOE","Fibroblast-BMP4","Endo-S100A8","Endo-SELP","Endo-DLL4","Endo-LYVE1","MSC/ASC-NT5E","Mesothelial-MSLN","Schwann-S100B")),tumorsite,fill=Enrichment),colour="white",size=1)+scale_fill_gradientn(colours = gplots::bluered(128),limit=c(-40,40),breaks=c(-40,0,40),labels=c("-40","0","40"),name="Enrichment \nmetastasis_tumor vs primary_tumor")+geom_point(data=meta8_2,aes(factor(xicluster,levels = c("Pericyte-RGS5","SMC-DES","SMC-RERGL","Fibroblast-INHBA","Fibroblast-ACTA2","Fibroblast-APOE","Fibroblast-BMP4","Endo-S100A8","Endo-SELP","Endo-DLL4","Endo-LYVE1","MSC/ASC-NT5E","Mesothelial-MSLN","Schwann-S100B")),tumorsite,size=Mean),shape=1)+scale_size_area(breaks=c(0.0001,0.001,0.01,0.05,0.1,0.2),labels=c(0.0001,0.001,0.01,0.05,0.1,">0.2"),name="among\nall Stromal cells")+
  labs(x="",y="")+scale_x_discrete(position = "top")+theme_bw()+theme(panel.grid.major = element_blank(),axis.text.x = element_text(angle = 90),legend.position = "bottom",legend.direction = "vertical",legend.title = element_text(angle=90),legend.title.align = 0.5,legend.box.just="left")

```

```{r}
p<-ggplot(data=meta8, aes(x=xicluster,y=Mean,fill=tumorsite))+geom_bar(stat="identity",position="identity",color="black",size=0.25)+theme(axis.text.x =  element_text(hjust = 1,angle = 90,vjust = .5))
```

```{r}
strmtvspt<-subset(meta5,subset=tumorsite!="metastasis_normal")
strmtvspt<-subset(strmtvspt,subset=tumorsite!="primary_normal")

d=list()
for (a in 3:9){
  
 var=colnames(strmtvspt)[a]
 p<-wilcox.test(strmtvspt[[var]]~tumorsite,data = strmtvspt,paired=T)
 d[[a-2]]=p$p.value
  
}
d=do.call(cbind,d)  
colnames(d)<-colnames(strmtvspt)[3:9]


value=round(p.adjust(d[1,],"BH"),4 )
d<-rbind(d,value)
```

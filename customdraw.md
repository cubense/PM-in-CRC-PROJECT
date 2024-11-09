```{r,cnmf plot}
usage <- read.delim("../cnmf/tumor2kgene/usage.txt")
rownames(usage)<-usage$X
usage<-usage[,-1]
metadata<-tumor@meta.data
identical(rownames(usage),rownames(metadata))



#minmaxscale<-function(val){return((val-min(val,na.rm = T))/(max(val,na.rm = T)-min(val,na.rm = T)))}
#usagescale<-as.data.frame(lapply(usage[1:42],minmaxscale))
#rownames(usagescale)<-rownames(usage)
meta<-cbind(metadata,usage)
tumor@meta.data<-meta
```

```{r,plot tumor program in umap}

FeaturePlot(tumor,features = c("program4","program25","program5","program38","program20","program42","program36","program28","program3","program19","program39","program35"),ncol = 4,min.cutoff = "q10",max.cutoff = "q90",cols = c("black","red"))
FeaturePlot(tumor,features = c("program4","program25","program5","program38","program20","program42","program36","program28","program3","program19","program39","program35"),ncol = 4,min.cutoff = "q10",max.cutoff = "q90",cols = c("yellow","red"))


FeaturePlot(tumor,features = c("program4","program25","program5","program38","program20","program42","program36","program28","program3","program19","program39","program35"),ncol = 4,min.cutoff = "q10",max.cutoff = "q90",cols = c("black","red"))

FeaturePlot(tumor,features = c("program4","program25","program5","program38","program20","program42","program36","program28","program3","program19","program39","program35"),ncol = 3,min.cutoff = "q10",max.cutoff = "q90",cols = c("yellow","red"))
```

```{r,plot cellcycle programs}
VlnPlot(tumor,features = c("program37"),group.by = "Phase",pt.size = 0,log = T,cols =c("#E64B357F", "#4DBBD57F", "#00A0877F"))+theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))
VlnPlot(tumor,features = c("program10"),group.by = "Phase",pt.size = 0,log = T,cols =c("#E64B357F", "#4DBBD57F", "#00A0877F"))+theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

```

```{r}
geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }

GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    
    draw_key = draw_key_polygon,
    
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    
    required_aes = c("x", "y")
  )

```


```{r,public data}
#cellepi
load("../cellmmrd/scoremeta.RData")

library(ggridges)
library(RColorBrewer)

ggplot(meta,aes(x="program11",y="epicluster",fill=..density..))+geom_density_ridges_gradient(scale=3,rel_min_height=0.00,size=0.3)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))

ggplot(meta,aes(x=program11,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-10,40))

ggplot(meta,aes(x=program15,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-0.8,1.5))

ggplot(meta,aes(x=program27,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-5,20))


meta$group<-gsub("T","Tumor",meta$SPECIMEN_TYPE)
meta$group<-gsub("N","Normal",meta$group)





#ngepi
load("D:/PROJECT/pm/ng/meta.RData")

ggplot(meta,aes(x=program11,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-5,40))

ggplot(meta,aes(x=program15,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-1,2))

ggplot(meta,aes(x=program27,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-3,7))

ggplot(meta,aes(x=program40,y=epicluster,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-5,40))



ggplot(data=meta,aes_string(x = "malignant",y ="program11",fill="malignant"))+geom_boxplot(data = meta,aes_string(x = "malignant",y = "program11",fill="malignant"))+geom_jitter(position=position_jitter(0.3),size=0.001)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line())+stat_compare_means(aes_string(group="malignant"),label = "p.signif")+scale_fill_manual(values = c("light blue","#FE7280"))
```

```{r}
load("../allep.RData")
meta<-allep@meta.data
meta$alltype<-factor(meta$alltype,levels = c("BEST4/OTOP2 cell","Colonocyte","Enteroendocrine cell","Goblet cell","Tuft cell","Undifferentiated cell","Tumor cell"))

names(meta)[names(meta)=="program151"]<-"program15"
names(meta)[names(meta)=="program271"]<-"program27"
names(meta)[names(meta)=="program401"]<-"program40"

ggplot(meta,aes(x=program11,y=alltype,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-0.5,2))
ggplot(meta,aes(x=program15,y=alltype,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())
ggplot(meta,aes(x=program27,y=alltype,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-0.2,0.8))

ggplot(meta,aes(x=program40,y=alltype,fill=..density..))+geom_density_ridges_gradient(scale=2.5,rel_min_height=0.00,size=0.6)+scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(11,"Spectral")))(32))+theme_light()+theme(panel.grid = element_blank())+coord_cartesian(xlim = c(-0.2,0.8))


```
```{r,radviz plot}
library(Radviz)
library(wesanderson)
das<-c("program11","program15","program27","program40")
S<-make.S(das)
scaled<-apply(meta[,das],2,do.L)
rv<-do.radviz(meta,S)
sim.mat<-cosine(scaled)
new<-do.optim(S,sim.mat,iter=100,n=1000)
new.S<-make.S(get.optim(new))
new.rv<-do.radviz(meta,new.S)
plot(new.rv)+geom_point(aes(color=alltype),alpha=0.05,size=2)+geom_density2d(aes(color=alltype),linemitre = 1000)+scale_discrete_manual(values = c(mypal[1:6],"grey"),aesthetics = "color")

```

```{r}
pop.cols<-setNames(c(wes_palette(n=7,name="Darjeeling1")),levels(meta$alltype))
bubbleRadviz(new.rv,bubble.color=mypal[as.integer(meta$alltype)],group = "alltype",bubble.fg="black",scale=0.05,decreasing=TRUE)

smoothRadviz(new.rv)
plot(new.rv)+geom_point(aes(color=group),alpha=0.2,size=0.2)
plot(new.rv)+geom_point(aes(color=group),alpha=0.2,size=0.2)+scale_discrete_manual(values = c("#FF0000","#00A08A","#F2AD00"),aesthetics = "color")


```

```{r}
mypal2<-c(ggsci::pal_jama("default",alpha = 0.5)(6),ggsci::pal_lancet("lanonc",alpha = 0.5)(3),ggsci::pal_simpsons("springfield",alpha = 0.5)(1),ggsci::pal_futurama(alpha = 0.6)(1),ggsci::pal_igv("default",alpha = 0.6)(25),ggsci::pal_uchicago("default",alpha = 0.5)(8),ggsci::pal_d3("category10",alpha = 0.8)(8))
mypal2<-unique(mypal2)
show_col(mypal2)

ggplot(filnoep@meta.data,aes(x=sample, fill=sample)) + 
    geom_bar() +
    theme_classic() +
    geom_text(stat='count',aes(label=..count..),position = position_dodge(0.9), vjust=0.5, color="black", size=3.5)+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("After QC")+scale_fill_manual(values = mypal2)

ggplot(filnoep@meta.data,aes(sample,fill=maincelltype)) + 
  	geom_bar(position="fill") +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 0)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +scale_fill_manual(values = mypal)+scale_y_continuous(labels = percent,name = "proportion")+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
```


```{r}
DimPlot(filnoep,split.by = "meth",group.by = "maincelltype")+scale_color_manual(values = mypal)


filnoep$tissue<-factor(filnoep$tissue,levels = c("crc","n","pm","pn","lm","ln"))
DimPlot(filnoep,split.by = "tissue",group.by = "maincelltype",ncol = 2)+scale_color_manual(values = mypal)

```







```{r}
data<-meta[,c("tnk_cluster","Cytotoxicity","Exhausition")]
cd8data<-subset(data,subset = tnk_cluster!="CD4_TH17")
cd8data<-subset(cd8data,subset = tnk_cluster!="CD4_TREG")
cd8data<-subset(cd8data,subset = tnk_cluster!="CD4_TNAIVE")
cd8data<-subset(cd8data,subset = tnk_cluster!="CD4_TFH")

ggplot(cd8data, aes(x=Cytotoxicity, col=tnk_cluster)) + 
  stat_ecdf(geom="smooth", se=F, size=1) + 
  theme_bw()  +
  theme(panel.grid = element_blank()) +
  labs(x="Cytotoxicity score", 
       y="ECDF",
       col="")  +
  scale_color_nejm()


ggplot(cd8data, aes(x=Exhausition, col=tnk_cluster)) + 
  stat_ecdf(geom="smooth", se=F, size=1) + 
  theme_bw()  +
  theme(panel.grid = element_blank()) +
  labs(x="Exhausition score", 
       y="ECDF",
       col="")  +
  scale_color_nejm()
```

```{r,public data}
load("../mmrd.RData")
tnk<-subset(mmrd,subset = clTopLevel=="TNKILC")
DotPlot(tnk,features = c("CD160","KIR2DL4","TMIGD2","ZBTB16"),group.by = "cl295v11SubFull",cols = "RdBu")+ theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

meta<-tnk@meta.data
meta$sig<-meta$cl295v11SubFull
meta$sig[which(str_detect(meta$cl295v11SubFull, "^cTNI18"))] <- "IEL_GDT"
meta$sig[which(str_detect(meta$cl295v11SubFull, "^cTNI17"))] <- "IEL_GDT"
meta$sig[which(str_detect(meta$cl295v11SubFull, "^cTNI19"))] <- "IEL_GDT"
meta$sig[which(str_detect(meta$cl295v11SubFull, "^cTNI20"))] <- "IEL_GDT"
meta$sig[which(str_detect(meta$cl295v11SubFull, "^cTNI21"))] <- "IEL_GDT"

fq<-prop.table(table(meta$sig,meta[,"sample"]),2)*100
df<-reshape2::melt(fq,value.name = "freq",varnames=c("sig","sample"))
group<-meta[!duplicated(meta$sample),c("sample","SPECIMEN_TYPE")]
df$tissuegroup<-group[match(df$sample, group$sample),2]
df$tissue<-"Normal"
df$tissue<-gsub("T","Tumor",df$tissuegroup)

df$tissue<-gsub("N","Normal",df$tissue)


data<-df[which(df$sig=="IEL_GDT"),]


ggplot(data=data,aes_string(x = "tissue",y ="freq",fill="tissue"))+labs(y="Proportion(%)")+geom_boxplot(data = data,aes_string(x = "tissue",y = "freq",fill="tissue"))+geom_jitter(position=position_jitter(0.3),size=0.001)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line())+stat_compare_means(aes_string(group="tissue"),label = "p.signif")+scale_fill_manual(values = c("light blue","#FE7280"))


Idents(mmrd)<-"clTopLevel"
b<-subset(mmrd,idents = c("B","Plasma"))
meta<-b@meta.data
fq<-prop.table(table(meta$cl295v11SubFull,meta[,"sample"]),2)*100
df<-reshape2::melt(fq,value.name = "freq",varnames=c("cl295v11SubFull","sample"))
group<-meta[!duplicated(meta$sample),c("sample","SPECIMEN_TYPE")]
df$tissuegroup<-group[match(df$sample, group$sample),2]
df$tissue<-"Normal"
df$tissue<-gsub("T","Tumor",df$tissuegroup)

df$tissue<-gsub("N","Normal",df$tissue)
data<-df[which(df$cl295v11SubFull=="cP1 (Plasma IgA)"),]

ggplot(data=data,aes_string(x = "tissue",y ="freq",fill="tissue"))+labs(y="Proportion(%)")+geom_boxplot(data = data,aes_string(x = "tissue",y = "freq",fill="tissue"))+geom_jitter(position=position_jitter(0.3),size=0.001)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line())+stat_compare_means(aes_string(group="tissue"),label = "p.signif")+scale_fill_manual(values = c("light blue","#FE7280"))


load("../ng.RData")
t<-subset(ng,subset = Cell_type=="T cells")
DotPlot(t,features = c("CD160","KIR2DL4","TMIGD2","ZBTB16"),group.by = "Cell_subtype",cols = "RdBu")+ theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

meta<-t@meta.data
fq<-prop.table(table(meta$Cell_subtype,meta[,"Sample"]),2)*100
df<-reshape2::melt(fq,value.name = "freq",varnames=c("Cell_subtype","Sample"))
group<-meta[!duplicated(meta$Sample),c("Sample","malignant")]
df$tissue<-group[match(df$Sample, group$Sample),2]
data<-df[which(df$Cell_subtype=="gamma delta T cells"),]

ggplot(data=data,aes_string(x = "tissue",y ="freq",fill="tissue"))+labs(y="Proportion(%)")+geom_boxplot(data = data,aes_string(x = "tissue",y = "freq",fill="tissue"))+geom_jitter(position=position_jitter(0.3),size=0.001)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line())+stat_compare_means(aes_string(group="tissue"),label = "p.signif")+scale_fill_manual(values = c("light blue","#FE7280"))

b<-subset(ng,subset = Cell_type=="B cells")


meta<-b@meta.data
fq<-prop.table(table(meta$Cell_subtype,meta[,"Sample"]),2)*100
df<-reshape2::melt(fq,value.name = "freq",varnames=c("Cell_subtype","Sample"))
group<-meta[!duplicated(meta$Sample),c("Sample","malignant")]
df$tissue<-group[match(df$Sample, group$Sample),2]
data<-df[which(df$Cell_subtype=="IgA+ Plasma"),]

ggplot(data=data,aes_string(x = "tissue",y ="freq",fill="tissue"))+labs(y="Proportion(%)")+geom_boxplot(data = data,aes_string(x = "tissue",y = "freq",fill="tissue"))+geom_jitter(position=position_jitter(0.3),size=0.001)+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line())+stat_compare_means(aes_string(group="tissue"),label = "p.signif")+scale_fill_manual(values = c("light blue","#FE7280"))


```







```{r,macrophage gsva}
#mfgsva

library(GSVA)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(DO.db)
library(BiocParallel)

expr<-tmpenv$expr
sets<-tmpenv$gene
es = gsva(expr,sets,method = "gsva",parallel.sz=12)
save(es,file="../exprgsva.RData")

data<-data.table::fread("../es.txt",data.table = F)
table(mlfil$ml_nocluster)
mf<-subset(mlfil,idents =c("Macrophage-FN1","Macrophage-IL1B","Macrophage-SELENOP","Macrophage-TREM2"))


meta<-as.data.frame(mf@meta.data[,c("orig.ident","ml_nocluster","tumorsite")])
meta<-meta %>% arrange(meta$ml_nocluster)
rownames(data)<-data[,1]
data<-data[,-1]
data<-data[,rownames(meta)]
identical(colnames(data),rownames(meta))

table(meta$ml_nocluster)
data$`Macrophage-FN1`<-apply(data[,1:3079],1,mean)
data$`Macrophage-IL1B`<-apply(data[,3080:5133],1,mean)
data$`Macrophage-SELENOP`<-apply(data[,5134:10677],1,mean)
data$`Macrophage-TREM2`<-apply(data[,10678:18776],1,mean)

test<-data[,c("Macrophage-FN1","Macrophage-IL1B","Macrophage-SELENOP","Macrophage-TREM2")]
result<-t(scale(t(test)))
result[result>2]=2
result[result<-2]=-2

list<-c("HALLMARK_COMPLEMENT","REACTOME_PD_1_SIGNALING","HALLMARK_OXIDATIVE_PHOSPHORYLATION","REACTOME_VEGF_LIGAND_RECEPTOR_INTERACTIONS","REACTOME_CELL_SURFACE_INTERACTIONS_AT_THE_VASCULAR_WALL","REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION","HALLMARK_IL6_JAK_STAT3_SIGNALING","REACTOME_LTC4_CYSLTR_MEDIATED_IL4_PRODUCTION","NABA_COLLAGENS","NABA_ECM_REGULATORS","REACTOME_TCR_SIGNALING","REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES","KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY","REACTOME_SIGNALING_BY_CSF3_G_CSF","REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS","REACTOME_UPTAKE_AND_ACTIONS_OF_BACTERIAL_TOXINS","REACTOME_ANTIVIRAL_MECHANISM_BY_IFN_STIMULATED_GENES","REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE","HALLMARK_TGF_BETA_SIGNALING","REACTOME_RUNX2_REGULATES_GENES_INVOLVED_IN_CELL_MIGRATION","REACTOME_REGULATION_OF_IFNA_SIGNALING")

result2<-result[list,]
result2<-as.data.frame(result2)

library(pheatmap)
p<-pheatmap(result2,cluster_rows=T,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = T,
                color =colorRampPalette(c("light green", "white","red"))(100),
                cellwidth = 15, cellheight = 10,
                fontsize = 9)
```



```{r}

nffil$nf_cluster<-factor(nffil$nf_cluster,levels = c("Neu-IL1B-1","Neu-IL1B-2","Neu-DUSP6","Neu-MMP9-1","Neu-MMP9-2","Neu-HSPA1A","Neu-ISG15-1","Neu-ISG15-2"))
 
 DimPlot(nffil,split.by = "tumorsite",group.by = "nf_cluster")+scale_color_manual(values = c('#E64B357F','#F39B7F7F','#4DBBD57F','#00A0877F','#8491B47F','#FED4397F','#3C54887F','#709AE17F'))
 
 
load("../nfresult.RData")

list<-c("HALLMARK_INFLAMMATORY_RESPONSE","REACTOME_TRAIL_SIGNALING","REACTOME_INTERFERON_SIGNALING","REACTOME_APOPTOSIS_INDUCED_DNA_FRAGMENTATION","REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION","REACTOME_ALTERNATIVE_COMPLEMENT_ACTIVATION","REACTOME_ADAPTIVE_IMMUNE_SYSTEM","REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION","REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION","REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES","REACTOME_PROGRAMMED_CELL_DEATH","REACTOME_APOPTOSIS","KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS","REACTOME_INNATE_IMMUNE_SYSTEM","REACTOME_NEUTROPHIL_DEGRANULATION","KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION","HALLMARK_HYPOXIA","REACTOME_ECM_PROTEOGLYCANS","REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES","REACTOME_COLLAGEN_CHAIN_TRIMERIZATION")
 result2<-result[list,]
result2<-as.data.frame(result2)
 
pheatmap(result2,cluster_rows=T,
                cluster_cols = F,
                show_rownames = T,
                show_colnames = T,
                color =colorRampPalette(c("light green", "white","red"))(100),
                cellwidth = 15, cellheight = 10,
                fontsize = 9)
```

```{r,metablism data plot}
library(EnhancedVolcano)
library(library)
library(readxl)
de <- read_excel("../metablism/de.xlsx")
data<-de[,c("Metabolites","Class","Sub Class","adj.P-value","log2(FC)")]

names(data)<-c("Metabolites","Class","Sub Class","pvalue","log2FoldChange")
EnhancedVolcano(data,lab = data$Metabolites,x="log2FoldChange",y="pvalue",pCutoff = 0.05,FCcutoff = 1,col = c("#695D73","#695D73","#695D73","#040772"),ylim = c(0,5),labSize = 5,selectLab = c("Dihydroxyacetone phosphate","Glucose pyruvate lactate","PE(22:0/18:1(12Z)-O(9S,10R))","PE-NMe(18:1(11Z)/16:0)","N-Acetylsphinganine","PC(12:0/20:1(11Z))","PC(20:4(6E,8Z,11Z,14Z)-OH(5S)/P-18:1(9Z))","PC(P-16:0/20:4(5Z,7E,11Z,14Z)-OH(9))","PC(16:0/18:2(9E,11E))","MG(0:0/22:4(7Z,10Z,13Z,16Z)/0:0)","PS(P-16:0/19:0)"),drawConnectors = TRUE,widthConnectors = 0.4,arrowheads = F,legendPosition = "top")


EnhancedVolcano(data,lab = data$Metabolites,x="log2FoldChange",y="pvalue",pCutoff = 0.05,FCcutoff = 5,col = c("#695D73","#695D73","#695D73","#040772"),ylim = c(0,5),labSize = 5,drawConnectors = TRUE,widthConnectors = 0.1,arrowheads = F,legendPosition = "top")

EnhancedVolcano(data,lab = data$Metabolites,x="log2FoldChange",y="pvalue",pCutoff = 0.05,FCcutoff = 1,col = c("#695D73","#695D73","#695D73","#040772"),ylim = c(0,5),labSize = 5,selectLab = c("Dihydroxyacetone phosphate","Glucose pyruvate lactate","PE(22:0/18:1(12Z)-O(9S,10R))","PE-NMe(18:1(11Z)/16:0)","N-Acetylsphinganine","PC(12:0/20:1(11Z))","PC(20:4(6E,8Z,11Z,14Z)-OH(5S)/P-18:1(9Z))","PC(P-16:0/20:4(5Z,7E,11Z,14Z)-OH(9))","PC(16:0/18:2(9E,11E))","MG(0:0/22:4(7Z,10Z,13Z,16Z)/0:0)","PS(P-16:0/19:0)","fructose-6-phosphate pyruvate","Methyl pyruvate","Indolepyruvate","fructose-6-phosphate lactate","Beta-D-Fructose 6-phosphate","beta-D-fructose 2,6-bisphosphate","1-Deoxy-1-morpholino-D-fructose","D-1-[(3-Carboxypropyl)amino]-1-deoxyfructose","cis-3-Hexenyl lactate","galactose lactate","Sodium lactate","Phosphocreatine lactate","3-(3,5-Diiodo-4-hydroxyphenyl)lactate","Glucose lactate ketone","Alanine lactate","D-Glucose, 6-deoxy-6-((7-nitro-4-benzofurazanyl)amino)-","D-Glucose, 2-deoxy-2-[[(methylnitrosoamino)carbonyl]amino]-"),drawConnectors = TRUE,widthConnectors = 0.4,arrowheads = F,legendPosition = "top")
```




library(Seurat)
library(hdf5r)

wt1 <- Read10X(data.dir = "liverWT1")
wt2 <- Read10X(data.dir = "liverWT2")
wt3 <- Read10X(data.dir = "liverWT3")
wt4 <- Read10X(data.dir = "liverWT4")

abc1 <- Read10X(data.dir = "liverABC1")
abc2 <- Read10X(data.dir = "liverABC2")
abc3 <- Read10X(data.dir = "liverABC3")
abc4 <- Read10X(data.dir = "liverABC4")



wt1obj <- CreateSeuratObject(counts = wt1, project = "liverWT", min.cells = 3, min.features = 300)
wt2obj <- CreateSeuratObject(counts = wt2, project = "liverWT", min.cells = 3, min.features = 300)
wt3obj <- CreateSeuratObject(counts = wt3, project = "liverWT", min.cells = 3, min.features = 300)
wt4obj <- CreateSeuratObject(counts = wt4, project = "liverWT", min.cells = 3, min.features = 300)
abc1obj <- CreateSeuratObject(counts = abc1, project = "liverABCC4KO", min.cells = 3, min.features = 300)
abc2obj <- CreateSeuratObject(counts = abc2, project = "liverABCC4KO", min.cells = 3, min.features = 300)
abc3obj <- CreateSeuratObject(counts = abc3, project = "liverABCC4KO", min.cells = 3, min.features = 300)
abc4obj <- CreateSeuratObject(counts = abc4, project = "liverABCC4KO", min.cells = 3, min.features = 300)


liver <- merge(wt1obj, y = c(wt2obj,wt3obj,wt4obj,abc1obj,abc2obj,abc3obj,abc4obj), add.cell.ids = c("wt1","wt2","wt3","wt4","abc1","abc2","abc3","abc4"), project = "liver")
head(liver[[]])

## description of the object with cell origins and normalization
list <- SplitObject(liver, split.by = "orig.ident")
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})



## find common anchors
anchors <- FindIntegrationAnchors(object.list = list, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"

save(combined,file="anchors.rda")


combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
ElbowPlot(combined,ndims = 50)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
library(pals)
DimPlot(combined, reduction = "umap",cols=cols25(),pt.size = 0.1)

combined <- RunTSNE(combined, reduction = "pca", dims = 1:30)
DimPlot(combined, reduction = "tsne",cols=cols25(),pt.size = 0.1)
FeaturePlot(combined, features = c("Cd24a", "Tnfrsf12a"), blend = T,split.by="orig.ident")+theme(text = element_text(size = 4))
FeatureScatter(combined, feature1 = "Cd24a", feature2 = "Tnfrsf12a",group.by = 'orig.ident',cols = cols25())
save(combined,file="clusters.rda")

print(combined[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(combined, dims = 1:2, reduction = "pca")
library(pals)
DimPlot(combined, reduction = "umap",cols=c("green","blue","red"))
DimPlot(combined, reduction = "tsne",cols=cols25())


combined <- FindNeighbors(combined, dims = 1:30)
combined <- FindClusters(combined, resolution = 0.2)
markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

FeaturePlot(combined, features = c("Tnfrsf12a"),min.cutoff = "q9",
            cols=c("lightgrey","darkblue"),split.by= "orig.ident",
            reduction = "umap",pt.size=0.3)
library(ggplot2)
VlnPlot(object = combined, features = 'Tnfrsf12a',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())

x<-combined[[]]
write.table(markers,file="ENUreg.txt")

## clusters
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5,algorithm = 2)

DimPlot(combined, reduction = "tsne",cols=cols25(),group.by = "seurat_clusters")

markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(markers,file="markers_clusters.txt")
FeaturePlot(combined, features = c("Vim"),min.cutoff = "q9",
            cols=c("lightgrey","darkblue"),split.by= "orig.ident",
            reduction = "umap",pt.size=0.3)

FeaturePlot(combined, features = c("Cyb5a"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),split.by= "orig.ident",
            reduction = "tsne",pt.size=1)

Idents(object = combined) <- 'orig.ident'
Idents(object = combined)
library(ggplot2)
library(pals)
VlnPlot(object = combined, features = 'RUNX1',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("RUNX1"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),split.by= "orig.ident",
            reduction = "tsne",pt.size=0.1)

##doheatmap
library(ggplot2)
library(dplyr)
markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(subset(combined, downsample = 50), features = top5$gene,size =2.5,group.bar = T) + NoLegend()+theme(text = element_text(size = 8))

downsample = 100




##KDR endo
RidgePlot(combined, features = "Epcam", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "KDR", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'KDR',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'KDR',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("KDR"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)
#MSX1 mesoderm
RidgePlot(combined, features = "MSX1", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "MSX1", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'MSX1',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'MSX1',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("MSX1"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)

#CD34
RidgePlot(combined, features = "CD34", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "CD34", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'MSX1',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'MSX1',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("CD34"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)

FeaturePlot(combined, features = c("CD34", "KDR"), blend = T,split.by="orig.ident")+theme(text = element_text(size = 4))
FeatureScatter(combined, feature1 = "Cd24a", feature2 = "Tnfrsf12a",group.by = 'orig.ident',cols = cols25())


#CD41
RidgePlot(combined, features = "ITGA2B", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "ITGA2B", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'ITGA2B',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'ITGA2B',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("ITGA2B"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)

#CD41
RidgePlot(combined, features = "ITGA2B", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "ITGA2B", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'ITGA2B',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'ITGA2B',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("ITGA2B"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)

#HBA1
RidgePlot(combined, features = "HBA1", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "HBA1", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'HBA1',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'HBG1',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("HBA1"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)

#MYC
RidgePlot(combined, features = "MYC", ncol = 1,group.by = "orig.ident")+scale_fill_manual(values = cols25())
RidgePlot(combined, features = "HBA1", ncol = 1,group.by = "seurat_clusters")+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'KRT19',  group.by = "seurat_clusters",pt.size=0.1)+scale_fill_manual(values = cols25())
VlnPlot(object = combined, features = 'KRT19',  group.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())
FeaturePlot(combined, features = c("HBA1"),min.cutoff = "q9",
            cols=c("beige","darkorchid4"),
            reduction = "tsne",pt.size=0.1)


# barplot


library(ggplot2)

df <- data.frame(x = x$orig.ident, z = x$seurat_clusters)
df <- as.data.frame(with(df, prop.table(table(x, z), margin = NULL)))



plot <- ggplot(data = df, aes(x = x, y = Freq, fill = z)) + 
  geom_bar(width = 0.9, position = "fill", stat = "identity") + 
  scale_fill_manual(values = cols25())+ 
  scale_y_continuous(expand = c(0.01, 0), labels = scales::percent_format()) + 
  xlab("patients") + 
  ylab("Percent") + 
  labs(fill = "seurat_clusters") + 
  theme_classic(base_size = 12, base_family = "sans") + 
  theme(legend.position = "right") + coord_flip()
plot

library(tidyr)
df %>%
  pivot_wider(names_from = z, values_from = Freq)->mat

write.csv2(mat,file="clusterprop.csv",row.names=F)



## fonction from toppfun enrichment
library(ggpubr)
ggdotchart(data, x = "Name", y = "NLOGFDR",
           group = "Name", color = "Name",
           rotate = TRUE,
           add = "segments",                             
           sorting = "ascending",dot.size = 6, label = round(data$NLOGFDR),
           font.label = list(color = "white", size = 9, 
                             vjust = 0.5),               
           ggtheme = theme_pubr(),
           y.text.col = TRUE )




## subsetting seurat 4.0
tail(combined[[]])
Idents(object = combined) <- 'orig.ident'
Idents(object = combined)

### subset graph
library(ggplot2)
VlnPlot(object = subset(combined, idents = "liverABCC4KO", subset = seurat_clusters == "0" | seurat_clusters ==  "5" | seurat_clusters == "11" , slot ="data"), 
        features = 'Epcam',  group.by = "seurat_clusters",split.by = "orig.ident",pt.size=0.1)+scale_fill_manual(values = cols25())


table<-subset(combined, idents = "liverABCC4KO", subset = seurat_clusters == "0" | seurat_clusters ==  "5" | seurat_clusters == "11")@assays[["RNA"]]@counts

meta<-subset(combined, idents = "liverABCC4KO", subset = seurat_clusters == "0" | seurat_clusters ==  "5" | seurat_clusters == "11")[[]]

write.table(table, file='sub_count.tsv', quote=FALSE, sep='\t', col.names = TRUE)
write.table(meta, file='sub_pheno.tsv', quote=FALSE, sep='\t', col.names = TRUE)

## downsampling
table<-as.matrix(table)
tmat_sub<-t(table)
tmat_sub<-as.data.frame(tmat_sub)

ident_sub<-as.data.frame(meta)
### merge count and identities
all<-merge(ident_sub,tmat_sub,by="row.names")

library(dplyr)
df <- all %>% sample_n(1316)
head(df[1:5,1:10])
table(df$integrated_snn_res.0.2)

row.names(df)<-df$Row.names
data<-df[,8:ncol(df)]
tmat<-t(data)

smallmeta<-df[,1:7]
dim(smallmeta)

write.table(tmat, file='sample_count.tsv', quote=FALSE, sep='\t', col.names = TRUE)
write.table(smallmeta, file='sample_pheno.tsv', quote=FALSE, sep='\t', col.names = TRUE)

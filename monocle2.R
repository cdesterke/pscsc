library(monocle)

data<-read.table("sample_count.tsv",h=T)

meta<-read.table("sample_pheno.tsv",h=T)

mat<-as.matrix(data)

meta$id<-colnames(mat)
meta$group<-meta[,7]
row.names(meta)<-meta$id
meta<-meta[,1:2]
colnames(meta)<-c("group","id")
pd <- new("AnnotatedDataFrame", data = meta)

genes<-as.data.frame(row.names(data))
colnames(genes)<-"gene_short_name"
row.names(genes)<-genes$gene_short_name

fd <- new("AnnotatedDataFrame", data = genes)

cds <- newCellDataSet(mat, phenoData = pd,featureData = fd,expressionFamily=negbinomial())



cth <- newCellTypeHierarchy()
Sox9_id <- row.names(subset(fData(cds), gene_short_name == "Sox9"))
Nr0b2_id <- row.names(subset(fData(cds), gene_short_name == "Nr0b2"))
cth <- addCellType(cth, "Sox9negNr0b2pos", classify_func = function(x) { x[Sox9_id,] < 1 & x[Nr0b2_id,] >= 1} )
cth <- addCellType(cth, "Sox9negNr0b2neg", classify_func = function(x) { x[Sox9_id,] < 1 & x[Nr0b2_id,] < 1} )
cth <- addCellType(cth, "Sox9posNr0b2neg", classify_func = function(x) { x[Sox9_id,] >= 1 & x[Nr0b2_id,] < 1} )
cth <- addCellType(cth, "Sox9posNr0b2pos", classify_func = function(x) { x[Sox9_id,] >= 1 & x[Nr0b2_id,] >= 1} )
cds <- classifyCells(cds, cth, 0.1)

table(pData(cds)$CellType)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

plot_pc_variance_explained(cds, return_all = FALSE)

cds <- reduceDimension(cds, max_components = 2, num_dim = 20, reduction_method = 'tSNE', verbose = TRUE)



cds <- clusterCells(cds)
table(pData(cds)$CellType)
#ggplot graph cell type
library(ggplot2)

pie <- ggplot(pData(cds),
              aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())


##plot_cell_clusters(cds, 1, 2, color = "Cluster") +  facet_wrap(~CellType)
plot_cell_clusters(cds, 1, 2, color = "group",cell_size=1)+ theme(legend.position = "none")

cds <- detectGenes(cds, min_expr = 0)
print(head(fData(cds)))
print(head(pData(cds)))

expressed_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
diff_test_res <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~CellType")

write.table(diff_test_res,file="cholangioDif.tsv",sep="\t")

sig_gene_names <- row.names(subset(diff_test_res, pval < 0.05))
sig_gene_names
save(cds,file="monoclehiearchyNr0b2.rda")

differential<-diff_test_res
write.table(differential,file="DEG.txt")

diff_test_res <- diff_test_res[-grep("Rpl", diff_test_res$gene_short_name), ]
diff_test_res <- diff_test_res[-grep("mt-", diff_test_res$gene_short_name), ]
diff_test_res <- diff_test_res[-grep("Rps", diff_test_res$gene_short_name), ]
sig_gene_names <- row.names(subset(diff_test_res, pval < 0.05))
my_ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]
my_ordering_genes

cds2 <- setOrderingFilter(cds, ordering_genes = my_ordering_genes)

#dimensional reduction and gene selection
cds2 <- reduceDimension(cds2, method = 'DDRTree')
gene_to_cluster <- row.names(diff_test_res)[order(diff_test_res$qval)][1:75] 



cds2 <- orderCells(cds2)

save(cds2,file = "cds2.rda")

#build graphs on trajectory 
plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="Sox9",cell_size=2 )

plot_cell_clusters(cds2, color_by = 'as.factor(CellType)',markers="Nr0b2",cell_size=2 )+  facet_wrap(~group)

plot_cell_trajectory(cds2, color_by = "group",cell_size=2)+ theme(legend.position = "none")
plot_cell_trajectory(cds2, color_by = "CellType",markers="Nr0b2",markers_linear = T,show_branch_points=T)+ theme(legend.position = "none")
plot_cell_trajectory(cds2, color_by = "group",markers="Nr0b2",markers_linear = F,show_branch_points=T)
plot_cell_trajectory(cds2, color_by = "Pseudotime",markers="Sox9",markers_linear = F) +  facet_wrap(~group)
plot_cell_trajectory(cds2, color_by = "Pseudotime",markers="Sox9",markers_linear = F) +  facet_wrap(~group)

plot_cell_trajectory(cds2, color_by = "CellType",cell_size=2 ) +  facet_wrap(~group)
plot_cell_trajectory(cds2, color_by = "CellType",cell_size=2 ) +  facet_wrap(~group)

my_pseudotime_cluster <- plot_pseudotime_heatmap(cds2[gene_to_cluster,],cores = 8,
                                                 show_rownames = TRUE,num_clusters=3,
                                                 return_heatmap = TRUE,cluster_rows = TRUE)


plot_genes_in_pseudotime(cds2[c("Nr0b2","Sox9","Id2","Gsta3","Alb","Tmem45a"),],cell_size = 1, color_by = "CellType",ncol = 1)+theme(legend.position = "none")

#library step----
library(biomaRt)
library(tidyverse)
library(tximport)
library(EnsDb.Hsapiens.v86)
library(ensembldb)



path <- file.path("mappedReads", targets$sample, "abundance.tsv")
targets=read_tsv("studydesign.txt",col_names=T)
file.exists(path)


tx=transcripts(EnsDb.Hsapiens.v86,columns=c("tx_id", "gene_name"))
tx = as_tibble(tx) %>%
  dplyr::rename(target_id = tx_id) %>%
  dplyr::select(target_id, gene_name)
txi_gene = tximport(
  path,
  type = "kallisto",
  tx2gene = tx,                
  txOut = FALSE, 
  countsFromAbundance = "lengthScaledTPM",
  ignoreTxVersion = TRUE
)



#data wrangling----
library(edgeR)
library(matrixStats)
library(cowplot)
library(ggplot2)
samplelabels=targets$sample
mycounts=txi_gene$counts
colSums(mycounts)
mytpm=txi_gene$abundance
mycounts_db=as_tibble(mycounts,rownames="gene_id")
mycounts_stats=(mycounts_db) %>%
  dplyr::mutate(sum=rowSums(mycounts),
                sd=rowSds(mycounts),
                md=rowMedians(mycounts))
ggplot(mycounts_stats)+
  aes(x=sd,y=md)+
  geom_point(shape=10,size=4)+
  geom_smooth()


mydgelist=DGEList(mycounts)
cpm=cpm(mydgelist)
cpm_log=cpm(mydgelist,log = TRUE)
cpm_log_df=as_tibble(cpm_log,rownames="geneID")
colnames(cpm_log_df) <- c("geneID", samplelabels)
cpm_log_df_pivot=pivot_longer(cpm_log_df,cols = HS01:CL13,
                              names_to = "samples",
                              values_to = "expression")
p1=ggplot(cpm_log_df_pivot)+
  aes(x=samples,y=expression,fill = samples)+
  geom_violin(trim = TRUE,show.legend = FALSE)+
  geom_smooth()+
  stat_summary(fun = "median",geom="point",shape=10,show.legend = FALSE)+
  theme_bw()+
  coord_flip()
keeps=rowSums(cpm>1)>=5
mydgelist_filtered=DGEList(cpm[keeps,])
mydgelist_filtered_norm=calcNormFactors(mydgelist_filtered,method = "TMM")
mydgelist_filtered_norm_cpm=cpm(mydgelist_filtered_norm,log=TRUE)
mydgelist_filtered_norm_cpm_df=as_tibble(mydgelist_filtered_norm_cpm,rownames="geneID")
colnames(mydgelist_filtered_norm_cpm_df)=c("geneID",samplelabels)
mydgelist_filtered_norm_cpm_df_pivot=pivot_longer(mydgelist_filtered_norm_cpm_df,
                                                  cols=HS01:CL13,
                                                  names_to = "samples",
                                                  values_to = "expression")
view(mydgelist_filtered_norm_cpm)
p2=ggplot(mydgelist_filtered_norm_cpm_df_pivot)+
  aes(x=samples,y=expression,fill = samples)+
  geom_violin(trim = FALSE,show.legend = FALSE)+
  geom_smooth()+
  stat_summary(fun = "median",geom="point",shape=10,show.legend = FALSE)+
  theme_bw()+
  coord_flip()
plot_grid(p1,p2,ncol=2,labels = c("raw","filtered_norm"))

#multivariate----
library(plotly)
group=targets$group
group=factor(group)

#denodogram plot
distance=dist(t(mydgelist_filtered_norm_cpm),method ="maximum" )
clusters=hclust(distance,method = "average")
plot(clusters,labels = samplelabels,main="Hierarchical clustering")

#pca analysis
pca.res=prcomp(t(mydgelist_filtered_norm_cpm),retx = TRUE,scale = FALSE)
screeplot(pca.res, type = "lines", main = "Scree Plot")
summary(pca.res)
pca.var=(pca.res$sdev^2)
pca.per=round(pca.var/sum(pca.var)*100,1)

pca.res.df = as_tibble(pca.res$x)

ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = samplelabels) + 
  geom_point(size = 5) +                        
  theme_bw() +
  xlab(paste0("PC1 (", pca.per[1], "%)")) +
  ylab(paste0("PC2 (", pca.per[2], "%)")) +
  ggtitle("PCA plot")

#small multiple charts
pca.res.df=as_tibble(pca.res$x[,1:4])
pca.res.df.pivot = pca.res.df %>%
  add_column(samples = samplelabels,
             type = group) %>%
  
  pivot_longer(cols = PC1:PC4,          
               names_to = "PC",
               values_to = "scores") 
ggplot(pca.res.df.pivot)+
  aes(x=samples,y=scores)+
  geom_bar(stat = "identity",aes(fill=type))+
  geom_smooth()+
  facet_wrap(~PC,scales = "free_y")

#logfc
mydata.df=mydgelist_filtered_norm_cpm_df %>%
  mutate(avghealthy=(HS01 + HS02 + HS03 + HS04 + HS05)/5,
         avgdisease=(CL08 + CL10 + CL11 + CL12 + CL13)/5,
         logfc=(avgdisease-avghealthy))%>%
  mutate_if(is.numeric,round,2)

myplot <- ggplot(mydata.df) + 
  aes(x=avghealthy, y=avgdisease,text=geneID) +
  geom_point(shape=16, size=1) +
  ggtitle("disease vs. healthy") +
  theme_bw()
ggplotly(myplot)

#diff genes----
library(gt)
group=factor(targets$group)
design_matrix=model.matrix(~0 +group)
colnames(design_matrix)=levels(group)


v.mydge.filtered.norm=voom(mydgelist_filtered_norm,design_matrix,plot = TRUE)
fit=lmFit(v.mydge.filtered.norm,design_matrix)

contrasts=makeContrasts(infection=disease-healthy,
                        levels = design_matrix)
fit2=contrasts.fit(fit,contrasts)
ebfit=eBayes(fit2,robust = TRUE)

top.genes=topTable(ebfit,adjust.method = "BH",coef = 1,sort.by = "logFC",number = Inf)

top.genes_df=as_tibble(top.genes,rownames = "geneID")
ggplot(top.genes_df,
       aes(x = logFC,
           y = -log10(adj.P.Val),
           text=geneID)) +
  geom_point(shape=16, size=1) +
  geom_smooth()+
  ggtitle("disease vs. healthy") +
  theme_bw()

results=decideTests(ebfit,method = "global",adjust.method = "BH",p.value = 0.01,lfc = 2)
class(results)

vennDiagram(results, include="up")

diff_genes=v.mydge.filtered.norm[results[,1]!=0,]
diff_gene_df=as_tibble(diff_genes$E,rownames = "geneID")
gt(diff_gene_df)

#heatmap----
library(pheatmap)
library(RColorBrewer)

sig_mat = diff_genes$E
rownames(sig_mat) = diff_gene_df$geneID

my_annotation = data.frame(Group = targets$group)
rownames(my_annotation) = samplelabels
colnames(sig_mat) = rownames(my_annotation)

pheatmap(sig_mat,
         scale = "row",         
         annotation_col = my_annotation,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         show_rownames = FALSE, 
         cluster_cols = TRUE,    
         main = "Significant Differentially Expressed Genes")

#functional enrichment----
library(clusterProfiler)
library(org.Hs.eg.db) 

sig_gene_names = diff_gene_df$geneID

go_results = enrichGO(gene          = sig_gene_names,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "SYMBOL", 
                      ont           = "BP",     
                      pAdjustMethod = "BH",     
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

go_plot = dotplot(go_results, showCategory = 15) +
  ggtitle("GO Enrichment Analysis: Biological Processes") +
  theme_bw()

print(go_plot)

# extracting----
library(readr)

up_regulated_genes = top.genes_df %>%
  filter(adj.P.Val < 0.01 & logFC >= 2)

down_regulated_genes = top.genes_df %>%
  filter(adj.P.Val < 0.01 & logFC <= -2)

write_tsv(top.genes_df, "DGE_all_genes_results.tsv") 
write_csv(up_regulated_genes, "up_regulated_genes.csv")
write_csv(down_regulated_genes, "down_regulated_genes.csv")

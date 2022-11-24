library("Rsubread")
library("pheatmap")
library("ggplot2")
library("DESeq2")
library("biomaRt")
library('edgeR')
library("dendextend")
library("matrixStats")
library("ggrepel")

bamfiles<-list.files(path ="/Volumes/GLFSSD/sequencing/RNA-seq/Mettl3/FASTQ_Generation_2022-01-11_12_04_04Z-511491981",
                     pattern = '*_duplicates.bam',full.names = T)
Rawcounts<-featureCounts(files=bamfiles,GTF.featureType="exon",
                         GTF.attrType="gene_id",useMetaFeatures=TRUE,
                         annot.inbuilt="mm10",
                         isPairedEnd=T,nthreads=32)

getwd()
write.csv(x = Rawcounts$counts,file = "202206_DuplicateRemoved_FeatureCounts_Rawcounts.csv")

Rawcounts<-read.csv("202206_DuplicateRemoved_FeatureCounts_Rawcounts.csv",header=T)
colnames(Rawcounts)<-c("gene_id","sample1","sample10","sample11","sample12","sample13",
                    "sample14","sample15","sample2","sample3","sample4",
                    "sample5","sample6","sample7","sample8","sample9")
rownames(Rawcounts)<-Rawcounts$gene_id 

cts<-Rawcounts[,c("gene_id","sample1","sample2","sample3","sample4","sample5",
                 "sample6","sample7","sample8","sample9","sample10",
               "sample11","sample12","sample13","sample14","sample15")]
rownames(cts)<-cts$gene_id 

experiment_design<-read.table(file = "MeRIP_seq_design.txt",header = T)
coldata<-experiment_design

####RPKM calculation####
counts<-cts
mart<- useDataset("mmusculus_gene_ensembl", 
                  useMart("ENSEMBL_MART_ENSEMBL"))
listAttributes(mart)
yIDs<-getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name",
                "entrezgene_id","strand","transcript_length"),
  values= rownames(counts),
  mart= mart)

dim(yIDs)
head(yIDs)
# remove duplicates
yIDs.matrix<-yIDs[-which(duplicated(yIDs$entrezgene_id)),]
dim(yIDs.matrix)
# select genes
rownames(yIDs.matrix)<-yIDs.matrix$entrezgene_id
sum(rownames(counts) %in% rownames(yIDs.matrix))
y.annot <- counts[rownames(counts) %in% rownames(yIDs.matrix),]
dim(y.annot)
# join the two tables
yID.annot <- merge(yIDs.matrix,y.annot,by='row.names',all=T)
# check the annotated table
head(yID.annot)
nrow(yID.annot)
rownames(yID.annot)<-yID.annot$entrezgene_id

yRpkm<-rpkm(y =counts[rownames(yID.annot),],gene.length = yID.annot$transcript_length,log = F )
yRpkm.annot <-merge(yID.annot,yRpkm,by='row.names',all=T)
colnames(yRpkm.annot)
yRpkm.merged<-yRpkm.annot[,-c(1,2,9,25)]
colnames(yRpkm.merged)
colnames(yRpkm.merged)<-c("ensembl_gene_id","mgi_symbol","chromosome_name","entrezgene_id",
                          "strand","transcript_length",paste("sample",c(1:15),"_Rawcount",sep=''),
                          paste("sample",c(1:15),"_RPKM",sep=''))

write.table(yRpkm.merged,file="RPKMAnno.txt",row.names=T,quote=F,sep="\t")
###pheatmap on RPKM###
yRpkm.merged.selected<-yRpkm.merged[apply(log2(yRpkm.merged[,27:36])>=1,1,sum)>=4,]
ColAnno<-data.frame("Fraction"=c(rep(c("Input","IgG-flowthrough","IgG-elution",
                                   "anti-m6A-flowthrough","anti-m6A-elution"),2)))
rownames(ColAnno)<-colnames(yRpkm.merged.selected[,27:36])
p1<-pheatmap(mat = yRpkm.merged.selected[,27:36],
         cluster_cols = T,treeheight_col = 40,
         cluster_rows = T,treeheight_row = 0,border_color = NA,
         show_rownames = F, show_colnames = F,
         scale = "row",display_numbers = F, annotation_col=ColAnno,
         color = colorRampPalette(c("blue", "white", "red"))(50))

col_dend<-p1[[2]]
col_dend <- rotate(col_dend, 
                   order = c(colnames(yRpkm.merged.selected[,27:36])
                             [c(3,8,1,6,2,7,4,9,5,10)]))
p1<-pheatmap(mat = yRpkm.merged.selected[,27:36],
             cluster_cols=as.hclust(col_dend),treeheight_col = 40,
             cluster_rows = T,treeheight_row = 0,border_color = NA,
             show_rownames = F, show_colnames = F,
             scale = "row",display_numbers = F, annotation_col=ColAnno,
             color = colorRampPalette(c("blue", "white", "red"))(50))
ggsave(filename = "pheatmap_BiologicalRepeat2&3.pdf",plot = p1,
       device = "pdf",width = 10,height = 8)


yRpkm.merged.selected1<-yRpkm.merged[apply(log2(yRpkm.merged[,c(22,25,26,27,30,31,32,35,36)])>=5,1,sum)>=4,]
ColAnno1<-data.frame("Fraction"=c(rep(c("Input",
                                   "anti-m6A-flowthrough","anti-m6A-elution"),3)))
rownames(ColAnno1)<-colnames(yRpkm.merged.selected1[,c(22,25,26,27,30,31,32,35,36)])
p3<-pheatmap(mat = yRpkm.merged.selected1[,c(22,25,26,27,30,31,32,35,36)],
             cluster_cols = T,treeheight_col = 40,
             cluster_rows = T,treeheight_row = 0,border_color = NA,
             show_rownames = F, show_colnames = F,
             scale = "row",display_numbers = F, annotation_col=ColAnno1,
             color = colorRampPalette(c("blue", "white", "red"))(50))

col_dend1<-p3[[2]]
col_dend1<- rotate(col_dend1, 
                   order = c(colnames(yRpkm.merged.selected1[,c(22,27,32,25,30,35,26,31,36)])))
p4<-pheatmap(mat = yRpkm.merged.selected1[,c(22,25,26,27,30,31,32,35,36)],
             cluster_cols=as.hclust(col_dend1),treeheight_col = 40,
             cluster_rows = T,treeheight_row = 0,border_color = NA,
             show_rownames = F, show_colnames = T,
             scale = "row",display_numbers = F, annotation_col=ColAnno1,
             color = colorRampPalette(c("blue", "white", "red"))(50))
ggsave(filename = "pheatmap_3Goups.pdf",plot = p4,
       device = "pdf",width = 10,height = 8)

####Data glimpse#####
Rawcountcounts<-cts[,-1]
# set the rownames to the sampleID to allow for ordering
rownames(experiment_design) <- paste("sample",c(1:15),sep='')
# order the design following the counts sample order
experiment_design.ord <- experiment_design[colnames(Rawcountcounts),]
# look at the design
experiment_design.ord
# list the ordered samples for future use
samples <- as.character(experiment_design.ord$sampleID)
# create factors for future plotting
group<-factor(experiment_design.ord$condition)
group

#Basic QC plots##
#Library size and distribution plots
y<-DGEList(Rawcountcounts)
y
barplot(y$samples$lib.size/1e06, names=colnames(y), 
        las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")

# density plot of raw read counts (log10)
logcounts <- log(Rawcountcounts[,1],10) 
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.3),xlab="Raw read counts per gene (log10)", ylab="Density")
for (s in 2:length(samples)){
  logcounts <- log(Rawcountcounts[,s],10) 
  d <- density(logcounts)
  lines(d)
}

## box plots of raw read counts (log10)
logcounts <- log((Rawcountcounts+1),10)
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
axis(2)
axis(1,at=c(1:length(samples)),labels=colnames(logcounts),las=2,cex.axis=0.8)
dev.off()

# select data for the 1000 most highly expressed genes
select <- order(rowMeans(Rawcountcounts[,]), decreasing=TRUE)[1:1000]
highexprgenes_counts <-Rawcountcounts[select,]
# heatmap with sample name on X-axis
pheatmap(mat = highexprgenes_counts,scale = "row",
         treeheight_row = 0)
# heatmap with condition group as labels
colnames(highexprgenes_counts)<- group
# plot
pheatmap(mat = highexprgenes_counts,scale = "row",
         treeheight_row = 0)

p0<-pheatmap(log2(highexprgenes_counts[order(highexprgenes_counts[,10]),]+1),
             scale="row",show_rownames = F,cellheight = 0.15,cellwidth = 15,
             cluster_rows = F,cluster_cols = T,
             border_color=NA,labels_col = group[6:15],
             color =colorRampPalette(c("blue","white","red"))(1000))
ggsave(filename = "Biological2&3_HighExpr_1000_heatmap.pdf",plot = p0,
       device = "pdf",width = 6,height = 6,path="/Volumes/GLFSSD/sequencing/RNA-seq/Mettl3/Mettl3")

# select data for the 1000 most highly variable expressed genes
logcounts <- cpm(Rawcountcounts[,],log=TRUE)
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:1000]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

p1<-pheatmap(highly_variable_lcpm[order(highly_variable_lcpm[,1]),],
             scale="row",show_rownames = F,cellheight = 0.2,cellwidth = 20,
             cluster_rows = F,cluster_cols = T,
             border_color=NA,labels_col = group,
             color =colorRampPalette(c("blue","white","red"))(1000))
ggsave(filename = "Biological2&3_HighVari_1000_heatmap.pdf",plot = p1,
       device = "pdf",width = 6,height = 6,path="/Volumes/GLFSSD/sequencing/RNA-seq/Mettl3/Mettl3")

####PCA############
# select data for the 100 most highly varible genes
select <- order(rowSds(as.matrix(Rawcountcounts[,6:15])), decreasing=TRUE)[1:100]
highexprgenes_counts <- Rawcountcounts[select,6:15]
# annotate the data with condition group as labels
colnames(highexprgenes_counts)<- group[6:15]
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
dim(data_for_PCA)
## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned

# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
pdf(file="p3_Biological2&3_High_Varible_genesPCA_PropExplainedVariance.pdf")
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
dev.off()
## calculate MDS
mds <- cmdscale(dist(data_for_PCA))
#Samples representation
pdf(file="p4_Biological2&3_High_varible_genesPCA.pdf")
plot(mds[,1], -mds[,2], type="p", xlab=paste("Dimension1_",eig_pc[1]), ylab=paste("Dimension2_",eig_pc[2]), main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
dev.off()

####DESEQ2######
coldata$sampleID<-paste("sample",c(1:15),sep='')
dds <- DESeqDataSetFromMatrix(countData = cts[,7:16],
                              colData = coldata[6:15,],
                              design = ~ condition)
dds$condition <- factor(dds$condition, 
                        levels = c("Input", 
                                   "IgG_flowthrough",
                                   "IgG_elution",
                                   "m6A_flowthrough",
                                   "m6A_elution" ))
dds$condition <- relevel(dds$condition, 
                         ref = "m6A_elution")
##DEGs
dds <- DESeq(dds)
resultsNames(dds) 

##Input_m6AElution##
Input_vs_m6A_elution_res <- results(dds,alpha =0.05,
                                    cooksCutoff = FALSE,
               independentFiltering = F,
               name = "condition_Input_vs_m6A_elution")
dim(Input_vs_m6A_elution_res)
dim(yRpkm.merged)
yRpkm.merged.input_vs_m6a<-yRpkm.merged[yRpkm.merged$entrezgene_id %in%
                                          rownames(Input_vs_m6A_elution_res) ,]
rownames(yRpkm.merged.input_vs_m6a)<-yRpkm.merged.input_vs_m6a$entrezgene_id
Input_vs_m6A_elution_res<-merge(x=as.data.frame(Input_vs_m6A_elution_res),
                                y=yRpkm.merged.input_vs_m6a,by=0)

write.csv(x=Input_vs_m6A_elution_res,file = "20220612_Input_vs_m6AElution.csv")
Input_m6AElution_Plotdata<-Input_vs_m6A_elution_res[!is.na(Input_vs_m6A_elution_res$log2FoldChange) 
                         & !is.na(Input_vs_m6A_elution_res$padj),]
# add a column of NAs
Input_m6AElution_Plotdata$diffexpressed <- "No"
# set "UP" 
Input_m6AElution_Plotdata$diffexpressed[Input_m6AElution_Plotdata$log2FoldChange > 1 & Input_m6AElution_Plotdata$padj < 0.05] <- "Down"
# set "DOWN"
Input_m6AElution_Plotdata$diffexpressed[Input_m6AElution_Plotdata$log2FoldChange < -1  & Input_m6AElution_Plotdata$padj < 0.05] <- "Up"
sum(Input_m6AElution_Plotdata$diffexpressed=="No")
sum(Input_m6AElution_Plotdata$diffexpressed=="Up")
sum(Input_m6AElution_Plotdata$diffexpressed=="Down")

Klf2x<--Input_m6AElution_Plotdata[Input_m6AElution_Plotdata$mgi_symbol=="Klf2",]$log2FoldChange
Klf2y<--log10(Input_m6AElution_Plotdata[Input_m6AElution_Plotdata$mgi_symbol=="Klf2",]$padj)

pvol_input_m6Aelution<-ggplot(Input_m6AElution_Plotdata,aes(x=-log2FoldChange,y=-log10(padj),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=0.5,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("darkorange3", "#9a9a9a", "darkorange3"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
 scale_x_continuous(limits = c(-12,12))+
 scale_y_continuous(limits = c(0,40))+
  geom_text_repel(
    data = subset(Input_m6AElution_Plotdata, mgi_symbol=="Klf2"),
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(Input_m6AElution_Plotdata$diffexpressed=="Up"),
                      "no",sum(Input_m6AElution_Plotdata$diffexpressed=="No"),
                      "down",sum(Input_m6AElution_Plotdata$diffexpressed=="Down")))
ggsave(filename = "20220612_pvol_input_m6Aelution.pdf",pvol_input_m6Aelution,"pdf",
       width=3, height = 4)

##m6Aflowthrough_m6AElution##
m6Aflowthrough_vs_m6A_elution_res <- results(dds,alpha =0.05,
                                    cooksCutoff = FALSE,
                                    independentFiltering = F,
                                    name = "condition_m6A_flowthrough_vs_m6A_elution")
dim(m6Aflowthrough_vs_m6A_elution_res)
dim(yRpkm.merged)
yRpkm.merged.m6Aflowthrough_vs_m6a<-yRpkm.merged[yRpkm.merged$entrezgene_id %in%
                                          rownames(m6Aflowthrough_vs_m6A_elution_res) ,]
rownames(yRpkm.merged.m6Aflowthrough_vs_m6a)<-yRpkm.merged.m6Aflowthrough_vs_m6a$entrezgene_id
m6Aflowthrough_vs_m6A_elution_res<-merge(x=as.data.frame(m6Aflowthrough_vs_m6A_elution_res),
                                y=yRpkm.merged.m6Aflowthrough_vs_m6a,by=0)

write.csv(x=m6Aflowthrough_vs_m6A_elution_res,file = "20220612_m6Aflowthrough_vs_m6AElution.csv")
m6Aflowthrough_m6AElution_Plotdata<-m6Aflowthrough_vs_m6A_elution_res[!is.na(m6Aflowthrough_vs_m6A_elution_res$log2FoldChange) 
                                                    & !is.na(m6Aflowthrough_vs_m6A_elution_res$padj),]
# add a column of NAs
m6Aflowthrough_m6AElution_Plotdata$diffexpressed <- "No"
# set "UP" 
m6Aflowthrough_m6AElution_Plotdata$diffexpressed[m6Aflowthrough_m6AElution_Plotdata$log2FoldChange > 1 & m6Aflowthrough_m6AElution_Plotdata$padj < 0.05] <- "Down"
# set "DOWN"
m6Aflowthrough_m6AElution_Plotdata$diffexpressed[m6Aflowthrough_m6AElution_Plotdata$log2FoldChange < -1  & m6Aflowthrough_m6AElution_Plotdata$padj < 0.05] <- "Up"
sum(m6Aflowthrough_m6AElution_Plotdata$diffexpressed=="No")
sum(m6Aflowthrough_m6AElution_Plotdata$diffexpressed=="Up")
sum(m6Aflowthrough_m6AElution_Plotdata$diffexpressed=="Down")

Klf2x<--m6Aflowthrough_m6AElution_Plotdata[m6Aflowthrough_m6AElution_Plotdata$mgi_symbol=="Klf2",]$log2FoldChange
Klf2y<--log10(m6Aflowthrough_m6AElution_Plotdata[m6Aflowthrough_m6AElution_Plotdata$mgi_symbol=="Klf2",]$padj)

pvol_m6Aflowthrough_m6Aelution<-ggplot(m6Aflowthrough_m6AElution_Plotdata,aes(x=-log2FoldChange,y=-log10(padj),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=0.5,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("darkorange3", "#9a9a9a", "darkorange3"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(0,20))+
  geom_text_repel(
    data = subset(m6Aflowthrough_m6AElution_Plotdata, mgi_symbol=="Klf2"),
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(m6Aflowthrough_m6AElution_Plotdata$diffexpressed=="Up"),
                      "no",sum(m6Aflowthrough_m6AElution_Plotdata$diffexpressed=="No"),
                      "down",sum(m6Aflowthrough_m6AElution_Plotdata$diffexpressed=="Down")))
ggsave(filename = "20220612_pvol_m6Aflowthrough_m6Aelution.pdf",
       pvol_m6Aflowthrough_m6Aelution,"pdf",
       width=3, height = 4)


####edgeR############
mycounts <- Rawcountcounts[,c(6:15)]
dim(mycounts)

#Keep genes with least 1 count-per-million reads (cpm) in at least 4 samples
isexpr <- rowSums(cpm(mycounts)>1) >= 2
table(isexpr)
mycounts <- mycounts[isexpr,]
genes <- rownames(mycounts)
dim(mycounts)

pheatmap(log2(mycounts+1),scale="row",
         cluster_rows = T,cluster_cols = T,
         border_color=NA,show_rownames = F,treeheight_row = 0,
         color =colorRampPalette(c("blue","white","red"))(50) )

# check your samples grouping
experiment_design.ord[colnames(mycounts),]$condition == group[c(6:15)]
# create design matrix for limma
design <- model.matrix(~0+group[c(6:15)])
# substitute "group" from the design column names
#colnames(design)<- gsub("group","",colnames(design))
colnames(design)<-levels(group)
# check your design matrix
design

# calculate normalization factors between libraries
nf <- calcNormFactors(mycounts)
#
# normalise the read counts with 'voom' function
yv <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
# extract the normalised read counts
counts.voom <- yv$E
# save normalised expression data into output dir
write.csv(counts.voom,file="edgeR_20220611_Biological2&3_counts.voom.csv")


##m6A_elution-Input##
# fit linear model for each gene given a series of libraries
fit <- lmFit(yv,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts(m6A_elution-Input,levels=design)
cont.matrix 
# compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit(fit, cont.matrix)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes(fit)
options(digits=3)
# check the output fit
dim(fit)
# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=2
# get the coefficient name for the comparison  of interest
colnames(fit$coefficients)
mycoef<-colnames(fit$coefficients)
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit,coef=mycoef)
# get the full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])
dim(limma.res)
write.csv(limma.res,"edgeR_20220611_Biological2&3_m6A_elution-Input_AllExpressedGene.csv")

limma.res$diffexpressed <- "No"
limma.res$diffexpressed[limma.res$logFC > 1 & limma.res$adj.P.Val < 0.05] <- "Up"
limma.res$diffexpressed[limma.res$logFC < -1& limma.res$adj.P.Val < 0.05] <- "Down"
sum(limma.res$diffexpressed=="No")
sum(limma.res$diffexpressed=="Up")
sum(limma.res$diffexpressed=="Down")

Klf2x<-limma.res["16598",]$logFC
Klf2y<--log10(limma.res["16598",]$adj.P.Val)

p3<-ggplot(limma.res,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=0.5,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("darkorange3", "#9a9a9a", "darkorange3"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,6))+
 geom_text_repel(
  data = limma.res["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res$diffexpressed=="Up"),
                      "no",sum(limma.res$diffexpressed=="No"),
                      "down",sum(limma.res$diffexpressed=="Down")))
ggsave(filename = "edgeR_20220612_pvol_m6Aelution_Input.pdf",
       p3,"pdf",
       width=3, height = 4)

##m6A_elution-m6Aflowthrough##
# fit linear model for each gene given a series of libraries
fit <- lmFit(yv,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts(m6A_elution-m6A_flowthrough,levels=design)
cont.matrix 
# compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit(fit, cont.matrix)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes(fit)
options(digits=3)
# check the output fit
dim(fit)
# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=2
# get the coefficient name for the comparison  of interest
colnames(fit$coefficients)
mycoef<-colnames(fit$coefficients)
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit,coef=mycoef)
# get the full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])
dim(limma.res)
write.csv(limma.res,"edgeR_20220611_Biological2&3_m6A_elution-m6Aflowthrough_AllExpressedGene.csv")

limma.res$diffexpressed <- "No"
limma.res$diffexpressed[limma.res$logFC > 1 & limma.res$adj.P.Val < 0.05] <- "Up"
limma.res$diffexpressed[limma.res$logFC < -1& limma.res$adj.P.Val < 0.05] <- "Down"
sum(limma.res$diffexpressed=="No")
sum(limma.res$diffexpressed=="Up")
sum(limma.res$diffexpressed=="Down")

Klf2x<-limma.res["16598",]$logFC
Klf2y<--log10(limma.res["16598",]$adj.P.Val)

p4<-ggplot(limma.res,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=0.5,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("darkorange3", "#9a9a9a", "darkorange3"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,5))+
  geom_text_repel(
    data = limma.res["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res$diffexpressed=="Up"),
                      "no",sum(limma.res$diffexpressed=="No"),
                      "down",sum(limma.res$diffexpressed=="Down")))
ggsave(filename = "edgeR_20220612_pvol_m6AElution_m6Aflowthrough.pdf",
       p4,"pdf",
       width=3, height = 4)


####20221008 useing 1/2*RPKM as filter####
data1<-yRpkm.merged
colnames(data1)
data1$average<-apply(data1[,c(22,27,32)],1,mean)
dim(data1)
length(data1$average[data1$average!=0])
median(data1$average[data1$average!=0])

data2<-data1[apply(data1[,c("sample1_RPKM" ,"sample6_RPKM","sample11_RPKM")],1,max)
                   >median(data1$average[data1$average!=0]),]
dim(data2)
colnames(data2)

###filter on DESeq2 result###
Input_vs_m6A_data2<-Input_vs_m6A_elution_res[Input_vs_m6A_elution_res$mgi_symbol %in% data2$mgi_symbol,]
dim(Input_vs_m6A_data2)
colnames(Input_vs_m6A_data2)
Input_vs_m6A_data2_noNA<-Input_vs_m6A_data2[!is.na(Input_vs_m6A_data2$log2FoldChange)&
                                              !is.na(Input_vs_m6A_data2$padj),]
Input_vs_m6A_data2_sig<-Input_vs_m6A_data2_noNA[abs(Input_vs_m6A_data2_noNA$log2FoldChange)>1 & 
                                                  Input_vs_m6A_data2_noNA$padj<0.05,]
dim(Input_vs_m6A_data2_sig)
colnames(Input_vs_m6A_data2_sig)
Input_vs_m6A_data2_sig$Klf2orNOT<-"F"
jj<-1
for (jj in (1:nrow(Input_vs_m6A_data2_sig))){
  if (Input_vs_m6A_data2_sig[jj,"mgi_symbol"]=="Klf2"){
    Input_vs_m6A_data2_sig$Klf2orNOT[jj]<-"T"
  }
}
Input_vs_m6A_data2_sig
pheatmap(mat =Input_vs_m6A_data2_sig[,c(34:43)],
        scale = 'row',cluster_cols = T,show_rownames =F )


limma.res_data2<-limma.res[rownames(limma.res) %in% data2$entrezgene_id,]
dim(limma.res_data2)
ggplot(limma.res_data2,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=0.5,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("darkorange3", "#9a9a9a", "darkorange3"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  scale_x_continuous(limits = c(-10,10))+
  scale_y_continuous(limits = c(0,5))+
  geom_text_repel(
    data = limma.res_data2["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res_data2$diffexpressed=="Up"),
                      "no",sum(limma.res_data2$diffexpressed=="No"),
                      "down",sum(limma.res_data2$diffexpressed=="Down")))


library("Rsubread")
library("pheatmap")
library("ggplot2")
library("DESeq2")
library("biomaRt")
library('edgeR')
library("dendextend")
library("matrixStats")
library("ggrepel")

getwd()
#load data
Rawcounts<-read.csv("./data/202206_DuplicateRemoved_FeatureCounts_Rawcounts.csv",header=T)
colnames(Rawcounts)<-c("gene_id",sub(pattern = "_marked_duplicates.bam",
                         replacement = "",x =colnames(Rawcounts)[-1]))
rownames(Rawcounts)<-Rawcounts$gene_id

cts<-Rawcounts[,c("gene_id","sample1","sample2","sample3","sample4","sample5",
                 "sample6","sample7","sample8","sample9","sample10",
               "sample11","sample12","sample13","sample14","sample15")]
rownames(cts)<-cts$gene_id 

experiment_design<-read.table(file = "./data/MeRIP_seq_design.txt",header = T)
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

write.csv(yRpkm.merged,file="./result/202210_RPKMAnno.csv")

###data filter based median value of all sample####
data1<-yRpkm.merged
colnames(data1)
data1$average<-apply(data1[,c(27:36)],1,median)
dim(data1)
length(data1$average[data1$average!=0])
median(data1$average[data1$average!=0])
kept<-apply(data1[,c(27,32)],1,mean)>median(data1$average[data1$average!=0]) |
  apply(data1[,c(28,33)],1,mean)>median(data1$average[data1$average!=0])|
 apply(data1[,c(29,34)],1,mean)>median(data1$average[data1$average!=0])|
  apply(data1[,c(30,35)],1,mean)>median(data1$average[data1$average!=0])|
  apply(data1[,c(31,36)],1,mean)>median(data1$average[data1$average!=0])    
#kept<-apply(data1[,c(27,30,31,32,35,36)],1,max)>median(data1$average[data1$average!=0])
data2<-data1[kept,]
dim(data2)
data2<-na.omit(data2)
colnames(data2)
rownames(data2)<-data2$entrezgene_id
data2[data2$mgi_symbol=="Klf2",]
#data2<-data2[data2$sample6_RPKM>0|data2$sample11_RPKM>0,]
#dim(data2)
###pheatmap on filtered RPKM###
rownames(yRpkm.merged)<-yRpkm.merged$entrezgene_id
Count_RPKM_MedianFiltered<-yRpkm.merged[rownames(data2),]
dim(Count_RPKM_MedianFiltered)
Count_RPKM_MedianFiltered<-na.omit(Count_RPKM_MedianFiltered)
colnames(Count_RPKM_MedianFiltered)
pHeatmap_FilteredRPKM<-pheatmap(Count_RPKM_MedianFiltered[,c(27:36)],
         scale = "row",show_rownames = F, 
         color =colorRampPalette(c("blue","white","red"))(1000))

ColAnno<-data.frame("sample"=colnames(Count_RPKM_MedianFiltered[,c(27:36)]))
ColAnno$fraction<-rep(c("Input",
                        "IgG_unbound","IgG_bound",
                    "antim6A_unbound","antim6A_bound"),2)
rownames(ColAnno)<-ColAnno$sample
ann_colors <- list(
  fraction=c("Input"="#00da97",
             "IgG_unbound"="#bdc000","IgG_bound"="#00caff",
             "antim6A_unbound"="#ff89ff",
             "antim6A_bound"="#ff9188"))

pHeatmap_FilteredRPKMcol_dend<-pHeatmap_FilteredRPKM[[2]]

pHeatmap_FilteredRPKMcol_dend <-dendextend::rotate(pHeatmap_FilteredRPKMcol_dend, 
                                                   order =c("sample6_RPKM","sample11_RPKM",
                                                           "sample8_RPKM","sample13_RPKM",
                                                           "sample7_RPKM","sample12_RPKM",
                                                            "sample9_RPKM" ,"sample14_RPKM",
                                                            "sample10_RPKM","sample15_RPKM"
                                                            ))
pHeatmap_FilteredRPKM2<-pheatmap(Count_RPKM_MedianFiltered[,c(27:36)],
                                scale = "row",show_rownames = F, 
                                annotation_col=ColAnno,
                                annotation_colors = ann_colors,
                                color =colorRampPalette(c("blue","white","red"))(1000))

ggsave(filename = "./result/20221031_Median_Filtred_RPKM_heatmap.pdf",device = "pdf",
       plot =pHeatmap_FilteredRPKM2,height=8,width=8)

##PCA
highexprgenes_counts <- Count_RPKM_MedianFiltered[,c(27:36)]
# annotate the data with condition group as labels
colnames(highexprgenes_counts)<- colnames(Count_RPKM_MedianFiltered)[c(27:36)]
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
as.data.frame(eig_pc)
PCA_FilteredRPKM_1<-ggplot(data =as.data.frame(eig_pc),
                           mapping = aes(x=c(1:10),y=eig_pc))+
  geom_bar(stat="identity")
## calculate MDS
mds <- cmdscale(dist(data_for_PCA))
#Samples representation
PCA_FilteredRPKM_2<-ggplot(data =as.data.frame(mds),
                           mapping = aes(x=as.data.frame(mds)[,1],
                                         y=as.data.frame(mds)[,2],
                                        label=rownames(as.data.frame(mds))))+
  geom_point(size=4)+
  geom_text(position = "nudge")+
  labs(x = element_text(paste("PCA_1_",round(eig_pc[1],2),sep="")),
       y= element_text(paste("PCA_2_",round(eig_pc[2],2),sep="")))+
  theme_bw()

ggsave(filename = "./result/20221031_Median_Filtred_RPKM_PCA.pdf",plot =PCA_FilteredRPKM_2,
       device = "pdf",height=8,width=8)
####DESEQ2######
coldata$sampleID<-paste("sample",c(1:15),sep='')
dds <- DESeqDataSetFromMatrix(countData = cts[,c(7:16)],
                              colData = coldata[c(6:15),],
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

#write.csv(x=Input_vs_m6A_elution_res,file = "20220612_Input_vs_m6AElution.csv")
Input_m6AElution_Plotdata<-Input_vs_m6A_elution_res[Input_vs_m6A_elution_res$entrezgene_id %in% data2$entrezgene_id,]
dim(Input_m6AElution_Plotdata)
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

Input_m6AElution_Plotdata

p_volDESeq2_input_m6Aelution<-ggplot(Input_m6AElution_Plotdata,aes(x=-log2FoldChange,y=-log10(padj),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=4)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
 scale_x_continuous(limits = c(-12,12))+
 scale_y_continuous(limits = c(0,20))+
  geom_text_repel(
    data = subset(Input_m6AElution_Plotdata, mgi_symbol=="Klf2"),
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(Input_m6AElution_Plotdata$diffexpressed=="Up"),
                      "no",sum(Input_m6AElution_Plotdata$diffexpressed=="No"),
                      "down",sum(Input_m6AElution_Plotdata$diffexpressed=="Down")))
ggsave(filename = "./result/20221031_DESeq2_pvol_input_m6Aelution.pdf",p_volDESeq2_input_m6Aelution,"pdf",
       width=8, height = 8)

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

#write.csv(x=m6Aflowthrough_vs_m6A_elution_res,file = "20220612_m6Aflowthrough_vs_m6AElution.csv")
m6Aflowthrough_m6AElution_Plotdata<-m6Aflowthrough_vs_m6A_elution_res[m6Aflowthrough_vs_m6A_elution_res$entrezgene_id %in% data2$entrezgene_id,]
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
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
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
ggsave(filename = "./result/20221031_DESeq2_pvol_m6Aflowthrough_m6Aelution.pdf",
       pvol_m6Aflowthrough_m6Aelution,"pdf",
       width=8, height = 8)

##IgG_elution_m6AElution##
IgG_elution_vs_m6A_elution_res <- results(dds,alpha =0.05,
                                          cooksCutoff = FALSE,
                                          independentFiltering = F,
                                          name = "condition_IgG_elution_vs_m6A_elution")
dim(IgG_elution_vs_m6A_elution_res)
dim(yRpkm.merged)
yRpkm.merged.IgG_elution_vs_m6a<-yRpkm.merged[yRpkm.merged$entrezgene_id %in%
                                                rownames(IgG_elution_vs_m6A_elution_res) ,]
rownames(yRpkm.merged.IgG_elution_vs_m6a)<-yRpkm.merged.IgG_elution_vs_m6a$entrezgene_id
IgG_elution_vs_m6A_elution_res<-merge(x=as.data.frame(IgG_elution_vs_m6A_elution_res),
                                      y=yRpkm.merged.IgG_elution_vs_m6a,by=0)

#write.csv(x=IgG_elution_vs_m6A_elution_res,file = "20220612_IgG_elution_vs_m6AElution.csv")
IgG_elution_m6AElution_Plotdata<-IgG_elution_vs_m6A_elution_res[IgG_elution_vs_m6A_elution_res$entrezgene_id %in% data2$entrezgene_id,]
# add a column of NAs
IgG_elution_m6AElution_Plotdata$diffexpressed <- "No"
# set "UP" 
IgG_elution_m6AElution_Plotdata$diffexpressed[IgG_elution_m6AElution_Plotdata$log2FoldChange > 1 & IgG_elution_m6AElution_Plotdata$padj < 0.05] <- "Down"
# set "DOWN"
IgG_elution_m6AElution_Plotdata$diffexpressed[IgG_elution_m6AElution_Plotdata$log2FoldChange < -1  & IgG_elution_m6AElution_Plotdata$padj < 0.05] <- "Up"
sum(IgG_elution_m6AElution_Plotdata$diffexpressed=="No")
sum(IgG_elution_m6AElution_Plotdata$diffexpressed=="Up")
sum(IgG_elution_m6AElution_Plotdata$diffexpressed=="Down")

Klf2x<--IgG_elution_m6AElution_Plotdata[IgG_elution_m6AElution_Plotdata$mgi_symbol=="Klf2",]$log2FoldChange
Klf2y<--log10(IgG_elution_m6AElution_Plotdata[IgG_elution_m6AElution_Plotdata$mgi_symbol=="Klf2",]$padj)

pvol_IgG_elution_m6Aelution<-ggplot(IgG_elution_m6AElution_Plotdata,aes(x=-log2FoldChange,y=-log10(padj),col=diffexpressed))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  scale_x_continuous(limits = c(-12,12))+
  scale_y_continuous(limits = c(0,20))+
  geom_text_repel(
    data = subset(IgG_elution_m6AElution_Plotdata, mgi_symbol=="Klf2"),
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(IgG_elution_m6AElution_Plotdata$diffexpressed=="Up"),
                      "no",sum(IgG_elution_m6AElution_Plotdata$diffexpressed=="No"),
                      "down",sum(IgG_elution_m6AElution_Plotdata$diffexpressed=="Down")))
ggsave(filename = "./result/20221031_DESeq2_pvol_IgG_elution_m6Aelution.pdf",
       pvol_IgG_elution_m6Aelution,"pdf",
       width=8, height = 8)

####edgeR############
mycounts <- cts[,c(7:16)]
dim(mycounts)

#Keep genes with least 1 count-per-million reads (cpm) in at least 4 samples
#isexpr <- rowSums(cpm(mycounts)>1) >= 2
#table(isexpr)
#mycounts <- mycounts[isexpr,]
genes <- rownames(mycounts)
dim(mycounts)
mycounts1<-na.omit(mycounts)
# pheatmap(mycounts1,scale="row",
#          cluster_rows = T,cluster_cols = F,
#          border_color=NA,show_rownames = F,treeheight_row = 0,
#          color =colorRampPalette(c("blue","white","red"))(50) )

# set the rownames to the sampleID to allow for ordering
rownames(experiment_design) <- paste("sample",c(1:15),sep='')
# order the design following the counts sample order
experiment_design.ord <- experiment_design[colnames(cts),]
# look at the design
experiment_design.ord<-experiment_design.ord[-1,]
experiment_design.ord
# list the ordered samples for future use
samples <- as.character(experiment_design.ord$sampleID)
# create factors for future plotting
group<-factor(experiment_design.ord$condition[c(6:15)])
group
# check your samples grouping
experiment_design.ord[colnames(mycounts),]$condition == group
# create design matrix for limma
design <- model.matrix(~0+group)
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
write.csv(counts.voom,file="./result/20221015_edgeR_Biological2&3_counts.voom.csv")


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
write.csv(limma.res,"./result/20221015_edgeR_Biological2&3_m6A_elution-Input_AllExpressedGene.csv")

limma.res.plotdata<-limma.res[rownames(limma.res) %in% data2$entrezgene_id,]
dim(limma.res.plotdata)
limma.res.plotdata$diffexpressed <- "No"
limma.res.plotdata$diffexpressed[limma.res.plotdata$logFC > 2 & limma.res.plotdata$adj.P.Val < 0.05] <- "Up"
limma.res.plotdata$diffexpressed[limma.res.plotdata$logFC < -2& limma.res.plotdata$adj.P.Val < 0.05] <- "Down"
sum(limma.res.plotdata$diffexpressed=="No")
sum(limma.res.plotdata$diffexpressed=="Up")
sum(limma.res.plotdata$diffexpressed=="Down")

Klf2x<-limma.res.plotdata["16598",]$logFC
Klf2y<--log10(limma.res.plotdata["16598",]$adj.P.Val)

p_edgeRm6A_elution_Input<-ggplot(limma.res.plotdata,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  # scale_x_continuous(limits = c(-10,10))+
  #scale_y_continuous(limits = c(0,10))+
 geom_text_repel(
  data = limma.res.plotdata["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res.plotdata$diffexpressed=="Up"),
                      "no",sum(limma.res.plotdata$diffexpressed=="No"),
                      "down",sum(limma.res.plotdata$diffexpressed=="Down")))
ggsave(filename = "./result/20221031_edgeR_pvol_m6Aelution_Input.pdf",
       p_edgeRm6A_elution_Input,"pdf", width=8, height = 8)

##m6A_elution-m6Aflowthrough##
# fit linear model for each gene given a series of libraries
fit1 <- lmFit(yv,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix1 <- makeContrasts(m6A_elution-m6A_flowthrough,levels=design)
cont.matrix1 
# compute estimated coefficients and standard errors for a given set of contrasts
fit1 <- contrasts.fit(fit1, cont.matrix1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit1 <- eBayes(fit1)
options(digits=3)
# check the output fit
dim(fit1)
# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=2
# get the coefficient name for the comparison  of interest
colnames(fit1$coefficients)
mycoef1<-colnames(fit1$coefficients)
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit1,coef=mycoef1)
# get the full table ("n = number of genes in the fit")
limma.res2 <- topTable(fit1,coef=mycoef1,n=dim(fit1)[1])
dim(limma.res2)
write.csv(limma.res2,"./result/20221015_edgeR_Biological2&3_m6A_elution-m6Aflowthrough_AllExpressedGene.csv")

limma.res2.Median.Filtered<-limma.res2[rownames(limma.res2) %in% data2$entrezgene_id,]
dim(limma.res2.Median.Filtered)
limma.res2.Median.Filtered$diffexpressed <- "No"
limma.res2.Median.Filtered$diffexpressed[limma.res2.Median.Filtered$logFC > 2 & limma.res2.Median.Filtered$adj.P.Val < 0.05] <- "Up"
limma.res2.Median.Filtered$diffexpressed[limma.res2.Median.Filtered$logFC < -2 & limma.res2.Median.Filtered$adj.P.Val < 0.05] <- "Down"
sum(limma.res2.Median.Filtered$diffexpressed=="No")
sum(limma.res2.Median.Filtered$diffexpressed=="Up")
sum(limma.res2.Median.Filtered$diffexpressed=="Down")

Klf2x<-limma.res2.Median.Filtered["16598",]$logFC
Klf2y<--log10(limma.res2.Median.Filtered["16598",]$adj.P.Val)

p4_edgeR_m6A_elution_m6Aflowthrough<-ggplot(limma.res2.Median.Filtered,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  #scale_x_continuous(limits = c(-10,10))+
  #scale_y_continuous(limits = c(0,5))+
  geom_text_repel(
    data = limma.res2.Median.Filtered["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res2.Median.Filtered$diffexpressed=="Up"),
                      "no",sum(limma.res2.Median.Filtered$diffexpressed=="No"),
                      "down",sum(limma.res2.Median.Filtered$diffexpressed=="Down")))
ggsave(filename = "./result/20221031_edgeR_pvol_m6AElution_m6Aflowthrough.pdf",
       p4_edgeR_m6A_elution_m6Aflowthrough,"pdf",
       width=8, height = 8)


##m6A_elution-IgG_elution##
# fit linear model for each gene given a series of libraries
fit3 <- lmFit(yv,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix3 <- makeContrasts(m6A_elution-IgG_elution,levels=design)
cont.matrix3 
# compute estimated coefficients and standard errors for a given set of contrasts
fit3 <- contrasts.fit(fit3, cont.matrix3)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit3 <- eBayes(fit3)
options(digits=3)
# check the output fit
dim(fit3)
# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=2
# get the coefficient name for the comparison  of interest
colnames(fit3$coefficients)
mycoef3<-colnames(fit3$coefficients)
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit3,coef=mycoef3)
# get the full table ("n = number of genes in the fit")
limma.res3 <- topTable(fit3,coef=mycoef3,n=dim(fit3)[1])
dim(limma.res3)

limma.res3.Median.Filtered<-limma.res3[rownames(limma.res3) %in% data2$entrezgene_id,]
dim(limma.res3.Median.Filtered)
limma.res3.Median.Filtered$diffexpressed <- "No"
limma.res3.Median.Filtered$diffexpressed[limma.res3.Median.Filtered$logFC > 2 & limma.res3.Median.Filtered$adj.P.Val < 0.05] <- "Up"
limma.res3.Median.Filtered$diffexpressed[limma.res3.Median.Filtered$logFC < -2& limma.res3.Median.Filtered$adj.P.Val < 0.05] <- "Down"
sum(limma.res3.Median.Filtered$diffexpressed=="No")
sum(limma.res3.Median.Filtered$diffexpressed=="Up")
sum(limma.res3.Median.Filtered$diffexpressed=="Down")

Klf2x<-limma.res3.Median.Filtered["16598",]$logFC
Klf2y<--log10(limma.res3.Median.Filtered["16598",]$adj.P.Val)

p5_m6A_elution_IgG_elution<-ggplot(limma.res3.Median.Filtered,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  # scale_x_continuous(limits = c(-10,10))+
  # scale_y_continuous(limits = c(0,5))+
  geom_text_repel(
    data = limma.res3.Median.Filtered["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res3.Median.Filtered$diffexpressed=="Up"),
                      "no",sum(limma.res3.Median.Filtered$diffexpressed=="No"),
                      "down",sum(limma.res3.Median.Filtered$diffexpressed=="Down")))
ggsave(filename = "./result/20221031_edgeR_pvol_m6AElution_IgG_elution.pdf",
       p5_m6A_elution_IgG_elution,"pdf",
       width=8, height = 8)

##m6A_elution-IgG_flowthrough##
# fit linear model for each gene given a series of libraries
fit4 <- lmFit(yv,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix4 <- makeContrasts(m6A_elution-IgG_flowthrough,levels=design)
cont.matrix4 
# compute estimated coefficients and standard errors for a given set of contrasts
fit4 <- contrasts.fit(fit4, cont.matrix4)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit4 <- eBayes(fit4)
options(digits=3)
# check the output fit
dim(fit4)
# set adjusted pvalue threshold and log fold change threshold
mypval=0.05
myfc=2
# get the coefficient name for the comparison  of interest
colnames(fit4$coefficients)
mycoef4<-colnames(fit4$coefficients)
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit4,coef=mycoef4)
# get the full table ("n = number of genes in the fit")
limma.res4 <- topTable(fit4,coef=mycoef4,n=dim(fit4)[1])
dim(limma.res4)

limma.res4.Median.Filtered<-limma.res4[rownames(limma.res4) %in% data2$entrezgene_id,]
dim(limma.res4.Median.Filtered)
limma.res4.Median.Filtered$diffexpressed <- "No"
limma.res4.Median.Filtered$diffexpressed[limma.res4.Median.Filtered$logFC > 2 & limma.res4.Median.Filtered$adj.P.Val < 0.05] <- "Up"
limma.res4.Median.Filtered$diffexpressed[limma.res4.Median.Filtered$logFC < -2& limma.res4.Median.Filtered$adj.P.Val < 0.05] <- "Down"
sum(limma.res4.Median.Filtered$diffexpressed=="No")
sum(limma.res4.Median.Filtered$diffexpressed=="Up")
sum(limma.res4.Median.Filtered$diffexpressed=="Down")

Klf2x<-limma.res4.Median.Filtered["16598",]$logFC
Klf2y<--log10(limma.res4.Median.Filtered["16598",]$adj.P.Val)

p6_m6A_elution_IgG_flowthrough<-ggplot(limma.res4.Median.Filtered,aes(x=logFC,y=-log10(adj.P.Val),col=diffexpressed))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5)+
  geom_point(alpha=12,size=1,pch=16) +
  geom_point(x=Klf2x,y=Klf2y,color="blue",size=1)+
  theme_bw(base_size = 12) + 
  scale_color_manual(values=c("#276acc", "#9a9a9a", "#f45864"))+ 
  theme(legend.position = "bottom",
        panel.grid=element_blank())+
  # scale_x_continuous(limits = c(-10,10))+
  # scale_y_continuous(limits = c(0,5))+
  geom_text_repel(
    data = limma.res4.Median.Filtered["16598",],
    aes(label = "Klf2"),
    size = 5,color="blue")+
  labs(subtitle=paste("up",sum(limma.res4.Median.Filtered$diffexpressed=="Up"),
                      "no",sum(limma.res4.Median.Filtered$diffexpressed=="No"),
                      "down",sum(limma.res4.Median.Filtered$diffexpressed=="Down")))
ggsave(filename = "./result/20221031_edgeR_pvol_m6AElution_IgG_flowthrough.pdf",
       p6_m6A_elution_IgG_flowthrough,"pdf",
       width=8, height = 8)


#####use edgeR_pvol_m6A_elution_Input#
#transcript enriched in m6A elution compared with Input, IgG elution, and m6A flowthrough
l1<-rownames(limma.res.plotdata[limma.res.plotdata$diffexpressed %in% c("Up","Down"),])
l2<-rownames(limma.res2.Median.Filtered[limma.res2.Median.Filtered$diffexpressed %in% c("Up","Down"),])
l3<-rownames(limma.res3.Median.Filtered[limma.res3.Median.Filtered$diffexpressed %in% c("Up","Down"),])

Diffgenes_id<-intersect(intersect(l1,l2),l3)
length(Diffgenes_id)

mart<- useDataset("mmusculus_gene_ensembl", 
                  useMart("ENSEMBL_MART_ENSEMBL"))
listAttributes(mart)
Diffgenes_id_names<-getBM(
  filters = "entrezgene_id",
  attributes= c("ensembl_gene_id", "mgi_symbol","chromosome_name",
                "entrezgene_id","strand","transcript_length"),
  values= Diffgenes_id,
  mart= mart)

dim(Diffgenes_id_names)
head(Diffgenes_id_names)

Diffgenes_id_names<-Diffgenes_id_names$mgi_symbol
Diffgenes_id_names<-Diffgenes_id_names[!duplicated(Diffgenes_id_names)]
length(Diffgenes_id_names)

####load mouse TFs from PMID 30204897##
AllTFs<-read.csv("./data/Mus_musculus_TF.txt",header = T,sep='\t')
TFnames<-AllTFs$Symbol
TFnames<-TFnames[!duplicated(TFnames)]

DiffTF_id_names<-intersect(Diffgenes_id_names,TFnames)


###Load diff genes between con and mut LepR cluster from scRNA-data
LepR_DiffGene<-read.csv(file = "./data/20210917_LepR_DiffGenes.csv")
LepR_DiffGene<-LepR_DiffGene[abs(LepR_DiffGene$avg_log2FC)>2 
                             & LepR_DiffGene$p_val_adj<0.01,]
LepR_DiffIF<-intersect(LepR_DiffGene$X,TFnames)

"Klf2" %in% intersect(DiffTF_id_names,LepR_DiffIF)

# Load library


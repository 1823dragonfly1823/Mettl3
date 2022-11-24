library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(sctransform)
library(multtest)
library(metap)
library(xlsx)
library(ggiraph)
library(ggrepel)
library(monocle3)
library(DoubletFinder)
library(SeuratWrappers)
library(SeuratDisk)
library(pheatmap)
library(rlang)
library(grid)
#########define function to DoMultiBarHeatmap########
DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          #x.divs <- pbuild$layout$panel_params[[1]]$x.major
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

#####ImportData############################
Dir_Ctrl<-'./data/LL001/analysis/210326_LEI_LONGFEI_2_MOUSE_10X-LL001-cellranger-count-default/LL001_cellranger_count_outs/raw_feature_bc_matrix'
Dir_Mut<-'./data/LL002/analysis/210326_LEI_LONGFEI_2_MOUSE_10X-LL002-cellranger-count-default/LL002_cellranger_count_outs/raw_feature_bc_matrix'

Ctrl.data <- Read10X(data.dir = Dir_Ctrl)
Ctrl <- CreateSeuratObject(counts = Ctrl.data, project = "Ctrl")
Ctrl

Mut.data <- Read10X(data.dir = Dir_Mut)
Mut <- CreateSeuratObject(counts = Mut.data, project = "Mut")
Mut

##Ctrl QC####
head(colnames(Ctrl))
tail(colnames(Ctrl))
table(Ctrl$orig.ident)
View(Ctrl@meta.data)
View(Ctrl@assays$RNA@counts@Dimnames[[1]])

# Add number of genes per UMI for each cell to metadata
Ctrl$log10GenesPerUMI <- log10(Ctrl$nFeature_RNA)/log10(Ctrl$nCount_RNA)
max(Ctrl$log10GenesPerUMI[Ctrl$log10GenesPerUMI!="NaN"])

##Extract mito genes from features list
mitogenes<-grep('^mt-',Ctrl@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)

# Compute percent mito ratio
Ctrl$mitoRatio <- PercentageFeatureSet(object = Ctrl, pattern = "^mt-")
Ctrl$mitoRatio <- Ctrl@meta.data$mitoRatio/100

#compute ribosomal gene ratio
Ribogenes<-grep('^Rp[sl][[:digit:]]',Ctrl@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
Ctrl$RiboRatio <- PercentageFeatureSet(object = Ctrl, pattern = "^Rp[sl][[:digit:]]")
Ctrl$RiboRatio <- Ctrl@meta.data$RiboRatio/100

# Create metadata dataframe
metadata <- Ctrl@meta.data
View(metadata)
# Rename columns
metadata<-dplyr::rename (metadata,seq_folder=orig.ident,
                         nUMI = nCount_RNA,
                         nGene = nFeature_RNA)
colnames(metadata)
View(metadata)
# Create sample column
metadata$sample <- NA
metadata$sample[metadata$seq_folder=='Ctrl'] <- "Control"
#metadata$sample[metadata$seq_folder=='Mut'] <- "Mutant"

# Add metadata back to Seurat object
Ctrl@meta.data <- metadata
View(Ctrl@meta.data)

CtrlnUMIplot<-ggplot(Ctrl@meta.data[Ctrl@meta.data$nUMI>100,],
                 aes(sample,nUMI))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=1000,color="red",size=1)+
  annotate(geom="text", x=0.5, y=1000, label="1000",
           size=6,color="blue")

CtrlnGeneplot<-ggplot(Ctrl@meta.data[Ctrl@meta.data$nUMI>100,],
                  aes(sample,nGene))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=1000,color="red",size=1)+
  annotate(geom="text", x=0.5, y=1000, label="1000",
           size=6,color="blue")

Ctrlnlog10GenesPerUMIplot<-ggplot(Ctrl@meta.data[Ctrl@meta.data$nUMI>100,],
                              aes(sample,log10GenesPerUMI))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.8,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.8, label="0.8",
           size=6,color="blue")

CtrlmitoRatioplot<-ggplot(Ctrl@meta.data[Ctrl@meta.data$nUMI>100,],
                      aes(sample,mitoRatio))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.1,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.10, label="0.1",
           size=6,color="blue")

CtrlRiboRatioplot<-ggplot(Ctrl@meta.data[Ctrl@meta.data$nUMI>100,],
                      aes(sample,RiboRatio))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.2,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.2, label="0.2",
           size=6,color="blue")

CtrlQCplot<-(CtrlnUMIplot+CtrlnGeneplot+Ctrlnlog10GenesPerUMIplot+
               CtrlmitoRatioplot+CtrlRiboRatioplot)
ggsave(filename = "./result/20221105_Ctrl_QCplot.pdf",plot=CtrlQCplot,
       device="pdf",width=8,height=6)

# Ctrl Filter####
##Filter out low quality reads using selected thresholds - these will change with experiment
Ctrl.filtered<-subset(x = Ctrl, 
                           subset= ((nUMI >= 1000) & 
                                      (nGene >= 1000)&
                                      (log10GenesPerUMI>0.8)&
                                      (mitoRatio < 0.1)&
                                      (RiboRatio<0.2)))
table(Ctrl.filtered@meta.data$sample)
##Scale and Normalize##
Ctrl.filtered<- CellCycleScoring(Ctrl.filtered, 
                                      g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                      s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                      set.ident=T)
Ctrl.filtered <- NormalizeData(Ctrl.filtered, normalization.method = "LogNormalize", scale.factor = 10000)
Ctrl.filtered <- FindVariableFeatures(Ctrl.filtered, selection.method = "vst", nfeatures = 2000)
Ctrl.filtered <- ScaleData(Ctrl.filtered, vars.to.regress = c("mitoRatio",
                                                              "S.Score", 
                                                              "G2M.Score",
                                                              "RiboRatio"))

Ctrl.filtered <- RunPCA(Ctrl.filtered, features = VariableFeatures(object = Ctrl.filtered))
DimPlot(Ctrl.filtered, reduction = "pca")
ElbowPlot(Ctrl.filtered)

Ctrl.filtered<-RunUMAP(Ctrl.filtered,dims = 1:20,reduction = "pca")
DimPlot(Ctrl.filtered,reduction = "umap")
Ctrl.filtered <- RunTSNE(Ctrl.filtered,dims = 1:20,reduction = "pca")
DimPlot(Ctrl.filtered,reduction = "tsne")

Ctrl.filtered <- FindNeighbors(object = Ctrl.filtered,dims = 1:20)
Ctrl.filtered <- FindClusters(object = Ctrl.filtered,resolution = c(0.1,0.2,0.3,0.4))
Idents(Ctrl.filtered)<-"RNA_snn_res.0.2"
##pK identify
sweep.res.list.Ctrl <- paramSweep_v3(Ctrl.filtered, PCs = 1:50, sct = FALSE)
sweep.stats.Ctrl <- summarizeSweep(sweep.res.list.Ctrl, GT = FALSE)
bcmvn.Ctrl <- find.pK(sweep.stats.Ctrl)
ggplot(bcmvn.Ctrl,mapping = aes(pK,BCmetric,group=1))+
  geom_point()+
  geom_line()
## Homotypic Doublet Proportion Estimate
homotypic.prop.Ctrl <- modelHomotypic(Ctrl.filtered@meta.data$RNA_snn_res.0.2)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi.Ctrl <- round(0.04*nrow(Ctrl.filtered@meta.data))  ## Assuming 4% doublet formation rate - tailor for your dataset
nExp_poi.Ctrl.adj <- round(nExp_poi.Ctrl*(1-homotypic.prop.Ctrl))

Ctrl.filtered <- doubletFinder_v3(Ctrl.filtered, PCs = 1:20, pN = 0.25, pK = 0.3, 
                                  nExp = nExp_poi.Ctrl.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(object = Ctrl.filtered,group.by = "DF.classifications_0.25_0.3_30")
VlnPlot(Ctrl.filtered, features = "nFeature_RNA", 
        group.by = "DF.classifications_0.25_0.3_30", pt.size = 0.1)
Ctrl.filtered.singlet<-Ctrl.filtered[, Ctrl.filtered@meta.data[, "DF.classifications_0.25_0.3_30"] == "Singlet"]
dim(Ctrl.filtered.singlet)

##Mut QC####
head(colnames(Mut))
tail(colnames(Mut))
table(Mut$orig.ident)
View(Mut@meta.data)
View(Mut@assays$RNA@counts@Dimnames[[1]])

# Add number of genes per UMI for each cell to metadata
Mut$log10GenesPerUMI <- log10(Mut$nFeature_RNA)/log10(Mut$nCount_RNA)
max(Mut$log10GenesPerUMI[Mut$log10GenesPerUMI!="NaN"])

##Extract mito genes from features list
mitogenes<-grep('^mt-',Mut@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)

# Compute percent mito ratio
Mut$mitoRatio <- PercentageFeatureSet(object = Mut, pattern = "^mt-")
Mut$mitoRatio <- Mut@meta.data$mitoRatio/100

#compute ribosomal gene ratio
Ribogenes<-grep('^Rp[sl][[:digit:]]',Mut@assays[["RNA"]]@counts@Dimnames[[1]],ignore.case = FALSE,value=T)
Mut$RiboRatio <- PercentageFeatureSet(object = Mut, pattern = "^Rp[sl][[:digit:]]")
Mut$RiboRatio <- Mut@meta.data$RiboRatio/100

# Create metadata dataframe
metadata <- Mut@meta.data
View(metadata)
# Rename columns
metadata<-dplyr::rename (metadata,seq_folder=orig.ident,
                         nUMI = nCount_RNA,
                         nGene = nFeature_RNA)
colnames(metadata)
View(metadata)
# Create sample column
metadata$sample <- NA
#metadata$sample[metadata$seq_folder=='Mut'] <- "Control"
metadata$sample[metadata$seq_folder=='Mut'] <- "Mutant"

# Add metadata back to Seurat object
Mut@meta.data <- metadata
View(Mut@meta.data)

MutnUMIplot<-ggplot(Mut@meta.data[Mut@meta.data$nUMI>100,],
                    aes(sample,nUMI))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=1000,color="red",size=1)+
  annotate(geom="text", x=0.5, y=1000, label="1000",
           size=6,color="blue")

MutnGeneplot<-ggplot(Mut@meta.data[Mut@meta.data$nUMI>100,],
                     aes(sample,nGene))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  scale_y_continuous(trans='log10')+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=1000,color="red",size=1)+
  annotate(geom="text", x=0.5, y=1000, label="1000",
           size=6,color="blue")

Mutnlog10GenesPerUMIplot<-ggplot(Mut@meta.data[Mut@meta.data$nUMI>100,],
                                 aes(sample,log10GenesPerUMI))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.8,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.8, label="0.8",
           size=6,color="blue")

MutmitoRatioplot<-ggplot(Mut@meta.data[Mut@meta.data$nUMI>100,],
                         aes(sample,mitoRatio))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.1,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.1, label="0.1",
           size=6,color="blue")

MutRiboRatioplot<-ggplot(Mut@meta.data[Mut@meta.data$nUMI>100,],
                         aes(sample,RiboRatio))+
  geom_jitter(size=0.1)+
  geom_violin(fill=NA,size=1,aes(color=sample))+
  theme(legend.position = "none")+
  theme(axis.line= element_line(colour = "black"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=0.25,color="red",size=1)+
  annotate(geom="text", x=0.5, y=0.25, label="0.25",
           size=6,color="blue")

MutQCplot<-(MutnUMIplot+MutnGeneplot+Mutnlog10GenesPerUMIplot+
              MutmitoRatioplot+MutRiboRatioplot)
ggsave(filename = "./result/20221105_Mut_QCplot.pdf",plot=MutQCplot,
       device="pdf",width=8,height=6)

# Mut Filter####
##Filter out low quality reads using selected thresholds - these will change with experiment
Mut.filtered<-subset(x = Mut, 
                     subset= ((nUMI >= 1000) & 
                                (nGene >= 1000)&
                                (log10GenesPerUMI>0.8)&
                                (mitoRatio < 0.1)&
                                (RiboRatio<0.25)))
table(Mut.filtered@meta.data$sample)
##Scale and Normalize##
Mut.filtered<- CellCycleScoring(Mut.filtered, 
                                g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                set.ident=T)
Mut.filtered <- NormalizeData(Mut.filtered, normalization.method = "LogNormalize", scale.factor = 10000)
Mut.filtered <- FindVariableFeatures(Mut.filtered, selection.method = "vst", nfeatures = 2000)
Mut.filtered <- ScaleData(Mut.filtered, vars.to.regress = c("mitoRatio",
                                                            "S.Score", 
                                                            "G2M.Score",
                                                            "RiboRatio"))

Mut.filtered <- RunPCA(Mut.filtered, features = VariableFeatures(object = Mut.filtered))
DimPlot(Mut.filtered, reduction = "pca")
ElbowPlot(Mut.filtered)

Mut.filtered<-RunUMAP(Mut.filtered,dims = 1:20,reduction = "pca")
DimPlot(Mut.filtered,reduction = "umap")
Mut.filtered <- RunTSNE(Mut.filtered,dims = 1:20,reduction = "pca")
DimPlot(Mut.filtered,reduction = "tsne")

Mut.filtered <- FindNeighbors(object = Mut.filtered,dims = 1:20)
Mut.filtered <- FindClusters(object = Mut.filtered,resolution = c(0.1,0.2,0.3,0.4))
Idents(Mut.filtered)<-"RNA_snn_res.0.2"
##pK identify
sweep.res.list.Mut <- paramSweep_v3(Mut.filtered, PCs = 1:50, sct = FALSE)
sweep.stats.Mut <- summarizeSweep(sweep.res.list.Mut, GT = FALSE)
bcmvn.Mut <- find.pK(sweep.stats.Mut)
ggplot(bcmvn.Mut,mapping = aes(pK,BCmetric,group=1))+
  geom_point()+
  geom_line()
## Homotypic Doublet Proportion Estimate 
homotypic.prop.Mut <- modelHomotypic(Mut.filtered@meta.data$RNA_snn_res.0.2)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi.Mut <- round(0.04*nrow(Mut.filtered@meta.data))  ## Assuming 4% doublet formation rate - tailor for your dataset
nExp_poi.Mut.adj <- round(nExp_poi.Mut*(1-homotypic.prop.Mut))

Mut.filtered <- doubletFinder_v3(Mut.filtered, PCs = 1:20, 
                                 pN = 0.25, pK =0.28 , 
                                 nExp = nExp_poi.Mut.adj, reuse.pANN = FALSE, sct = FALSE)
DimPlot(object = Mut.filtered,group.by = "DF.classifications_0.25_0.28_32")
VlnPlot(Mut.filtered, features = "nFeature_RNA", 
        group.by = "DF.classifications_0.25_0.28_32", pt.size = 0.1)
Mut.filtered.singlet<-Mut.filtered[, Mut.filtered@meta.data[, "DF.classifications_0.25_0.28_32"] == "Singlet"]
dim(Mut.filtered.singlet)


###Combine Ctrl and Mut Cells####
filtered_combined<-merge(x=Ctrl.filtered.singlet,y=Mut.filtered.singlet)
View(filtered_combined@meta.data)
#SCTransform
options(future.globals.maxSize = 4000 * 1024^2)
##perform the cell cycle scoring and sctransform on all samples
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(filtered_combined, split.by = "sample")
split_seurat <- split_seurat[c("Control", "Mutant")]
i<-1
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], 
                                        g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)), 
                                        s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),
                                        set.ident=T)
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], 
                                   vars.to.regress = c("mitoRatio",
                                                       "S.Score", 
                                                       "G2M.Score",
                                                       "RiboRatio"),
                                   variable.features.n=3000)
}

###Integrate samples using shared highly variable genes
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features) 
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

###UMAP visualization
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated,
                            ndims.print = 1:50,nfeatures.print = 5)

# Plot PCA
PCAPlot(seurat_integrated) 
ElbowPlot(seurat_integrated,ndims = 50)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:50,
                             reduction = "pca")
# Plot UMAP                             
DimPlot(seurat_integrated,reduction = "umap",group.by = "sample")


# Run TSNE
seurat_integrated <- RunTSNE(seurat_integrated, 
                             dims = 1:50,
                             reduction = "pca")
# Plot TSNE                             
DimPlot(seurat_integrated,reduction = "tsne",group.by = "sample")

# Determine the K-nearest neighbor graph
seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                   dims = 1:50)
# saveRDS(seurat_integrated,"seurat_integrated.RDS")
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                                  resolution = c(0.1,0.2,0.4))
View(seurat_integrated@meta.data)  

ggplot(data = seurat_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,nUMI))+
  geom_boxplot()+
  ggplot(data = seurat_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,nGene))+
  geom_boxplot()+
  ggplot(data = seurat_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,mitoRatio))+
  geom_boxplot()+
  ggplot(data = seurat_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,RiboRatio))+
  geom_boxplot()
  

# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.2"

#Find markers for each cluster
DefaultAssay(seurat_integrated) <- "RNA"
j<-0
for (j in (0:9)){
  conserved_markers <- FindConservedMarkers(seurat_integrated, 
                                            ident.1 = j,
                                            grouping.var = "sample",
                                            verbose = T,only.pos=T) 
  write.csv(x=conserved_markers,file = paste("./result/cluster_",j,".csv",sep=''))
}

seurat_integrated <- RenameIdents(object = seurat_integrated, 
                                  "0" = "LepR+_1",
                                  "1" = "LepR+_2",
                                  "2" = "Osteoblast",
                                  "3" = "Cycling",
                                  "4" = "Osteocyte",
                                  "5" = "Myloid",
                                  "6" = "Endothelial cell",
                                  "7" = "Megakaryocyte",
                                  "8" = "Fibroblast",
                                  "9" = "B cell")


# Plot the UMAP
ALLUMAP<-DimPlot(seurat_integrated,
                 reduction = "umap",
                 label = T,repel = T,
                 label.size = 5,
                 pt.size = 0.5)
ggsave(filename ="./result/20221105_ALLUMAP.pdf",
       plot = ALLUMAP,device = "pdf", 
       width =10,height = 8)

AllTSNE<-DimPlot(seurat_integrated,
                 reduction = "tsne",
                 label = T,repel = T,
                 label.size = 5,
                 pt.size = 0.5)

ggsave(filename = "./result/20221105_AllTSNE.pdf",
       plot = AllTSNE,device = "pdf",
       width =10,height = 8)
# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(seurat_integrated, 
                     vars = c("ident", "seq_folder")) %>%
  dplyr::count(ident, seq_folder) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)
write.csv(x = n_cells,file = "./result/n_cells_All.csv")
# UMAP of cells in each cluster by sample
DimPlot(seurat_integrated, 
        label = TRUE, 
        split.by = "sample")  + NoLegend()
# Explore whether clusters segregate by cell cycle phase
DimPlot(seurat_integrated,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio","RiboRatio")
FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

DefaultAssay(seurat_integrated) <- "RNA"
AllVln1<-VlnPlot(object = seurat_integrated, 
                 features = c("Cxcl12","Kitl","Lepr","Bglap",
                              "Bglap2","Mki67","Dmp1","Lyz2",
                              #"Ngp",
                              "Cdh5","Pf4","Col3a1","Cd79a"
                             # ,"Igkc"
                              ),
                 log=F,ncol = 4)
ggsave(filename ="./result/20221105_AllVln.pdf",plot = AllVln1,
       device = "pdf", 
       width =18,height = 10)

##heatmap
DefaultAssay(seurat_integrated) <- "RNA"
HeatMapGenes<-read.xlsx(file = './data/20221106_HeatamapGenes.xlsx',sheetIndex = 1,header = T)
HMap<-DoHeatmap(subset(ScaleData(seurat_integrated),downsample =100) ,
                features = HeatMapGenes$Genes, size = 3)
ggsave(filename ="./result/20221105_heatmap_ScaleRNA.pdf",
       plot = HMap,device = "pdf", 
       width =10,height = 10)

DefaultAssay(seurat_integrated) <- "integrated"
HeatMapGenes<-read.xlsx(file = './data/20221106_HeatamapGenes.xlsx',sheetIndex = 1,header = T)
HMap1<-DoHeatmap(subset(seurat_integrated,downsample =100) ,
                features = HeatMapGenes$Genes, size = 3)
ggsave(filename ="./result/20221105_heatmap_integrated.pdf",
       plot = HMap1,device = "pdf", 
       width =10,height = 16)


ggplot(seurat_integrated@meta.data,
       aes(integrated_snn_res.0.2,G2M.Score)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("")+
  geom_hline(yintercept = 0.3,color="red",size=1)+
  annotate(geom = "text", x =0, y = 0.3, 
           label = "0.3", hjust="left",color="blue")

ggplot(seurat_integrated@meta.data,
       aes(log10GenesPerUMI,RiboRatio,color=mitoRatio)) + 
  geom_point() + 
  theme_classic()+
  geom_vline(xintercept = 0.8,color="red")+
  geom_hline(yintercept = 0.3,color="red")+
  facet_grid(~integrated_snn_res.0.2)

#####Stromal cells#########################################
stromal<-subset(seurat_integrated, ident = c("LepR+_1","LepR+_2",
                                             "Osteoblast","Osteocyte",
                                             "Fibroblast","Endothelial cell"))
Stromal_seurat <- SplitObject(stromal, split.by = "sample")
Stromal_seurat <- Stromal_seurat[c("Control", "Mutant")]
###Extract stromal cell barcode
CtrlStromalBarcode<-colnames(Stromal_seurat[[1]])
write.table(x =CtrlStromalBarcode,file= "./result/20221106_BARCODE_CtrlStromal_Barcode.csv",
            quote = F,row.names = F,col.names = F)

CtrlStromalBarcode<-colnames(Stromal_seurat[[2]])
write.table(x =MutStromalBarcode,file= "./result/20221106_BARCODE_MutStromal_Barcode.csv",
            quote = F,row.names = F,col.names = F)
length(CtrlStromalBarcode)
length(MutStromalBarcode)
####run Velocyto Command line on HPC####
####load loom files from VElOCYTO###
Ctrl_Velocyto_loom<-ReadVelocity("./data/VelocytoResult/LL001_possorted_genome_bam_filtered_LR11A.loom")
Ctrl_Velocyto_loom_seurat<-as.Seurat(Ctrl_Velocyto_loom)
colnames(Ctrl_Velocyto_loom_seurat)
Ctrl_Velocyto_loom_seurat<-RenameCells(Ctrl_Velocyto_loom_seurat, 
                                       new.names =paste(sub(pattern = "LL001_possorted_genome_bam_filtered_LR11A:",
                                                            replacement = "",
                                                            x = colnames(Ctrl_Velocyto_loom_seurat)),"-1",sep=""))

Mut_Velocyto_loom<-ReadVelocity("./data/VelocytoResult/LL002_possorted_genome_bam_filtered_EKDN1.loom")
Mut_Velocyto_loom_seurat<-as.Seurat(Mut_Velocyto_loom)
colnames(Mut_Velocyto_loom_seurat)
Mut_Velocyto_loom_seurat<-RenameCells(Mut_Velocyto_loom_seurat, 
                                      new.names =paste(sub(pattern = "LL002_possorted_genome_bam_filtered_EKDN1:",
                                                           replacement = "",
                                                           x = colnames(Mut_Velocyto_loom_seurat)),"-1",sep=""))

Stromal_seurat[[1]][["spliced"]]<-Ctrl_Velocyto_loom_seurat[["spliced"]]
Stromal_seurat[[1]][["unspliced"]]<-Ctrl_Velocyto_loom_seurat[["unspliced"]]
Stromal_seurat[[1]][["ambiguous"]]<-Ctrl_Velocyto_loom_seurat[["ambiguous"]]

Stromal_seurat[[2]][["spliced"]]<-Mut_Velocyto_loom_seurat[["spliced"]]
Stromal_seurat[[2]][["unspliced"]]<-Mut_Velocyto_loom_seurat[["unspliced"]]
Stromal_seurat[[2]][["ambiguous"]]<-Mut_Velocyto_loom_seurat[["ambiguous"]]

i<-1
for (i in 1:length(Stromal_seurat)) {
  Stromal_seurat[[i]] <- CellCycleScoring(Stromal_seurat[[i]],
                                          g2m.features=stringr::str_to_title(tolower(cc.genes$g2m.genes)),
                                          s.features=stringr::str_to_title(tolower(cc.genes$s.genes)),set.ident=T)
  Stromal_seurat[[i]] <- SCTransform(Stromal_seurat[[i]],
                                     vars.to.regress = c("mitoRatio","S.Score", 
                                                         "G2M.Score","nUMI","RiboRatio"),
                                     variable.features.n=3000)
}

View(Stromal_seurat)
###Integrate samples using shared highly variable genes
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = Stromal_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
Stromal_seurat <- PrepSCTIntegration(object.list = Stromal_seurat, 
                                     anchor.features = integ_features) 
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = Stromal_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
Stromal_integrated <- IntegrateData(anchorset = integ_anchors, 
                                    normalization.method = "SCT")

###UMAP visualization
# Run PCA
Stromal_integrated <- RunPCA(object = Stromal_integrated,
                             ndims.print = 1:10,nfeatures.print = 5)

# Plot PCA
PCAPlot(Stromal_integrated) 
ElbowPlot(Stromal_integrated,ndims = 50)

# Run UMAP & TSNE
Stromal_integrated <- RunUMAP(object = Stromal_integrated, 
                              reduction="pca",
                              dims = 1:50)
Stromal_integrated <- RunTSNE(Stromal_integrated, 
                              dims = 1:50,
                              reduction = "pca")
# Plot UMAP & TSNE                            
stromalumap0<-DimPlot(Stromal_integrated,reduction = "umap",ncol=1,
                      group.by = "sample",split.by = "sample",cols = c("#12777f","#f35e5a"))
ggsave(filename ="./result/20221106_stromalumap.pdf",plot =stromalumap0,device = "pdf", 
        width =6,height = 10)
DimPlot(Stromal_integrated,reduction = "tsne",group.by = "sample")
# Determine the K-nearest neighbor graph
Stromal_integrated <- FindNeighbors(object = Stromal_integrated, 
                                    dims = 1:50)
# saveRDS(seurat_integrated,"seurat_integrated.RDS")
# Determine the clusters for various resolutions                                
Stromal_integrated <- FindClusters(object = Stromal_integrated,
                                   resolution = c(0.1,0.2,0.3,0.4,0.5,0.6))
View(Stromal_integrated@meta.data)  

# Assign identity of clusters
Idents(object = Stromal_integrated) <- "integrated_snn_res.0.2"

ggplot(data = Stromal_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,nUMI))+
  geom_boxplot()+
  ggplot(data = Stromal_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,nGene))+
  geom_boxplot()+
  ggplot(data = Stromal_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,mitoRatio))+
  geom_boxplot()+
  ggplot(data = Stromal_integrated@meta.data,mapping = aes(integrated_snn_res.0.2,RiboRatio))+
  geom_boxplot()+
  ggplot(data = Stromal_integrated@meta.data,mapping = aes(integrated_snn_res.0.1,log10GenesPerUMI))+
  geom_boxplot()

##Find markers
DefaultAssay(Stromal_integrated) <- "RNA"
n<-0
for (n in 0:5){
  conserved_markers <- FindConservedMarkers(Stromal_integrated, 
                                            ident.1 = n,
                                            grouping.var = "sample",
                                            verbose = T,only.pos=T) 
  write.csv(x=conserved_markers,file = paste("./result/Stromalcluster_",n,".csv",sep=''))
}

Stromal_integrated <- RenameIdents(object = Stromal_integrated, 
                                   "0" = "LepR+",
                                   "1" = "LepR+",
                                   "2" = "Osteoblast",
                                   "3" = "Osteocyte",
                                   "4" = "Endothelial cell",
                                   "5" = "Fibroblast")
my_levels <- c( "LepR+",
                "Osteoblast", "Osteocyte",
               "Endothelial cell",
                "Fibroblast")
Stromal_integrated@active.ident <- factor(x = Stromal_integrated@active.ident, levels = my_levels)

# Plot the UMAP
stromalumap<-DimPlot(Stromal_integrated,
                     reduction = "umap",
                     label = T,repel = T,
                     label.size = 5,ncol=1,
                     pt.size = 2)
ggsave(filename ="./result/20221106_stromalumap_0.pdf",plot =stromalumap,device = "pdf", 
       width =10,height = 10)

stromalVln<-VlnPlot(Stromal_integrated,ncol=4,
                    features=c("Lepr","Bglap","Bglap2","Spp1",
                               "Dmp1","Cdh5","Col3a1"),
                    split.by = "sample",split.plot = T,
                    cols = c("#12777f","#f35e5a"))

ggsave(filename = "./result/20221106_stromalVln.pdf",
       plot = stromalVln,device = "pdf",
       width =16,height = 8)

##heatmap
DefaultAssay(Stromal_integrated) <- "RNA"
Stromal_integrated$celltype<-Idents(Stromal_integrated)
HeatMapGenes<-read.xlsx(file = './data/20221106_HeatamapGenes copy.xlsx',sheetIndex = 1,header = T)
# StromalHMap<-DoHeatmap(subset(Stromal_integrated,downsample =100) ,assay = "SCT",
#                        features = HeatMapGenes$Genes, size = 3)

StromalHMap<-DoHeatmap(ScaleData(Stromal_integrated),
                features = HeatMapGenes$Genes, size = 3)
ggsave(filename ="./result/20221106_heatmapStromal.pdf",
       plot = StromalHMap,device = "pdf", 
       width =10,height = 10)

pDoMultiBarHeatmap<-DoMultiBarHeatmap(object = ScaleData(Stromal_integrated),
                  features=HeatMapGenes$Genes,group.by = "celltype",
                  additional.group.by = "sample",
                  label=F, ##required, otherwise error will pump out#
                  additional.group.sort.by = "sample",lines.width = 2)

ggsave(filename ="./result/20221106_heatmapStromal_BYPheatmap.pdf",
       pDoMultiBarHeatmap,device = "pdf", 
        width =10,height = 10)

###use scillus package to generate heatmap
#from https://scillus.netlify.app/
library(Scillus)
markers <- FindAllMarkers(Stromal_integrated, logfc.threshold = 0.1, min.pct = 0, only.pos = T)
pheatmap_Scillus<-plot_heatmap(dataset = ScaleData(Stromal_integrated,use_raster=F),
             n = 12,
             markers =HeatMapGenes$Genes,
             sort_var = c("celltype","sample"),
             anno_var = c("celltype","sample"),
             anno_colors = list(c("#f8766d","#a3a500","#00bf7d","#00b0f6","#e76bf3"),
                                c("#12777f","#f35e5a")),
             hm_limit = c(-1,0,1),
             hm_colors = c("purple","black","yellow")
             )

ggsave(filename ="./result/20221106_pheatmap_BYScillus.pdf",
       pheatmap_Scillus,device = "pdf", 
       width =10,height = 10)

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
stromaln_cells <- FetchData(Stromal_integrated, 
                            vars = c("ident", "seq_folder")) %>%
  dplyr::count(ident, seq_folder) %>%
  tidyr::spread(ident, n)

# View table
View(data.frame(t(stromaln_cells)))
write.csv(x =stromaln_cells,file = "./result/n_cells_Stromal.csv" )
table(Stromal_integrated@meta.data$orig.ident)

# UMAP of cells in each cluster by sample
p1<-DimPlot(Stromal_integrated, 
        label = TRUE, pt.size = 2,
        split.by = "sample",repel=T)  + NoLegend()
ggsave(filename ="./result/20221106_stromalumap_3.pdf",plot =p1,device = "pdf", 
       width =20,height = 10)
# Explore whether clusters segregate by cell cycle phase
DimPlot(Stromal_integrated,
        label = TRUE, 
        split.by = "Phase",repel = T)  + NoLegend()

# Determine metrics to plot present in seurat_integrated@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score", "mitoRatio","RiboRatio")
FeaturePlot(Stromal_integrated, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

DefaultAssay(Stromal_integrated)<-"RNA"
Nichefactorplot<-FeaturePlot(Stromal_integrated, 
                             reduction = "umap", 
                             features = c(#"Lepr","Pdgfra",
                                          "Cxcl12","Kitl"), 
                             label = TRUE,cols=c('gray90','red'),ncol=1,
                             #  min.cutoff = 'q1',
                             max.cutoff = 'q90',
                             label.size = 5,
                             pt.size = 2)
ggsave("./result/20221106_Nichefactorplot_3.pdf",Nichefactorplot,"pdf",
       width=20,height=40)

DefaultAssay(Stromal_integrated) <- "RNA"
LepRLinVln<-VlnPlot(object = Stromal_integrated, idents = "LepR+",
                    features =c("Lepr","Pdgfra","Cxcl12","Kitl",
                                "Adipoq","Lpl","Cebpa","Pparg",
                                "Col1a1","Bglap","Bglap2","Spp1",
                                "Sp7","Alpl","Acan","Col2a1",
                                "Dmp1",
                                "Sost",
                                "Sox9","Runx2","Klf2","Mettl3"),
                    split.by = "sample",ncol = 4,
                    log=F,cols = c("#12777f","#f35e5a"))
ggsave(filename ="./result/20221106_LepRLinVln.pdf",plot = LepRLinVln,device = "pdf", 
       width =8,height = 18)

Chondro_Vln<-VlnPlot(object = Stromal_integrated, idents = "LepR+",
                     features =c("Sox9","Col2a1","Acan"),assay = "RNA",
                     split.by = "sample",
                     log=T,cols = c("#12777f","#f35e5a"))
ggsave(filename ="./result/20221106_Chondro_Vln.pdf",plot = Chondro_Vln,device = "pdf", 
       width =12,height = 4)

##save h5ad file for scVelo use
#####1_1.save stromal1 with cluster infromation to be used for scVelo in Python
##http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html
DefaultAssay(Stromal_integrated) <- "RNA"
SaveH5Seurat(Stromal_integrated, filename = "./result/Velocyto/20221109_Stromal_postVelocyto.h5Seurat")
Convert("./result/Velocyto/20221109_Stromal_postVelocyto.h5Seurat", dest = "h5ad")

Stromal_integrated_Ctrl<-subset(Stromal_integrated,sample=="Control")
DefaultAssay(Stromal_integrated_Ctrl) <- "RNA"
SaveH5Seurat(Stromal_integrated_Ctrl, filename = "./result/Velocyto/20221109_Stromal_Ctrl_postVelocyto.h5Seurat")
Convert("./result/Velocyto/20221109_Stromal_Ctrl_postVelocyto.h5Seurat", dest = "h5ad")

Stromal_integrated_Mut<-subset(Stromal_integrated,sample=="Mutant")
DefaultAssay(Stromal_integrated_Mut) <- "RNA"
SaveH5Seurat(Stromal_integrated_Mut, filename = "./result/Velocyto/20221109_Stromal_Mut_postVelocyto.h5Seurat")
Convert("./result/Velocyto/20221109_Stromal_Mut_postVelocyto.h5Seurat", dest = "h5ad")


stromalctrlumap<-DimPlot(Stromal_integrated_Ctrl,
                     reduction = "umap",
                     label = T,repel = T,
                     label.size = 5,ncol=1,
                     pt.size = 2)
ggsave(filename ="./result/20221106_stromalCtrlumap_0.pdf",
       plot =stromalctrlumap,device = "pdf", 
       width =10,height = 10)

stromalmutumap<-DimPlot(Stromal_integrated_Mut,
                         reduction = "umap",
                         label = T,repel = T,
                         label.size = 5,ncol=1,
                         pt.size = 2)
ggsave(filename ="./result/20221106_stromalMutumap_0.pdf",
       plot =stromalmutumap,device = "pdf", 
       width =10,height = 10)

#####monocleTrajectory##############################
library("monocle3")
library("SeuratWrappers")

# a helper function to identify the root principal points:
LepR<-subset(Stromal_integrated, ident = c("LepR+","Osteoblast","Osteocyte"))
LepRumap<-DimPlot(LepR,
                  reduction = "umap",
                  label = T,repel = T,
                  label.size = 5,
                  pt.size = 0.5,split.by="sample")
LepRtsne<-DimPlot(LepR,
                  reduction = "tsne",
                  label = T,repel = T,
                  label.size = 5,
                  pt.size = 0.5)
LepR$celltype<-Idents(LepR)

LepR_integrated.cds <- as.cell_data_set(LepR)
LepR_integrated.cds <- cluster_cells(cds =LepR_integrated.cds, reduction_method = "UMAP")
plot_cells(LepR_integrated.cds, color_cells_by="celltype", 
           reduction_method = "UMAP")

LepR_integrated.cds <- learn_graph(LepR_integrated.cds, use_partition = TRUE)

LepR_integrated.cds<-order_cells(LepR_integrated.cds)

LepRMonoplot1<-plot_cells(LepR_integrated.cds,
                          color_cells_by = "pseudotime",
                          label_groups_by_cluster=FALSE,
                          label_leaves=FALSE,
                          label_branch_points=FALSE,cell_size = 1)
ggsave(filename ="./result/Velocyto/20221109_Monocle3_LepR.pdf",plot = LepRMonoplot1,device = "pdf", 
       width =10,height = 8)

LepR<-AddMetaData(
  object = LepR,
  metadata = LepR_integrated.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "pseudotime")

LepRmono<-FeaturePlot(LepR,"pseudotime",pt.size = 1,label = T,repel=T) & scale_color_viridis_c()
ggsave(filename ="./result/Velocyto/20221109_Monocle3_LepR_Pseudotime.pdf",plot = LepRmono,device = "pdf", 
       width =10,height = 8)
##LepRctrl
LepRctrl<-subset(LepR,sample=="Control")

LepRctrl_integrated.cds <- as.cell_data_set(LepRctrl)
LepRctrl_integrated.cds <- cluster_cells(cds =LepRctrl_integrated.cds, reduction_method = "UMAP")
plot_cells(LepRctrl_integrated.cds, color_cells_by="celltype", 
           reduction_method = "UMAP")

LepRctrl_integrated.cds <- learn_graph(LepRctrl_integrated.cds, use_partition = TRUE)

LepRctrl_integrated.cds<-order_cells(LepRctrl_integrated.cds)


LepRctrlMonoplot1<-plot_cells(LepRctrl_integrated.cds,
                              color_cells_by = "pseudotime",
                              label_groups_by_cluster=FALSE,
                              label_leaves=FALSE,
                              label_branch_points=FALSE,cell_size = 1)
ggsave(filename ="./result/Velocyto/20221109_Monocle3_LepRctrl.pdf",plot = LepRctrlMonoplot1,device = "pdf", 
       width =10,height = 8)

LepR<-AddMetaData(
  object = LepR,
  metadata = LepRctrl_integrated.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "ctrlpseudotime")
LepRctrl<-AddMetaData(
  object = LepRctrl,
  metadata = LepRctrl_integrated.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "ctrlpseudotime")

LepRctrlmono<-FeaturePlot(LepRctrl,"ctrlpseudotime",pt.size = 1,label = T,repel=T) & scale_color_viridis_c()
ggsave(filename ="./result/Velocyto/20221109_Monocle3_LepRctrl_Pseudotime.pdf",plot = LepRctrlmono,device = "pdf", 
       width =10,height = 8)

##LepRmut
LepRmut<-subset(LepR,sample=="Mutant")

LepRmut_integrated.cds <- as.cell_data_set(LepRmut)

LepRmut_integrated.cds <- cluster_cells(cds =LepRmut_integrated.cds, reduction_method = "UMAP")
plot_cells(LepRmut_integrated.cds, color_cells_by="celltype", 
           reduction_method = "UMAP")

LepRmut_integrated.cds <- learn_graph(LepRmut_integrated.cds, use_partition = TRUE)

LepRmut_integrated.cds<-order_cells(LepRmut_integrated.cds,reduction_method = "UMAP")

LepRmutMonoplot1<-plot_cells(LepRmut_integrated.cds,
                             color_cells_by = "pseudotime",
                             label_groups_by_cluster=FALSE,
                             label_leaves=FALSE,
                             label_branch_points=FALSE,cell_size = 1)
ggsave(filename ="./result/Velocyto/20221109_Monocle3_LepRmut.pdf.pdf",plot = LepRmutMonoplot1,device = "pdf", 
       width =10,height = 8)

LepR<-AddMetaData(
  object = LepR,
  metadata = LepRmut_integrated.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "mutpseudotime")

LepRmut<-AddMetaData(
  object = LepRmut,
  metadata = LepRmut_integrated.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "mutpseudotime")

LepRmutmono<-FeaturePlot(LepRmut,"mutpseudotime",pt.size = 1,label = T,repel=T) & scale_color_viridis_c()
ggsave(filename ="./result/Velocyto/20221109_Monocle3_LepRmut_Pseudotime.pdf",
       plot = LepRmutmono,device = "pdf", 
       width =10,height = 8)


DefaultAssay(LepR) <- "RNA"
SaveH5Seurat(LepR, filename = "./result/Velocyto/20221109_LepR_postVelocyto.h5Seurat")
Convert("./result/Velocyto/20221109_LepR_postVelocyto.h5Seurat", dest = "h5ad")

DefaultAssay(LepRctrl) <- "RNA"
SaveH5Seurat(LepRctrl, filename = "./result/Velocyto/20221109_LepRctrl_postVelocyto.h5Seurat")
Convert("./result/Velocyto/20221109_LepRctrl_postVelocyto.h5Seurat", dest = "h5ad")

DefaultAssay(LepRmut) <- "RNA"
SaveH5Seurat(LepRmut, filename = "./result/Velocyto/20221109_LepRmut_postVelocyto.h5Seurat")
Convert("./result/Velocyto/20221109_LepRmut_postVelocyto.h5Seurat", dest = "h5ad")

# load required packages

req_packages = c("Biobase", "cluster", "clusterProfiler", "cowplot",
                 "data.table", "edgeR", "ggpubr", 
                 "ggrepel", 
                 "GO.db", "goseq", "grid", "gridExtra", "lattice", 
                 "pheatmap", "qvalue", "RColorBrewer", 
                 "Rmisc", "RUVSeq","splitstackshape", "statmod", "stringr", "tidyverse",
                 "VennDiagram", "viridis")
## load them, quietly
invisible(suppressWarnings(suppressMessages(
    lapply(req_packages, require, character.only = TRUE)
)))

## The Cowplot package changes the default themes of ggplot2. Set the bw theme with larger font sizes like so:
theme_set(theme_bw(base_size = 16))
## ... or set the default theme
# theme_set(theme_gray())

## suppress excessive VennDiagram log files
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# load featureCounts_file

cm <- read.csv("featureCounts_gene_counts_matrix_mod.txt", header=T, row.names=1, sep='\t', check.names = T)
head(cm)

colnames(cm) %>% 
    as.data.frame %>% 
    rename(sample = ".") %>% 
    separate(sample, into = c("cross", "timepoint", "replicate"), remove = FALSE) %>% 
    mutate(treatment = if_else(cross == "VA", "heterospecific", if_else(cross == "VV", "conspecific", "unmated")),
          sample_name = paste(cross, timepoint, sep = "_")) -> sample_info

nrow(cm)

cpm_cm <- cpm(cm)
thresh_cm <- cpm_cm > 2
keep_cm <- rowSums(thresh_cm) >= 3
countsKeep_cm <- cm[keep_cm,]
# countsKeep_qM <- subset(countsKeep_qM, rownames(countsKeep_qM) %!in% wierd_genes)
table(keep_cm)

groups = factor(sample_info$sample_name)
design = model.matrix(~0+groups)
colnames(design) <- levels(groups)
rownames(design) <- sample_info$sample


#create DGEList
dgeList <- DGEList(counts = countsKeep_cm, group = groups)
dgeList <- calcNormFactors(dgeList)
dgeList <- estimateDisp(dgeList, design)
# dgeList <- estimateGLMTagwiseDisp(dgeList, design)
dgeList_fit <- glmQLFit(dgeList, design)
summary(dgeList$tagwise.dispersion)

options(repr.plot.width = 12, repr.plot.height = 9)
par(mfrow=c(2,2))
# Biological coefficient of variation
plotBCV(dgeList)
# mean-variance trend
voom = voom(dgeList, design, plot=TRUE)
# QQ-plot
g.v <- gof(dgeList_fit, plot = T, pcutoff = 0.05, adjust = "holm")
# z.v <- zscoreGamma(g.v
df/2,scale=2)
# qqnorm(z.v); qqline(z.v, col = 4,lwd=1,lty=1)
# log2 transformed and normalize boxplot of counts across samples
boxplot(voom$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(voom$E),col="blue")

rm(voom)




## Plot sample correlation
data = log2(RUVrNormalizedCounts+1)
# colnames(data) = gsub("Female_", "", colnames(data))
data = as.matrix(data)
sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')

options(repr.plot.width = 14, repr.plot.height = 18)
pheatmap(
  mat               = sample_cor,
  color             = inferno(50),
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
#   filename          = "Figures/sample_correlations_all.pdf",
#   width             = 6, 
#   height            = 5,
  fontsize          = 12    
)
rm(data)
rm(sample_cor)
     

plotMDS(dgeList)

mdsObj_Ae = merge(mdsObj_Ae, sample.info_Amaro_Ae_pre, by.x = "replicate", by.y = "replicate")
mdsObj_Ae$treat = gsub("none", "Non-injected", mdsObj_Ae
treat = gsub("MAG", "MAG extract", mdsObj_Ae$treat)
mdsObj_Ae %<>% separate(replicate, c("rep_id", "samp_id")) %>% mutate(rep_id = gsub("Ae", "", rep_id))

options(repr.plot.width = 4.5, repr.plot.height = 2.5)
(ggscatter(filter(mdsObj_Ae, timepoint != "6hpi"),
              x = "dim1", 
              y = "dim2",
              color = "treat",
#               shape = "timepoint",
              size = 4,
              alpha = 0.8, 
#               ellipse = T,
              ggtheme = theme_bw(),
              repel = "timepoint",
         ) + 
                theme(axis.text = element_text(size = 10), legend.title = element_blank(), axis.title = element_text(size = 12), legend.text = element_text(size = 12)) +
#                 theme_monokai_full() +
                labs ( x = "Dimension 1", y = "Dimension 2") +
                geom_text_repel(aes(label = rep_id), box.padding = 0.5) +
#                 facet_wrap(~source, scale = "free") +
                scale_colour_manual(values=c("#007d30", "black", "#ffad3d")) -> preBatch_mds)
     ##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

# TPM function
tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
}

##------------------------------------------------------------------------##
##------------------------------------------------------------------------##

## Miscellaneous operators
'%!in%' <- function(x,y)!('%in%'(x,y))

##------------------------------------------------------------------------##
##---------------------

                       set <- newSeqExpressionSet(as.matrix(countsKeep_cm), phenoData = data.frame(groups, row.names = colnames(countsKeep_cm)))
set <- betweenLaneNormalization(set, which="upper")

y <- DGEList(counts=counts(set), group=groups)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateDisp(y, design, robust = T)
fit <- glmQLFit(y, design, dispersion = y$tagwise.dispersion, robust = T)
res <- residuals(fit, type="deviance")

batch_ruv_res = RUVr(set,rownames(countsKeep_cm),k=10,res)
RUVrNormalizedCounts = normCounts(object = batch_ruv_res)
rownames(RUVrNormalizedCounts) = rownames(countsKeep_cm)

# RUVrNormalizedTPM = tpm(counts = RUVrNormalizedCounts, lengths = 2922)

# tmp.tpmMatrix<-RUVrNormalizedTPM
# tmp.tpmMatrix.m <- as.data.frame(melt(as.matrix(tmp.tpmMatrix)))
# colnames(tmp.tpmMatrix.m) <- c("gene_id", "replicate", "TPM")
# Amaro_Ae.geneTPM.table.rc <- merge(sample_info, tmp.tpmMatrix.m, by.x = "replicate", by.y = "sample", all.y = T)
# rm(tmp.tpmMatrix, tmp.tpmMatrix.m)
     

nonRUVrNormalizedTPM = tpm(counts = cm, lengths = 2922)
     

design_2 <- model.matrix(~ 0 + groups + W_1 + W_2 + W_3 + W_4 + W_5 + W_6 + W_7 + W_8 + W_9 + W_10, data=pData(batch_ruv_res))
colnames(design_2) <- gsub("groups", "", colnames(design_2))

options(repr.plot.width = 5, repr.plot.height = 3)
plot(batch_ruv_res$W_9,pch=19,main="RUVseq residuals")
     

     
     

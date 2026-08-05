## Differential expression analysis of TEs ##
## Takes raw count tables from quantification, one for each sample, named after each sample, with .tsv file extension ##

# Load packages
library(tidyverse)
library(DESeq2)
library("pheatmap")
library(apeglm)
library(scales)

## Load in the count tables, filter, and make into matrix ##

# Check directory
getwd()
#setwd("./project")

# Read in the sample metadata
metadata <- read.csv("../data/metadata/mouse_TE_metadata.csv")

# Inspect
metadata
dim(metadata)
str(metadata)

# Read in count data
folder_path <- "../data/counts/"
files_path <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = T) 
files_list <- lapply(files_path, read_tsv)

# Make into df (joined by transcript then matrix)
str(files_list)
count_data_allsamples <- files_list %>% purrr::reduce(full_join, by='transcript') 

# Remove short repeats/low complexity sequences (have parentheses in them)- before maing gene id the row id so there is a column to filter on
# Also have to filter strings with "rich" as well as these are low complexity 
str(count_data_allsamples)
dim(count_data_allsamples)
head(count_data_allsamples)

# how many have ()?
count_data_allsamples %>%
  filter(str_detect(transcript, "^\\(")) %>%
  nrow() #148686 are short repeats

count_data_allsamples %>%
  filter(str_detect(transcript, "rich")) %>%
  nrow() #24929 are low complexity

# Together
count_data_allsamples %>%
  filter(str_detect(transcript, "^\\(|rich")) %>%
  nrow() # 173615 is the total of the two above

# Make filtered dataframe by negating match
count_data_allsamples_filtered <- count_data_allsamples %>%
  filter(!(str_detect(transcript, "^\\(|rich")))

# Check
head(count_data_allsamples_filtered)
dim(count_data_allsamples_filtered) # 1381748 TEs remaining

# make the Geneid the row id 
count_data_allsamples_filtered <- count_data_allsamples_filtered %>% remove_rownames %>% column_to_rownames(var="transcript")
head(count_data_allsamples_filtered)
dim(count_data_allsamples_filtered)

# Make sample names the row names
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="sample") 
all(colnames(count_data_allsamples_filtered) %in% rownames(metadata)) # Check they both have all samples

# Make rownames/colnames in same order
count_data_allsamples_filtered <- count_data_allsamples_filtered[, rownames(metadata)]
all(rownames(metadata) == colnames(count_data_allsamples_filtered)) # True

# Replace all NA values with 0 in your raw count matrix
count_data_allsamples_filtered[is.na(count_data_allsamples_filtered)] <- 0

# Create deseq object
dds <- DESeqDataSetFromMatrix(countData = count_data_allsamples_filtered,
                              colData = metadata,
                              design = ~ sex + treatment + sex:treatment) # add interaction term

# Prefilter based on smallest group size
smallestGroupSize <- 5 # sex/treatment
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize # Suggested low count
dds <- dds[keep,] # 50333 genes remaining

# Relevel just in case, Sham will be the reference
dds$treatment <- relevel(dds$treatment, ref = 'SHAM')

# Normalise 
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = T)

# Write normalised counts out to file
#write.table(normalized_counts, file="../data/normalised_counts_allsamples.tsv", sep = '\t', quote=F, col.names = NA)

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
summary(res) # with the interaction term there are 361 genes upreg and 547 downreg

# write out results
#write.csv(as.data.frame(res), file="../data/DE_TEs_plus.csv")

# Apply LFC shrinkage 
resultsNames(dds) # we are interested in sham vs sni which is coef 3
resNorm <- lfcShrink(dds, coef=3, type="apeglm") # apeglm is recommended
class(resNorm)
summary(resNorm, alpha = 0.05) # 188 upregulated, 143 downregulated, 7 outliers

# MA plot
plotMA(resNorm, main="apeglm")

# Convert the shrunken/interaction term exclusive results to df
# Visualising results
DE_TEs_int <- data.frame(resNorm)
head(DE_TEs_int)
dim(DE_TEs_int)

DE_TEs_int %>%
  filter(padj < 0.05) %>%
  nrow() # 331

# Write these out as a dataframe
DE_TEs_int_signif_all <- DE_TEs_int %>% 
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))

# Saved again 27.07 after filtering low complexity as well
#write.csv(DE_TEs_int_signif_all, "../results_data/signif_DE_TEs.csv")

## Visualisation
# rlog transform for purposes of PCA transformation
rld <- rlog(dds)
pca_dat <- plotPCA(rld, intgroup=c("treatment", "sex"), returnData = T) # Get PCA data for plot

# Get % for each PCA
plotPCA(rld, intgroup=c("treatment", "sex"))
#PC1 is 15% and PC2 is 9%

pca_plot <- ggplot(pca_dat, aes(x = PC1, y = PC2, colour = treatment, shape = sex)) + geom_point(size = 4) + theme_classic() + scale_shape_discrete(name = "Sex") + scale_colour_manual(name = "Treatment", values = c("orangered1", "skyblue2"), labels = c("Control","SNI")) + xlab("PC1: 15% variance") + ylab("PC2: 9% variance") + labs(title = "PCA of TE expression")

# Saved 28/07
#ggsave("../figures/TE_PCA_plot.png", pca_plot, dpi=300)

# Dispersion plot
plotDispEsts(dds)

# Volcano plot

# Make label vector
# Choosing a lfc threshold of 1 for plotting
#28.7- change to 0.5 as per Ana

DE_TEs_int$signif <- ifelse(DE_TEs_int$padj < 0.05 & abs(DE_TEs_int$log2FoldChange) > 0.5, ifelse(DE_TEs_int$log2FoldChange > 0, "Up", "Down"), "NotSig")

# Count how many there are:
DE_TEs_int %>%
  group_by(signif) %>%
  summarise(n())
# 1 Down      81
#2 NotSig 22042
#3 Up       168


# Omit NAs, they are low counts and outliers
DE_TEs_int %>% filter(if_any(everything(), is.na)) %>%
  nrow() # 26483
DE_TEs_int <- na.omit(DE_TEs_int)

# Plot
volcanoplot <- ggplot(DE_TEs_int, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(size = 1)+labs(x = "Log2 Fold Change", y = "-Log10 p-value",title = "Differentially Expressed TEs") + scale_colour_manual(labels = c("Downregulated: 81", "Not significant: 22042", "Upregulated: 168"), values = c("Down" = "blue", "NotSig" = "darkgrey", "Up" = "red")) + labs(color = stringr::str_wrap("Expression: p < 0.05 and Log2FC > 0.5", 20)) + theme_classic() + theme(legend.position = "right", legend.text = element_text(size = 8), legend.title = element_text(size = 9)) + scale_x_continuous(limits = c(-3, 2)) + scale_y_continuous(limits = c(0, 6)) + theme(legend.position = "none")

# Save plot, saved 30/07
#ggsave("../figures/TE_volcano_plot.png", volcanoplot, dpi=300)

# Heatmap
# Top 20 TEs for clarity
select <- order(rowMeans(counts(dds,normalized=T)),
                decreasing=T)[1:20]

df <- as.data.frame(colData(dds)[,c("treatment","sex")]) # To label

# Transform for plotting
ntd <- normTransform(dds)

# needs work on colours, etc
heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=T, show_rownames=F, show_colnames = T, cluster_cols=T, annotation_col=df, fontsize_row = 6, fontsize_col = 6)

# Try without clustering rows just to see
pheatmap(assay(ntd)[select,], cluster_rows=F, show_rownames=F, show_colnames = T, cluster_cols=F, annotation_col=df, fontsize_row = 6, fontsize_col = 6)

# Given the separation of sex on PCA- redo analysis on each sex

# Read in metadata again for dividing by sex

metadata <- read.csv("../data/metadata/mouse_TE_metadata.csv")

females <- metadata %>% 
  filter(sex=="Female") %>%
  remove_rownames %>% 
  column_to_rownames(var="sample") 

males <- metadata %>% 
  filter(sex=="Male") %>%
  remove_rownames %>% 
  column_to_rownames(var="sample") 

females
males

# Make count matrix for males and females
counts_filtered_female <- count_data_allsamples_filtered[, rownames(females)]
counts_filtered_male <- count_data_allsamples_filtered[, rownames(males)]

#Check
all(rownames(males) == colnames(counts_filtered_male)) # True
all(rownames(females) == colnames(counts_filtered_female)) # True

# Replace all NA values with 0 in your raw count matrix
counts_filtered_female[is.na(counts_filtered_female)] <- 0
counts_filtered_male[is.na(counts_filtered_male)] <- 0


# Create deseq object
ddsfemale <- DESeqDataSetFromMatrix(countData = counts_filtered_female,
                                    colData = females,
                                    design = ~ treatment) # Only treatment now

ddsmale <- DESeqDataSetFromMatrix(countData = counts_filtered_male,
                                  colData = males,
                                  design = ~ treatment) # Only treatment now

# Prefilter based on smallest group size
smallestGroupSize <- 5 # treatment
keep <- rowSums(counts(ddsmale) >= 10) >= smallestGroupSize
keep <- rowSums(counts(ddsfemale) >= 10) >= smallestGroupSize# Suggested low count
ddsmale <- ddsmale[keep,] # 34881
ddsfemale <- ddsfemale[keep,] # 34881 


# Relevel just in case, Sham will be the reference (control)
ddsmale$treatment <- relevel(ddsmale$treatment, ref = 'SHAM')
ddsfemale$treatment <- relevel(ddsfemale$treatment, ref = 'SHAM')

# Normalise 
ddsmale <- estimateSizeFactors(ddsmale)
normalised_countsmale <- counts(ddsmale, normalized = T)
ddsfemale <- estimateSizeFactors(ddsfemale)
normalised_countsfemale <- counts(ddsfemale, normalized = T)

# Write normalised counts out to file
#write.table(normalized_counts, file="../data/normalised_counts_allsamples.tsv", sep = '\t', quote=F, col.names = NA)

# Differential expression analysis
ddsmale <- DESeq(ddsmale)
ddsfemale <- DESeq(ddsfemale)
resmale <- results(ddsmale, alpha = 0.05)
resfemale <- results(ddsfemale, alpha = 0.05)
summary(resmale) # 1 upregulated, 4 downregulated
summary(resfemale) # 320 upregulated, 151 downregulated

# write out results
#write.csv(as.data.frame(res), file="../data/DE_TEs_plus.csv")

# Apply LFC shrinkage 
resultsNames(ddsmale) # we are interested in sham vs sni which is coef 2
resultsNames(ddsfemale)
resNormmale <- lfcShrink(ddsmale, coef=2, type="apeglm") # apeglm is recommended
resNormfemale <- lfcShrink(ddsfemale, coef=2, type="apeglm") # apeglm is recommended

summary(resNormmale, alpha = 0.05) # 0 upregulated, 4 downregulated, 22 outliers
summary(resNormfemale, alpha = 0.05) # 302 upregulated, 161 downregulated, 9 outliers

# MA plot
plotMA(resNormmale, main="apeglm")
plotMA(resNormfemale, main="apeglm")

# Convert the shrunken/interaction term exclusive results to df
# Visualising results
DE_TEs_male <- data.frame(resNormmale)
DE_TEs_female <- data.frame(resNormfemale)

head(DE_TEs_male)
head(DE_TEs_female)
dim(DE_TEs_male)
dim(DE_TEs_female)

DE_TEs_male %>%
  filter(padj < 0.05) %>%
  nrow() # 4

DE_TEs_female %>%
  filter(padj < 0.05) %>%
  nrow() # 463

# Write these out as a dataframe
female_TEs_signif <- DE_TEs_female %>% 
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))

male_TEs_signif <- DE_TEs_male %>% 
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))

# Saved again 28.07 
#write.csv(female_TEs_signif, "../results_data/female_TEs_signif.csv")
#write.csv(male_TEs_signif, "../results_data/male_TEs_signif.csv")

## Visualisation
# rlog transform for purposes of PCA transformation
rldmale <- rlog(ddsmale)
pca_datmale <- plotPCA(rldmale, intgroup=("treatment"), returnData = T)

rldfemale <- rlog(ddsfemale)
pca_datfemale <- plotPCA(rldfemale, intgroup=("treatment"), returnData = T)

# Get % for each PCA
plotPCA(rldfemale, intgroup="treatment")
#PC1 is 19% and PC2 is 16%
plotPCA(rldmale, intgroup="treatment")
#PC1 is 20% and PC2 is 14%

pca_plotmale <- ggplot(pca_datmale, aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(size = 4) + theme_classic() + scale_colour_manual(name = "Treatment", values = c("orangered1", "skyblue2")) + xlab("PC1: 20% variance") + ylab("PC2: 14% variance") + labs(title = "PCA of TE expression in males") 

pca_plotfemale <- ggplot(pca_datfemale, aes(x = PC1, y = PC2, colour = treatment)) +
  geom_point(size = 4) + theme_classic() + scale_colour_manual(name = "Treatment", values = c("orangered1", "skyblue2")) + xlab("PC1: 19% variance") + ylab("PC2: 16% variance") + labs(title = "PCA of TE expression in females") 

# Saved 28/07
#ggsave("../figures/TE_PCA_plotmale.png", pca_plotmale, dpi=300)
#ggsave("../figures/TE_PCA_plotfemale.png", pca_plotfemale, dpi=300)

# Dispersion plot
plotDispEsts(ddsmale)
plotDispEsts(ddsfemale)

# Volcano plot

# Make label vector
# Choosing a lfc threshold of 0.5
DE_TEs_female$signif <- ifelse(DE_TEs_female$padj < 0.05 & abs(DE_TEs_female$log2FoldChange) > 0.5, ifelse(DE_TEs_female$log2FoldChange > 0, "Up", "Down"), "NotSig")
DE_TEs_male$signif <- ifelse(DE_TEs_male$padj < 0.05 & abs(DE_TEs_male$log2FoldChange) > 0.5, ifelse(DE_TEs_male$log2FoldChange > 0, "Up", "Down"), "NotSig")

# Count how many there are:
DE_TEs_male %>%
  group_by(signif) %>%
  summarise(n())
#1 Down       3
#2 NotSig 34797
#3 NA        81

DE_TEs_female %>%
  group_by(signif) %>%
  summarise(n())
#Down      72
#2 NotSig 34198
#3 Up       255
#4 NA       356

# Omit NAs, they are low counts and outliers
DE_TEs_male  %>% filter(if_any(everything(), is.na)) %>%
  nrow() # 28412
DE_TEs_female %>% filter(if_any(everything(), is.na)) %>%
  nrow() #16916

DE_TEs_male <- na.omit(DE_TEs_male)
DE_TEs_female <- na.omit(DE_TEs_female)

# write out the significant TEs as a column

TE_names_female <- DE_TEs_female %>%
  rownames_to_column(var = "TE_names") %>%
  filter(!signif=="NotSig") %>%
  dplyr::select(TE_names) %>%
  as.data.frame()

TE_names_male <- DE_TEs_male %>%
  rownames_to_column(var = "TE_names") %>%
  filter(!signif=="NotSig") %>%
  dplyr::select(TE_names) %>%
  as.data.frame()

# write out, 4/8
#write.table(TE_names_female, "../results_data/signifTEs_names/TE_names_female.txt", row.names = F, col.names = F, quote = F)
#write.table(TE_names_male, "../results_data/signifTEs_names/TE_names_male.txt", row.names = F, col.names = F, quote = F)

# Plot
femalevolcanoplot <- ggplot(DE_TEs_female, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(size = 1)+labs(x = "Log2 Fold Change", y = "-Log10 p-value",title = "Differentially Expressed TEs in Females") + scale_colour_manual(labels = c("Downregulated: 72", "Not significant: 34198", "Upregulated: 255"), values = c("Down" = "blue", "NotSig" = "darkgrey", "Up" = "red")) + theme_classic() + theme(legend.position = "top", legend.text = element_text(size = 7), legend.title = element_text(size = 9))+ scale_x_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2)) + scale_y_continuous(limits = c(0, 7)) + labs(title = NULL, colour = "Expression: p < 0.05 and Log2FC > 0.5") + guides(color = guide_legend(title.position = "top")) 

malevolcanoplot <- ggplot(DE_TEs_male, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(size = 1)+labs(x = "Log2 Fold Change", y = "-Log10 p-value",title = "Differentially Expressed TEs in Males") + scale_colour_manual(labels = c("Downregulated: 3", "Not significant: 34797"), values = c("Down" = "blue", "NotSig" = "darkgrey", "Up" = "red")) + theme_classic() + theme(legend.position = "top", legend.text = element_text(size = 7), legend.title = element_text(size = 9)) + scale_x_continuous(limits = c(-3, 2)) + scale_y_continuous(limits = c(0, 6)) + labs(title = NULL, colour = "Expression: p < 0.05 and Log2FC > 0.5") + guides(color = guide_legend(title.position = "top"))

# Save plot, saved 31/07
#ggsave("../figures/male_volcano_plot.png", malevolcanoplot, dpi=300)
#ggsave("../figures/female_volcano_plot.png", femalevolcanoplot, dpi=300)

## End ##
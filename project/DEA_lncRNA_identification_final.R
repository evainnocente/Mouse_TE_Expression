# New script #

# Load packages
library(tidyverse)
library(DESeq2)
library("pheatmap")
#install.packages("BiocManager")
#BiocManager::install("biomaRt")
library(biomaRt)

# Check directory
getwd()
setwd("./project")

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

pca_plot <- ggplot(pca_dat, aes(x = PC1, y = PC2, colour = treatment, shape = sex)) +
  geom_point(size = 4) + theme_classic() + scale_shape_discrete(name = "Sex") + scale_colour_manual(name = "Treatment", values = c("orangered1", "skyblue2")) + xlab("PC1: 15% variance") + ylab("PC2: 9% variance") + labs(title = "PCA of TE expression") 

# Saved 27/07
#ggsave("../figures/TE_PCA_plot.png", pca_plot, dpi=300)

# Dispersion plot
plotDispEsts(dds)

# Volcano plot

# Make label vector
# Choosing a lfc threshold of 1 for plotting

DE_TEs_int$signif <- ifelse(DE_TEs_int$padj < 0.05 & abs(DE_TEs_int$log2FoldChange) > 1, ifelse(DE_TEs_int$log2FoldChange > 0, "Up", "Down"), "NotSig")

# Omit NAs, they are low counts and outliers
DE_TEs_int %>% filter(if_any(everything(), is.na)) %>%
  nrow() # 26483
DE_TEs_int <- na.omit(DE_TEs_int)

# Plot
volcanoplot <- ggplot(DE_TEs_int, aes(x = log2FoldChange, y = -log10(padj), color = signif)) +
  geom_point(size = 1)+labs(x = "Log2 Fold Change", y = "-Log10 p-value",title = "Differentially Expressed TEs") +theme(legend.position = "right")+ scale_colour_manual(name = "Expression", labels = c("Downregulated: 143", "Not significant", "Upregulated: 188"), values = c("Down" = "blue", "NotSig" = "darkgrey", "Up" = "red")) + theme_classic()

# Save plot, saved 26/07
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

# Try withoout clustering rows just to see
pheatmap(assay(ntd)[select,], cluster_rows=F, show_rownames=F, show_colnames = T, cluster_cols=F, annotation_col=df, fontsize_row = 6, fontsize_col = 6)

# Save plot, rough draft save 23/07, not saving again
#ggsave("../figures/TE_heatmap.png", heatmap, dpi=300)

# Mapping significantly differentially expressed TEs to see if any are lncRNAs
  
# remove low complexity/simple repeats: starting with a ()
# No longer needed as I did it upstream
#DE_TEs_int_signif_df %>%
  #filter(str_detect(TE_names, "^\\(")) %>%
  #nrow() # 221 out of 706

#DE_signif_TEs_only_ranked_filtered <- DE_TEs_int_signif_df %>%
  #filter(!(str_detect(TE_names, "^\\("))) # 485

# Searching for genes
listEnsembl() 

# Make my biomart object
ensembl <- useEnsembl(biomart = "genes")

# Search mouse genes only
searchDatasets(mart = ensembl, pattern = "mmusculus")

#mmusculus_gene_ensembl	is the dataset I want
 
# Make the biomart object
ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", mart = ensembl)

# Get TE names from written out datatset
# load back in the dataframe and turn rownames to TE_names
DE_TEs_int_signif_all
DE_TEs_int_signif_df <- DE_TEs_int_signif_all %>%
  rownames_to_column(var = "TE_names")

# remove everything after first underscore inclusive as this is the chromosome name and position left over from earlier- potentially could impact the information
TE_names_only <- DE_TEs_int_signif_df %>%
  mutate(TEs = str_remove(TE_names, "_.*")) %>%
  dplyr::select(TEs)

# Search the database using my list of TEs as keys (external gene name), return the transcript or gene biotype which should tell you whether it is a lncRNA
# Need to return chromosome coordinates too so I know which Bc1 it is

lncRNA <- getBM(attributes = c("external_gene_name", "transcript_biotype", "gene_biotype", "chromosome_name", "start_position", "end_position"), filters = c("external_gene_name"), values = TE_names_only, mart = ensembl)

# one lncRNA was returned: Bc1
lncRNA

# Save results, saved 26/7
#write.csv(lncRNA, "../results_data/lncRNAs.csv")

## End ##

## Still in progress to get correct output ##

# output the TE names with chr etc attached still
TE_names_for_matching <- DE_TEs_int_signif_df %>%
  dplyr::select(TE_names) %>%
  as.data.frame()
#export
#write.table(TE_names_for_matching, "../results_data/TE_names_formatching.txt", row.names = F, col.names = F, quote = F)

# Load necessary libraries
library(affy)
library(hgu133a.db)
library(AnnotationDbi)
library(limma)
library(vsn)

## Data from E-GEOD-12288 - Gene expression patterns in peripheral blood correlate with the
extent of coronary artery disease

##AIM: To use statistical methods to find differential expression, and rank them based on p-value.

# Set your working directory to the folder containing the CEL files 
setwd("C:/Users/Priyamvadha/Documents/Analytics/Assignment-2/Data")

# List the CEL files in the directory
cel_files <- list.files(pattern = ".CEL$", full.names = TRUE)
# Read the CEL files and create an AffyBatch object
data <-read.affybatch(dir(patt="CEL"))
clinical <- read.csv("experiment.csv", header=T)
# Perform RMA normalization
data_rma <- affy::rma(data)
# View the normalized data
exprs_data <- exprs(data_rma)
head(exprs_data)  # Display the first few rows of the expression matrix

##t-test
t2 <- vector()
pval.t2 <- vector()
group <- clinical$group
#iteration
for(j in 1:nrow(exprs_data))
{
  temp <- exprs_data[j,]
  res <- t.test(temp[group==1], temp[group==0], var.equal=T)
  t2[j] <- res$stat
  pval.t2[j] <- res$p.val
  }
#Adjusting p value
adj.pval.t2 <- p.adjust(pval.t2, "BH") ##Benjamini-Hochberg correction
result.table2 <- data.frame(ID=rownames(exprs_data), t.stat=t2,
                           pvalue=pval.t2, fdr.pvalue=adj.pval.t2)
result.table2.sorted <- result.table2[order(adj.pval.t2),]
result.table2.sorted[1:10,] # listing the top 10 genes
top10_result <- result.table2.sorted[1:10,]
probeIDs <- top10_result$ID ## separating probe IDs
##Mapping Probe Id to gene
mappedProbeToGene <- select(hgu133a.db, keys=probeIDs, columns = c("PROBEID", "ENSEMBL", "SYMBOL"))
##Final tabular output
final_ttest <- dplyr::left_join(top10_result, mappedProbeToGene, by= c("ID" = "PROBEID"))
print(final_ttest)

#write.csv(mappedProbeToGene, "top50.csv", row.names =  FALSE)

####ATTEMPT USING eBAYES
mod_clinical <- clinical
rownames(mod_clinical) <- mod_clinical$Array.Data.File
mod_clinical$Array.Data.File <- NULL
mod_clinical

fold_change_threshold <- 1.5
fit1 <- lmFit(exprs_data, mod_clinical)
fit1 <- eBayes(fit1)
results_eBayes <- decideTests(fit1, method = "global", adjust.method = "fdr", p.value = 0.5, lfc = 2)
adj_pval <- p.adjust(fit1$p.value, method = "BH")
significant_genes <- rownames(exprs_data)[adj_pval <- 0.5]
significant_genes

data_vsn <- vsn2(data)

mean_data <- apply(data_vsn, 1, mean)
sd_data <- apply(data_vsn, 1, sd) 

# Create a data frame with mean and standard deviation
mean_sd_df <- data.frame(Mean = mean_data, SD = sd_data)

meanSdPlot(mean_sd_df)

adj.pval.t2_bonferroni <- p.adjust(pval.t2, "bonferroni")
adj.pval.t2_holm <- p.adjust(pval.t2, "holm")

#bonferroni
result.table2_bon = data.frame(ID=rownames(exprs_data), t.stat=t2,
                               pvalue=pval.t2, fdr.pvalue=adj.pval.t2_bonferroni)
result.table2.sorted_bon = result.table2_bon[order(adj.pval.t2),]
result.table2.sorted_bon[1:10,]

#holm
result.table2_holm = data.frame(ID=rownames(exprs_data), t.stat=t2,
                                pvalue=pval.t2, fdr.pvalue=adj.pval.t2_holm)
result.table2.sorted_holm = result.table2_holm[order(adj.pval.t2),]
result.table2.sorted_holm[1:10,]

result_holm <- result.table2_holm[result.table2_holm$fdr.pvalue < 0.5, ]

# Sort the filtered data based on adjusted p-values
result_holm <- result_holm[order(result_holm$fdr.pvalue), ]



exprs_data <- as.data.frame(exprs_data)
##trial on split and mean
merged_data <- exprs_data
merged_data$ID <- NA
merged_data <- merge(exprs_data, clinical, by.x=rownames(exprs_data), by.y = "Array.Data.File")



exprs_transposed <- t(exprs_data)
transpose_sorted <- exprs_transposed[order(rownames(exprs_transposed)), ]
clinical <- clinical[order(clinical$Array.Data.File),]
transpose_sorted$RowNames <- row.names(transpose_sorted)
merged_data <- merge(transpose_sorted, clinical, by.x= rownames(transpose_sorted), by.y = "Array.Data.File")

exprs_data <- as.data.frame(exprs_data)

t_clinical <- clinical
row.names(t_clinical) <- t_clinical$Array.Data.File
t_clinical$Array.Data.File <- NULL
t_clinical <- as.data.frame(t(t_clinical))

exp_row <- t_clinical[1,]
sample_exprs <- exprs_data
sample_exprs <- rbind(sample_exprs, exp_row)

sample_exprs <- as.data.frame(t(sample_exprs))
sample_group0 <- sample_exprs[sample_exprs$group ==0, ]
sample_group1 <- sample_exprs[sample_exprs$group ==1, ]
sample_group0$group <- NULL
sample_group1$group <- NULL
sample_group0 <- as.data.frame(t(sample_group0))
sample_group1 <- as.data.frame(t(sample_group1))




meanrow_group0 <- rowMeans(sample_group0)
meanrow_group1 <- rowMeans(sample_group1)

logfoldchange <- log2(meanrow_group0/meanrow_group1)
logfoldchange <- as.data.frame(logfoldchange)

logfoldchange <- logfoldchange[order(logfoldchange$logfoldchange), ]

print(as.data.frame(logfoldchange))


#####WILCOXON
# Perform Wilcoxon rank sum test for each gene
wilcox_result <- apply(exprs_data, 1, function(x) wilcox.test(x ~ group))
# Extract p-values from the test results
p_values <- sapply(wilcox_result, function(x) x$p.value)
# Adjust p-values for multiple testing using Benjamini-Hochberg correction
adjusted_p_values <- p.adjust(p_values, method = "BH")
# Identify significantly differentially expressed genes based on adjusted p-values
significant_genes <- rownames(exprs_data)[adjusted_p_values < 0.59]
# Create a data frame with gene names and adjusted p-values
significant_genes_df <- data.frame(
  Gene = rownames(exprs_data),
  Adjusted_P_Value = adjusted_p_values)
# Filter the data frame to include only significant genes (adjusted p-value < 0.05)
significant_genes_df <- significant_genes_df[significant_genes_df$Adjusted_P_Value < 0.6, ]
significant_genes_df <- significant_genes_df[order(significant_genes_df$Adjusted_P_Value), ]
top10_wilcoxon <- significant_genes_df[1:10, ] ##top 10 genes
##Mapping probe ID to gene
mappedProbeToGene_2 <- select(hgu133a.db, keys=top10_wilcoxon$Gene, columns = c("PROBEID", "ENSEMBL", "SYMBOL"))
##Final tabular output
final_wilcoxon <- dplyr::left_join(top10_wilcoxon, mappedProbeToGene_2, by= c("Gene" = "PROBEID"))
print(final_wilcoxon)



#####KRUSKAL-WALLIS
# Perform Kruskal-Wallis test for each gene
kw_results <- apply(exprs_data, 1, function(x) kruskal.test(x ~ group))

# Extract p-values from the test results
p_values <- sapply(kw_results, function(res) res$p.value)

# Adjust p-values if necessary (e.g., using the Benjamini-Hochberg method)
# p_values_adj <- p.adjust(p_values, method = "BH")
p_values <- as.data.frame(p_values)
p_values_adj <- p.adjust(p_values$p_values, method = "BH")
p_values_adj
p_values <- p_values[, order(p_values$p_values)]
p_values$probes <- rownames(p_values)

top10_kruskal <- p_values[1:10, ]
mappedProbeToGene_kruskal <- mapIds(hgu133a.db, keys = top10_kruskal$probes, column = "SYMBOL", keytype = "PROBEID")
mappedProbeToGene_kruskal

###RMA2 normalization?
library(oligo)
data_rma2 <- oligo::rma(data, target = "core")

###GCRMA ---- Does not suit this data
library(gcrma)
data_gcrma <- gcrma(data)

exprs_data_gcrma <- affy::exprs(data_gcrma)
head(exprs_data_gcrma)  # Display the first few rows of the expression matrix

t3 = vector()
pval.t3 = vector()

#iteration
for(j in 1:nrow(exprs_data_gcrma))
{
  temp1 <- exprs_data_gcrma[j,]
  res1 <- t.test(temp1[group==1], temp1[group==0], var.equal=T)
  t3[j] <- res1$stat
  pval.t3[j] <- res1$p.val
}
str(res1)
adj.pval.t3 <- p.adjust(pval.t3, "BH")
result.table3 <- data.frame(ID=rownames(exprs_data_gcrma), t.stat=t3,
                            pvalue=pval.t3, fdr.pvalue=adj.pval.t3)
head(result.table3)
result.table3.sorted <- result.table3[order(adj.pval.t3),]
result.table3.sorted[1:10,] # listing the top 10 genes

top10_result_gcrma <- result.table3.sorted[1:10,]
probeIDs <- top10_result$ID

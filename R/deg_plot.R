# Load required libraries
library(dplyr)
library(Seurat)
library(ggplot2)
library(data.table)
library(ggrepel)

load("data/pseudo_seurat_111_20240121.RData")

# Define a list of cell types to be used
ctlist <- c("B-cell", "DC", "Epithelial-cells", "Fibroblasts", "Keratinocytes", 
            "Macrophage", "Monocyte", "Neutrophils", "NK-cell", "T-cells")

# aggregate celltypes of interest and conditions of interest
test <- c(paste(ctlist, "HS", sep = "_"), paste(ctlist, "PSO", sep = "_"))

# Generate a violin plot for IL1B expression across the specified cell types
VlnPlot(pseudo_seurat, features = c("IL1B"), idents = test, group.by = "celltype.stim", pt.size = 1)

# List of genes to analyze
genes <- c("CAMP", "CCL20", "CD83", "CD86", "CLCN2", "CLCN3", "CLCN5", "CSF3", 
           "CXCL1", "CXCL10", "CXCL2", "CXCL8", "CXCR2", "DEFB103B", "DEFB4A", 
           "HLA-DRA", "IFNA1", "IFNGR1", "IL1B", "IL23A", "IL36A", "IL36B", 
           "IL36G", "IL36RN", "IL17C", "IL6", "S100A7", "TNF", "TNFRSF1A")

# Extract the relevant data from the Seurat object and convert it to a data.table
data <- as.data.table(FetchData(pseudo_seurat, vars = c(genes, "celltype.stim", "sample", "disease")))
data[, ident := Idents(pseudo_seurat)]  # Add the 'ident' column to the data table

# Filter the data to include only the celltypes and conditions of interest
filtered_data <- data[ident %in% test]

# Reshape data from wide format to long format
long_data <- melt(filtered_data, 
                  id.vars = c("celltype.stim", "sample", "disease", "ident"),  # Identify these columns
                  variable.name = "gene",  # New column for gene names
                  value.name = "expression")  # New column for expression values

# Identify outliers based on the IQR rule for each gene and cell type
long_data[, is_outlier := expression > (quantile(expression, 0.75) + 1.5 * IQR(expression)) | 
            expression < (quantile(expression, 0.25) - 1.5 * IQR(expression)), 
          by = .(gene, celltype.stim)]

# Define the gene of interest
gene_of_interest <- "CCL20"  # Replace with any gene of interest

# Filter the data for the selected gene
filtered_data_gene <- long_data[gene == gene_of_interest,]

# Add columns for zero-expression checks within each cell type
filtered_data_gene[, all_zero := all(expression == 0), by = celltype.stim]
filtered_data_gene[, majority_above_zero := mean(expression > 0) > 0.5, by = celltype.stim]

# Extract the cell type information by removing condition labels (HS/PSO/PPP)
filtered_data_gene <- filtered_data_gene %>%
  mutate(celltype = sub("(_HS|_PSO|_PPP)$", "", celltype.stim))

# Determine y-axis limits based on the expression range for the gene
y_min <- min(filtered_data_gene$expression, na.rm = TRUE) - 0.001
y_max <- max(filtered_data_gene$expression, na.rm = TRUE) + 0.001
print(c(y_min, y_max))  # Print the calculated y-axis limits for reference

# Define color scheme for the disease conditions
current_colors <- c("PSO" = "#F8766D", "HS" = "#00BFC4", "PPP" = "#875BA6")

# Calculate the overall standard deviation for expression
overall_sd <- sd(filtered_data_gene$expression, na.rm = TRUE)

# Function to assign bandwidth based on overall standard deviation
assign_bw <- function(sd_val) {
  if (sd_val < 1e-6) {
    return(0.0001)  # Small bandwidth for low variability
  } else if (sd_val < 0.01) {
    return(0.001)  # Medium bandwidth for moderate variability
  } else {
    return(0.1)  # Larger bandwidth for high variability
  }
}

# Calculate the bandwidth value based on the overall standard deviation
bw_value <- assign_bw(overall_sd)

####################
# Perform t-test for each cell type comparing expression by disease
df <- filtered_data_gene[gene == gene_of_interest,]
t_test_results <- df %>%
  group_by(celltype) %>%
  filter(n_distinct(disease) == 2) %>%  # Ensure comparison is between two disease conditions
  summarise(
    t_test_p_value = t.test(expression ~ disease, var.equal = TRUE)$p.value  # Perform t-test
  )

# View the t-test results
t_test_results

####################

# Filter significant p-values and prepare annotations
pvalue_annotations <- t_test_results %>%
  mutate(label = ifelse(is.nan(t_test_p_value), "p = NaN", paste0("p = ", round(t_test_p_value, 3)))) %>%
  filter(t_test_p_value < 0.05) %>%  # Only keep p-values below 0.05
  left_join(df, by = "celltype") %>%
  group_by(celltype) %>%
  summarise(
    label = dplyr::first(label),  # Get the first label for each cell type
    celltype.stim = dplyr::first(celltype.stim),  # Position the label
    disease = dplyr::first(disease),  # Position based on disease
    y_position = ifelse(all(is.na(expression)), NA, max(expression, na.rm = TRUE) + 0.1)  # Handle all NA values
  ) %>%
  ungroup()


# Plot the gene expression for the selected gene
ggplot(filtered_data_gene, aes(x = celltype.stim, y = expression, fill = disease)) +
  # Violin plot for non-zero groups
  geom_violin(trim = FALSE, bw = bw_value, draw_quantiles = TRUE, data = filtered_data_gene[all_zero == FALSE]) +
  # Box plot for groups where all values are zero
  geom_boxplot(width = 0.2, alpha = 0.7, data = filtered_data_gene[all_zero == TRUE]) +
  # Scatter plot for individual data points
  geom_jitter(size = 1, alpha = 0.6, aes(color = "black"), position = position_jitter(width = 0.2, height = 0)) +
  # Label outliers using ggrepel for non-overlapping text
  geom_text_repel(aes(label = ifelse(is_outlier, as.character(sample), "")),
                  size = 3, box.padding = 0.35, point.padding = 0.5,
                  segment.color = 'grey50', max.overlaps = Inf) + 
  # Set y-axis limits
  scale_y_continuous(limits = c(0, y_max)) +  
  # Set manual fill and color scales for disease groups
  scale_fill_manual(values = current_colors) +  
  scale_color_manual(values = current_colors) +  
  # Titles and axis labels
  labs(title = paste(gene_of_interest, " (PSO vs HS)"),
       x = "",
       y = "Expression") +
  # Add p-value annotations
  geom_text(data = pvalue_annotations, 
            aes(x = celltype.stim, y = y_max - 0.001, label = label), 
            vjust = -0.5, 
            angle = 45) +
  # Theme adjustments
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
        plot.title = element_text(hjust = 0.5))


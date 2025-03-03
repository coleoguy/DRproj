### Megan Copeland
### July 9th, 2024
### Plot gene and repeat counts

######## Load and filter data ############################

# Load required libraries
library(ggplot2)      
library(rtracklayer)

# Create data frame for chromosome lengths
chrom_lengths <- data.frame(
  Chromosome = c("NC_035107.1", "NC_035108.1", "NC_035109.1"),
  Length = c(310827022, 474425716, 409777670))

# Read in GFF data
gffData <- as.data.frame(import("data/genomic.gff"))

# Subset to keep protein coding exons for main chromosomes
exonData <- subset(gffData, 
                   seqnames %in% chrom_lengths$Chromosome & 
                   type == "exon" & gbkey == "mRNA")

# Remove dup exons belonging to multiple transcripts
exonData <- exonData[!duplicated(exonData[, c("seqnames", "start", "end")]), ]

# Calculate midpoint for each exon
exonData$Midpoint <- (exonData$start + exonData$end) / 2

# Read in direct repeat data
drs <- read.csv("results/condensed_results.csv")

# Filter out data for unplaced scafs
drs <- subset(drs, Chromosome %in% chrom_lengths$Chromosome)

# Get midpoint between repeats using midpoint of the repeat and its copy
drs$Midpoint <- ((drs$Start_Position + drs$End_Position) / 2 + (drs$Match_Position + drs$Match_End_Position) / 2) / 2

####### Get repeat and gene counts in windows ###############

# Define sliding window size and step
window_size <- 200000  
step_size <- 200000    

# Initialize list to store results
results <- list()

# Get chromosome names from repeat data
chr_names <- unique(drs$Chromosome)

# Get counts in 200kb windows
for (i in 1:length(chr_names)) {
  # Get length of current chromosome from chrom_lengths
  chr_length <- chrom_lengths$Length[chrom_lengths$Chromosome == chr_names[i]]
  # Gets midpoints of genes and repeats for the current chromosome
  gene_positions <- exonData$Midpoint[exonData$seqnames == chr_names[i]]
  repeat_positions <- drs$Midpoint[drs$Chromosome == chr_names[i]]
  # Create sliding windows
  windows <- seq(1, chr_length, by = step_size)
  # Count genes in each window
  gene_counts <- sapply(windows, function(w) {
    sum(gene_positions >= w & gene_positions < (w + window_size))
  })
  # Count repeats in each window
  repeat_counts <- sapply(windows, function(w) {
    sum(repeat_positions >= w & repeat_positions < (w + window_size))
  })
  # Store results for each chromosome
  results[[chr_names[i]]] <- data.frame(
    Window_Start = windows,
    Gene_Count = gene_counts,
    Repeat_Count = repeat_counts,
    Chromosome = chr_names[i]
  )
}

# Combine results for each chromosome into a single df
results_df <- do.call(rbind, results)

# Get log of repeat counts
results_df$Log_Repeat_Count <- log10(results_df$Repeat_Count + 1)
max(results_df$Repeat_Count)

######## Plot results & fit to linear model #################

ggplot(results_df, aes(x = Gene_Count, y = Log_Repeat_Count)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(x = "Gene count",
       y = "Repeat count") +
  scale_y_continuous(
    breaks = log10(c(1, 10, 100, 1000)), # set the breaks for log scale
    labels = c(1, 10, 100, 1000) # manually set the labels to reflect unlogged values
  ) +
  coord_cartesian(ylim = c(log10(1), log10(1033))) + # Set the y-axis limits
  theme_minimal() +
  theme(
    axis.title.x = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Fit the linear model
lm_model <- lm(results_df$Repeat_Count ~ Gene_Count, data = results_df)

# Display the summary statistics for the model
summary(lm_model)


### Name: Megan Copeland
### Date: Sep 14th, 2024
### Purpose: Create null dist of exons flanked by drs

################ Loading libraries ###########################################
# load the required packages
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(viridis)

############### Loading function for monte carlo sim ########################

run_simulation <- function(chrom_filter = c("NC_035109.1"), num_simulations = 100) {
 
  # read in data files: gff data, repeat data, and chromosome lengths
  gffData <- import("data/genomic.gff")
  drs <- read.csv("results/condensed_results.csv")
  chrom_lengths <- read.csv("data/chrom_lengths.csv")
  
  # Filter for selected chromosomes in gff and repeat data
  exonData <- as.data.table(subset(gffData, 
                                   seqnames %in% chrom_filter & 
                                     type == "exon" & gbkey == "mRNA"))
  
  # Remove duplicate exons
  duplicates <- duplicated(exonData[, .(seqnames, start, end)])
  exonData <- exonData[!duplicates]
  
  # Remove unplaced scaffolds from repeat data
  drs <- subset(drs, Chromosome %in% chrom_filter)
  
  # Create an empty vector to store the results
  null_dist <- numeric(num_simulations)
  
  # Simulate num_simulations times
  for (i in 1:num_simulations) {
    print(paste("Running simulation:", i))
    
    # Create an empty data frame to store random repeat locations
    random_drs <- data.frame(
      Chromosome = drs$Chromosome,
      New_Start_Position = numeric(nrow(drs)),
      New_End_Position = numeric(nrow(drs)),
      New_Match_Position = numeric(nrow(drs)),
      New_Match_End_Position = numeric(nrow(drs))
    )
    
    # Loop through each row in drs to assign new start positions
    for (j in 1:nrow(drs)) {
      row <- drs[j,]
      
      # Is this the first row or is the previous position more than 10bp away?
      if (j == 1 || (row$Start_Position - drs$End_Position[j - 1]) > 10) {
        # If so, find the corresponding chrom length and sample new start position
        chrom_length <- chrom_lengths$End[chrom_lengths$Chromosome == row$Chromosome]
        new_start <- sample(1:(chrom_length - (row$End_Position - row$Start_Position)), 1)
      } else {
        # Set new_start relative to the previous repeat
        prev_new_end <- random_drs$New_End_Position[j - 1]
        new_start <- prev_new_end + (row$Start_Position - drs$End_Position[j - 1])
      }
      
      # Calculate the spacing according to the original data
      spacing <- row$End_Position - row$Start_Position
      match_spacing <- row$Match_End_Position - row$Match_Position
      
      # Apply spacing to get new positions
      random_drs$New_Start_Position[j] <- new_start
      random_drs$New_End_Position[j] <- new_start + spacing
      random_drs$New_Match_Position[j] <- new_start + (row$Match_Position - row$Start_Position)
      random_drs$New_Match_End_Position[j] <- random_drs$New_Match_Position[j] + match_spacing
    }
    
    # Convert random_drs to a data table
    random_drs_dt <- as.data.table(random_drs)
    
    # Perform a non-equi join with exonData
    flanked_data <- random_drs_dt[exonData, 
                                  .(ExonID = i.ID, 
                                    Chromosome = x.Chromosome, 
                                    New_Start_Position = x.New_Start_Position, 
                                    New_Match_End_Position = x.New_Match_End_Position), 
                                  on = .(Chromosome = seqnames, 
                                         New_Start_Position < end, 
                                         New_Match_End_Position > start),
                                  nomatch = 0L]
    
    # Combine IDs of flanked exons into a single string
    combined_dat <- flanked_data[, .(FlankedExonIDs = paste(ExonID, collapse = "; ")), 
                                 by = .(Chromosome, New_Start_Position, New_Match_End_Position)]
    
    # Merge combined data with original repeat data and reorder columns
    drs_final_dt <- merge(random_drs_dt, combined_dat, 
                          by = c("Chromosome", "New_Start_Position", "New_Match_End_Position"), 
                          all.x = TRUE)
    columns <- c("Chromosome", "New_Start_Position", "New_End_Position", 
                 "New_Match_Position", "New_Match_End_Position", "FlankedExonIDs")
    drs_final_selected <- drs_final_dt[, ..columns]
    
    # Store the count of unique flanked exon IDs
    null_dist[i] <- uniqueN(unlist(strsplit(na.omit(drs_final_selected$FlankedExonIDs), "; ")))
  }
  
  return(null_dist)
}


############### Run sim for whole genome & individual chroms ########################

# Example of how users can run the simulation for a specific chromosome filter
whole_genome_results <- run_simulation(c("NC_035107.1","NC_035108.1","NC_035109.1"))
chr1_results <- run_simulation(c("NC_035107.1"))
chr2_results <- run_simulation(c("NC_035108.1"))
chr3_results <- run_simulation(c("NC_035109.1"))

############# Plot obs vs. exp for entire genome ###################################

colors <- plasma(100)

# Create the genome-wide  plot
plot(density(whole_genome_results), 
     main = "",
     xlab = "Unique Flanked Exons", 
     ylab = "Density", 
     col = colors[45], 
     lwd = 2, 
     xlim = c(1000, max(c(whole_genome_results) + 500)))

# Add polygon for shading
polygon(density(whole_genome_results), col = adjustcolor(colors[50], alpha.f = 0.75))

# Add vertical line for observed value 
# Observed value are obtained from the flanked_exons script
abline(v = 5850, col = adjustcolor(colors[50]), lwd = 2)

# Save plot as PDF (optional)
dev.copy(pdf, file = "figures/wg_nulldist.pdf", width = 4, height = 6)
dev.off()

############ Plot obs vs. exp for each chromosome ##############################

# Get colors
colors <- plasma(12)

# Create the plot for first chromosome
plot(density(chr1_results), 
     main = "",
     xlab = "Unique Flanked Exons", 
     ylab = "Density", 
     col = colors[1], 
     lwd = 2, 
     xlim = c(0, max(c(chr1_results, chr2_results, chr3_results)) + 5000))

# Add density lines for other chromosomes
lines(density(chr2_results), col = colors[6], lwd = 2)
lines(density(chr3_results), col = colors[11], lwd = 2)

# Fill in null distributions
polygon(density(chr1_results), col = adjustcolor(colors[1], alpha.f = 0.3))
polygon(density(chr2_results), col = adjustcolor(colors[6], alpha.f = 0.5))
polygon(density(chr3_results), col = adjustcolor(colors[11], alpha.f = 0.6))

# Add vertical lines for observed values 
# Observed values are obtained from the flanked_exons script
abline(v = 1274, col = adjustcolor(colors[1]), lwd = 2)
abline(v = 2443, col = adjustcolor(colors[6]), lwd = 2.5)
abline(v = 2122, col = adjustcolor(colors[11]), lwd = 2.5)

# Add legend
legend(x = 12000, y = 0.0042, # Specify coordinates
       legend = c("Chr1", "Chr2", "Chr3"),
       col = colors[c(1, 6, 11)],
       lwd = 2,
       bg = adjustcolor("white"),
       cex = .6, 
       bty = "n")

# Save plot as PDF (optional)
dev.copy(pdf, file = "figures/chromosome_nulldist.pdf", width = 4, height = 6)
dev.off()


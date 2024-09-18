### Name: Megan Copeland
### Date: Sep 14, 2024
### Purpose: ID exons flanked by DRs
### Contact: mcc146@tamu.edu

########################### Load and filter data #############################

library(data.table)     
library(rtracklayer)     
library(GenomicRanges)

# Import GFF file and direct repeat (DR) data
gffData <- import("data/genomic.gff")
drs <- fread("results/condensed_results.csv")

# Filter out unplaced scaffolds & exons for non-protein coding genes
exonData <- as.data.table(subset(gffData,
                                 seqnames %in% c("NC_035107.1", "NC_035108.1", "NC_035109.1") &
                                   type == "exon" & gbkey == "mRNA"))

# Remove duplicate exons entries for exons belonging to multiple transcripts
exonData <- exonData[!duplicated(exonData[, .(seqnames, start, end)])]

# Filter out unplaced scaffold data from drs data
drs <- subset(drs, Chromosome %in% c("NC_035107.1", "NC_035108.1", "NC_035109.1"))

########################### Identify flanked exons ########################

# Perform a non-equi join to find exons flanked by DRs
flanked_data <- drs[exonData,
                    .(ExonID = i.ID,                  # From exonData
                      Chromosome = x.Chromosome,      # From drs data
                      Start_Position = x.Start_Position,
                      Match_End_Position = x.Match_End_Position),
                    on = .(Chromosome = seqnames,     # Join conditions
                           Start_Position < end, 
                           Match_End_Position > start),
                    nomatch = 0L]                     # Exclude rows without matches

# Combine IDs of flanked exons into a single string per repeat
combined_dat <- flanked_data[, .(FlankedExonIDs = paste(ExonID, collapse = "; ")),
                             by = .(Chromosome, Start_Position, Match_End_Position)]

# Merge combined data with the original repeat data and reorder columns
drs_final <- merge(drs, combined_dat, 
                   by = c("Chromosome", "Start_Position", "Match_End_Position"), 
                   all.x = TRUE)

# Reorder columns for clarity
column_order <- c("Chromosome", "Start_Position", "End_Position", 
                  "Match_Position", "Match_End_Position", "FlankedExonIDs")
drs_final <- drs_final[, ..column_order]

########################### Get flanked exon counts ###########################

# Remove rows where FlankedExonIDs is NA
condensed_drs_final <- na.omit(drs_final, cols = "FlankedExonIDs")

# Calculate the total number of unique flanked exons
num_unique_exonIDs <- uniqueN(unlist(strsplit(na.omit(drs_final$FlankedExonIDs), "; ")))
cat("Total number of unique flanked exons:", num_unique_exonIDs, "\n")

# Calculate the number of flanked exons per chromosome
exons_flanked_per_chromosome <- flanked_data[, .(Num_Exons_Flanked = uniqueN(ExonID)), by = .(Chromosome)]
cat("Number of flanked exons per chromosome:\n")
print(exons_flanked_per_chromosome)

########################### Export results ################################

# Write final results to CSV
write.csv(condensed_drs_final, "results/flanked_exons_data.csv", row.names = FALSE)

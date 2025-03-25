### Name: Megan Copeland
### Date: January 22, 2024
### Purpose: ID exons flanked by drs
### Contact: mcc146@tamu.edu

# Load required packages
library(data.table)
library(rtracklayer)
library(GenomicRanges)

# Import GFF data
gffData <- import("../data/genomic.gff")

# Import direct repeats data
drs <- fread("../results/aedes_merged.csv")

# Extract everything before the first space
drs$Chromosome <- sapply(strsplit(drs$Chromosome, " "), `[`, 1)

# Filter exon entries for specific chromosomes and feature types
exonData_gr <- subset(gffData,
                      seqnames %in% c("NC_035107.1", "NC_035108.1", "NC_035109.1") &
                        type == "exon" &
                        (gbkey == "mRNA"))

# Remove redundant exons
exonData_gr_unique <- unique(exonData_gr)

# Convert GRanges object to data.table
exonData_df <- as.data.frame(exonData_gr_unique)
exonData <- setDT(exonData_df)

# Find exons flanked by DRs using non-equi join
flanked_data <- drs[exonData,
                    .(ExonID = i.ID,
                      Chromosome = x.Chromosome,
                      Start_Position = x.Start_Position,
                      Match_End_Position = x.Match_End_Position),
                    on = .(Chromosome = seqnames,
                           Start_Position < end,
                           Match_End_Position > start),
                    nomatch = 0L]

# Combine flanked exon IDs for each DR
combined_dat <- flanked_data[, .(FlankedExonIDs = paste(ExonID, collapse = "; ")),
                             by = .(Chromosome, Start_Position, Match_End_Position)]

# Merge combined data with original DR data and reorder columns
drs_final <- merge(drs, combined_dat, by = c("Chromosome", "Start_Position", "Match_End_Position"), all.x = TRUE)
column_order <- c("Chromosome", "Start_Position", "End_Position", "Match_Position", "Match_End_Position", "FlankedExonIDs")
drs_final <- drs_final[, ..column_order]

# Remove rows where FlankedExonIDs is NA and save to CSV
condensed_drs_final <- na.omit(drs_final, cols = "FlankedExonIDs")
write.csv(condensed_drs_final, "../results/flanked_exons_data.csv", row.names = F)

# Get the total number of unique exons flanked by DRs
num_unique_exonIDs <- uniqueN(unlist(strsplit(na.omit(drs_final$FlankedExonIDs), "; ")))
print(num_unique_exonIDs)

# Get the number of flanked exons per chromosome
exons_flanked_per_chromosome <- flanked_data[, .(Num_Exons_Flanked = uniqueN(ExonID)), by = .(Chromosome)]
print(exons_flanked_per_chromosome)

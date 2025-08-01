### Name: Megan Copeland
### Date: January 22, 2024
### Purpose: ID exons flanked by drs
### Contact: mcc146@tamu.edu

#############################
##  1. LOAD PACKAGES      ##
#############################
library(data.table)
library(rtracklayer)
library(GenomicRanges)

#############################
##  2. IMPORT DATA        ##
#############################
# GFF annotation file
gffData <- import("../data/genomic.gff")
# Direct repeats results
drs     <- fread("../results/aedes_merged.csv")

#############################
##  3. CLEAN & PREPARE    ##
#############################
# Keep only the accession (before any space)
drs$Chromosome <- sapply(strsplit(drs$Chromosome, " "), `[`, 1)

#############################
##  4. FILTER EXONS       ##
#############################
# Specify chromosomes of interest
keep_ids <- c("NC_035107.1", "NC_035108.1", "NC_035109.1")
# Subset to exons on those chromosomes with gbkey "mRNA"
exonData_gr <- subset(
  gffData,
  seqnames %in% keep_ids &
    type == "exon" &
    gbkey == "mRNA"
)
# Remove any duplicate exon entries
exonData_gr_unique <- unique(exonData_gr)

#############################
##  5. CONVERT TO DATA.TABLE ##
#############################
# Turn GRanges into a data.table for joining
exonData_df <- as.data.frame(exonData_gr_unique)
exonData    <- setDT(exonData_df)

#############################
##  6. IDENTIFY FLANKED EXONS ##
#############################
# Non-equi join: DRs that overlap exon regions
flanked_data <- drs[exonData,
  .(
    ExonID             = i.ID,
    Chromosome         = x.Chromosome,
    Start_Position     = x.Start_Position,
    Match_End_Position = x.Match_End_Position
  ),
  on = .(
    Chromosome         = seqnames,
    Start_Position < end,
    Match_End_Position > start
  ),
  nomatch = 0L
]

#############################
##  7. COMBINE & MERGE    ##
#############################
# Aggregate ExonIDs per DR
combined_dat <- flanked_data[, .(FlankedExonIDs = paste(ExonID, collapse = "; ")), 
                             by = .(Chromosome, Start_Position, Match_End_Position)]

# Merge back onto original DR table, preserving all DRs
drs_final <- merge(drs, combined_dat,
  by = c("Chromosome", "Start_Position", "Match_End_Position"),
  all.x = TRUE
)

# Reorder columns for clarity
column_order <- c("Chromosome", "Start_Position", "End_Position", "Match_Position", "Match_End_Position", "FlankedExonIDs")
drs_final <- drs_final[, ..column_order]

#############################
##  8. WRITE OUTPUT       ##
#############################
# Remove DRs with no flanked exons and save results
condensed_drs_final <- na.omit(drs_final, cols = "FlankedExonIDs")
write.csv(
  condensed_drs_final,
  "../results/flanked_exons_data.csv",
  row.names = FALSE
)

#############################
##  9. SUMMARY STATISTICS ##
#############################
# Total unique exons flanked by DRs
num_unique_exonIDs <- uniqueN(unlist(strsplit(na.omit(drs_final$FlankedExonIDs), "; ")))
print(num_unique_exonIDs)

# Number of flanked exons per chromosome
exons_flanked_per_chromosome <- flanked_data[, .(Num_Exons_Flanked = uniqueN(ExonID)), by = .(Chromosome)]
print(exons_flanked_per_chromosome)

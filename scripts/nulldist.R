#############################
##  Observed vs. Expected ##
#############################

# Name: Megan Copeland
# Date: January 24th, 2024
# Purpose: Create null distribution of exons flanked by direct repeats (DRs) with faster logic,
#          returning both whole-genome and per-chromosome counts.

# Load required packages
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(parallel)    # for parallelization
library(viridis)     # for color palettes (plots)

######################
##  1. DATA IMPORT  ##
######################

# Import GFF data
gffData <- import("../data/genomic.gff")

# Import direct repeats data
drs <- fread("../results/aedes_merged.csv")

# Extract everything before the first space
drs$Chromosome <- sapply(strsplit(drs$Chromosome, " "), `[`, 1)
drs <- drs[Chromosome %in% c("NC_035107.1", "NC_035108.1", "NC_035109.1")]

# Import chromosome lengths
chrom_lengths <- fread("../data/chrom_lengths.csv")  # Adjust read.csv/fread as needed

########################
##  2. DATA FILTERING ##
########################

# Filter exons for specific chromosomes and relevant feature types
exonData_gr <- subset(
  gffData,
  seqnames %in% c("NC_035107.1", "NC_035108.1", "NC_035109.1") &
    type == "exon" &
    (gbkey == "mRNA" | gbkey == "misc_RNA" | gbkey == "lncRNA")
)

# Remove redundant exons
exonData_gr_unique <- unique(exonData_gr)

###################################
##  3. PRECOMPUTE HELPER COLUMNS ##
###################################

drs <- drs[, .(
  Chromosome,
  Start_Position,
  End_Position,
  Match_Position,
  Match_End_Position,
  spacing       = End_Position - Start_Position,
  match_spacing = Match_End_Position - Match_Position,
  offset        = Match_Position - Start_Position
)]

##########################################
##  4. SPLIT REPEATS BY CHROMOSOME     ##
##     & RANDOMIZATION HELPER FUNCTION ##
##########################################

drs_split <- split(drs, by = "Chromosome")

# Create a lookup for chromosome lengths
chrom_len_map <- setNames(chrom_lengths$End, chrom_lengths$Chromosome)

randomize_chromosome <- function(drs_chr, chrom_length) {
  n <- nrow(drs_chr)
  new_starts <- numeric(n)
  
  # 1) First repeat: random start
  new_starts[1] <- sample(
    seq_len(chrom_length - drs_chr$spacing[1]),
    size = 1
  )
  
  # 2) Loop over repeats in this chromosome
  for (k in 2:n) {
    dist_from_prev <- drs_chr$Start_Position[k] - drs_chr$End_Position[k - 1]
    
    if (dist_from_prev > 10) {
      new_starts[k] <- sample(
        seq_len(chrom_length - drs_chr$spacing[k]),
        size = 1
      )
    } else {
      prev_new_end <- new_starts[k - 1] + drs_chr$spacing[k - 1]
      new_starts[k] <- prev_new_end + dist_from_prev
    }
  }
  
  drs_chr[, New_Start_Position := new_starts]
  drs_chr[, New_End_Position   := new_starts + spacing]
  drs_chr[, New_Match_Position := new_starts + offset]
  drs_chr[, New_Match_End_Position := New_Match_Position + match_spacing]
  
  return(drs_chr)
}

############################################
##  5. SINGLE SIMULATION FUNCTION (PAR)   ##
############################################

do_one_simulation <- function(iter, drs_split, chrom_len_map, exonData_gr_unique) {
  # 1) Randomize each chromosome separately
  random_drs_list <- lapply(names(drs_split), function(chrom) {
    drs_chr <- drs_split[[chrom]]
    chrom_length <- chrom_len_map[chrom]
    randomize_chromosome(drs_chr, chrom_length)
  })
  
  # 2) Combine into a single data.table for the entire genome
  random_drs_dt <- rbindlist(random_drs_list, use.names = TRUE, fill = TRUE)
  
  # 3) Build a single GRanges object from the combined data
  random_drs_gr <- GRanges(
    seqnames = random_drs_dt$Chromosome,
    ranges   = IRanges(
      start  = random_drs_dt$New_Start_Position,
      end    = random_drs_dt$New_Match_End_Position
    )
  )
  
  # 4) Overlap the entire random set with all exons
  hits <- findOverlaps(random_drs_gr, exonData_gr_unique)
  
  # 5) Identify the *unique exon indices* that were overlapped
  flanked_exons <- unique(subjectHits(hits))
  
  # 6) Count how many unique exons per chromosome and overall
  #    (We need the chromosome for each overlapped exon)
  if (length(flanked_exons) == 0) {
    # If no exons were overlapped, return 0s
    #   for each chromosome + total
    out_list <- as.list(setNames(rep(0, length(drs_split)), names(drs_split)))
    out_list$total <- 0
  } else {
    # Find the chromosome for each overlapped exon
    overlapped_exon_chroms <- as.character(seqnames(exonData_gr_unique))[flanked_exons]
    
    # Create a data.table of unique exons + their chromosome
    dt_exons <- data.table(
      exon_index = flanked_exons,
      chrom      = overlapped_exon_chroms
    )
    dt_exons <- unique(dt_exons)
    
    # Count how many exons per chromosome
    chr_counts <- dt_exons[, .N, by = chrom]  # .N = number of rows
    
    # Turn this into a named list
    # Initialize 0 for each chromosome of interest
    out_list <- as.list(setNames(rep(0, length(drs_split)), names(drs_split)))
    
    # Fill in counts for any chromosome that has flanked exons
    for (i in seq_len(nrow(chr_counts))) {
      out_list[[chr_counts$chrom[i]]] <- chr_counts$N[i]
    }
    
    # Add total
    out_list$total <- length(unique(dt_exons$exon_index))
  }
  
  # Return as a data.table row
  return(as.data.table(out_list))
}

##############################################
##  6. RUN THE MONTE CARLO USING MULTICORE  ##
##############################################

# Number of iterations
n_iter <- 100

# Detect number of cores, or set as desired
n_cores <- max(1, detectCores() - 1)

# Run in parallel, collecting a data.table row per iteration
null_dist_chr_list <- mclapply(
  X         = seq_len(n_iter),
  FUN       = do_one_simulation,
  drs_split = drs_split,
  chrom_len_map = chrom_len_map,
  exonData_gr_unique = exonData_gr_unique,
  mc.cores  = n_cores
)

# Combine all rows
null_dist_chr_dt <- rbindlist(null_dist_chr_list, fill = TRUE)

# Inspect columns: each chromosome + 'total'
# e.g., 'NC_035107.1', 'NC_035108.1', 'NC_035109.1', 'total'
# head(null_dist_chr_dt)

# Separate out columns for convenience
null_dist_chr1 <- null_dist_chr_dt[["NC_035107.1"]]
null_dist_chr2 <- null_dist_chr_dt[["NC_035108.1"]]
null_dist_chr3 <- null_dist_chr_dt[["NC_035109.1"]]
null_dist_wg   <- null_dist_chr_dt[["total"]]

############################
##  7. PLOTTING RESULTS   ##
############################

# For demonstration, we’ll do two plots:
# 1) Whole-genome distribution
# 2) Per-chromosome distribution (3 chromosomes)

colors <- plasma(12)

## -- 1) WHOLE GENOME PLOT --

plot(
  density(null_dist_wg), 
  xlab = "Number of Flanked Exons (Unique)", 
  ylab = "Density", 
  col  = colors[5], 
  lwd  = 2,
  xlim = c(0, max(null_dist_wg) + 500)
)

polygon(
  density(null_dist_wg),
  col    = adjustcolor(colors[6], alpha.f = 0.5),
  border = NA
)

# Example observed line at 5882 (adjust as needed)
abline(v = 5882, col = adjustcolor(colors[6]), lwd = 2)


## -- 2) PER-CHROMOSOME PLOT --

# We'll use the three vectors:
# null_dist_chr1, null_dist_chr2, null_dist_chr3
plot(density(null_dist_chr1),
  xlab = "Number of Flanked Exons (Unique)", 
  ylab = "Density",
  col  = colors[1],
  lwd  = 2,
  xlim = c(0, max(c(null_dist_chr1, null_dist_chr2, null_dist_chr3)) + 500))

lines(density(null_dist_chr2), col = colors[6],  lwd = 2)
lines(density(null_dist_chr3), col = colors[11], lwd = 2)

# Polygons for shading
polygon(density(null_dist_chr1),
  col    = adjustcolor(colors[1], alpha.f = 0.2),
  border = NA)

polygon(density(null_dist_chr2),
  col    = adjustcolor(colors[6], alpha.f = 0.2),
  border = NA)

polygon(density(null_dist_chr3),
  col    = adjustcolor(colors[11], alpha.f = 0.2),
  border = NA)

# Example observed lines per chromosome (adjust values as needed)
abline(v = 1270, col = colors[1],  lwd = 2)   # Chr1
abline(v = 2453, col = colors[6],  lwd = 2)   # Chr2
abline(v = 2159, col = colors[11], lwd = 2)   # Chr3

legend("topright",
  legend = c("Chr1", "Chr2", "Chr3"),
  col    = colors[c(1, 6, 11)],
  lwd    = 2,
  inset  = 0.02,
  bg     = adjustcolor("white", alpha.f = 0.8),
  cex    = 0.5,
  bty    = "n")

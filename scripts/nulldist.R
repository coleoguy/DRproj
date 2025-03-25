#############################
##  Observed vs. Expected ##
#############################

# Name: Megan Copeland
# Date: March 13th, 2024
# Purpose: Create null distribution of exons flanked by direct repeats (DRs) with faster logic,
#          returning both whole-genome and per-chromosome counts.

# Load required packages
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(parallel)    # for parallelization
library(beeswarm)
library(plotrix)

######################
##  1. DATA IMPORT  ##
######################

# Import GFF data
gffData <- import("../data/genomic.gff")

# Import direct repeats data
drs <- fread("../results/aedes_merged.csv")

# Extract chrom names
chrom_names <- unique(drs$Chromosome)

# Import chromosome lengths
chrom_lengths <- fread("../data/chrom_lengths.csv")  # Adjust read.csv/fread as needed

########################
##  2. DATA FILTERING ##
########################

# Filter exons for specific chromosomes and relevant feature types
exonData_gr <- subset(
  gffData,
  seqnames %in% chrom_names &
    type == "exon" &
    (gbkey == "mRNA")
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
    
    if (dist_from_prev > 26) {
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
##  7) Visualization ##
############################

# Define the observed (empirical) counts from the actual data for each chromosome and genome-wide.
obs_chr1 <- 1247
obs_chr2 <- 2427
obs_chr3 <- 2108
obs_wg   <- 5782

# Set color values for plotting the chromosome-level data and the whole-genome data.
chrom_color <- "steelblue"
wg_color    <- "tomato"

y_chr_range <- range(null_dist_chr1, null_dist_chr2, null_dist_chr3,
                     obs_chr1, obs_chr2, obs_chr3)
y_wg_range  <- range(null_dist_wg, obs_wg)

# Forward transform: WG -> chromosome scale
trans_wg <- function(x) {
  (x - y_wg_range[1]) / (y_wg_range[2] - y_wg_range[1]) *
    (y_chr_range[2] - y_chr_range[1]) + y_chr_range[1]
}

# We'll place Chr1 at x=1, Chr2 at x=4, Chr3 at x=8, and Genome at x=12
x_positions <- c(1, 4, 8, 12)

#############################
##  7a) Setup Plot & Axes   ##
#############################
op <- par(mar = c(7, 5, 4, 7) + 0.1, xpd = TRUE)  # more bottom & right margin

plot(NA, type = "n",
     xlim = c(0, 13),
     ylim = y_chr_range,
     xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")


# Make left axis ticks horizontal (las=1) & bold, in steelblue
axis(2, las=1, col.axis = chrom_color, font.axis = 2, cex.axis = 0.8)
mtext("Number of flanked exons", side = 2, line = 4, col=chrom_color, font=2)

#######################################
##  7b) Bottom axis with a "break"   ##
#######################################
# 1) Draw ticks for Chr1, Chr2, Chr3
axis(1, at = c(1,4,8), labels = c("Chr1","Chr2","Chr3"), las = 1)

# 2) Insert a slash break around x=9.5 (between 8 and 12)
#    purely cosmetic, doesn't "shrink" the scale.
axis.break(axis = 1, breakpos = 9.5, style = "slash", brw = 0.02)

# 3) Manually add a tick/label for "Genome" at x=12
axis(1, at = 12, labels = "Genome", las = 1)

#########################
##  7c) Beeswarm Points ##
#########################
beeswarm(null_dist_chr1, at = x_positions[1], add = TRUE,
         pch = 16, col = adjustcolor(chrom_color, alpha.f = 0.65), method = "hex",
         corral = "wrap", spacing = .6, cex = 0.75)

beeswarm(null_dist_chr2, at = x_positions[2], add = TRUE,
         pch = 16, col = adjustcolor(chrom_color, alpha.f = 0.65), method = "hex",
         corral = "wrap", spacing = .6, cex = 0.75)

beeswarm(null_dist_chr3, at = x_positions[3], add = TRUE,
         pch = 16, col = adjustcolor(chrom_color, alpha.f = 0.65), method = "hex",
         corral = "wrap", spacing = .6, cex = 0.75)

# Transform WG onto chromosome scale
wg_points_trans <- trans_wg(null_dist_wg)
beeswarm(wg_points_trans, at = x_positions[4], add = TRUE,
         pch = 16, col = adjustcolor(wg_color, alpha.f = 0.65), method = "hex",
         corral = "wrap", spacing = .6, cex = 0.75)

##############################
##  7d) Observed (Empirical) ##
##############################
points(x_positions[1], obs_chr1,
       pch = 21, bg = chrom_color, col = "black",
       cex = 1, lwd = 2)
points(x_positions[2], obs_chr2,
       pch = 21, bg = chrom_color, col = "black",
       cex = 1, lwd = 2)
points(x_positions[3], obs_chr3,
       pch = 21, bg = chrom_color, col = "black",
       cex = 1, lwd = 2)
points(x_positions[4], trans_wg(obs_wg),
       pch = 21, bg = wg_color, col = "black",
       cex = 1, lwd = 2)

##############################
##  7e) Right-Side Axis (WG) ##
##############################
old_xpd <- par("xpd")
par(xpd = FALSE)  # clip to plot region

wg_ticks <- pretty(y_wg_range, n = 4)
axis_tick_positions <- trans_wg(wg_ticks)

# Make tick labels horizontal (las=1), bold, tomato color
axis(4, at = axis_tick_positions, labels = wg_ticks, line = 0, 
     las=1, cex.axis = 0.8, col.axis = wg_color, font.axis = 2)
mtext("Number of flanked exons", side = 4, line = 4, col=wg_color, font=2)

par(xpd = old_xpd)

############################
##  7f) Legend & Cleanup   ##
############################
legend("topleft",
       legend = c("Null", "Empirical"),
       pch    = c(16, 21),
       pt.bg  = c("grey", "grey"),
       col    = c("grey", "black"),
       pt.cex = 1.2,
       bty    = "n",
       yjust  = 0.5)

par(op)  # restore par settings

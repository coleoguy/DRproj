# Name: Megan Copeland
# Date: April 25th, 2025
# Purpose: Download Aedes aegypti genome (GCF_002204515.2_AaegL5.0),
#          filter to keep only NC_035107.1, NC_035108.1, NC_035109.1,
#          then run DirectRepeateR on the filtered FASTA.

#############################
##  1. SETUP & PREP      ##
#############################

# Define FTP URL and expected file size (bytes)
ftp_url       <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/",
  "GCF_002204515.2_AaegL5.0/",
  "GCF_002204515.2_AaegL5.0_genomic.fna.gz"
)
gz_file       <- basename(ftp_url)    # "GCF_002204515.2_AaegL5.0_genomic.fna.gz"
expected_size <- 410300578L           # compressed size

# Define output and filtering parameters
out_fa        <- "data/aedes_filt.fasta"
keep_ids      <- c("NC_035107.1", "NC_035108.1", "NC_035109.1")
header_pattern <- paste0("^>(", paste(keep_ids, collapse="|"), ")\\b")

#############################
##  2. DOWNLOAD GENOME    ##
#############################

download_with_curl <- function(url, dest) {
  message("Downloading via curl: ", url)
  res <- system2("curl", args = c("-L", "-o", shQuote(dest), shQuote(url)))
  if (res != 0) stop("curl failed (exit status ", res, ")")
}

if (!file.exists(gz_file) || file.size(gz_file) < expected_size) {
  if (file.exists(gz_file)) {
    message("→ Incomplete download (", file.size(gz_file),
            " < ", expected_size, "); removing and retrying.")
    file.remove(gz_file)
  }
  download_with_curl(ftp_url, gz_file)
} else {
  message("Found complete file: ", gz_file)
}

#############################
##  3. FILTER FASTA       ##
#############################

message("Filtering for: ", paste(keep_ids, collapse = ", "))
con_in  <- gzfile(gz_file, open = "rt")
con_out <- file(out_fa, open = "w")

keep_block <- FALSE
repeat {
  lines <- readLines(con_in, n = 1000, warn = FALSE)
  if (length(lines) == 0) break
  for (ln in lines) {
    if (startsWith(ln, ">")) {
      keep_block <- grepl(header_pattern, ln)
      if (keep_block) writeLines(ln, con = con_out)
    } else if (keep_block) {
      writeLines(ln, con = con_out)
    }
  }
}

close(con_in)
close(con_out)
message("Filtered FASTA written to: ", out_fa)

#############################
##  4. INSTALL & LOAD     ##
#############################

# Install devtools if needed, then DirectRepeateR
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
library(devtools)
install_github('coleoguy/DirectRepeateR', build_vignettes = TRUE)
library(DirectRepeateR)

#############################
##  5. RUN DirectRepeateR ##
#############################

# Arguments: fasta, min repeat length, max gap, step size, output CSV
GetRepeats(out_fa, 25, 20000, 50, "results/aedes_merged.csv")

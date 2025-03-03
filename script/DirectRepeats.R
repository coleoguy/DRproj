### Megan Copeland
### June 16th, 2024
### Direct Repeats
### mcc146@tamu.edu

# Load libraries
library(seqinr) # used to read in fasta file
library(data.table) # used for rbindlist

# Define parameters with default values
query_length <- 25
maxdist <- 20000
minlength <- 25

# Read in fasta, obtain chromosome lengths and names
genome <- read.fasta("GCF_000001405.40_GRCh38.p14_genomic.fna")
SeqLength <- sapply(genome, length)
chromosome_names <- names(genome)

# Function to create an empty results data.table
create_empty_results <- function(n = 1000000) {
  data.frame(Chromosome = character(n), 
             Start_Position = integer(n),
             Match_Position = integer(n),
             End_Position = integer(n),
             Match_End_Position = integer(n))
}

# Create a folder for temp files if it doesn't exist
if (!dir.exists("temp")) {
  dir.create("temp")
}

# Create a folder for temp files if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
}

# Condensing function
condense_results <- function(data) {
  chr <- unique(data$Chromosome)
  condensed_list <- list() 
  
  for (i in seq_along(chr)) {
    chromosome_data <- data[data$Chromosome == chr[i],] # subset data to the current chromosome
    chromosome_condensed <- list()  # store condensed results for the current chromosome
    result_index <- 1  # keep track of the index in the list
    
    j <- 1 
    n <- nrow(chromosome_data)
    while (j <= n) {
      start_pos <- chromosome_data$Start_Position[j]
      match_start_pos <- chromosome_data$Match_Position[j]
      end_pos <- chromosome_data$End_Position[j]
      match_end_pos <- chromosome_data$Match_End_Position[j]
      
      curset <- j
      
      k <- j + 1
      while (k <= n && chromosome_data$Start_Position[k] <= (end_pos + 1)) {
        if (chromosome_data$Start_Position[k] == (end_pos + 1) && 
            chromosome_data$Match_Position[k] == (match_end_pos + 1)) {
          end_pos <- chromosome_data$End_Position[k]
          match_end_pos <- chromosome_data$Match_End_Position[k]
          curset <- c(curset, k)
        }
        k <- k + 1
      }
      
      condensed_set <- data.table(Chromosome = chr[i], 
                                  Start_Position = start_pos, 
                                  End_Position = end_pos, 
                                  Match_Position = match_start_pos, 
                                  Match_End_Position = match_end_pos)
      chromosome_condensed[[result_index]] <- condensed_set
      result_index <- result_index + 1
      
      j <- curset[length(curset)] + 1
    }
    
    condensed_list[[chr[i]]] <- rbindlist(chromosome_condensed)
  }
  
  return(rbindlist(condensed_list))
}

# Initialize variables
raw_results <- create_empty_results() # store raw results
insert_idx <- 1 # keep track of where next set of results should be put in results df
file_count <- 1 # keep count of temp files

# Loop through genome and segment sequences into start positions according to query length
for (i in seq_along(genome)) {
  sequence <- paste(genome[[i]], collapse = "") # Get i-th sequence and collapse into single string
  seq_length <- SeqLength[i] # Get sequence length
  starts <- seq(from = 1, by = query_length, length.out = floor(seq_length / query_length)) # Get positions in increments of query length
  
  # Loop through start positions, define pattern and field, search for matches
  for (j in seq_along(starts)) {
    start_pos <- starts[j] # Get the j-th start position
    upper_bound <- min(start_pos + query_length - 1 + maxdist, seq_length) # define upper bound of field
    pattern <- substr(sequence, start_pos, start_pos + query_length - 1) # define pattern
    search_field <- substr(sequence, start_pos + query_length, upper_bound) # define field
    
    # Search for matches
    start_matches <- gregexpr(pattern = pattern, text = search_field, fixed = TRUE)[[1]]
    
    # If there is a match, proceed
    if (any(start_matches > 0)) {
      curr_match_pos <- start_pos + query_length - 1 + start_matches # calculate match position
      num_matches <- length(curr_match_pos) # find number of matches
      
      # Check if raw_results data frame is full, if so, proceed
      if ((insert_idx + num_matches - 1) > nrow(raw_results)) {
        # Remove rows with empty Chromosome column
        temp_results <- raw_results[1:(insert_idx - 1), ]
        temp_results <- temp_results[temp_results$Chromosome != "", ]
        
        # Apply condensing function before writing to temp file
        temp_results <- condense_results(temp_results)
        
        # Write current sets of results to a temp file
        fwrite(temp_results, paste0("temp/temp.", file_count, ".csv"))
        file_count <- file_count + 1 # increment file count
        raw_results <- create_empty_results() # initialize new, empty results df
        insert_idx <- 1 # reset insert index
      }
      
      # If data frame is not full, add results
      raw_results[insert_idx:(insert_idx + num_matches - 1), ] <- data.frame(
        Chromosome = rep(chromosome_names[i], num_matches), # rep used because of multiple matches
        Start_Position = rep(start_pos, num_matches),
        Match_Position = curr_match_pos,
        End_Position = rep(start_pos + query_length - 1, num_matches), # Calculate End_Position here
        Match_End_Position = curr_match_pos + query_length - 1 # Calculate Match_End_Position here
      )
      
      # Update insert index
      insert_idx <- insert_idx + num_matches
    }
  }
}

# Write any remaining results to a file
if (insert_idx > 1) {
  temp_results <- raw_results[1:(insert_idx - 1), ]
  temp_results <- temp_results[temp_results$Chromosome != "", ]
  temp_results <- condense_results(temp_results)
  fwrite(temp_results, paste0("temp/temp.", file_count, ".csv"))
}

# Combine all temp files
temp_files <- list.files("temp", pattern = "temp.*\\.csv", full.names = TRUE)
combined_results <- rbindlist(lapply(temp_files, fread))

# Remove temporary files
file.remove(temp_files)

# Apply the condensing function on the combined results
condensed_results <- condense_results(combined_results)

# Calculate repeat length
condensed_results[, Repeat_Length := End_Position - Start_Position]

# Filter rows where Repeat_Length is greater than minlength
filtered_results <- condensed_results[Repeat_Length > minlength]

# Drop the Repeat_Length column before saving
final_output <- filtered_results[, .(Chromosome, Start_Position, End_Position, Match_Position, Match_End_Position)]

# Save the final results
fwrite(final_output, "results/condensed_results.csv")



PlotRepeatsRollingAverage <- function(
    data,
    window_size    = 200000,  # Size of each bin (bp)
    step_size      = 200000,  # Step between bin starts (bp)
    rolling_window = 5        # Number of bins in the rolling average
) {
  # Load libraries (or ensure they are loaded in your session)
  library(ggplot2)
  library(zoo)  # for rollmean()
  
  # Get the list of unique chromosomes
  chr_names <- unique(data$Chromosome)
  
  # Loop through each chromosome
  for (chr in chr_names) {
    # Subset for the current chromosome
    chr_data <- subset(data, Chromosome == chr)
    
    # Determine the chromosome length
    chr_length <- max(c(chr_data$End_Position, chr_data$Match_End_Position))
    
    # Calculate midpoint of each repeat
    repeat_midpoints <- ((chr_data$Start_Position + chr_data$End_Position) / 2 +
                           (chr_data$Match_Position + chr_data$Match_End_Position) / 2) / 2
    
    # Define bin starts
    windows <- seq(1, chr_length, by = step_size)
    
    # Count repeats in each bin
    repeat_counts <- sapply(windows, function(w) {
      sum(repeat_midpoints >= w & repeat_midpoints < (w + window_size))
    })
    
    # Create a data frame with the binned results
    results_df <- data.frame(
      Window_Start = windows,
      Repeat_Count = repeat_counts,
      Chromosome   = chr
    )
    
    # Compute rolling average (centered, ignoring edges)
    # If rolling_window = 5, each point becomes the average of ±2 neighbors
    results_df$Rolling_Avg <- zoo::rollmean(
      results_df$Repeat_Count,
      k     = rolling_window,
      fill  = NA,   # or "extend" if you want to avoid NA at edges
      align = "center"
    )
    
    # Build the plot
    p <- ggplot(results_df, aes(x = Window_Start, y = Rolling_Avg)) +
      geom_line(color = "blue", size = 1) +
      geom_ribbon(
        aes(ymin = 0, ymax = Rolling_Avg),
        fill  = "lightblue",
        alpha = 0.4
      ) +
      labs(
        title = paste("Chromosome:", chr),
        x     = "Genomic Position (binned)",
        y     = "Rolling Average of Repeats"
      ) +
      theme_minimal()
    
    # Display the plot for the current chromosome
    print(p)
  }
}
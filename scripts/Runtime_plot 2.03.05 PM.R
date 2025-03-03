library(ggplot2)
library(ggrepel)

# Organism names
species <- c("C. elegans", "A. gambiae", "V. vinifera", 
             "O. latipes", "G. gallus", "A. aegypti", 
             "M. musculus", "H. sapiens")

# Genome sizes
genome_sizes <- c(100.3, 264.5, 494.9, 734, 1100, 1300, 2700, 3100)

# Updated Times in minutes
times_minutes <- c(4.15, 10, 20.28, 29.71, 
                   43.57, 48.16, 2.56 * 60, 3.61 * 60)

# Create data frame
data <- data.frame(Species = species, Genome_Size = genome_sizes, Time = times_minutes)

# Plot using geom_text_repel
ggplot(data, aes(x = Genome_Size, y = Time)) +
  geom_point(color = "black") +
  geom_text_repel(aes(label = Species),
                  size = 3.5,
                  fontface = "italic",
                  max.overlaps = Inf, hjust = .75) +  # allow as many repels as necessary
  labs(x = "Genome Size (Mb)",
       y = "Time (minutes)") +  
  xlim(min(genome_sizes) * 0.9, max(genome_sizes) * 1.1) +  
  theme(
    axis.text.x = element_text(hjust = .5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

### Megan Copeland
### August 24th, 2024
### DirectRepeateR Run Time Plot

# Install and load ggplot2
library(ggplot2)

# Organism names
species <- c("C. elegans", "A. gambiae", "V. vinifera", 
             "O. latipes", "G. gallus", "A. aegypti", 
             "M. musculus", "H. sapiens")
# Genome sizes
genome_sizes <- c(100.3, 264.5, 494.9, 734, 1100, 1300, 2700, 3100)

# Times in hours
times_hours <- c((11 + 50/60) / 60, 1.72, 2.23, 1.68, 
                 2.26, 3.32, 8.5, 17.72) 

# Create with names, genome sizes, and run time
data <- data.frame(Species = species, Genome_Size = genome_sizes, Time = times_hours)

# Plot
ggplot(data, aes(x = Genome_Size, y = Time, label = Species)) +
  geom_point(stat = "identity", fill = "black") +
  geom_text(aes(label = Species), vjust = -0.75, hjust = 0.45, size = 2.5, fontface = "italic") + 
  labs(x = "Genome Size (Mb)",
       y = "Time (hours)") +  
  xlim(min(genome_sizes) * 0.9, max(genome_sizes) * 1.1) +  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

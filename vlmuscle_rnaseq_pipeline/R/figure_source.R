##### Figure source file ###########






# theme function 
pl.theme <- function() {
  
  theme_bw() +
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(), 
          axis.line = element_line(size = line_size), 
          axis.ticks = element_line(size = line_size), 
          axis.text = element_text(color = "black", size = 8), 
          axis.title = element_text(color = "black", size = 8),
          legend.title = element_blank(), 
          legend.background = element_rect(fill = "white"),
          legend.margin = margin(t = 0, r = 1, b = 1, l = 1, unit = "pt"),
          legend.key = element_rect(fill = "white"),
          legend.position = c(0.85, 0.9)) 
  
  
  
}



# Color scale 

multiple.color <- "#636363"
single.color <- "#bdbdbd"

# line sizes
line_size <- 0.3

# text sizes
label.size <- 10

text.size <- 8


sequential.colors <- c("#fee5d9",  "#fcae91",  "#fb6a4a",  "#de2d26",  "#a50f15")

# diff.colors <- c("#e66101", "#fdb863", "#b2abd2", "#5e3c99")

diff.colors <- c("#7fc97f", "#beaed4",  "#fdc086")

# Colors used in volcano plots 
volcano.colors <- c("gray30", "red")




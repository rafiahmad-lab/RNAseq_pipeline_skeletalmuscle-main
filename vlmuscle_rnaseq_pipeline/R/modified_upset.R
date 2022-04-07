
library(UpSetR) # needed for the fromList function

# A function for upset plots using our three sets of models. 
# Be aware, the function will probably generalize bad to other situations

list.up <- venn_list_up[1:3]
list.down <- venn_list_down[1:3]




modified_upset <- function(list.up, # lists of genes/go terms that are up regulated
                           list.down, # list of genes/go terms down regulated
                           set.order = c("to_w2", "els_w2", "naive_w2"),
                           set.names = c("Tissue-offset", 
                                         "Effective library size", 
                                         "Naïve"), # Set names, order of display in set size fig, 
                           # should correspond to order
                           legend.labels = c("Week 2 > Week 0", "Week 2 < Week 0"),
                           legend.title = "Gene regulation",
                           line.size = 0.3, # size of lines in figures (may be changed later)
                           regulation.col = c("green", "red"), # Regulation coloring (up vs. down) 
                           bar.width = 0.4, # Width of bars in bar plots
                           text.size = 7, # text size in all figures
                           point.size = 4) { # point size in intersection matrix 
  
  # Calculate total count of each set (up regulated...)
  up <- data.frame(set.size = colSums(fromList(list.up)),
             set.name = colnames(fromList(list.up)), 
             regulation = "up",
             row.names = NULL)
  # and down regulated
  down <- data.frame(set.size = colSums(fromList(list.down)),
                     set.name = colnames(fromList(list.down)),
                     regulation = "down",
                     row.names = NULL)
    
  
# Calculate intersection sizes, i.e. genes/go terms that are regulated Up/down in in each intersection.
# This corresponds to distinct version https://jokergoo.github.io/ComplexHeatmap-reference/book/08-upset_files/figure-html/unnamed-chunk-7-1.png


  
  
intersections.up <- data.frame(fromList(list.up)) %>%
  dplyr::select(set1 = set.order[1], 
                set2 = set.order[2], 
                set3 = set.order[3]) %>%
    mutate(i1 = if_else(set1 == 1 & set2 == 0 & set3 == 0, 1, 0), 
           i2 = if_else(set1 == 0 & set2 == 1 & set3 == 0, 1, 0),
           i3 = if_else(set1 == 0 & set2 == 0 & set3 == 1, 1, 0), 
           i4 = if_else(set1 == 1 & set2 == 1 & set3 == 0, 1, 0),
           i5 = if_else(set1 == 0 & set2 == 1 & set3 == 1, 1, 0),
           i6 = if_else(set1 == 1 & set2 == 0 & set3 == 1, 1, 0),
           i7 = if_else(set1 == 1 & set2 == 1 & set3 == 1, 1, 0)) %>%
    dplyr::select(i1:i7) %>%
    summarize_all(list( ~ sum(.))) %>%
    
    pivot_longer(names_to = "intersection", values_to = "count", cols = i1:i7) %>%
    mutate(regulation = "up")
 
intersections.down <- data.frame(fromList(list.down)) %>%
  dplyr::select(set1 = set.order[1], 
                set2 = set.order[2], 
                set3 = set.order[3]) %>%
  mutate(i1 = if_else(set1 == 1 & set2 == 0 & set3 == 0, 1, 0), 
         i2 = if_else(set1 == 0 & set2 == 1 & set3 == 0, 1, 0),
         i3 = if_else(set1 == 0 & set2 == 0 & set3 == 1, 1, 0), 
         i4 = if_else(set1 == 1 & set2 == 1 & set3 == 0, 1, 0),
         i5 = if_else(set1 == 0 & set2 == 1 & set3 == 1, 1, 0),
         i6 = if_else(set1 == 1 & set2 == 0 & set3 == 1, 1, 0),
         i7 = if_else(set1 == 1 & set2 == 1 & set3 == 1, 1, 0)) %>%
  dplyr::select(i1:i7) %>%
  summarize_all(list( ~ sum(.))) %>%
  
  pivot_longer(names_to = "intersection", values_to = "count", cols = i1:i7) %>%
  mutate(regulation = "down")

intersections <- rbind(intersections.up,intersections.down)


max.intersections <- intersections %>%
  group_by(intersection) %>%
  summarise(sum = sum(count)) %>%
  ungroup() %>%
  summarise(max = max(sum)) %>%
  pull(max)




## Round up function 



roundUpNice <- function(x) {
  
  
  if(x > 10000) {return(round(x + 50, -2))}
  if(x < 10000 & x >= 1000) {return(round(x + 5, -1))}
  if(x < 1000) {return(round(x + 2, -1))}
  
  }




  ##### Intersection figs


intersection.fig <- intersections %>%
  
  mutate(regulation = factor(regulation, levels = c("up", "down"), labels = legend.labels)) %>%
  ggplot(aes(intersection, count, fill = regulation)) + 
  
  geom_bar(stat = "identity", width = bar.width) +
  
  scale_fill_manual(values = regulation.col)  +
  
  labs(fill = legend.title) + 
  

  
  theme(legend.title = element_text(size = 7), 
        legend.text = element_text(size = 7), 
        legend.key.size = unit(0.3, "cm"))
  
  
  
legend <- get_legend(intersection.fig)

intersection.fig <- intersection.fig +
  scale_y_continuous(limits = c(0, roundUpNice(max.intersections)), 
                     breaks = c(0, 
                                (roundUpNice(max.intersections)/2), 
                                roundUpNice(max.intersections)), 
                     expand = c(0,0)) +
  
  theme_classic() +
  
  theme(axis.title = element_blank(), 
        
        plot.margin = unit(c(0.2, 0.2, -0.2, 0.2), "cm"),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = text.size, color = "black"),
        
     
        axis.line.x = element_blank(),
        axis.line.y = element_line(size = line.size),
        axis.ticks.y = element_line(size = line.size),

        axis.ticks.x = element_blank())
  


  
  
  
  #### Set sizes fig 
    
   set.size.dat <- rbind(up, down)



  max.set <- set.size.dat %>%
    group_by(set.name) %>%
    summarise(sum = sum(set.size)) %>%
    ungroup() %>%
    summarise(max = max(sum)) %>%
    pull(max)
     
  max.set <- roundUpNice(max.set)
  
  set.size.fig <- set.size.dat %>%
    mutate(set.name = factor(set.name, levels = rev(set.order), labels = rev(set.names))) %>%

    
    
    ggplot(aes(-set.size, set.name, fill = regulation)) + 
    
    geom_bar(stat = "identity", width = set.column.width) +
    
    scale_x_continuous(limits = c(-max.set, 0), 
                       breaks = c(-max.set, -(max.set/2), 0), 
                       labels = c(max.set,round((max.set/2),0), 0), 
                       expand = c(0,0)) +
    
    scale_y_discrete(position = "right") +
    
    scale_fill_manual(values = regulation.col) +
    
    theme_classic() +
    
    theme(axis.title = element_blank(), 
          
          plot.margin = unit(c(0, 0, 0,0), "cm"),
          
          legend.position = "none",
          axis.text.x = element_text(hjust = c(0, 0.5, 0.5),
                                     margin =  margin(10, 0, 0, 0), 
                                     color = "black", 
                                     size = text.size),
          axis.text.y = element_text(size = text.size, color = "black"),
          
          axis.line.x = element_line(size = line.size),
          axis.line.y = element_blank(),
          axis.ticks.length.x=unit(-0.1, "cm"), 
          axis.ticks.x = element_line(size = line.size),
          axis.ticks.y = element_blank())
    
  #### Make intersection matrix #####
  
  
  intersect.fig <- ggplot(data.frame(x = c(1, 7), 
                                     y = c(1, 3)), 
                          aes(x = x, y = y)) + 

    annotate("rect", 
             xmin = c(0.5, 0.5), 
             xmax = c(7.5, 7.5), 
             ymin = c(0.7, 2.7), 
             ymax = c(1.3, 3.3), 
             color = "gray90", 
             fill = "gray90") +

    annotate("segment", 
             x =    c(4, 5, 6, 7), 
             xend = c(4, 5, 6, 7), 
             y =    c(2, 1, 1, 1), 
             yend = c(3, 2, 3, 3)) +

    annotate("point", 
             x = c(1, 1, 1,
                   2, 2, 2,
                   3, 3, 3,
                   4, 4, 4,
                   5, 5, 5,
                   6, 6, 6,
                   7, 7, 7), 
             y = c(3, 2, 1, 
                   3, 2, 1,
                   3, 2, 1,
                   3, 2, 1,
                   3, 2, 1,
                   3, 2, 1,
                   3, 2, 1),
             size = point.size, 
             color = "grey70") + 
    
    
    annotate("point", x = c(1, 2, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7), 
             y = c(3, 2, 1, 3, 2, 2, 1, 3, 1,3,2,1),
             size = point.size, color = "black") + 
    
    scale_x_continuous(limits = c(0.5, 7.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 3.5), expand = c(0, 0)) +
    
    theme_classic() +
    theme(axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))
    
    

  ## create legend for inset plot 
  

  
  # Return the figure parts as a list 
  return(list(set.size.fig = set.size.fig, 
              intersection.fig = intersection.fig, 
              intersection.matrix = intersect.fig, 
              legend = legend))
  

}






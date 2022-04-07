#### Figure 5 --- DE analysis of different normalization methods TIME EFFECTS #####


source("./R/lib_fun.R")
source("./R/figure_source.R")
source("./R/modified_upset.R")


#### Assess the effect of Poisson vs. negative binomial model -- 
# Is the extra parameter needed?



## Time effects
mm2_2 <- readRDS(file = "./data/derivedData/DE/mixedmodel2_results2.RDS")


mm2_1 <- readRDS(file = "./data/derivedData/DE/mixedmodel2_results1.RDS")

### Sets comparisons #####

  volcano_data <- mm2_2 %>%
  filter(coef %in% c("timew2pre", "timew12")) %>%
  mutate(coef = factor(coef, levels = c("timew2pre", "timew12"),
                       labels = c("Week 2 - Week 0", 
                                  "Week 12 - Week 0")),
         model = factor(model, levels = c("naive", "lib_size_normalized", "tissue_offset_lib_size_normalized"), 
                        labels = c("No normalization", 
                                   "Effective library size", 
                                   "Per tissue mass +\nEffective library size"))) %>%
  
  group_by(coef, model) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"), 
         sig = if_else(fcthreshold == "s" & pthreshold == "s", "s", "ns"), 
         regulation = if_else(sig == "s" & log2fc > 0.5, "up", 
                              if_else(sig == "s" & log2fc < -0.5, "down", "ns")), 
         regulation = factor(regulation, levels = c("ns", "up", "down"))) %>%
  print()





##### Settings specifically for figure 3 ##################

## Color themes in figure source 
# using volcano.regulation.color


## Sizes volcano plots

volcano_size <- c(0.8, 1.1, 1.1)
volcano_alpha <- c(0.2, 0.8, 0.8)



############## Week 2 Volcano plots #####################################

week2_naive <- volcano_data %>%
  filter(model == "No normalization" & coef == "Week 2 - Week 0") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  # Annotate plot
  annotate("text", x = -2.5, y = 76, label = "Naïve", hjust = 0, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2,  0, 2, 4, 6, 8), limits = c(-3, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80), expand = c(0, 0)) +
  
  labs(x = "", 
       y = "") +
  
  pl.theme() +
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),  
        # panel.border = element_rect(colour = "black"),
        legend.position = "none")  

week2_libsize <- volcano_data %>%
  filter(model == "Effective library size" & coef == "Week 2 - Week 0") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -2.5, y = 79, label = "Effective library-\nsize", hjust = 0,vjust = 1, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2,  0, 2, 4, 6, 8), limits = c(-3, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "") +
  
  pl.theme() +
  
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),  
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week2_tissue <- volcano_data %>%
  filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 2 - Week 0") %>%
  
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -2.5, y = 76, label = "Tissue offset", hjust = 0, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2,  0, 2, 4, 6, 8), limits = c(-3, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80), limits = c(0, 80), expand = c(0, 0)) +
  labs(x = "", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),  
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 

#### Upset plots for time-effects ###############



venn_list_up <- list() 
venn_list_down <- list()

Coef <- c("timew2pre", "timew2pre", "timew2pre",
          "timew12", "timew12", "timew12")
Model <- c("naive", "lib_size_normalized", "tissue_offset_lib_size_normalized", 
           "naive", "lib_size_normalized", "tissue_offset_lib_size_normalized")


for(i in 1:6){
  
  temp_up <- mm2_2 %>%
    group_by(model, coef) %>%
    mutate(adj.p = p.adjust(p.val, method = "fdr")) %>% 
    mutate(pt = if_else(adj.p  < 0.05, "s", "ns"), 
           sig.gene = if_else(pt == "s", gene, "ns")) %>%
    
    filter(sig.gene != "ns", 
           estimate > 0) %>%
    filter(coef == Coef[i] & model == Model[i]) %>%
    
    data.frame()
  
  
  temp_down <- mm2_2 %>%
    group_by(model, coef) %>%
    mutate(adj.p = p.adjust(p.val, method = "fdr")) %>% 
    mutate(pt = if_else(adj.p  < 0.05, "s", "ns"), 
           sig.gene = if_else(pt == "s", gene, "ns")) %>%
    
    filter(sig.gene != "ns", 
           estimate < 0) %>%
    filter(coef == Coef[i] & model == Model[i]) %>%
    
    data.frame()
  
  
  venn_list_up[[i]] <- temp_up$sig.gene
  venn_list_down[[i]] <- temp_down$sig.gene
  
}


names(venn_list_up) <- c("naive_w2", "els_w2", "to_w2", "naive_w12", "els_w12", "to_w12")
names(venn_list_down) <- c("naive_w2", "els_w2", "to_w2", "naive_w12", "els_w12", "to_w12")



upset.list.w2 <- modified_upset(list.up = venn_list_up[1:3], 
                         list.down = venn_list_down[1:3], 
                         set.order = c("to_w2", "els_w2", "naive_w2"),
                         set.names = c("Tissue-offset", 
                                       "Effective\nlibrary size", 
                                       "Naïve"), # Set names, order of display in set size fig, 
                         # should correspond to order
                         
                         line.size = 0.4, # size of lines in figures (may be changed later)
                         regulation.col = c("#2ca25f","#f03b20"), # Regulation coloring (up vs. down) 
                         bar.width = 0.4, # Width of bars in bar plots
                         text.size = 7, # text size in all figures
                         point.size = 3)  # point size in intersection matrix )


upset.list.w12 <- modified_upset(list.up = venn_list_up[4:6], 
                                list.down = venn_list_down[4:6], 
                                set.order = c("to_w12", "els_w12", "naive_w12"),
                                set.names = c("Tissue-offset", 
                                              "Effective\nlibrary size", 
                                              "Naïve"), # Set names, order of display in set size fig, 
                                # should correspond to order
                                legend.labels = c("Week 12 > Week 0", "Week 12 < Week 0"),
                                legend.title = "Gene regulation",
                                line.size = 0.4, # size of lines in figures (may be changed later)
                                regulation.col = c("#2ca25f","#f03b20"), # Regulation coloring (up vs. down) 
                                bar.width = 0.4, # Width of bars in bar plots
                                text.size = 7, # text size in all figures
                                point.size = 3)  # point size in intersection matrix )





upset_w2pre <- ggdraw(plot_grid(NULL,
                         plot_grid(upset.list.w2$set.size.fig, upset.list.w2$intersection.matrix, 
                                   ncol = 2, align = "h"), nrow = 2)) +
  draw_plot(upset.list.w2$intersection.fig, x = 0.42, y = 0.47, width = 0.59, height = 0.55) +
  draw_plot(upset.list.w2$legend, x = 0.05, y = 0.55, width = 0.3, height = 0.3)
  




upset_w12 <- ggdraw(plot_grid(NULL,
                       plot_grid(upset.list.w12$set.size.fig, upset.list.w12$intersection.matrix, 
                                 ncol = 2, align = "h"), nrow = 2)) +
  draw_plot(upset.list.w12$intersection.fig, x = 0.42, y = 0.47, width = 0.59, height = 0.55) +
  draw_plot(upset.list.w12$legend, x = 0.05, y = 0.55, width = 0.3, height = 0.3)




## Change parameters in the individual plots







####### Week 12 volcano plots ##################################

week12_naive <- volcano_data %>%
  filter(model == "No normalization" & coef == "Week 12 - Week 0") %>%

  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2,  0, 2, 4, 6), limits = c(-3, 6), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0,10, 20,30, 40), limits = c(0, 40), expand = c(0, 0)) +
  
  annotate("text", x = -2.5, y = 38, label = "Naïve", hjust = 0, size = 2.2) +
  
  labs(x = "", 
       y = "") +
  
  pl.theme() +
  
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),  
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 

week12_libsize <- volcano_data %>%
  filter(model == "Effective library size" & coef == "Week 12 - Week 0") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -2.5, y = 39, label = "Effective library-\nsize", hjust = 0,vjust = 1, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2,  0, 2, 4, 6), limits = c(-3, 6), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0,10, 20,30, 40), limits = c(0, 40), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "") +
  
  pl.theme() +
  
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),  
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week12_tissue <- volcano_data %>%
  filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 12 - Week 0") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -2.5, y = 38, label = "Tissue-offset", hjust = 0, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2,  0, 2, 4, 6), limits = c(-3, 6), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0,10, 20,30, 40), limits = c(0, 40), expand = c(0, 0)) +
  labs(x = "", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),  
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 



week2_grid <-  plot_grid(plot_grid(week2_tissue, week2_libsize, week2_naive, ncol = 3),
                         
                         plot_grid(NULL, upset_w2pre, NULL, ncol = 3, rel_widths = c(0.01, 0.6, 0.01)), 
                         rel_heights  = c(0.4, 0.6), nrow = 2)




week12_grid <-  plot_grid(plot_grid(week12_tissue, week12_libsize, week12_naive,   ncol = 3),
                          
                          plot_grid(NULL, upset_w12, NULL,  ncol = 3, rel_widths = c(0.01, 0.6, 0.01)), 
                          rel_heights = c(0.4, 0.6), nrow = 2)





figure5 <- plot_grid(week2_grid, 
                     
                     week12_grid,
                     
                     rel_heights = c(0.5, 0.5),
 
                     ncol = 1) +
  
  
  draw_plot_label(label=c("A",  "B"),
                  x =   c(0.02, 0.02), 
                  y =   c(0.99, 0.49),
                  hjust=.5, vjust=.5, size = label.size)




# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure5.pdf", plot = figure5, width = 8.5, height = 18, 
       dpi = 600,
       units = "cm", device=cairo_pdf)






#### Figure 4 --- DE analysis of different normalization methods #####


source("./R/lib_fun.R")
source("./R/figure_source.R")



#### Assess the effect of Poisson vs. negative binomial model -- 
# Is the extra parameter needed?

## Interaction model 
mm_interaction <- readRDS("./data/derivedData/DE/mixedmodel2_results1.RDS")

## Time model
mm_time <- readRDS("./data/derivedData/DE/mixedmodel2_results2.RDS")



######## Fig a: Muscle weight figure ############



# Muscle weights in cDNA synthesis are loaded ##
mw <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc * elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  inner_join(read_csv2("./data/oneThreeSetLeg.csv") %>%
               gather("sets", "leg", multiple:single)) %>%
  filter(rnaseq_include == "incl") %>%
  filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  filter(!(subject == "FP15" & time == "w2pre")) %>%
  print()



m0 <- lme(log(tissue) ~ time + time:sets, random = list(subject = ~1), data = mw, 
          control = lmeControl(msMaxIter = 100, msVerbose = TRUE, 
                               opt = "optim"))

m1 <- lme(log(tissue) ~ time + time:sets, random = list(subject = ~1 + time), data = mw, 
          control = lmeControl(msMaxIter = 100, msVerbose = TRUE, 
                               opt = "optim"))

anova(m0, m1)
plot(m1)


mw_stats <- data.frame(broom::tidy(m1, effects = "fixed")) %>%
  rowwise() %>%
  mutate(p.flag = if_else(term %in% c("timew12", "timew2pre"), publR::pval(p.value, flag = TRUE, sign = "\U2020"), 
                          if_else(term %in% c("timew0:setssingle", 
                                              "timew2pre:setssingle", 
                                              "timew12:setssingle"), 
                                  publR::pval(p.value, flag = TRUE), "Itercept"))) %>%
  data.frame() %>%
  print()




muscle_weight_fig <- emmeans(m1, specs = ~"time|sets") %>%
  data.frame() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
         sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
  ggplot(aes(time, exp(emmean), fill = sets)) + 
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                position = position_dodge(width = 0.35), 
                width = 0, 
                size = line_size) +
  geom_point(position = position_dodge(width = 0.35), 
             shape = 21, 
             size = 1.5) +
  
  scale_fill_manual(values = c("#bdbdbd", "#636363")) +
  scale_y_continuous(limits = c(2, 4), 
                     breaks = c(2, 2.5, 3,3.5, 4), 
                     labels = c(2,  "", 3, "", 4),
                     expand = c(0, 0)) +
  
  
  ggplot2::annotate("text", x = c(2,3), y = c(2.9, 3.05), 
                    label = c(mw_stats[6,6], mw_stats[5,6]), 
                    size = 5) +
  
  
  ggplot2::annotate("text", x = c(2,3), y = c(3.2, 3.35), label = c(mw_stats[3,6], mw_stats[2,6]), 
                    size = 2.5) +
  
  
  
  labs(fill = "", 
       x = "Time-point", 
       y = "Muscle biopsy mass in\n cDNA synthesis (mg)") +
  pl.theme() + 
  theme(axis.title.x = element_blank(), 
        
        legend.position = c(0.53, 0.94), 
        legend.text = element_text(size = 7, margin = margin(t = 0)), 
        legend.spacing.x = unit(0.1, 'cm'),
        legend.key.size = unit(0.3, 'lines'),
        legend.margin = margin(t = -0.2, r = -0.2, b = -0.2, l = -0.2),
        axis.title.y = element_text(size = 7, margin = margin(t = 0, r = -0.2, b = 0, l = 0))) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2))




### ################### Total mRNA counts (library size) per time-point ##########



# read DGEList

dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")



# using rsem for this example #

selected.method <- "rsem"


# Muscle weight in cDNA synthesis

# Muscle weights in cDNA synthesis are loaded ##
lib.size.data <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc * elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  inner_join(read_csv2("./data/oneThreeSetLeg.csv") %>%
               gather("sets", "leg", multiple:single)) %>%
  filter(rnaseq_include == "incl") %>%
  filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  filter(!(subject == "FP15" & time == "w2pre")) %>% 
  # This filters out one subject with bad estimation of RNA to mg tissue estimate
  # See script in training study for details
  inner_join(dge_lists[[selected.method]]$samples %>% 
               mutate(eff.lib = lib.size * norm.factors)) %>%
  mutate(eff.lib = as.integer(round(eff.lib, 0)), 
         eff.lib.scaled = eff.lib/10^6) %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12"))) %>%
  print()



# Modeling effective library size 
lib.m0 <- lme(log(eff.lib.scaled) ~ time + time:sets, random = list(subject = ~ 1), 
              data = lib.size.data)

lib.m1 <- lme(log(eff.lib.scaled) ~ time + time:sets, random = list(subject = ~ 1),
              weights = varIdent(form = ~ 1|time), 
              data = lib.size.data)


anova(lib.m0, lib.m1)

plot(lib.m0, main = "m0"); plot(lib.m1, main = "m1")

# Evidence for accounting for heteroscedasticity, using model lib.m1

# Including random slope per time

lib.m2 <- lme(log(eff.lib.scaled) ~ time + time:sets, random = list(subject = ~ 1 + time),
              weights = varIdent(form = ~ 1|time), 
              data = lib.size.data)

anova(lib.m1, lib.m2) # not called for 




# Modeling library size per tissue weight

# negative binomial model

# library per tissue
lt.m0 <- lme(log(eff.lib.scaled/tissue) ~ time + time:sets, random = list(subject = ~1), 
             data = lib.size.data)


lt.m1 <- lme(log(eff.lib.scaled/tissue) ~ time + time:sets, random = list(subject = ~1), 
             data = lib.size.data, weights = varIdent(form = ~ 1|time))



# Accounting for heteroscedasticity? 

plot(lt.m0, main = "Library size per tissue, m0"); plot(lt.m1, main = "Library size per tissue, m1") 

anova(lt.m0, lt.m1) 

# Clear evidence of better models using weights

# Random slopes? 

lt.m2 <- lme(log(eff.lib.scaled/tissue) ~ time + time:sets, random = list(subject = ~1 + time), 
             data = lib.size.data, weights = varIdent(form = ~ 1|time),
             control = lmeControl(msMaxIter = 100, msVerbose = TRUE, 
                                  opt = "optim"))

anova(lt.m1, lt.m2) # not called for



# Extract statistics
lt.m1.stats <- data.frame(broom::tidy(lt.m1, effects = "fixed")) %>%
  rowwise() %>%
  mutate(p.flag = if_else(term %in% c("timew12", "timew2pre"), publR::pval(p.value, flag = TRUE, sign = "\U2020"), 
                          if_else(term %in% c("timew0:setssingle", 
                                              "timew2pre:setssingle", 
                                              "timew12:setssingle"), 
                                  publR::pval(p.value, flag = TRUE), "Itercept"))) %>%
  data.frame() %>%
  print()


lib.m1.stats <- data.frame(broom::tidy(lib.m1, effects = "fixed")) %>%
  rowwise() %>%
  mutate(p.flag = if_else(term %in% c("timew12", "timew2pre"), publR::pval(p.value, flag = TRUE, sign = "\U2020"), 
                          if_else(term %in% c("timew0:setssingle", 
                                              "timew2pre:setssingle", 
                                              "timew12:setssingle"), 
                                  publR::pval(p.value, flag = TRUE), "Itercept"))) %>%
  data.frame() %>%
  print()



# Plotting the data 
per.rna <- emmeans(lib.m1, specs = ~ "time|sets") %>%
  data.frame() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
         sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
  
  ggplot(aes(time, exp(emmean), fill = sets )) +
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                position = position_dodge(width = 0.35), 
                width = 0, 
                size = line_size) +
  geom_point(position = position_dodge(width = 0.35), 
             shape = 21, 
             size = 1.5) +
  
  scale_fill_manual(values = c("#bdbdbd", "#636363")) +
  scale_y_continuous(limits = c(8, 16), 
                     breaks = c(8, 10, 12, 14, 16), 
                     labels = c(8, "", 12, "", 16),
                     expand = c(0, 0)) +
  
  # Sets effects stats segment
  
  # geom_segment(data = data.frame(x1 = c(1.9, 2.9), 
  #                                x2 = c(2.1, 3.1), 
  #                                y1 = c(2.92, 3.07), 
  #                                y2 = c(2.92, 3.07)), 
  #              aes(x = x1, xend = x2, y = y1, yend = y2), 
  #              inherit.aes = FALSE) +
  ggplot2::annotate("text", x = c(2,3), y = c(14, 15.3), label = c(lib.m1.stats[2,6], lib.m1.stats[3,6]), 
                    size = 2.5) +
  
  
  ggplot2::annotate("text", x = c(2,3), y = c(13.8, 15.2), label = c(lib.m1.stats[5,6], lib.m1.stats[6,6]), 
                    size = 5) +
  
  labs(fill = "", 
       x = "Time-point", 
       y = bquote("Effective library size<br>(million counts)")) +
  pl.theme() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_markdown(size = 7, margin = margin(t = 0, r = -0.2, b = 0, l = 0))) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2))



per.tissue <- emmeans(lt.m1, specs = ~ "time|sets") %>%
  data.frame() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
         sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
  
  ggplot(aes(time, exp(emmean), fill = sets )) +
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                position = position_dodge(width = 0.35), 
                width = 0, 
                size = line_size) +
  geom_point(position = position_dodge(width = 0.35), 
             shape = 21, 
             size = 1.5) +
  
  scale_fill_manual(values = c("#bdbdbd", "#636363")) +
  scale_y_continuous(limits = c(2, 6), 
                     breaks = c(2, 3, 4, 5, 6), 
                     labels = c(2, "", 4, "", 6),
                     expand = c(0, 0)) +
  
  # Sets effects stats segment
  
  # geom_segment(data = data.frame(x1 = c(1.9, 2.9), 
  #                                x2 = c(2.1, 3.1), 
  #                                y1 = c(2.92, 3.07), 
  #                                y2 = c(2.92, 3.07)), 
  #              aes(x = x1, xend = x2, y = y1, yend = y2), 
  #              inherit.aes = FALSE) +
  annotate("text", x = c(2,3), y = c(5.4, 5.9), label = c(lt.m1.stats[2,6], lt.m1.stats[3,6]), 
           size = 2.5) +
  
  
  labs(fill = "", 
       x = "Time-point", 
       y = "Effective library size<br>(million counts mg^-1)") +
  pl.theme() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        axis.title.y = element_markdown()) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))







mm2_2 <- readRDS(file = "./data/derivedData/DE/mixedmodel2_results2.RDS")


mm2_1 <- readRDS(file = "./data/derivedData/DE/mixedmodel2_results1.RDS")

### Sets comparisons #####

volcano_data <- mm2_1 %>%
  filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                     "timew12:setsmultiple")) %>%
  mutate(coef = factor(coef, levels = c("timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple"),
                       labels = c("Week 0: Multiple - Single-set", 
                                  "Week 2: Multiple - Single-set", 
                                  "Week 12: Multiple - Single-set")),
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
  filter(model == "No normalization" & coef == "Week 2: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  # Annotate plot
  annotate("text", x = -3.5, y = 7.8, label = "Naïve", hjust = 0, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
  
  labs(x = "", 
       y = "") +
  
  pl.theme() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none")  

week2_libsize <- volcano_data %>%
  filter(model == "Effective library size" & coef == "Week 2: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -3.5, y = 7.8, label = "Effective library-\nsize", hjust = 0,vjust = 1, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
  
  
  labs(x = "Log2 fold change", 
       y = "") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week2_tissue <- volcano_data %>%
  filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 2: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -3.5, y = 7.8, label = "Tissue offset", hjust = 0, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
  
  labs(x = "", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 

#### Upset plots for sets-effects ###############


venn_list_up <- list() 
venn_list_down <- list()

Coef <- c("timew2pre:setsmultiple", "timew2pre:setsmultiple", "timew2pre:setsmultiple",
          "timew12:setsmultiple", "timew12:setsmultiple", "timew12:setsmultiple")
Model <- c("naive", "lib_size_normalized", "tissue_offset_lib_size_normalized", 
           "naive", "lib_size_normalized", "tissue_offset_lib_size_normalized")



for(i in 1:6){
  
  temp_up <- mm2_1 %>%
    group_by(model, coef) %>%
    mutate(adj.p = p.adjust(p.val, method = "fdr")) %>% 
    mutate(pt = if_else(adj.p  < 0.05, "s", "ns"), 
           sig.gene = if_else(pt == "s", gene, "ns")) %>%
    
    filter(sig.gene != "ns", 
           estimate > 0) %>%
    filter(coef == Coef[i] & model == Model[i]) %>%
    
    data.frame()
  
  
  temp_down <- mm2_1 %>%
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
                                legend.labels = c("Multiple-set > Single-set", "Multiple-set < Single-set"),
                                legend.title = "Gene regulation Week 2",
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
                                 legend.labels = c("Multiple-set > Single-set", "Multiple-set < Single-set"),
                                 legend.title = "Gene regulation Week 12",
                                 line.size = 0.4, # size of lines in figures (may be changed later)
                                 regulation.col = c("#2ca25f","#f03b20"), # Regulation coloring (up vs. down) 
                                 bar.width = 0.4, # Width of bars in bar plots
                                 text.size = 7, # text size in all figures
                                 point.size = 3)  # point size in intersection matrix )





upset_w2pre <- ggdraw(plot_grid(NULL,
                                plot_grid(upset.list.w2$set.size.fig, upset.list.w2$intersection.matrix, 
                                          ncol = 2, align = "h"), nrow = 2)) +
  draw_plot(upset.list.w2$intersection.fig, x = 0.45, y = 0.47, width = 0.56, height = 0.55) +
  draw_plot(upset.list.w2$legend, x = 0.05, y = 0.55, width = 0.3, height = 0.3)





upset_w12 <- ggdraw(plot_grid(NULL,
                              plot_grid(upset.list.w12$set.size.fig, upset.list.w12$intersection.matrix, 
                                        ncol = 2, align = "h"), nrow = 2)) +
  draw_plot(upset.list.w12$intersection.fig, x = 0.45, y = 0.47, width = 0.56, height = 0.55) +
  draw_plot(upset.list.w12$legend, x = 0.05, y = 0.55, width = 0.3, height = 0.3)







####### Week 12 volcano plots ##################################

week12_naive <- volcano_data %>%
  filter(model == "No normalization" & coef == "Week 12: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
  
  annotate("text", x = -3.5, y = 7.8, label = "Naïve", hjust = 0, size = 2.2) +
  
  labs(x = "", 
       y = "") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 

week12_libsize <- volcano_data %>%
  filter(model == "Effective library size" & coef == "Week 12: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -3.5, y = 7.8, label = "Effective library-\nsize", hjust = 0,vjust = 1, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
  
  
  labs(x = "Log2 fold change", 
       y = "") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week12_tissue <- volcano_data %>%
  filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 12: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = regulation, alpha = regulation, size = regulation)) + 
  geom_point(shape = 21) +
  
  annotate("text", x = -3.5, y = 7.8, label = "Tissue-offset", hjust = 0, size = 2.2) +
  
  # Scales
  scale_color_manual(values = volcano.regulation.color) + 
  scale_fill_manual(values = volcano.regulation.color) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8), limits = c(0, 8), expand = c(0,0)) +
  
  
  labs(x = "", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 




######### Volcano and upset configuration 

week2_grid <-  plot_grid(plot_grid(week2_tissue, week2_libsize, week2_naive, ncol = 3),
                         
                         plot_grid(NULL, upset_w2pre, NULL, ncol = 3, rel_widths = c(0.01, 0.6, 0.01)), 
                         rel_heights  = c(0.4, 0.6), nrow = 2)




week12_grid <-  plot_grid(plot_grid(week12_tissue, week12_libsize, week12_naive,   ncol = 3),
                          
                          plot_grid(NULL, upset_w12, NULL,  ncol = 3, rel_widths = c(0.01, 0.6, 0.01)), 
                          rel_heights = c(0.4, 0.6), nrow = 2)








############## Complete figure ###################



figure4 <- plot_grid(plot_grid(muscle_weight_fig, per.rna, per.tissue, ncol = 3, 
                               rel_widths = c(0.3, 0.3, 0.3), 
                               align = "h"), 
                     
                     week2_grid, 
                     
                     week12_grid,
                     
                     rel_heights = c(0.18, 0.4, 0.4),
                     
                     #timeweek2_grid, timeweek12_grid, 
                     ncol = 1) +
  
  
  draw_plot_label(label=c("A",  "B",  "C",  "D",  "E",   "F", "G"),
                  x =   c(0.02, 0.37, 0.71, 0.02,  0.12, 0.02, 0.12), 
                  y =   c(0.99, 0.99, 0.99, 0.81,  0.62, 0.4,  0.22),
                  hjust=.5, vjust=.5, size = label.size)




# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure4.pdf", plot = figure4, width = 8.5, height = 23, 
       dpi = 600,
       units = "cm", device=cairo_pdf)


# GO overlap alternatively gene-ovelap between methods?



### Time figures ###
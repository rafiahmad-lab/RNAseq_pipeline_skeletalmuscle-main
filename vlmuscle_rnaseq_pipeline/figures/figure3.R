#### Figure 3 --- DE analysis of different normalization methods #####


source("./R/lib_fun.R")
source("./R/figure_source.R")



#### Assess the effect of Poisson vs. negative binomial model -- 
# Is the extra parameter needed?

## Interaction model 
mm_interaction <- readRDS("./data/derivedData/DE/mixedmodel2_results1.RDS")

## Time model
mm_time <- readRDS("./data/derivedData/DE/mixedmodel2_results2.RDS")



mm_time %>%
  filter(coef == "timew2pre", 
         model == "tissue_offset_lib_size_normalized") %>%
  ggplot(aes(poisson.test)) + geom_histogram(bins = 100) 













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
                position = position_dodge(width = 0.25), 
                width = 0.2, 
                size = line_size) +
  geom_point(position = position_dodge(width = 0.25), 
             shape = 21, 
             size = 2) +
  
  scale_fill_manual(values = c("#bdbdbd", "#636363")) +
  scale_y_continuous(limits = c(2.25, 3.25), 
                     breaks = c(2.25, 2.5, 2.75, 3, 3.25), 
                     labels = c("", 2.50, "", "3.0", ""),
                     expand = c(0, 0)) +
  
  # Sets effects stats segment
  
  # geom_segment(data = data.frame(x1 = c(1.9, 2.9), 
  #                                x2 = c(2.1, 3.1), 
  #                                y1 = c(2.92, 3.07), 
  #                                y2 = c(2.92, 3.07)), 
  #              aes(x = x1, xend = x2, y = y1, yend = y2), 
  #              inherit.aes = FALSE) +
  annotate("text", x = c(2,3), y = c(2.9, 3.05), label = c(mw_stats[6,6], mw_stats[5,6]), 
           size = 5) +
  
  
  annotate("text", x = c(2,3), y = c(2.98, 3.12), label = c(mw_stats[3,6], mw_stats[2,6]), 
           size = 2.5) +
  
  
  
  labs(fill = "", 
       x = "Time-point", 
       y = "Muscle biopsy mass in\n cDNA synthesis (mg)") +
  pl.theme() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none")




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
                position = position_dodge(width = 0.25), 
                width = 0.1, 
                size = line_size) +
  geom_point(position = position_dodge(width = 0.25), 
             shape = 21, 
             size = 2) +
  
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
  annotate("text", x = c(2,3), y = c(14, 15.3), label = c(lib.m1.stats[2,6], lib.m1.stats[3,6]), 
           size = 2.5) +
  
  
  annotate("text", x = c(2,3), y = c(13.8, 15.1), label = c(lib.m1.stats[5,6], lib.m1.stats[6,6]), 
           size = 5) +
  
  labs(fill = "", 
       x = "Time-point", 
       y = bquote("Effective library size<br>(million counts)")) +
  pl.theme() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_markdown())





per.tissue <- emmeans(lt.m1, specs = ~ "time|sets") %>%
  data.frame() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12"), labels = c("Week 0", "Week 2", "Week 12")), 
         sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
  
  ggplot(aes(time, exp(emmean), fill = sets )) +
  geom_errorbar(aes(ymin = exp(lower.CL), ymax = exp(upper.CL)), 
                position = position_dodge(width = 0.25), 
                width = 0.1, 
                size = line_size) +
  geom_point(position = position_dodge(width = 0.25), 
             shape = 21, 
             size = 2) +
  
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
  annotate("text", x = c(2,3), y = c(5.3, 5.8), label = c(lt.m1.stats[2,6], lt.m1.stats[3,6]), 
           size = 2.5) +
  
  
  labs(fill = "", 
       x = "Time-point", 
       y = "Effective library size<br>(million counts mg^-1^)") +
  pl.theme() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        axis.title.y = element_markdown())







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
         sig = if_else(fcthreshold == "s" & pthreshold == "s", "s", "ns")) %>%
  print()



##### Settings specifically for figure 3 ##################

## Color themes 
volcano.color.naive <- c("gray30", diff.colors[1])
volcano.color.ls <- c("gray30", diff.colors[2])
volcano.color.tissue <- c("gray30", diff.colors[3])


## Sizes volcano plots

volcano_size <- c(1, 1.4)
volcano_alpha <- c(0.2, 0.8)



############## Week 2 comparisons #####################################






week2_naive <- volcano_data %>%
  filter(model == "No normalization" & coef == "Week 2: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.naive) + 
  scale_fill_manual(values = volcano.color.naive) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none")  

week2_libsize <- volcano_data %>%
  filter(model == "Effective library size" & coef == "Week 2: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.ls) + 
  scale_fill_manual(values = volcano.color.ls) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week2_tissue <- volcano_data %>%
  filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 2: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.tissue) + 
  scale_fill_manual(values = volcano.color.tissue) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 
#### Venn diagrams for sets-effects ###############


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


names(venn_list_up) <- c("Naïve", "Effective library-size", "Tissue offset", "Naïve", "Effective library-size", "Tissue offset")
names(venn_list_down) <- c("Naïve", "Effective library-size", "Tissue offset", "Naïve", "Effective library-size", "Tissue offset")




week2_venn_up <- ggvenn(venn_list_up[1:3], 
                        show_elements = FALSE, 
                        value_type = "count",
                        fill_color = diff.colors, 
                        stroke_color = "gray90", 
                        set_name_size = 2, 
                        text_size = 2, 
                        stroke_size = 0.4) +
  labs(title = "Multiple-set > Single-set") + 
  theme(plot.title = element_text(size = 7))

week2_venn_down <- ggvenn(venn_list_down[1:3], 
                     show_elements = FALSE, 
                     value_type = "count",
                     fill_color = diff.colors, 
                     stroke_color = "gray90", 
                     set_name_size = 2, 
                     text_size = 2, 
                     stroke_size = 0.4) +
  labs(title = "Multiple-set < Single-set ")+ 
  theme(plot.title = element_text(size = 7))




week12_venn_up <- ggvenn(venn_list_up[4:6], 
                         show_elements = FALSE, 
                         value_type = "count",
                         fill_color = diff.colors, 
                         stroke_color = "gray90", 
                         set_name_size = 2, 
                         text_size = 2, 
                         stroke_size = 0.4) +
  labs(title = "Multiple-set > Single-set")+ 
  theme(plot.title = element_text(size = 7))

week12_venn_down <- ggvenn(venn_list_down[4:6], 
                           show_elements = FALSE, 
                           value_type = "count",
                           fill_color = diff.colors, 
                           stroke_color = "gray90", 
                           set_name_size = 2, 
                           text_size = 2, 
                           stroke_size = 0.4) +
  labs(title = "Multiple-set < Single-set ")+ 
  theme(plot.title = element_text(size = 7))





week2_grid <-  plot_grid(week2_venn_down, 
                         plot_grid(week2_naive, week2_libsize, week2_tissue, ncol = 1), 
            week2_venn_up, ncol = 3)




####### Week 12 comparisons ##################################

week12_naive <- volcano_data %>%
  filter(model == "No normalization" & coef == "Week 12: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.naive) + 
  scale_fill_manual(values = volcano.color.naive) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 

week12_libsize <- volcano_data %>%
  filter(model == "Effective library size" & coef == "Week 12: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.ls) + 
  scale_fill_manual(values = volcano.color.ls) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week12_tissue <- volcano_data %>%
  filter(model == "Per tissue mass +\nEffective library size" & coef == "Week 12: Multiple - Single-set") %>%
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.tissue) + 
  scale_fill_manual(values = volcano.color.tissue) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), limits = c(-4, 4), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 







week12_grid <- plot_grid(week12_venn_down, 
                         plot_grid(week12_naive, week12_libsize, week12_tissue, ncol = 1), 
                          week12_venn_up, ncol = 3)







############################ Week 2 time effect ###############

#### Venn diagrams for sets-effects ###############


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


names(venn_list_up) <- c("Naïve", "Effective library-size", "Tissue offset", "Naïve", "Effective library-size", "Tissue offset")
names(venn_list_down) <- c("Naïve", "Effective library-size", "Tissue offset", "Naïve", "Effective library-size", "Tissue offset")





timeweek2_venn_up <- ggvenn(venn_list_up[1:3], 
                        show_elements = FALSE, 
                        value_type = "count",
                        fill_color = diff.colors, 
                        stroke_color = "gray90", 
                        set_name_size = 2, 
                        text_size = 2, 
                        stroke_size = 0.4) +
  labs(title = "Week 2 > Week 0") + 
  theme(plot.title = element_text(size = 7))

timeweek2_venn_down <- ggvenn(venn_list_down[1:3], 
                          show_elements = FALSE, 
                          value_type = "count",
                          fill_color = diff.colors, 
                          stroke_color = "gray90", 
                          set_name_size = 2, 
                          text_size = 2, 
                          stroke_size = 0.4) +
  labs(title = "Week 2 < Week 0")+ 
  theme(plot.title = element_text(size = 7))




timeweek12_venn_up <- ggvenn(venn_list_up[4:6], 
                         show_elements = FALSE, 
                         value_type = "count",
                         fill_color = diff.colors, 
                         stroke_color = "gray90", 
                         set_name_size = 2, 
                         text_size = 2, 
                         stroke_size = 0.4) +
  labs(title = "Multiple-set > Single-set")+ 
  theme(plot.title = element_text(size = 7))

timeweek12_venn_down <- ggvenn(venn_list_down[4:6], 
                           show_elements = FALSE, 
                           value_type = "count",
                           fill_color = diff.colors, 
                           stroke_color = "gray90", 
                           set_name_size = 2, 
                           text_size = 2, 
                           stroke_size = 0.4) +
  labs(title = "Multiple-set < Single-set ")+ 
  theme(plot.title = element_text(size = 7))









volcano_time_data <- mm2_2 %>%
  group_by(coef, model) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns"), 
         sig = if_else(fcthreshold == "s" & pthreshold == "s", "s", "ns")) %>%
  print()


week2_time_naive <- volcano_time_data %>%
  filter(model == "naive" & coef == "timew2pre") %>% 
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.naive) + 
  scale_fill_manual(values = volcano.color.naive) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), limits = c(-2, 8), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 




week2_time_libsize <- volcano_time_data %>%
  filter(model == "lib_size_normalized" & coef == "timew2pre") %>% 
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.ls) + 
  scale_fill_manual(values = volcano.color.ls) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), limits = c(-2, 8), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week2_time_tissue <- volcano_time_data %>%
  filter(model == "tissue_offset_lib_size_normalized" & coef == "timew2pre") %>% 
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.tissue) + 
  scale_fill_manual(values = volcano.color.tissue) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), limits = c(-2, 8), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 








########### w12 time effect #############################

week12_time_naive <- volcano_time_data %>%
  
  filter(model == "naive" & coef == "timew12") %>% 
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.naive) + 
  scale_fill_manual(values = volcano.color.naive) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), limits = c(-2, 8), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 



week12_time_libsize <- volcano_time_data %>%
  filter(model == "lib_size_normalized" & coef == "timew12") %>% 
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.ls) + 
  scale_fill_manual(values = volcano.color.ls) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), limits = c(-2, 8), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 


week12_time_tissue <- volcano_time_data %>%
  filter(model == "tissue_offset_lib_size_normalized" & coef == "timew2pre") %>% 
  
  ggplot(aes(log2fc, -log10(p.val), fill = sig, alpha = sig, size = sig)) + 
  geom_point(shape = 21) +
  
  # Scales
  scale_color_manual(values = volcano.color.tissue) + 
  scale_fill_manual(values = volcano.color.tissue) +
  scale_size_manual(values = volcano_size) +
  scale_alpha_manual(values = volcano_alpha) +
  
  scale_x_continuous(breaks = c(-2, 0, 2, 4, 6, 8), limits = c(-2, 8), expand = c(0, 0)) +
  
  labs(x = "Log2 fold change", 
       y = "-Log10(P-value)") +
  
  pl.theme() +
  
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 7),
        # panel.border = element_rect(colour = "black"),
        legend.position = "none") 



timeweek2_grid <-  plot_grid(plot_grid(week2_time_naive, week2_time_libsize, week2_time_tissue, ncol = 1), 
                              plot_grid(timeweek2_venn_up, timeweek2_venn_down, ncol = 1))


timeweek12_grid <-  plot_grid(plot_grid(week12_time_naive, week12_time_libsize, week12_time_tissue, ncol = 1), 
                          plot_grid(timeweek12_venn_up, timeweek12_venn_down, ncol = 1))








############## Complete figure ###################



figure3 <- plot_grid(plot_grid(muscle_weight_fig, per.rna, per.tissue, ncol = 3, rel_widths = c(0.3, 0.3, 0.36), 
                               align = "h"), 
                     
                     week2_grid, 
                     
                     week12_grid,
                     
                     rel_heights = c(0.2, 0.4, 0.4),
                     
                     #timeweek2_grid, timeweek12_grid, 
                     ncol = 1)




#          draw_plot_label(label=c("A", "B"),
#                          x = c(0.02, 0.02), 
#                          y = c(0.97, 0.47),
#                          hjust=.5, vjust=.5, size = label.size)
          
          
          
          
          # Width of figure = 1x columns 8.9 cm
          # height of figure = full page = 23 cm
          
          
          ggsave("figures/figure3.pdf", plot = figure3, width = 8.9 * 2, height = 23, 
                 dpi = 600,
                 units = "cm", device=cairo_pdf)
          

# GO overlap alternatively gene-ovelap between methods?



### Time figures ###










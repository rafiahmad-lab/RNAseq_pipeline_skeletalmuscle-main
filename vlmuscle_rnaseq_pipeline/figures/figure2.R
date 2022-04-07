

### Figure 2. Quality control of RNA-seq experiment and check of assumptions: 
# 1. training induced changes
# 2. Constant cell/volume amount in cDNA synthesis
# 

source("./R/lib_fun.R")
source("./R/figure_source.R")


##### Quality control of RNA-seq experiment ######

######## Gene count density per method ##############



# Retrieve logCPM values for plotting
logCPM_data <- readRDS(file = "./data/derivedData/prel_analysis/logCPM_data.RDS")



# annotation data frame
ann_text <- data.frame(method = factor("hisat", levels = c("hisat", "star","salmon","kallisto", "rsem"), 
                                       labels = c("HISAT2", "STAR", "Salmon", "kallisto", "RSEM")), 
                       filtering = rep("pre", 4),
                       x1 = c(-1.8, 0), 
                       x2 = c(3, 3.1), 
                       y1 = c(0.5, 0.4),
                       y2 = c(0.65, 0.5), 
                       lab = rep(c("Pre-\nfiltering", "Post-\nfiltering"), 2))


# Density plot prior to filtering  
abundance_density_fig <- logCPM_data %>%
  mutate(method = factor(method, levels = c("hisat", "star","salmon","kallisto", "rsem"), 
                         labels = c("HISAT2", "STAR", "Salmon", "kallisto", "RSEM"))) %>%
  ggplot(aes(x = logCPM,
           #  y=..scaled..,
             fill = filtering, group = filtering)) + 
  geom_density(alpha = 0.2, stat = "density", size = line_size) + 
  scale_x_continuous(limits = c(-5, 15), expand = c(0,0), 
                     breaks = c(-5, 0, 5, 10, 15), 
                     labels = c("", 0, 5, 10, "")) +
  
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c(0, "", 0.4, "", 0.8, ""), 
                     limits = c(0, 2.5), 
                     expand = c(0,0)) +
  
  coord_cartesian(ylim=c(0,0.7)) +
  
  facet_grid(. ~ method) +
  
  # "Custom" legend pointing to data 
 geom_segment(data = ann_text, aes(x = x1, xend = x2, y = y1, yend = y2), size = line_size) + 
  geom_text(data = ann_text, aes(x = x2 + 0.05, y = y2, label = lab), size = 2, hjust = 0) +
  
  labs(y = "Gene count density", 
       x = "log2-Counts per million") +
  pl.theme() +
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = 'black', size = 7), 
        legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.y = element_blank())



######### Quality filtering figure #############

## Load data 

qfdata <- list()

qfdata$beforeqc <- read.delim("./data/qfdata/beforefilteringmqc_fastqc_per_base_sequence_quality_plot_1.txt", header = TRUE, row.names = 1)
qfdata$trimmomatc <- read.delim("./data/qfdata/trimmomaticmqc_fastqc_per_base_sequence_quality_plot_1.txt", header = TRUE, row.names = 1)
qfdata$trimgalore <- read.delim("./data/qfdata/trimgaloremqc_fastqc_per_base_sequence_quality_plot_1.txt", header = TRUE, row.names = 1)


# Combine long data sets and add method 
for(i in 1:length(qfdata)) {
  
  qfdata[[i]] <- qfdata[[i]] %>%
    mutate(method = names(qfdata)[i], 
           sample = rownames(.)) %>%
    pivot_longer(names_to = "bp", 
                 values_to = "quality", cols = X1:X150) %>%
    print()
  
  
}


qfdata <- bind_rows(qfdata) # Combining the list


# quality filtering figure
# Change from 10 th and 90 th percentile if needed.


qffig <- qfdata %>%
  mutate(bp = as.numeric(gsub("X", "", bp))) %>%
  group_by(method, bp) %>%
  summarise(m = mean(quality), 
            med = median(quality), 
            p10 = quantile(quality, prob = 0.1),  
            p90 = quantile(quality, prob = 0.9), 
            p25 = quantile(quality, prob = 0.25), 
            p75 = quantile(quality, prob = 0.75)) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = c("trimmomatc", "trimgalore", "beforeqc"), 
                         labels = c("Trimmomatic", 
                                    "Trim Galore", 
                                    "Prior to QF"))) %>%
  ggplot(aes(bp, med, group = method, fill = method)) + 
  geom_line(size = 0.4, aes(color = method)) + 
  geom_ribbon(aes(ymin = p10, ymax = p90, color = NULL), alpha = 0.2) + 
  
  pl.theme() + 
  scale_fill_manual(values = diff.colors) +
  scale_color_manual(values = diff.colors) +
  scale_y_continuous(breaks = c(26, 31, 36, 41), 
                     labels = c(26, 31, 36, 41), 
                     expand = c(0,0), 
                     limits = c(26, 41)) +
  labs(x = "Position in read (bp)", 
       y = "Quality score") +
  
  coord_cartesian(ylim=c(31,41)) +
  
  guides(color = guide_legend(override.aes = list(size = 0.2))) +
  theme(legend.position = c(0.45, 0.25), 
        legend.text = element_text(size = 7), 
        legend.key.size = unit(0.5,"line"), 
        axis.title.x = element_text(margin = margin(t =0 , r =0 , b =0 , l=0)))


qffig


######### Paired variation per method ###############


# Load data


paired_sd_data <- readRDS(file = "./data/derivedData/paired_sd_data.RDS")



## Plot results from analysis ##

sd_a <- paired_sd_data %>%
  filter(hk == TRUE) %>%
  ggplot(aes(A, SD, color = method, fill = method)) + 
  
  geom_smooth(se = FALSE, size = line_size) +
  pl.theme() + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5)) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20), 
                     labels = c(0, "", 10, "", 20)) + 
  coord_cartesian(ylim=c(0, 1.5),
                  xlim = c(-1, 15), 
                  expand = 0) +
  
  scale_color_manual(values = sequential.colors) +
  scale_fill_manual(values = sequential.colors) +
  theme(legend.position = "none") +
  labs(y = "Average log2-difference\nbetween replicates", 
       x = "Average log2-count (A)")





ncount_medium <- paired_sd_data %>%
  filter(hk == TRUE) %>%
  mutate(abundance.strata = if_else(A <= 1, "low", 
                                    if_else(A>1 & A <6, "medium", 
                                            if_else(A >= 6, "high", "undetermined")))) %>%
  mutate(abundance.strata = factor(abundance.strata, levels = c("low", "medium", "high"), 
                                   labels = c("Low abundace\n(A<1)",
                                              "Medium abundace\n(1 \u2264 A \u2265 6)",
                                              "High abundace\n(A > 6)"))) %>%
  group_by(method, abundance.strata) %>%
  filter(abundance.strata == "Medium abundace\n(1 \u2264 A \u2265 6)") %>%
  summarise(n = n()) %>%
  print()



ann_text <- data.frame(abundance.strata = factor( "Medium count\n(1 \u2264 A \u2265 6)", 
                                                  levels = c("Low abundace\n(A<1)",
                                                             "Medium count\n(1 \u2264 A \u2265 6)",
                                                             "High abundace\n(A > 6)")), 
                       
                       x = c(1, 2, 3, 4, 5), 
                       method = c("rsem", "kallisto", "salmon", "star", "hisat"),
                       n = ncount_medium$n + 30,
                       lab = c("RSEM", 
                                   "kallisto", 
                                   "Salmon", 
                                   "STAR", 
                                   "HISAT2"))



n_counts_fig <- paired_sd_data %>%
  filter(hk == TRUE) %>%
  mutate(abundance.strata = if_else(A <= 1, "low", 
                                    if_else(A>1 & A <6, "medium", 
                                            if_else(A >= 6, "high", "undetermined")))) %>%
  mutate(abundance.strata = factor(abundance.strata, levels = c("low", "medium", "high"), 
                                   labels = c("Low count\n(A<1)",
                                              "Medium count\n(1 \u2264 A \u2265 6)",
                                              "High count\n(A > 6)"))) %>%
  group_by(method, abundance.strata) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  ggplot(aes(method, n, fill = method)) + 
  geom_bar(stat = "identity") +
  
  pl.theme() +
  theme(strip.text.y = element_text(angle = 0, size = 7), 
        strip.background = element_blank(), 
        legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  scale_fill_manual(values = sequential.colors) +
  scale_y_continuous(breaks = c(0, 3000, 6000), 
                     labels = c(0, " ", 6000), 
                     limits = c(0, 6000), 
                     expand = c(0, 0), 
                     name = "n gene counts") +
    facet_grid(abundance.strata ~ .) + 
  geom_text(data = ann_text, aes(label = lab), hjust = 0, size = 2) +
  coord_flip()
  
############## RNA integrity scores and selection of samples ################

rna_integrity <- read_excel("./data/RNAintegrity.xlsx", sheet = 2)

temp <- data.frame(str_split_fixed(rna_integrity$`Sample Name`, " ", 3))

temp$subject <- paste("FP", temp$X1, sep="")
colnames(temp) <- c("subjnr", "leg", "timepoint", "subject")


RNA_integrity <- cbind(temp[, c(2:4)], rna_integrity[,c(3, 4, 5, 6, 7)])

colnames(RNA_integrity) <- c("leg", "timepoint", "subject","RNA_area",  "RNA_conc", "RNA_ratio", "RQI", "Classification")

rna_subjinfo <- read_excel("./data/RNAintegrity.xlsx", sheet = 1)

### Figure for insets into larger figure
dxa_sets_classification <- read_excel("./data/bodycomp_DXA.xlsx") %>%
  dplyr::select(subject, timepoint, L = lean.left_leg, R = lean.right_leg) %>%
  pivot_longer(cols = L:R, names_to = "leg", values_to = "lean.mass") %>%
  inner_join(read_csv2("./data/oneThreeSetLeg.csv") %>%
               gather("sets", "leg", multiple:single)) %>%
  filter(include == "incl") %>%
  spread(timepoint, lean.mass) %>%
  mutate(percentage.change = ((post/pre)-1)*100, 
         rnaseq_include = factor(rnaseq_include, levels = c("incl", "excl"))) %>%
  dplyr::select(subject:sets, percentage.change, -leg) %>%
  spread(sets, percentage.change) %>%
  
  
  ggplot(aes(single, multiple, fill = rnaseq_include)) + 
  scale_fill_manual(values = diff.colors) +
  
  geom_abline(slope = 1, intercept = 0, color = "gray40", lty = 2, alpha = 0.2) +
  
  geom_point(shape = 21, size = 1.5) +
  
  scale_x_continuous(breaks = c(-5, 0, 5, 10, 15), 
                     labels = c("", 0, "", 10, ""), 
                     limits = c(-10, 15), 
                     expand = c(0, 0)) +
  
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15), 
                     labels = c("", 0, "", 10, ""), 
                     limits = c(-10, 15), 
                     expand = c(0, 0)) +
  
  annotate("text", x = 5, y = -6, label = "Single-set\nLBM %-change", size = 2, hjust = 0) +
  annotate("text", x = -9, y = 13, label = "Multiple-set\nLBM %-change", size = 2, hjust = 0) +
  pl.theme() +
  theme(legend.position = "none", 
        axis.title = element_blank(), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        axis.ticks.length = unit(-0.075, "cm"),
        axis.text.x = element_text(vjust = 0.2), 
        axis.text.y = element_text(hjust = 1, margin = margin(t = 0, r = 4, b = 0, l = 0, unit = "pt"))) 
  
## Figure for side-by-side plotting

dxa_sets_classification2 <- read_excel("./data/bodycomp_DXA.xlsx") %>%
  dplyr::select(subject, timepoint, L = lean.left_leg, R = lean.right_leg) %>%
  pivot_longer(cols = L:R, names_to = "leg", values_to = "lean.mass") %>%
  inner_join(read_csv2("./data/oneThreeSetLeg.csv") %>%
               gather("sets", "leg", multiple:single)) %>%
  filter(include == "incl") %>%
  spread(timepoint, lean.mass) %>%
  mutate(percentage.change = ((post/pre)-1)*100, 
         rnaseq_include = factor(rnaseq_include, levels = c("incl", "excl"))) %>%
  dplyr::select(subject:sets, percentage.change, -leg) %>%
  spread(sets, percentage.change) %>%
  
  
  ggplot(aes(single, multiple, fill = rnaseq_include)) + 
  scale_fill_manual(values = diff.colors) +
  
  geom_abline(slope = 1, intercept = 0, color = "gray40", lty = 2, alpha = 0.2) +
  
  geom_point(shape = 21, size = 1.5) +
  
  scale_x_continuous(breaks = c(-5, 0, 5, 10, 15), 
                     labels = c("", 0, "", 10, ""), 
                     limits = c(-10, 15), 
                     expand = c(0, 0)) +
  
  scale_y_continuous(breaks = c(-5, 0, 5, 10, 15), 
                     labels = c("", 0, "", 10, ""), 
                     limits = c(-10, 15), 
                     expand = c(0, 0)) +
  
  labs(x = "Single-set LBM %-change", 
       y = "Multiple-set LBM %-change") +
  

  pl.theme() +
  theme(legend.position = "none", 
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) 





RNA <- RNA_integrity%>%
  inner_join(rna_subjinfo[ , c(1,3,4,5,6,7,8, 9)])%>%
  mutate(RQI = as.numeric(RQI))

rna_integrity_fig <- RNA %>% 
  group_by(subject) %>%
  mutate(min.rna.quality = min(RQI), 
         score = if_else(RQI > 7, "OK", "notOK")) %>%
  mutate(score = factor(score, levels = c("OK", "notOK")), 
         rna.yield = (elution.volume * conc)/prot.mrna1) %>%

  
  ggplot(aes(rna.yield, RQI, fill = score)) + 
  
  geom_hline(yintercept = 7, lty = 2, color = "gray40") +
  
  geom_point(shape=21, size = 1.5, stroke = 0.2) +
  scale_fill_manual(values = diff.colors)+
  xlab(bquote('RNA yield '~(ng~mg^-1))) + 
  scale_y_continuous(breaks = c(0, 2.5, 5, 7.5, 10), 
                     labels = c(0, "", 5, "", 10), 
                     limits = c(0, 10), 
                     expand = c(0, 0)) + 
  scale_x_continuous(breaks = c(100, 200, 300, 400, 500, 600, 700, 800), 
                     labels = c("", 200, "",400, "",600,"", 800),
                     limits = c(100, 800), 
                     expand = c(0, 0)) +
  pl.theme() + 
  theme(legend.position = "none", 
        plot.margin = margin(t = 0, r = 7, b = 0, l = 7))


# rna_integrity_fig <- ggdraw(rna_integrity_fig) + 
#   draw_plot(dxa_sets_classification, 0.55, 0.20, 0.42, 0.45)

####### Training induced changes ############

### DXA results 

dxa <- read_excel("./data/bodycomp_DXA.xlsx") %>%
  dplyr::select(subject, timepoint, L = lean.left_leg, R = lean.right_leg) %>%
  pivot_longer(cols = L:R, names_to = "leg", values_to = "lean.mass") %>%
  inner_join(read_csv2("./data/oneThreeSetLeg.csv") %>%
               gather("sets", "leg", multiple:single)) %>%
  filter(rnaseq_include == "incl") %>%
  print()


dxa_m1 <- dxa %>%
  pivot_wider(names_from = timepoint, values_from = lean.mass) %>%
  mutate(change = post-pre) %>%
  group_by(sex) %>%
  mutate(pre = pre - mean(pre)) %>%
  lme(change ~  sex + pre + sets, random = list(subject = ~ 1), data =.) 


dxa_pval <- publR::pval(data.frame(broom::tidy(dxa_m1, effects = "fixed"))[4, 5], 
                        flag = TRUE)




dxa_results_fig <- emmeans(dxa_m1, specs = ~ "sets") %>%
  data.frame() %>%
  mutate(sets = factor(sets, levels = c("single", "multiple"), labels = c("Single-set", "Multiple-set"))) %>%
  ggplot(aes(sets, emmean, fill = sets)) + 
  geom_hline(yintercept = 0, size = line_size, color = "gray40", lty = 2) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) + 
  geom_point(shape = 21, size = 3) + 
  scale_fill_manual(values = c("#bdbdbd", "#636363")) +
  scale_y_continuous(limits = c(-100, 600), 
                     breaks = c(-100, 0, 100, 200, 300, 400, 500, 600), 
                     labels = c("", 0, "", 200, "", 400, "", 600),
                     expand = c(0, 0)) +
  # Segments for p-value 
 # geom_segment(data = data.frame(x1 = 1, x2 = 2, y1 = 500, y2 = 500), 
 #              aes(x = x1, xend = x2, y = y1, yend = y2), 
 #              inherit.aes = FALSE) +
  annotate("text", x = 1.5, y = 550, label = dxa_pval, 
           size = 5) +
  
  
  labs(fill = "", 
       y = "Change in lower \n extremity lean-mass (g)") +
  pl.theme() +
  theme(axis.title.x = element_blank(), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank())








##### FIGURE 2 COMPOSITION #####



  
  figure2 <- plot_grid(
    
    # Row 1 
    plot_grid(rna_integrity_fig, dxa_sets_classification2, rel_widths = c(0.5, 0.5), ncol = 2, align = "h"), 
    
    # Row 2
    plot_grid(dxa_results_fig, 
                                 NULL,
                                 qffig, 
                          
                                 rel_widths = c(0.35,0.05, 0.60), 
                                 align = "h", ncol = 3), 
                       
    # Row 3
                        abundance_density_fig, 
                     
    # Row 4
                       plot_grid(
                         # Row 4, column 1
                         plot_grid(sd_a, n_counts_fig, rel_heights = c(0.5, 0.5), ncol = 1),
                       
                         NULL, # adding white space
                         
                       # Total library sizes and muscle weight
                       # Row 4, column 2
                         plot_grid(muscle_weight_fig, per.rna, per.tissue, ncol = 1, rel_heights = c(0.3, 0.3, 0.36), 
                                 align = "v"), ncol = 3, rel_widths = c(0.6, 0.05, 0.4)), 
                       
                       
                       
                       
                      
                       ncol = 1, rel_heights = c(0.2, 0.2, 0.2, 0.4)) +
    
    draw_plot_label(label=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
                    x = c(0.02, 0.50, 0.02, 0.45, 0.02, 0.02, 0.02, 0.62, 0.62, 0.62), 
                    y = c(0.99, 0.99, 0.8,  0.8,  0.58,  0.4,  0.20, 0.4, 0.25, 0.15),
                    hjust=.5, vjust=.5, size = label.size)
  
  
  

# Width of figure = 1x columns 8.9 cm
# height of figure = full page = 23 cm


ggsave("figures/figure2.pdf", plot = figure2, width = 8.9, height = 23, 
       dpi = 600,
       units = "cm", device=cairo_pdf)












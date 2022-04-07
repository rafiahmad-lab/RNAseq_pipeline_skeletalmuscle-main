
#### Suppl Figure 2 ##############

# External reference gene
# RNA to muscle tissue weight



source("./R/lib_fun.R")
source("./R/figure_source.R")





# Load leg-condition data
legs_condition <- read.csv("./data/oneThreeSetLeg.csv", sep = ";") %>%
  pivot_longer(names_to = "sets", values_to = "leg", multiple:single) %>%
  print()






### Read lambda data
counts_vs_cq <- read_delim("./data/lambda_abundance.txt", delim = ",") %>%

  
  pivot_longer(names_to = "sample", 
               values_to = "count", 
               cols =  "100-FP23w12L":"9-FP13w0L" ) %>%
  filter(X1 %in% c("gene_35", "gene_41")) %>%
  rowwise() %>%
  mutate(sample = as.character(sample), 
         subject = str_split_fixed(sample, "w", 2)[1], 
         leg = substrRight(sample, 1), 
         time = gsub(subject, "", sample), 
         time = gsub(leg, "", time)) %>%
  separate(subject, into = c("nr", "subject")) %>%
  dplyr::select(- nr) %>%
  inner_join(legs_condition) %>%
  dplyr::select(sample, subject,  leg, timepoint = time, sets, gene = X1, count) %>%
  group_by(subject, leg, timepoint, gene) %>%
  summarise(count = mean(count)) %>%

  
  inner_join(read_feather("./data/qpcrdat_replicates") %>%
               filter(target == "Lambda Kit") %>%
               group_by(subject, timepoint, leg) %>%
               summarise(cq = mean(cq, na.rm = TRUE))) %>%

  ggplot(aes(log(2^-cq), log(count), color = gene)) + geom_point() + 
  
  geom_smooth(method = "lm") + 
  
  labs(x = "log(2^-cq) [qPCR]", 
       y = "log(Counts) [RNA-seq]", 
       color = "Gene") + 
  
pl.theme() +
  theme(legend.position = c(0.9, 0.2))
  


### Alternative counts vs. cq where counts are averaged

### Read lambda data
counts_vs_cq2 <- read_delim("./data/lambda_abundance.txt", delim = ",") %>%
  
  
  pivot_longer(names_to = "sample", 
               values_to = "count", 
               cols =  "100-FP23w12L":"9-FP13w0L" ) %>%
  filter(X1 %in% c("gene_35", "gene_41")) %>%
  rowwise() %>%
  mutate(sample = as.character(sample), 
         subject = str_split_fixed(sample, "w", 2)[1], 
         leg = substrRight(sample, 1), 
         time = gsub(subject, "", sample), 
         time = gsub(leg, "", time)) %>%
  separate(subject, into = c("nr", "subject")) %>%
  dplyr::select(- nr) %>%
  inner_join(legs_condition) %>%
  dplyr::select(sample, subject,  leg, timepoint = time, sets, gene = X1, count) %>%
  group_by(subject, leg, timepoint, gene) %>%
  summarise(count = mean(count)) %>%
  
  
  inner_join(read_feather("./data/qpcrdat_replicates") %>%
               filter(target == "Lambda Kit") %>%
               group_by(subject, timepoint, leg) %>%
               summarise(cq = mean(cq, na.rm = TRUE))) %>%
  
  filter(timepoint != "w2post") %>%

  
  ggplot(aes(log(2^-cq), log(count), fill = gene)) + 
  
  geom_point(shape = 21, size = 0.9, alpha = 0.5) + 
  scale_shape_manual(values = c(21, 22)) +
  geom_smooth(method = "lm", color = "gray40") + 
  
  labs(x = "log(2^-cq) [qPCR]", 
       y = "log(Counts) [RNA-seq]", 
       color = "Gene", 
       title = "qPCR relative abundance vs. RNA-seq counts") + 
  
  pl.theme() +
  theme(legend.position = c(0.85, 0.2), 
        plot.title = element_text(size = 10))



counts_vs_cq2

###### Muscle weight vs counts ##############


### Muscle weight data 
mw <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, leg, weight = prot.mrna1, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  print()



counts_vs_weight <- read_delim("./data/lambda_abundance.txt", delim = ",") %>%
  
  
  pivot_longer(names_to = "sample", 
               values_to = "count", 
               cols =  "100-FP23w12L":"9-FP13w0L" ) %>%
  filter(X1 %in% c("gene_35", "gene_41")) %>%
  rowwise() %>%
  mutate(sample = as.character(sample), 
         subject = str_split_fixed(sample, "w", 2)[1], 
         leg = substrRight(sample, 1), 
         time = gsub(subject, "", sample), 
         time = gsub(leg, "", time)) %>%
  separate(subject, into = c("nr", "subject")) %>%
  dplyr::select(- nr) %>%
  inner_join(legs_condition) %>%
  dplyr::select(sample, subject,  leg,  time, sets, gene = X1, count) %>%
  inner_join(mw) %>%
  
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12")))  %>%
  
  ggplot(aes(log(weight), log(count), fill = gene)) + 
  geom_point(shape = 21, size = 0.9, alpha = 0.5) + 
  
  geom_smooth(method = "lm", color = "gray40") + 
  
  labs(x = "log(muscle weight) ", 
       y = "log(Counts) [RNA-seq]", 
       color = "Gene", 
       title = "Muscle weight vs. RNA-seq counts") + 
  
  pl.theme() +
  theme(legend.position = c(0.2, 0.2), 
        plot.title = element_text(size = 10))
  
  

figureS2 <- plot_grid(counts_vs_cq2, counts_vs_weight, 
                      ncol = 1, rel_heights = c(0.5, 0.5))


ggsave("figures/figureS2.pdf", plot = figureS2, width = 8.9, height = 16, 
       dpi = 600,
       units = "cm", device=cairo_pdf)
|


  


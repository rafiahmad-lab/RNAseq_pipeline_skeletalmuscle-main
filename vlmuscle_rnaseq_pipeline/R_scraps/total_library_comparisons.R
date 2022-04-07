


dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")

mgcv_results <- bind_rows(readRDS(file = "./data/derivedData/DE/mgcv_results.RDS"))
estimates_mixed <- readRDS(file = "./data/derivedData/DE/glmmTMB_results.RDS")

# get genes information 

genes <- list(data.frame(gene = rownames(dge_lists[[1]])),
              data.frame(gene = rownames(dge_lists[[2]])),
              data.frame(gene = rownames(dge_lists[[3]])),
              data.frame(gene = rownames(dge_lists[[4]])),
              data.frame(gene = rownames(dge_lists[[5]]))) 
                         

genes <- bind_rows(genes)


cbind(genes, mgcv_results, data.frame(model = "mgcv")) %>%
  rbind(cbind(genes, mgcv_results, data.frame(model = "glmmTMB"))) %>%
  dplyr::select(model, gene, method, coef, estimate, p.val) %>%
  filter(method == "kallisto", 
         !(coef %in% c("(Intercept)", "nf", "timew2pre", "timew2post", "timew12"))) %>%
  group_by(model, coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%

  ggplot(aes(estimate, -log10(p.val), color = fcthreshold, fill = pthreshold)) + 
  geom_point(shape = 21, alpha = 0.4) +
  scale_color_manual(values = c("blue", "red")) + 
  scale_fill_manual(values = c("blue", "red")) +
  facet_grid(model ~ coef)
  


?glmmTMB

mgcv_results %>%
  print()
?gam

estimates_mixed %>%
  print()



##### Average counts per milligram tissue ###############



dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")

kallisto <- dge_lists[[2]]

mw <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%

  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  print()


samples <- kallisto$samples %>%
  mutate(sample = rownames(.), 
         eff.lib = lib.size * norm.factors) %>%
  inner_join(mw) %>%
  dplyr::select(sample, subject, sex, time, sets, eff.lib, tissue) %>%
  print()



total.counts <- kallisto$counts %>%
  data.frame(gene = rownames(.)) %>%
  
  pivot_longer(names_to = "sample", values_to = "count", cols = FP11w0L:FP9w2preR) %>%

  inner_join(samples) %>%
  mutate(count = as.integer(round(count, 0))) %>%
  print()


tc_subset <- total.counts %>%
  filter(gene %in% sample(unique(gene), 100)) %>%
  data.frame() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w12")), 
         sets = factor(sets, levels = c("single", "multiple")), 
         subject = factor(subject), 
         sample = factor(sample)) %>%
  print()




### Total model ! ####
 m1 <- glmer.nb(count ~ time + time:sets + (1|subject) + (1|sample), data = tc_subset, 
                offset = log(tissue), 
                control = glmerControl(calc.derivs = FALSE))



m1  <- gam(count ~  time * sets + s(subject, bs = "re")+ s(sample, bs = "re"), 
    offset = log(tissue),
    data = tc_subset,
    family = nb, method = "ML")

summary(m1)

exp(0.37)


m1  <- gam(eff.lib ~  time * sets + s(subject, bs = "re")+ s(sample, bs = "re"), 
           offset = log(tissue),
           data = tc_subset,
           family = nb, method = "ML")
?lme4





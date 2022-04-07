








# read DGEList

dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")



## Download gene family information from HUGO

# Myosin heavy chain gene family 

url <- "https://www.genenames.org/cgi-bin/genegroup/download?id=1098&type=node"


myosin_family <- read_delim(url, delim = "\t") %>%
  dplyr::select(symbol = 'Approved symbol', gene = 'Ensembl gene ID') %>%
  print()



results <- list()

for(i in 1:length(dge_lists)) {
  

  
  temp  <- dge_lists[[i]][rownames(dge_lists[[i]]) %in% myosin_family$gene,]
  
myhc_all  <- temp$counts %>%
    data.frame(gene = rownames(.)) %>%
    gather(sample, count, FP11w0L:FP9w2preR) %>%
    inner_join(temp$samples %>%
                 mutate(sample = rownames(.))) %>%
    inner_join(myosin_family) %>%
  filter(symbol %in% c("MYH1", "MYH2", "MYH7")) %>%
    dplyr::select(subject, time, sets, gene, symbol, count) %>%
    group_by(subject, time, sets) %>%
    mutate(sum_count = sum(count)) %>%
    ungroup() %>%
    mutate(relative_count = (count/sum_count)*100, 
           method = names(dge_lists[i])) 
    
  
results[[i]] <- myhc_all


}


bind_rows(results) %>%
  dplyr::select(subject, time, sets, symbol, method, relative_count) %>%
  pivot_wider(names_from = method, values_from = relative_count) %>%
  filter(symbol %in% c("MYH1", "MYH2", "MYH7")) %>%
  ggplot(aes(kallisto, star, color = symbol)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)  +
  geom_smooth(method = "lm")
  print()




myhc_all %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w2post", "w12")),
         sets = factor(sets, levels = c("single", "multiple")),
         symbol = fct_relevel(symbol, "MYH7", "MYH2", "MYH1")) %>%
  ggplot(aes(time, percentage.median, fill = sets)) + 
  geom_point(shape = 21) + 
  facet_grid(.~ symbol)




temp$counts %>%
  data.frame(gene = rownames(.)) %>%
  gather(sample, count, FP11w0L:FP9w2preR) %>%
  inner_join(temp$samples %>%
               mutate(sample = rownames(.))) %>%
  inner_join(myosin_family) %>%
  dplyr::select(subject, time, sets, gene, symbol, count) %>%
  group_by(subject, time, sets) %>%
  mutate(sum_count = sum(count)) %>%
  ungroup() %>%
  mutate(relative_count = (count/sum_count)*100) %>%
  ungroup() %>%
  mutate(time = factor(time, levels = c("w0", "w2pre", "w2post", "w12")),
         sets = factor(sets, levels = c("single", "multiple")),
         symbol = fct_relevel(symbol, "MYH7", "MYH2", "MYH1")) %>%
  
  
  
  ggplot(aes(time, relative_count, fill = symbol)) + 
  geom_bar(stat = "identity") + 
  facet_grid(sets~subject)












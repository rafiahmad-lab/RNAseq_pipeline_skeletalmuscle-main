
########## Go analysis of DE genes #################

# DE genes have been identified in three models (results from mixed_models2.R)
# This script performs filtering of genes (by fold-change and p-value
# criteria), reduces go terms and compile summary data frames.


# Load libraries
source("./R/lib_fun.R")


# Load data 

## Interaction model 
mm_interaction <- readRDS("./data/derivedData/DE/mixedmodel2_results1.RDS")

## Time model
mm_time <- readRDS("./data/derivedData/DE/mixedmodel2_results2.RDS")


## Combine data for covenint handling ##

sig.genes <- mm_interaction %>%
  dplyr::select(gene, method, model, coef, estimate, p.val) %>%
  mutate(interaction = TRUE) %>%
  rbind(mm_time %>%
          dplyr::select(gene, method, model, coef, estimate, p.val) %>%
          mutate(interaction = FALSE)) %>% 
  group_by(model, interaction, coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"),  # Uses coef-wise fdr
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.5, "s", "ns")) %>% # ~ 1.4 fold change
  
  # Select coefficients of interest
  
  filter((coef %in% c("timew2pre", "timew12") & interaction == FALSE) ||
           coef %in% c("timew0:setsmultiple",  
                       "timew2pre:setsmultiple", 
                       "time12:setsmultiple")) %>%
  
  # Filter 'significant' genes 
  
  filter(adj.p < 0.05 & fcthreshold == "s") %>%
  
  print()


sig.genes %>%
 filter(      coef == "timew12:setsmultiple") %>%
  print()


# Extract gene sets 



models <- c("tissue_offset_lib_size_normalized", 
            "lib_size_normalized", 
            "naive")

coefs <- c("timew2pre", 
           "timew12", 
           "timew2pre:setsmultiple", 
           "timew12:setsmultiple")


# Backrgound gene vector

background <- bitr(unique(mm_interaction$gene), fromType = "ENSEMBL",
                   toType = c( "ENTREZID", "SYMBOL"),
                   OrgDb = org.Hs.eg.db)


# Storage 


i <- 2
j <- 3

# All gene ontology terms (reduced)
go_terms <- list()

#


for(i in 1:3) {
  
  temp <- sig.genes %>%
    filter(model == models[i])
   
  
  model_go <- list()  
  
  for(j in 1:4) {
    
    goi_up <- temp %>%
      filter(coef == coefs[j]) %>%
      filter(estimate > 0) %>%
      ungroup() %>%
      dplyr::select(gene) %>%
      data.frame()
    
    goi_down <- temp %>%
      filter(coef == coefs[j]) %>%
      filter(estimate < 0) %>%
      ungroup() %>%
      dplyr::select(gene) %>%
      data.frame()
    
    
    if(nrow(goi_down) > 5) {
      
      entrezid_down <- bitr(goi_down[,1], fromType = "ENSEMBL",
                      toType = c( "ENTREZID", "SYMBOL"),
                      OrgDb = org.Hs.eg.db)
      
      
      bp_down <- enrichGO(entrezid_down[,2], 
                        OrgDb = 'org.Hs.eg.db', 
                        ont = "BP", 
                        universe = background[,2], 
                        readable = TRUE)
    
      
      cc_down <-   enrichGO(entrezid_down[,2], 
                          OrgDb = 'org.Hs.eg.db', 
                          ont = "CC", 
                          universe = background[,2], 
                          readable = TRUE)
      
      mf_down <-  enrichGO(entrezid_down[,2], 
                         OrgDb = 'org.Hs.eg.db', 
                         ont = "MF", 
                         universe = background[,2], 
                         readable = TRUE)
      
      
      bp_down2 <- simplify(bp_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
      cc_down2 <- simplify(cc_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
      mf_down2 <- simplify(mf_down, cutoff = 0.7, by = "p.adjust", select_fun = min)
      
      
      simpl.df.down <- data.frame(bp_down2) %>%
        mutate(ont = "BP") %>%
        rbind(data.frame(cc_down2) %>%
                mutate(ont = "CC")) %>%
        rbind(data.frame(mf_down2) %>%
                mutate(ont = "MF")) %>%
        mutate(change = "down", 
               model = models[i], 
               coef = coefs[j])
      
      
      
      
    } else {
      
      simpl.df.down <- data.frame(ID = character(), 
                 Description = character(),
                 GeneRatio = character(),
                 BgRatio = character(),
                 pvalue = numeric(), 
                 p.adjust = numeric(), 
                 qvalue = numeric(), 
                 geneID = character(), 
                 Count = integer(), 
                 ont = character(), 
                 change = character(),
                 model = character(),
                 coef = character())

      
      
      
    } 
      
    
    
    if(nrow(goi_up) > 5) {
      
      entrezid_up <- bitr(goi_up[,1], fromType = "ENSEMBL",
                            toType = c( "ENTREZID", "SYMBOL"),
                            OrgDb = org.Hs.eg.db)
      
      bp_up <- enrichGO(entrezid_up[,2], 
                        OrgDb = 'org.Hs.eg.db', 
                        ont = "BP", 
                        universe = background[,2], 
                        readable = TRUE)
      
      
      cc_up <-   enrichGO(entrezid_up[,2], 
                                  OrgDb = 'org.Hs.eg.db', 
                                  ont = "CC", 
                                  universe = background[,2], 
                                  readable = TRUE)
      
      mf_up <-  enrichGO(entrezid_up[,2], 
                                  OrgDb = 'org.Hs.eg.db', 
                                  ont = "MF", 
                                  universe = background[,2], 
                                  readable = TRUE)
      

      bp_up2 <- simplify(bp_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
      cc_up2 <- simplify(cc_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
      mf_up2 <- simplify(mf_up, cutoff = 0.7, by = "p.adjust", select_fun = min)
      
      
      simpl.df.up <- data.frame(bp_up2) %>%
        mutate(ont = "BP") %>%
        rbind(data.frame(cc_up2) %>%
                mutate(ont = "CC")) %>%
        rbind(data.frame(mf_up2) %>%
                mutate(ont = "MF")) %>%
        mutate(change = "up", 
               model = models[i], 
               coef = coefs[j])
        


      
    } else {
      
      simpl.df.up <- data.frame(ID = character(), 
                                Description = character(),
                                GeneRatio = character(),
                                BgRatio = character(),
                                pvalue = numeric(), 
                                p.adjust = numeric(), 
                                qvalue = numeric(), 
                                geneID = character(), 
                                Count = integer(), 
                                ont = character(), 
                                change = character(),
                                model = character(),
                                coef = character())

    } 
    
    model_go[[j]] <- rbind(simpl.df.up, simpl.df.down)

  }

  go_terms[[i]] <- bind_rows(model_go)

}







go_terms <- bind_rows(go_terms)

colnames(go_terms)

models <- c("tissue_offset_lib_size_normalized", 
            "lib_size_normalized", 
            "naive")


geneRatio.plots <- list()

for(i in 1:3){
 
  
  goid <- go_terms %>%
    filter(coef %in% c("timew2pre", "timew12")) %>%
    filter(model == models[1], 
           ont == "BP") %>%
    group_by(coef) %>%
    top_n(10, -log10(p.adjust)) %>%
    ungroup() %>%
    separate(GeneRatio, into = c("k", "n"), sep = "/", convert = TRUE) %>%
    mutate(geneRatio = k/n, 
           Description = factor(as.character(Description))) %>% 
    mutate(Description = fct_reorder(Description, geneRatio)) %>%
    print()
  
  

  
  
go_terms %>%
  filter(coef %in% c("timew2pre", "timew12")) %>%
  filter(ID %in% goid$ID) %>%
  filter(ont == "BP") %>%
    separate(GeneRatio, into = c("k", "n"), sep = "/", convert = TRUE) %>%
    ungroup() %>%
    mutate(geneRatio = k/n, 
           Description = factor(as.character(Description), 
                                levels =   levels(goid$Description))) %>% 
    mutate(coef = factor(coef, levels = c("timew2pre", "timew12"), 
                         labels = c("Week 2 - Week 0", 
                                    "Week 12 - Week 0"))) %>%
    
    
    ggplot(aes(Description, geneRatio)) + 
    geom_point(aes(size = Count, fill = -log10(p.adjust), shape = model), 
               position = position_dodge(width = 0.2)) + 
  scale_shape_manual(values =  c(21, 23, 24)) +
    facet_grid(. ~ coef) + 
    labs(size = "N Genes in category", 
         fill = "-Log10(Adjusted p-value)", 
         y = "Gene ratio") +
    coord_flip() + 
    pl.theme() +
    theme(strip.background = element_blank(),
          strip.text = element_text(size = 7),
          axis.title.y = element_blank(), 
          legend.position = "bottom")
  
    
    
  names(geneRatio.plots)[i] <- models[i]


}


plot_grid(geneRatio.plots$tissue_offset_lib_size_normalized, 
          geneRatio.plots$lib_size_normalized, 
          geneRatio.plots$naive, ncol = 3)











go_unique <- go_terms %>%
  filter(coef %in% c("timew2pre", "timew12")) %>%
  dplyr::select(-geneID) %>%
  separate(GeneRatio, into = c("k", "n"), sep = "/", convert = TRUE) %>%
  mutate(geneRatio = k/n) %>% 
  
  filter(Count  > 10) %>%
  
  group_by(coef, model, change) %>%
  slice(tail(row_number(), 10)) %>%
  ungroup() %>%

  dplyr::select(ID) %>%
  print()
  
go_terms %>%
  separate(GeneRatio, into = c("k", "n"), sep = "/", convert = TRUE) %>%
  mutate(geneRatio = k/n) %>% 
  filter(ID %in% unique(go_unique$ID)) %>%
    mutate(ID = fct_reorder(ID, geneRatio)) %>%
  
  
  ggplot(aes(ID, geneRatio)) + 
  geom_point(aes(size = Count, fill = p.adjust, shape = model), 
             position = position_dodge(width = 0.2)) + 
  scale_shape_manual(values = c(21, 23, 24)) +
  facet_grid(change ~ coef) + 
  coord_flip()
  
  
  print()
  
  










selection <- entrezid
background <- as.character(mapIds(org.Hs.eg.db, keys = unique(de_interact$gene), 
                                  keytype = "ENSEMBL", 
                                  column = 'ENTREZID'))







# GOSemSim clustering of genes

# Gene similarity
hsGO_bp <- godata('org.Hs.eg.db', ont="BP")


sim <- mgeneSim(genes = selection, semData = hsGO, 
                measure = "Wang", verbose = FALSE)

# Extract clustering information for plotting using ggdendrogram

model <- hclust(as.dist(1-sim), "ward.D2")
dhc <- as.dendrogram(model)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")


dendrogram <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + 
  geom_text(data = ddata$labels, aes(x = x, y = y - 0, label = label), hjust = 0)


### Go terms



# Term similarity






go <- goana(de = selection)


# Go terms 
setOntology("BP")
temp <- getGOInfo(selection)
str(temp)

View(getGOInfo)



go_terms <- data.frame(got) %>%
  mutate(go_id = rownames(.)) %>%
  mutate(p.adj = p.adjust(P.DE)) %>%
  filter(p.adj < 0.05) %>%
  arrange(p.adj) %>%
  print()


termSim <- getTermSim(go_terms[,6])

model <- hclust(as.dist(1-termSim), "ward.D2")
dhc <- as.dendrogram(model)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")


dendrogram <- ggplot(segment(ddata)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0.2, 0)) + 
  geom_text(data = ddata$labels, aes(x = x, y = y - 0, label = label), hjust = 0)





library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'entrezgene_id', 'go_id'),
                   filters = "ensembl_gene_id", values = goi[,1], mart = ensembl)




goid_annotate <-  gene.data %>%


  dplyr::select(-hgnc_symbol) %>%
  table() %>%
  data.frame() %>%

  spread(key = go_id, value = Freq) %>%
  dplyr::select(-V1) %>%
  gather(go_id, Freq, "GO:0000139":"GO:2001202") %>%

 
  mutate(entrezgene_id = factor(entrezgene_id, levels = ddata$labels$label)) %>% 
#          go_id = factor(go_id, levels = go_terms$go_id)) %>%
  ggplot(aes(go_id, entrezgene_id, fill = Freq)) + geom_tile() + 

  theme(axis.text.x = element_text(angle = 90)) 
  


plot_grid(dendrogram, goid_annotate, ncol = 2)







topGO(go, ontology = "BP")






?compute_SS_distances


keg <- kegga(entrezid)  

topKEGG(keg)


bind_rows(mgcv_results5) %>%
  filter(coef %in% c("sexmale",  "timew2pre:sexmale", "timew12:sexmale")) %>%
  mutate(interact.sig = if_else(gene %in% sex_interaction_genes[,1], "s", "ns")) %>%
  
  group_by(coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
  ggplot(aes(log2fc, -log10(p.val), color = fcthreshold, fill = interact.sig)) + 
  geom_point(shape = 21, alpha = 0.4) +
  scale_color_manual(values = c("blue", "red")) + 
  scale_fill_manual(values = c("blue", "red")) +
  facet_grid(. ~ coef)


filter(coef == "(Intercept)") %>%
  
  
  
  filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                     "timew12:setsmultiple")) %>%
  mutate(coef = factor(coef, levels = c("timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple"))) %>%
  
  group_by(model, coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
  
  ggplot(aes(log2fc, -log10(p.val), color = fcthreshold, fill = pthreshold)) + 
  geom_point(shape = 21, alpha = 0.4) +
  scale_color_manual(values = c("blue", "red")) + 
  scale_fill_manual(values = c("blue", "red")) +
  facet_grid(model ~ coef)


de_genes <- bind_rows(mgcv_results4) %>%
  filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                     "timew12:setsmultiple")) %>%
  mutate(coef = factor(coef, levels = c("timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple"))) %>%
  
  group_by(model, coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
  filter(pthreshold == "s") %>%
  ungroup() %>%
  filter(model == "no_sex_covariate") %>%
  data.frame() %>%
  dplyr::select(gene) 

print()

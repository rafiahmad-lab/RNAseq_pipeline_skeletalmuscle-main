

######## Mixed model approach for DE ###########


# Preliminaries: dge_lists are store in ./data/derivedData/dge_lists/dge_lists.RDS


# RNA-seq data analysis using itearitive fitting of a Poisson- or negative binomial
# generalized mixed model. 
# The main purpose of the analysis is to identify differentially expressed genes between 
# training volume-conditions. This implementation is inspired by:

# Cui S, Ji T, Li J, Cheng J, Qiu J. 
# What if we ignore the random effects when analyzing RNA-seq data in a multifactor experiment. 
# Stat Appl Genet Mol Biol. 2016;15(2):87–105. doi:10.1515/sagmb-2015-0011


source("./R/lib_fun.R")


# FULL MODEL -- all time points
# Contructing the general loop for model based normalization.
# Effective library size (lib size * norm factor) is included in the model


## Setting up parallelization
## This may need to be changed



############ mgcv solution  #########################################

# mgcv has considerable advantage in terms of speed compared to glmmTMB
# a description of random effects in negative binomial gams:
# https://www.fromthebottomoftheheap.net/2017/05/04/compare-mgcv-with-glmmTMB/




library(mgcv)

# read DGEList

dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")



## Save results in list
mgcv_results <- list()



for(j in 1:length(dge_lists)) {
  
  genes <- rownames(dge_lists[[j]])
  
  
  # set up parallel processing 
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) ### Change this when live-looping!
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  ### foreach loop 
  results <- foreach(i = 1:length(genes), 
                     .packages = c("mgcv", "dplyr"), 
                     .options.snow = opts) %dopar% {
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           
                           which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
                             if (! na.rm & any (is.na (x)))
                               return (rep (NA_integer_, length (probs)))
                             
                             o <- order (x)
                             n <- sum (! is.na (x))
                             o <- o [seq_len (n)]
                             
                             nppm <- n * probs - 0.5
                             j <- floor(nppm)
                             h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
                             j <- j + h
                             
                             j [j == 0] <- 1
                             o[j]
                           }
                           
                           
                           
                           ## Calculate reference library (this is set as the median library)
                           reference.lib <- dge_lists[[j]]$samples[which.quantile(dge_lists[[j]]$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                           # Note:
                           # The normalization factor used in Cui et al:
                           
                           # "The library scaling factor adjusts for the differential expression 
                           # caused by both the differential sequencing depths and the differential RNA 
                           # compositions of different RNA samples. It can be obtained by dividing the 
                           # effective library size of each library to that of a reference library. 
                           # Here the effective library size refers to the product of the original 
                           # library size, the total number of read counts per library, and a normalization 
                           # factor to adjust for the  RNA  composition  effect.  Throughout  the  paper,  
                           # we  use  the  trimmed  mean  method  (TMM)  of  Robinson  and Oshlack (2010) 
                           # to calculate the normalization factor for RNA composition effect, which uses a 
                           # weighted trimmed  mean  of  the  log  expression  ratios  across  genes  to  
                           # estimate  the  global  fold  change  of  two  samples  to adjust for the RNA 
                           # composition effect, assuming that majority of genes are non-differentially expressed."
                           
                           
                           ## Extract data for each sample and calculate normalization factor
                           
                           dat <- dge_lists[[j]]$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = dge_lists[[j]]$counts[genes[i],], 
                                             sample = colnames(dge_lists[[j]]$counts)) %>%
                             inner_join(dat) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject)) 
                           
                           
                           ## The mgcv negative binomial model  ##
                           
                           m1 <- gam(counts ~ nf + time + time:sets + s(subject, bs = "re"), data = dat,
                                     family = nb, method = "ML")
                           
                           
                           # save results
                           results <- data.frame(summary(m1)$p.table) %>%
                             mutate(coef = rownames(.), 
                                    gene = genes[1], 
                                    method = names(dge_lists[j])) %>%
                             dplyr::select(gene, 
                                           method, 
                                           coef, 
                                           estimate = Estimate, 
                                           se = Std..Error, 
                                           z.val = z.value, 
                                           p.val = Pr...z..) 
                           
                           results
                           
                         },
                         error = function(e){
                           message('** ERR at ', Sys.time(), " **")
                           print(e)
                           
                           
                         })
                       
                       
                       
                     }
  
  close(pb)
  stopCluster(cl)
  
  mgcv_results[[j]] <- bind_rows(results)
  
  
  
}


saveRDS(mgcv_results, file = "./data/derivedData/DE/mgcv_results.RDS")



### Models with offset term and sex as covariate #### 



# read DGEList

dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")


## Save results in list
mgcv_results2 <- list()



for(j in 1:length(dge_lists)) {
  
  genes <- rownames(dge_lists[[j]])
  
  
  # set up parallel processing 
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) ### Change this when live-looping!
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  ### foreach loop 
  results <- foreach(i = 1:length(genes),
                     .packages = c("mgcv", "dplyr"), 
                     .options.snow = opts) %dopar% {
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           
                           which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
                             if (! na.rm & any (is.na (x)))
                               return (rep (NA_integer_, length (probs)))
                             
                             o <- order (x)
                             n <- sum (! is.na (x))
                             o <- o [seq_len (n)]
                             
                             nppm <- n * probs - 0.5
                             j <- floor(nppm)
                             h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
                             j <- j + h
                             
                             j [j == 0] <- 1
                             o[j]
                           }
                           
                           
                           
                           ## Calculate reference library (this is set as the median library)
                           reference.lib <- dge_lists[[j]]$samples[which.quantile(dge_lists[[j]]$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                           # Note:
                           # The normalization factor used in Cui et al:
                           
                           # "The library scaling factor adjusts for the differential expression 
                           # caused by both the differential sequencing depths and the differential RNA 
                           # compositions of different RNA samples. It can be obtained by dividing the 
                           # effective library size of each library to that of a reference library. 
                           # Here the effective library size refers to the product of the original 
                           # library size, the total number of read counts per library, and a normalization 
                           # factor to adjust for the  RNA  composition  effect.  Throughout  the  paper,  
                           # we  use  the  trimmed  mean  method  (TMM)  of  Robinson  and Oshlack (2010) 
                           # to calculate the normalization factor for RNA composition effect, which uses a 
                           # weighted trimmed  mean  of  the  log  expression  ratios  across  genes  to  
                           # estimate  the  global  fold  change  of  two  samples  to adjust for the RNA 
                           # composition effect, assuming that majority of genes are non-differentially expressed."
                           
                           
                           ## Extract data for each sample and calculate normalization factor
                           
                           dat <- dge_lists[[j]]$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = dge_lists[[j]]$counts[genes[i],], 
                                             sample = colnames(dge_lists[[j]]$counts)) %>%
                             inner_join(dat) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject), 
                                    eff.lib = lib.size * norm.factors) 
                           
                           
                           ## The mgcv negative binomial model  ##
                           
                           m1 <- gam(counts ~  time + time:sets + s(subject, bs = "re") , 
                                   offset = log(eff.lib),
                                     data = dat,
                                     family = nb, method = "REML")

                          m3 <- gam(counts ~  nf + time + time:sets + s(subject, bs = "re") , 
                                     data = dat,
                                     family = nb, method = "REML")
                           
                           # naive model (no normalization)
                           m4 <- gam(counts ~  time + time:sets + s(subject, bs = "re") , 
                                     data = dat,
                                     family = nb, method = "REML")
                           
                           
                           # save results
                           
                        results  <-  rbind(data.frame(summary(m1)$p.table) %>%
                                 mutate(coef = rownames(.), 
                                        model = "eff.lib_offset", 
                                    gene = genes[i], 
                                    method = names(dge_lists[j])) , 
                                 
                                 data.frame(summary(m2)$p.table) %>%
                                   mutate(coef = rownames(.), 
                                          model = "relative.nf_offset", 
                                          gene = genes[i], 
                                          method = names(dge_lists[j])) , 
                                 
                                 data.frame(summary(m3)$p.table) %>%
                                   mutate(coef = rownames(.), 
                                          model = "nf_covariate", 
                                          gene = genes[i], 
                                          method = names(dge_lists[j])) , 
                           
                                 data.frame(summary(m4)$p.table) %>%
                                   mutate(coef = rownames(.), 
                                          model = "naive", 
                                          gene = genes[i], 
                                          method = names(dge_lists[j]))) %>%
                             dplyr::select(gene, 
                                           method, 
                                           model,
                                           coef, 
                                           estimate = Estimate, 
                                           se = Std..Error, 
                                           z.val = z.value, 
                                           p.val = Pr...z..) 
                                 
                      # return results
                           results
                           
                         },
                         error = function(e){
                           message('** ERR at ', Sys.time(), " **")
                           print(e)
                           
                           
                         })
                       
                       
                       
                     }
  
  close(pb)
  stopCluster(cl)
  
  mgcv_results2[[j]] <- bind_rows(results)
  
  
  
}



bind_rows(mgcv_results2) %>%
  filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                     "timew12:setsmultiple")) %>%
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



bind_rows(mgcv_results2) %>%
  filter(coef == "nf") %>%
  ggplot(aes(estimate)) + geom_density()
 

bind_rows(mgcv_results2) %>%
  filter(coef == "nf") %>%
  ggplot(aes(-log10(p.val))) + geom_density()




# Notes on interpretation of models (fitted above).

# Assumptions in models: The "naive" model makes assumptions about constant mRNA/cell and library size.
# And, does not control for technical factors
# The effective library size offset model, and nf as a covariate makes assumptions about mRNA/cell.




# The offset model converts counts to a rate per offset, when the offset is the effective library size, 
# the model explains counts per library size. The offset parameter could potentially be set to muscle weight.
# As in previous analyses (qPCR), weight measures in w2post are not reliable. A resting-samples analysis 
# is therefore the only alternative. 








# read DGEList

dge_lists <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")

# Muscle weights in cDNA synthesis are loaded ##
mw <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  print()


## Save results in list
mgcv_results3 <- list()


for(j in 1:length(dge_lists)) {
  
  genes <- rownames(dge_lists[[j]])
  
  
  # set up parallel processing 
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) ### Change this when live-looping!
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  ### foreach loop 
  results <- foreach(i = 1:length(genes), 
                     .packages = c("mgcv", "dplyr"), 
                     .options.snow = opts) %dopar% {
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           
                           which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
                             if (! na.rm & any (is.na (x)))
                               return (rep (NA_integer_, length (probs)))
                             
                             o <- order (x)
                             n <- sum (! is.na (x))
                             o <- o [seq_len (n)]
                             
                             nppm <- n * probs - 0.5
                             j <- floor(nppm)
                             h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
                             j <- j + h
                             
                             j [j == 0] <- 1
                             o[j]
                           }
                           
                           
                           
                           ## Calculate reference library (this is set as the median library)
                           reference.lib <- dge_lists[[j]]$samples[which.quantile(dge_lists[[j]]$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                           # Note:
                           # The normalization factor used in Cui et al:
                           
                           # "The library scaling factor adjusts for the differential expression 
                           # caused by both the differential sequencing depths and the differential RNA 
                           # compositions of different RNA samples. It can be obtained by dividing the 
                           # effective library size of each library to that of a reference library. 
                           # Here the effective library size refers to the product of the original 
                           # library size, the total number of read counts per library, and a normalization 
                           # factor to adjust for the  RNA  composition  effect.  Throughout  the  paper,  
                           # we  use  the  trimmed  mean  method  (TMM)  of  Robinson  and Oshlack (2010) 
                           # to calculate the normalization factor for RNA composition effect, which uses a 
                           # weighted trimmed  mean  of  the  log  expression  ratios  across  genes  to  
                           # estimate  the  global  fold  change  of  two  samples  to adjust for the RNA 
                           # composition effect, assuming that majority of genes are non-differentially expressed."
                           
                           
                           ## Extract data for each sample and calculate normalization factor
                           
                           # Subset the dge-list
                           dge <- dge_lists[[j]][, dge_lists[[j]]$samples$time != "w2post",]
                           
                           
                           
                           dat <- dge$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = dge_lists[[j]]$counts[genes[i],], 
                                             sample = colnames(dge_lists[[j]]$counts)) %>%
                             inner_join(dat) %>%
                             inner_join(mw) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject), 
                                    sets = factor(sets, levels = c("single", "multiple")), 
                                    time = factor(time, levels = c("w0", "w2pre", "w12")),
                                    eff.lib = lib.size * norm.factors) 
                           
                           
                        
                           
                           ## Possibly, calcNormFactors should be use again?
                           
                           
                           ## The mgcv negative binomial model  ##
                           
                           # Offset log(tissue)
                           m1 <- gam(counts ~  time + time:sets + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, method = "REML")
                           
                          # Offset log(tissue) and normalization factor in model
                           m2 <- gam(counts ~  nf + time + time:sets + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, method = "REML")
                           
                       
                           
                           # naive model (no normalization)
                           m3 <- gam(counts ~  time + time:sets + s(subject, bs = "re") , 
                                     data = dat,
                                     family = nb, method = "REML")
                           
                           
                         
                           # save results
                           
                           results  <-  rbind(data.frame(summary(m1)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "tissue_offset", 
                                                       gene = genes[i], 
                                                       method = names(dge_lists[j])) , 
                                              
                                              data.frame(summary(m2)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "tissue_offset_nf_covariate", 
                                                       gene = genes[i], 
                                                       method = names(dge_lists[j])) , 
                                         
                                              
                                              data.frame(summary(m3)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "naive", 
                                                       gene = genes[i], 
                                                       method = names(dge_lists[j]))) %>%
                             dplyr::select(gene, 
                                           method, 
                                           model,
                                           coef, 
                                           estimate = Estimate, 
                                           se = Std..Error, 
                                           z.val = z.value, 
                                           p.val = Pr...z..) 
                           
                           # return results
                           results
                           
                         },
                         error = function(e){
                           message('** ERR at ', Sys.time(), " **")
                           print(e)
                           
                           
                         })
                       
                       
                       
                     }
  
  close(pb)
  stopCluster(cl)
  
  mgcv_results3[[j]] <- bind_rows(results)
  
  
  
}


saveRDS(mgcv_results3, file = "./data/derivedData/DE/mgcv_results_weightCorrected.RDS")





bind_rows(mgcv_results3) %>%
  filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                     "timew12:setsmultiple")) %>%
  mutate(coef = factor(coef, levels = c("timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple"))) %>%
  filter(method == "kallisto") %>%
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


bind_rows(mgcv_results3) %>%
  filter(coef %in% c("timew2pre", "timew12")) %>%
  mutate(coef = factor(coef, levels = c("timew2pre", "timew12"))) %>%
  filter(method == "kallisto") %>%
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





bind_rows(mgcv_results3) %>%
  filter(coef %in% c("timew2pre:setsmultiple")) %>%
  filter(method == "kallisto", 
         -log10(p.val) > 6) %>%
  print()


########### mgcv single method foreach #####

# read DGEList

temp <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")

dge_lists <- list()


dge_lists[[1]] <- temp[[2]]

names(dge_lists) <- c("kallisto")
#

# Muscle weights in cDNA synthesis are loaded ##
mw <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  print()


## Save results in list
mgcv_results4 <- list()



for(j in 1:length(dge_lists)) {
  
  genes <- rownames(dge_lists[[j]])
  
  
  # set up parallel processing 
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) ### Change this when live-looping!
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  ### foreach loop 
  results <- foreach(i = 1:length(genes), 
                     .packages = c("mgcv", "dplyr"), 
                     .options.snow = opts) %dopar% {
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           
                           which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
                             if (! na.rm & any (is.na (x)))
                               return (rep (NA_integer_, length (probs)))
                             
                             o <- order (x)
                             n <- sum (! is.na (x))
                             o <- o [seq_len (n)]
                             
                             nppm <- n * probs - 0.5
                             j <- floor(nppm)
                             h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
                             j <- j + h
                             
                             j [j == 0] <- 1
                             o[j]
                           }
                           
                           
                           
                           ## Calculate reference library (this is set as the median library)
                           reference.lib <- dge_lists[[j]]$samples[which.quantile(dge_lists[[j]]$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                           # Note:
                           # The normalization factor used in Cui et al:
                           
                           # "The library scaling factor adjusts for the differential expression 
                           # caused by both the differential sequencing depths and the differential RNA 
                           # compositions of different RNA samples. It can be obtained by dividing the 
                           # effective library size of each library to that of a reference library. 
                           # Here the effective library size refers to the product of the original 
                           # library size, the total number of read counts per library, and a normalization 
                           # factor to adjust for the  RNA  composition  effect.  Throughout  the  paper,  
                           # we  use  the  trimmed  mean  method  (TMM)  of  Robinson  and Oshlack (2010) 
                           # to calculate the normalization factor for RNA composition effect, which uses a 
                           # weighted trimmed  mean  of  the  log  expression  ratios  across  genes  to  
                           # estimate  the  global  fold  change  of  two  samples  to adjust for the RNA 
                           # composition effect, assuming that majority of genes are non-differentially expressed."
                           
                           
                           ## Extract data for each sample and calculate normalization factor
                           
                           # Subset the dge-list
                           dge <- dge_lists[[j]][, dge_lists[[j]]$samples$time != "w2post",]
                           
                           # Possibly re-do calcNormFactors
                           
                           
                           dat <- dge$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = dge_lists[[j]]$counts[genes[i],], 
                                             sample = colnames(dge_lists[[j]]$counts)) %>%
                             inner_join(dat) %>%
                             inner_join(mw) %>%
                             # Remove subject/time with bad muscle weight
                             filter(!(subject == "FP15" & time == "w2pre")) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject), 
                                    sets = factor(sets, levels = c("single", "multiple")), 
                                    time = factor(time, levels = c("w0", "w2pre", "w12")),
                                    eff.lib = lib.size * norm.factors)  
                           
                           
                           
                           
                           ## Possibly, calcNormFactors should be use again?
                           
                           
                           ## The mgcv negative binomial model  ##
                           
                           # Offset log(tissue)
                           
                           
                           # Offset log(tissue) and normalization factor in model
                           m1 <- gam(counts ~  nf + time + time:sets + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, method = "ML")
                           
                           
                           
                           # naive model (no normalization)
                           m2 <- gam(counts ~  sex + nf + time + time:sets + s(subject, bs = "re") , 
                                    offset = log(tissue),
                                    data = dat,
                                    family = nb, method = "ML")
                           
                           
                           
                           # save results
                           
                           results  <-  rbind(data.frame(summary(m1)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "no_sex_covariate", 
                                                       gene = genes[i], 
                                                       method = names(dge_lists[j])) , 
                                              
                                              data.frame(summary(m2)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "sex_covariate", 
                                                       gene = genes[i], 
                                                       method = names(dge_lists[j]))) %>%
                             dplyr::select(gene, 
                                           method, 
                                           model,
                                           coef, 
                                           estimate = Estimate, 
                                           se = Std..Error, 
                                           z.val = z.value, 
                                           p.val = Pr...z..) 
                           
                           # return results
                           results
                           
                         },
                         error = function(e){
                           message('** ERR at ', Sys.time(), " **")
                           print(e)
                           
                           
                         })
                       
                       
                       
                     }
  
  close(pb)
  stopCluster(cl)
  
  mgcv_results4[[j]] <- bind_rows(results)
  
  
  
}


bind_rows(mgcv_results4) %>%
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


de_genes<- bind_rows(mgcv_results4) %>%
  filter(coef == "timew2pre", 
         model == "sex_covariate") %>%
  group_by(model, coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 1, "s", "ns")) %>%
  ungroup() %>%
  filter(fcthreshold == "s", pthreshold == "s") %>%
  filter(log2fc < -1) %>%
  dplyr::select(gene) %>%
  data.frame()
  print()




de_genes <- bind_rows(mgcv_results4) %>%
  filter(coef %in% c("timew0:setsmultiple",  "timew2pre:setsmultiple", "timew2post:setsmultiple", 
                     "timew12:setsmultiple")) %>%
  mutate(coef = factor(coef, levels = c("timew0:setsmultiple", "timew2pre:setsmultiple", "timew12:setsmultiple"))) %>%
  
  group_by(model, coef) %>%
  mutate(adj.p = p.adjust(p.val, method = "fdr"), 
         pthreshold = if_else(adj.p > 0.05, "ns", "s"), 
         log2fc = estimate/log(2), 
         fcthreshold = if_else(abs(log2fc) > 0.585, "s", "ns")) %>%
  filter(pthreshold == "s" & fcthreshold == "s") %>%
  ungroup() %>%
  filter(model == "no_sex_covariate") %>%
  data.frame() %>%
  dplyr::select(gene) 

  print()
  
  
  
  ?goana()
  
  
  library(org.Hs.eg.db)


# Get entrez gene IDs that are mapped to Ensembl ID
x <- org.Hs.egENSEMBL

mapped_genes <- mappedkeys(x)

xx <- as.list(x[mapped_genes])

entrezid <- as.character(mapIds(org.Hs.eg.db, keys = de_genes[,1], keytype = "ENSEMBL", column = 'ENTREZID'))

?topGO
go <- goana(de = entrezid)


topGO(go)

keg <- kegga(entrezid)  

topKEGG(keg)

data.frame(topGO(go)) %>%
  filter(Ont == "MF") %>%
  print()




####### Sex-time interactions ##################

########### mgcv single method foreach #####

# read DGEList

temp <- readRDS("./data/derivedData/dge_lists/dge_list.RDS")

dge_lists <- list()


dge_lists[[1]] <- temp[[2]]

names(dge_lists) <- c("kallisto")
#

# Muscle weights in cDNA synthesis are loaded ##
mw <- read_excel("./data/RNAamount.xlsx") %>%
  mutate(RNA.per.mg = (conc*elution.volume) / prot.mrna1) %>% 
  filter(include == "incl") %>%
  filter(timepoint != "w2post") %>%
  dplyr::select(subject, sets, time = timepoint, RNA.per.mg) %>%
  mutate(tissue = 1000 / RNA.per.mg, 
         tl = log(tissue))  %>%
  print()


## Save results in list
mgcv_results5 <- list()



for(j in 1:length(dge_lists)) {
  
  genes <- rownames(dge_lists[[j]])
  
  
  # set up parallel processing 
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) ### Change this when live-looping!
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  ### foreach loop 
  results <- foreach(i = 1:length(genes), 
                     .packages = c("mgcv", "dplyr"), 
                     .options.snow = opts) %dopar% {
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           
                           which.quantile <- function (x, probs = 0.5, na.rm = FALSE){
                             if (! na.rm & any (is.na (x)))
                               return (rep (NA_integer_, length (probs)))
                             
                             o <- order (x)
                             n <- sum (! is.na (x))
                             o <- o [seq_len (n)]
                             
                             nppm <- n * probs - 0.5
                             j <- floor(nppm)
                             h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
                             j <- j + h
                             
                             j [j == 0] <- 1
                             o[j]
                           }
                           
                           
                           
                           ## Calculate reference library (this is set as the median library)
                           reference.lib <- dge_lists[[j]]$samples[which.quantile(dge_lists[[j]]$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                           # Note:
                           # The normalization factor used in Cui et al:
                           
                           # "The library scaling factor adjusts for the differential expression 
                           # caused by both the differential sequencing depths and the differential RNA 
                           # compositions of different RNA samples. It can be obtained by dividing the 
                           # effective library size of each library to that of a reference library. 
                           # Here the effective library size refers to the product of the original 
                           # library size, the total number of read counts per library, and a normalization 
                           # factor to adjust for the  RNA  composition  effect.  Throughout  the  paper,  
                           # we  use  the  trimmed  mean  method  (TMM)  of  Robinson  and Oshlack (2010) 
                           # to calculate the normalization factor for RNA composition effect, which uses a 
                           # weighted trimmed  mean  of  the  log  expression  ratios  across  genes  to  
                           # estimate  the  global  fold  change  of  two  samples  to adjust for the RNA 
                           # composition effect, assuming that majority of genes are non-differentially expressed."
                           
                           
                           ## Extract data for each sample and calculate normalization factor
                           
                           # Subset the dge-list
                           dge <- dge_lists[[j]][, dge_lists[[j]]$samples$time != "w2post",]
                           
                           # Possibly re-do calcNormFactors
                           
                           
                           dat <- dge$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = dge_lists[[j]]$counts[genes[i],], 
                                             sample = colnames(dge_lists[[j]]$counts)) %>%
                             inner_join(dat) %>%
                             inner_join(mw) %>%
                             # Remove subject/time with bad muscle weight
                             filter(!(subject == "FP15" & time == "w2pre")) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject), 
                                    sets = factor(sets, levels = c("single", "multiple")), 
                                    time = factor(time, levels = c("w0", "w2pre", "w12")),
                                    eff.lib = lib.size * norm.factors)  
                           
                           
                           
                           
                           ## Possibly, calcNormFactors should be use again?
                           
                           
                           ## The mgcv negative binomial model  ##
                           
                           # Offset log(tissue)
                           
                           
                           # Offset log(tissue) and normalization factor in model
                           m1 <- gam(counts ~  nf + time * sex + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, method = "ML")
                           
                           
                           
                           # naive model (no normalization)
                           m2 <- gam(counts ~   nf + time + sex + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, method = "ML")
                           
                        global.p  <-  anova(m2, m1, test = "Chisq")[2, 5]
                           
                           
                           
                           # save results
                           
                           results  <-  data.frame(summary(m1)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "no_sex_covariate", 
                                                       global.p = global.p,
                                                       gene = genes[i], 
                                                       method = names(dge_lists[j])) %>%
                             dplyr::select(gene, 
                                           method, 
                                           model,
                                           coef, 
                                           estimate = Estimate, 
                                           se = Std..Error, 
                                           z.val = z.value, 
                                           p.val = Pr...z.., 
                                           sex.interaction.p = global.p) 
                           
                           # return results
                           results
                           
                         },
                         error = function(e){
                           message('** ERR at ', Sys.time(), " **")
                           print(e)
                           
                           
                         })
                       
                       
                       
                     }
  
  close(pb)
  stopCluster(cl)
  
  mgcv_results5[[j]] <- bind_rows(results)
  
  
  
}





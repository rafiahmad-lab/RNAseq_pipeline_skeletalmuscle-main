##### Mixed models 2 ############################

# This is an update of the mixed model approach used in mixed_models.R



source("./R/lib_fun.R")


# Preliminaries: 
# dge_lists are store in ./data/derivedData/dge_lists/dge_lists.RDS
# dge_lists are used for gene-wise iterative modeling. 
# The loop is used for one selected list...
# See below for settings on fitting algorithm 
# and diagnostics.



# Notes: 

# RNA-seq data analysis using itearitive fitting of a Poisson- or negative binomial
# generalized mixed model. 
# The main purpose of the analysis is to identify differentially expressed genes between 
# training volume-conditions and time-point regardless of condition. This implementation is inspired by:

# Cui S, Ji T, Li J, Cheng J, Qiu J. 
# What if we ignore the random effects when analyzing RNA-seq data in a multifactor experiment. 
# Stat Appl Genet Mol Biol. 2016;15(2):87–105. doi:10.1515/sagmb-2015-0011

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

# Diagnostics (Under construction!)
# If diagnostics are used, the DHARMa package will test for overdispersion/underdispersion
# and zero inflation 


# 

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


# See that these correspond
selected.dge <- dge_lists$rsem
selected.method <- "rsem"

# A list of genes (available in the dge-list)
  genes <- rownames(selected.dge)
  
  
  # set up parallel processing using the foreach package 
  cores <- detectCores()
  cl <- makeCluster(cores[1]) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) 
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  # foreach loop 
  results <- foreach(i = 1:length(genes), 
                     .packages = c("mgcv", "glmmTMB", "dplyr", "DHARMa"),  # Packages used in fitting/data wrangling
                     .options.snow = opts) %dopar% {
                       
                       
                       # settings
                       # fitting
                       
                       fit.mgcv <- FALSE
                       fit.glmmTMB <- TRUE
                       diagnostics <- TRUE
                       
                       # glmmTMB allows for diagnostics
                    
                       
                       
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           # the function is used for selecting the median sized library
                           # which in turn is used as reference library as suggested in Cui
                           
                           
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
                           reference.lib <- selected.dge$samples[which.quantile(selected.dge$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                          # Extract data for each sample and calculate normalization factor
                           
                           # Subset the dge-list
                           dge <- selected.dge[, selected.dge$samples$time != "w2post",]
                           
                          
                           dat <- dge$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = selected.dge$counts[genes[i],], 
                                             sample = colnames(selected.dge$counts)) %>%
                             inner_join(dat) %>%
                             inner_join(mw) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject), 
                                    sets = factor(sets, levels = c("single", "multiple")), 
                                    time = factor(time, levels = c("w0", "w2pre", "w12")),
                                    eff.lib = lib.size * norm.factors) %>%
                          # One subject/timepoint/condition deviated from the expected 
                          # RNA to tissue relationship. This observation is filtered out. 
                            filter(!(subject == "FP15" & time == "w2pre" & sets == "multiple"))

                           
                           
                           # Possibly, calcNormFactors should be use again?
                           #TODO check the effect of norm factor calculation on subset
                           
                           
                           # A list of results
                           results.models <- list()
                           
                          
                           ## The mgcv negative binomial model  ##
                           
                           if(fit.mgcv) { 
                             
                             # The naive model (no normalization)
                             
                               mgcv.m1 <- gam(counts ~  time + time:sets + s(subject, bs = "re") , 
                                     data = dat,
                                     family = nb, 
                                     method = "REML")
                               
                    
                             # The normalized model (eff library size)
                             
                               mgcv.m2 <- gam(counts ~  nf + time + time:sets + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, 
                                     method = "REML")
                             
                               # The muscle tissue corrected model (offset and eff.lib)
                             
                               mgcv.m3 <- gam(counts ~  time + time:sets + s(subject, bs = "re") , 
                                     offset = log(tissue),
                                     data = dat,
                                     family = nb, 
                                     method = "REML")
                             
                               # mgcv results 
                               
                             mgcv.results  <-  rbind(data.frame(summary(mgcv.m1)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "naive", 
                                                       gene = genes[i], 
                                                       method = selected.method, 
                                                       theta = mgcv.m1$family$getTheta(TRUE)) , 
                                              
                                              data.frame(summary(mgcv.m2)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "lib_size_normalized", 
                                                       gene = genes[i], 
                                                       method = selected.method, 
                                                       theta = mgcv.m2$family$getTheta(TRUE)) , 
                                              
                                              
                                              data.frame(summary(mgcv.m3)$p.table) %>%
                                                mutate(coef = rownames(.), 
                                                       model = "tissue_offset_lib_size_normalized", 
                                                       gene = genes[i], 
                                                       method = selected.method, 
                                                       theta = mgcv.m3$family$getTheta(TRUE))) %>%
                               
                               mutate( dispersion.statistic = NA, 
                                       dispersion.p        = NA,
                                       uniform.statistic   = NA,
                                       uniform.p      = NA) %>%
                               
                               dplyr::select(gene, 
                                             method, 
                                             model,
                                             coef, 
                                             estimate = Estimate, 
                                             se = Std..Error, 
                                             z.val = z.value, 
                                             p.val = Pr...z.., 
                                             theta,
                                             dispersion.statistic,
                                             dispersion.p,        
                                             uniform.statistic,   
                                             uniform.p) %>%
                               mutate(fitting = "mgcv")
                             
                             
                             results.models[[1]] <- mgcv.results
                             
                             
                             }
                           
                           if(fit.glmmTMB) {
                             
                             # The naive model (no normalization)
                             
                             glmmTMB.m1 <- glmmTMB(counts ~  time + time:sets + (1|subject) , 
                                            data = dat,
                                            family = "nbinom2")
                             
                             
                             # The normalized model (eff library size)
                             
                             glmmTMB.m2 <- glmmTMB(counts ~  nf + time + time:sets + (1|subject) , 
                                               data = dat,
                                               family = "nbinom2")
                             
                             # The muscle tissue corrected model (offset and eff.lib)
                             
                             glmmTMB.m3 <- glmmTMB(counts ~  nf + time + time:sets + (1|subject) , 
                                               offset = log(tissue),
                                               data = dat,
                                               family = "nbinom2")
                             
                             if(diagnostics) {
                               
                              res.glmmTMB.m1 <-  simulateResiduals(glmmTMB.m1, n = 1000)
                              res.glmmTMB.m2 <-  simulateResiduals(glmmTMB.m2, n = 1000)
                              res.glmmTMB.m3 <-  simulateResiduals(glmmTMB.m3, n = 1000)
                               
                              disp.m1    <- testDispersion(res.glmmTMB.m1, plot = FALSE)
                              disp.m2    <- testDispersion(res.glmmTMB.m2, plot = FALSE)
                              disp.m3    <- testDispersion(res.glmmTMB.m3, plot = FALSE)
                               
                              uniform.m1 <- testUniformity(res.glmmTMB.m1, plot = FALSE)
                              uniform.m2 <- testUniformity(res.glmmTMB.m2, plot = FALSE)
                              uniform.m3 <- testUniformity(res.glmmTMB.m3, plot = FALSE)
                              
                              
                              
                             } else {
                               
                               disp.m1 <- list(statistic = NA, p.value = NA)  
                               disp.m2   <- list(statistic = NA, p.value = NA) 
                               disp.m3   <- list(statistic = NA, p.value = NA) 
                               
                               uniform.m1 <- list(statistic = NA, p.value = NA) 
                               uniform.m2 <- list(statistic = NA, p.value = NA) 
                               uniform.m3 <- list(statistic = NA, p.value = NA) 
                               
                               
                               
                             }
                             
                             
                             glmmTMB.results  <-  rbind(data.frame(summary(glmmTMB.m1)$coef$cond) %>%
                                                       mutate(coef = rownames(.), 
                                                              model = "naive", 
                                                              gene = genes[i], 
                                                              method = selected.method, 
                                                              theta = sigma(glmmTMB.m1), 
                                                              dispersion.statistic = disp.m1$statistic, 
                                                              dispersion.p = disp.m1$p.value, 
                                                              uniform.statistic = uniform.m1$statistic,
                                                              uniform.p = uniform.m1$p.value) , 
                                                     
                                                       data.frame(summary(glmmTMB.m2)$coef$cond) %>%
                                                         mutate(coef = rownames(.), 
                                                                model = "lib_size_normalized", 
                                                                gene = genes[i], 
                                                                method = selected.method, 
                                                                theta = sigma(glmmTMB.m2), 
                                                                dispersion.statistic = disp.m2$statistic, 
                                                                dispersion.p = disp.m2$p.value, 
                                                                uniform.statistic = uniform.m2$statistic,
                                                                uniform.p = uniform.m2$p.value)  , 
                                                     
                                                     
                                                       data.frame(summary(glmmTMB.m3)$coef$cond) %>%
                                                         mutate(coef = rownames(.), 
                                                                model = "tissue_offset_lib_size_normalized", 
                                                                gene = genes[i], 
                                                                method = selected.method, 
                                                                theta = sigma(glmmTMB.m3), 
                                                                dispersion.statistic = disp.m3$statistic, 
                                                                dispersion.p         = disp.m3$p.value, 
                                                                uniform.statistic    = uniform.m3$statistic,
                                                                uniform.p            = uniform.m3$p.value)) %>%
                               dplyr::select(gene, 
                                             method, 
                                             model,
                                             coef, 
                                             estimate = Estimate, 
                                             se = Std..Error, 
                                             z.val = z.value, 
                                             p.val = Pr...z.., 
                                             theta,
                                             dispersion.statistic,
                                             dispersion.p,        
                                             uniform.statistic,   
                                             uniform.p) %>%
                               mutate(fitting = "glmmTMB")
                             
                             
                             results.models[[2]] <- glmmTMB.results
                             
                           } 
                           
                           # save results
                           results <- bind_rows(results.models)
                        
                        
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
  
  results <- bind_rows(results)
  
  
### This has been saved!  
  saveRDS(results, file = "./data/derivedData/DE/mixedmodel2_results1.RDS")
  
  

 

  
  #### Time-only model ######
  
  # Notes:
  
  # Similarly to what is done above, different normalization strategies are used
  # but the data are combined to reflect the average increase from time-point
  # w2pre.
  
  
  
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
  
  
  # See that these correspond
  selected.dge <- dge_lists$rsem
  selected.method <- "rsem"
  
  # A list of genes (available in the dge-list)
  genes <- rownames(selected.dge)
  
  
  # set up parallel processing using the foreach package 
  cores <- detectCores()
  cl <- makeCluster(cores[1]-4) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- length(genes) 
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  # foreach loop 
  results_timeonly <- foreach(i = 1:length(genes), 
                     .packages = c("mgcv", "glmmTMB", "dplyr", "DHARMa"),  # Packages used in fitting/data wrangling
                     .options.snow = opts) %dopar% {
                       
                       
                       # settings
                       # fitting
                       
                       fit.mgcv <- FALSE
                       fit.glmmTMB <- TRUE
                       diagnostics <- TRUE
                       
                       # glmmTMB allows for diagnostics
                       
                       
                       
                       
                       
                       tryCatch(
                         expr = {
                           
                           
                           # Define which quantile function 
                           # the function is used for selecting the median sized library
                           # which in turn is used as reference library as suggested in Cui
                           
                           
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
                           reference.lib <- selected.dge$samples[which.quantile(selected.dge$samples$lib.size, na.rm = TRUE),c(2,3)]
                           
                           ## Calculate effective library size for the reference library
                           reference.lib <- reference.lib[1,1] * reference.lib[1,2]
                           
                           # Extract data for each sample and calculate normalization factor
                           
                           # Subset the dge-list
                           dge <- selected.dge[, selected.dge$samples$time != "w2post",]
                           
                           
                           dat <- dge$samples %>%
                             mutate(sample = rownames(.), 
                                    nf = (lib.size * norm.factors) / reference.lib) 
                           
                           ## Extract gene counts and put together in the data frame
                           dat <- data.frame(counts = selected.dge$counts[genes[i],], 
                                             sample = colnames(selected.dge$counts)) %>%
                             inner_join(dat) %>%
                             inner_join(mw) %>%
                             mutate(counts = as.integer(round(counts, 0)), 
                                    subject = factor(subject), 
                                    sets = factor(sets, levels = c("single", "multiple")), 
                                    time = factor(time, levels = c("w0", "w2pre", "w12")),
                                    eff.lib = lib.size * norm.factors) %>%
                             # One subject/timepoint/condition deviated from the expected 
                             # RNA to tissue relationship. This observation is filtered out. 
                             filter(!(subject == "FP15" & time == "w2pre" & sets == "multiple"))
                           
                           
                           
                           # Possibly, calcNormFactors should be use again?
                           #TODO check the effect of norm factor calculation on subset
                           
                           
                           # A list of results
                           results.models <- list()
                           
                           
                           ## The mgcv negative binomial model  ##
                           
                           if(fit.mgcv) { 
                             
                             # The naive model (no normalization)
                             
                             mgcv.m1 <- gam(counts ~  time + s(subject, bs = "re") , 
                                            data = dat,
                                            family = nb, 
                                            method = "REML")
                             
                             
                             # The normalized model (eff library size)
                             
                             mgcv.m2 <- gam(counts ~  nf + time + s(subject, bs = "re") , 
                                            offset = log(tissue),
                                            data = dat,
                                            family = nb, 
                                            method = "REML")
                             
                             # The muscle tissue corrected model (offset and eff.lib)
                             
                             mgcv.m3 <- gam(counts ~  time + s(subject, bs = "re") , 
                                            offset = log(tissue),
                                            data = dat,
                                            family = nb, 
                                            method = "REML")
                             
                             # mgcv results 
                             
                             mgcv.results  <-  rbind(data.frame(summary(mgcv.m1)$p.table) %>%
                                                       mutate(coef = rownames(.), 
                                                              model = "naive", 
                                                              gene = genes[i], 
                                                              method = selected.method, 
                                                              theta = mgcv.m1$family$getTheta(TRUE)) , 
                                                     
                                                     data.frame(summary(mgcv.m2)$p.table) %>%
                                                       mutate(coef = rownames(.), 
                                                              model = "lib_size_normalized", 
                                                              gene = genes[i], 
                                                              method = selected.method, 
                                                              theta = mgcv.m2$family$getTheta(TRUE)) , 
                                                     
                                                     
                                                     data.frame(summary(mgcv.m3)$p.table) %>%
                                                       mutate(coef = rownames(.), 
                                                              model = "tissue_offset_lib_size_normalized", 
                                                              gene = genes[i], 
                                                              method = selected.method, 
                                                              theta = mgcv.m3$family$getTheta(TRUE))) %>%
                               
                               mutate( dispersion.statistic = NA, 
                                       dispersion.p        = NA,
                                       uniform.statistic   = NA,
                                       uniform.p      = NA) %>%
                               
                               dplyr::select(gene, 
                                             method, 
                                             model,
                                             coef, 
                                             estimate = Estimate, 
                                             se = Std..Error, 
                                             z.val = z.value, 
                                             p.val = Pr...z.., 
                                             theta,
                                             dispersion.statistic,
                                             dispersion.p,        
                                             uniform.statistic,   
                                             uniform.p) %>%
                               mutate(fitting = "mgcv")
                             
                             
                             results.models[[1]] <- mgcv.results
                             
                             
                           }
                           
                           if(fit.glmmTMB) {
                             
                             # The naive model (no normalization)
                             
                             glmmTMB.m1 <- glmmTMB(counts ~  time + (1|subject) , 
                                                   data = dat,
                                                   family = "nbinom2")
                             
                             
                             # The normalized model (eff library size)
                             
                             glmmTMB.m2 <- glmmTMB(counts ~  nf + time + (1|subject) , 
                                                   data = dat,
                                                   family = "nbinom2")
                             
                             # The muscle tissue corrected model (offset and eff.lib)
                             
                             glmmTMB.m3 <- glmmTMB(counts ~  nf + time + (1|subject) , 
                                                   offset = log(tissue),
                                                   data = dat,
                                                   family = "nbinom2")
                             
                             
                            
                             
                             
                          
                             
                             if(diagnostics) {
                               
                               res.glmmTMB.m1 <-  simulateResiduals(glmmTMB.m1, n = 1000)
                               res.glmmTMB.m2 <-  simulateResiduals(glmmTMB.m2, n = 1000)
                               res.glmmTMB.m3 <-  simulateResiduals(glmmTMB.m3, n = 1000)
                               
                               disp.m1    <- testDispersion(res.glmmTMB.m1, plot = FALSE)
                               disp.m2    <- testDispersion(res.glmmTMB.m2, plot = FALSE)
                               disp.m3    <- testDispersion(res.glmmTMB.m3, plot = FALSE)
                               
                               uniform.m1 <- testUniformity(res.glmmTMB.m1, plot = FALSE)
                               uniform.m2 <- testUniformity(res.glmmTMB.m2, plot = FALSE)
                               uniform.m3 <- testUniformity(res.glmmTMB.m3, plot = FALSE)
                               
                                  # Test if the extra dispersion parameter is called for (i.e. compare nbinom with poisson).
                             glmmTMB.m4 <- glmmTMB(counts ~  nf + time + (1|subject) , 
                                                   offset = log(tissue),
                                                   data = dat,
                                                   family = "poisson") 
                             
                             poisson.test <- anova(glmmTMB.m3, glmmTMB.m4)
                               
                             
                             poisson.test <- poisson.test[2, 8]
                             
                             } else {
                               
                               disp.m1 <- list(statistic = NA, p.value = NA)  
                               disp.m2   <- list(statistic = NA, p.value = NA) 
                               disp.m3   <- list(statistic = NA, p.value = NA) 
                               
                               uniform.m1 <- list(statistic = NA, p.value = NA) 
                               uniform.m2 <- list(statistic = NA, p.value = NA) 
                               uniform.m3 <- list(statistic = NA, p.value = NA) 
                               
                               poisson.test <- NA
                               
                             }
                             
                             
                             glmmTMB.results  <-  rbind(data.frame(summary(glmmTMB.m1)$coef$cond) %>%
                                                          mutate(coef = rownames(.), 
                                                                 model = "naive", 
                                                                 gene = genes[i], 
                                                                 method = selected.method, 
                                                                 theta = sigma(glmmTMB.m1), 
                                                                 dispersion.statistic = disp.m1$statistic, 
                                                                 dispersion.p = disp.m1$p.value, 
                                                                 uniform.statistic = uniform.m1$statistic,
                                                                 uniform.p = uniform.m1$p.value, 
                                                                 poisson.test = NA) , 
                                                        
                                                        data.frame(summary(glmmTMB.m2)$coef$cond) %>%
                                                          mutate(coef = rownames(.), 
                                                                 model = "lib_size_normalized", 
                                                                 gene = genes[i], 
                                                                 method = selected.method, 
                                                                 theta = sigma(glmmTMB.m2), 
                                                                 dispersion.statistic = disp.m2$statistic, 
                                                                 dispersion.p = disp.m2$p.value, 
                                                                 uniform.statistic = uniform.m2$statistic,
                                                                 uniform.p = uniform.m2$p.value, 
                                                                 poisson.test = NA)  , 
                                                        
                                                        
                                                        data.frame(summary(glmmTMB.m3)$coef$cond) %>%
                                                          mutate(coef = rownames(.), 
                                                                 model = "tissue_offset_lib_size_normalized", 
                                                                 gene = genes[i], 
                                                                 method = selected.method, 
                                                                 theta = sigma(glmmTMB.m3), 
                                                                 dispersion.statistic = disp.m3$statistic, 
                                                                 dispersion.p         = disp.m3$p.value, 
                                                                 uniform.statistic    = uniform.m3$statistic,
                                                                 uniform.p            = uniform.m3$p.value, 
                                                                 poisson.test = poisson.test)) %>%
                               dplyr::select(gene, 
                                             method, 
                                             model,
                                             coef, 
                                             estimate = Estimate, 
                                             se = Std..Error, 
                                             z.val = z.value, 
                                             p.val = Pr...z.., 
                                             theta,
                                             dispersion.statistic,
                                             dispersion.p,        
                                             uniform.statistic,   
                                             uniform.p, 
                                             poisson.test) %>%
                               mutate(fitting = "glmmTMB")
                             
                             
                             results.models[[2]] <- glmmTMB.results
                             
                           } 
                           
                           # save results
                           results <- bind_rows(results.models)
                           
                           
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
  
  results <- bind_rows(results_timeonly)
  
  ### This has been saved!  
  saveRDS(results, file = "./data/derivedData/DE/mixedmodel2_results2.RDS")
  
  
  ##### Acute exercise model ###############
  
  # Notes:
  
  # In the acute phase, we do not have valid measurements of muscle weight.
  # As an effect, we have to assume equal input of material in each reaction. 
  # We also have to assume equal mRNA to RNA ratios over the two time-points (Week 2)
  
  # Model-wise it is the interaction term that is interesting.
  
  
  
  
  
  
  
  

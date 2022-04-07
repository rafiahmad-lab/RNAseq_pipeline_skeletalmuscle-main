######### Differential expression analysis on DGELists ############
#
# Content:
# 1. voom/limma approach with duplicate correlation (full data set)
# 2. edgeR glmQLFit with subjects as "fixed effects" on subsets of the 
# dge_lists to simplify comparisons (correlations btw subj)
# 3. dream/variancePartition
# 

# Todo:
# Control for sex, look at p-value distributions with/without sex as a covariate
# Do acute phase diff-expr with baseline correction ...

# source libraries and functions

source("./R/lib_fun.R")


## Load DGELists 


dge_lists <- readRDS(file = "./data/derivedData/dge_lists/dge_list.RDS")


######### lmFit approach ####################


## Estimating dispersion 
# Note: this adds (modifies) the design matrix 




for(i in 1:length(dge_lists)) {
  
  dge_lists[[i]]$samples <- dge_lists[[i]]$samples %>%
    mutate(condition = factor(paste0(time, "_", sets), 
                              levels = c("w0_single", 
                                         "w0_multiple", 
                                         "w2pre_single", 
                                         "w2pre_multiple",
                                         "w2post_single",
                                         "w2post_multiple",
                                         "w12_single", 
                                         "w12_multiple"))) 
  
  
  dge_lists[[i]]$design <- model.matrix(~ 0 + condition, data = dge_lists[[i]]$samples)

  dge_lists[[i]] <- estimateDisp(dge_lists[[i]], design = dge_lists[[i]]$design, robust = TRUE)

}

### Explorative plots of dispersion parameters #### 

par(mfrow = c(2, 3))

for(i in 1:length(dge_lists)) {
  
  plotBCV(dge_lists[[i]], main = names(dge_lists[i]))
  
}
par(mfrow = c(1,1))

## Dispersion parameters per method ##

data.frame(method = names(dge_lists), 
           dispersion = c(dge_lists[[1]]$common.dispersion, 
                          dge_lists[[2]]$common.dispersion,
                          dge_lists[[3]]$common.dispersion,
                          dge_lists[[4]]$common.dispersion,
                          dge_lists[[5]]$common.dispersion)) %>%
  mutate(BVC = sqrt(dispersion)) %>%
  print()




#### Store voom/limma results in list 

vl_results<- list()  

# NOTES: 
# The design matrix (stored in the DGELists), is ~ 0 + condition.
# This is a simplified model to retrieve custom contrasts (i.e. mean of w0 vs. mean of w2pre).
# Within time-point comparisons can also be made.


# Make custom contrasts


contr <- makeContrasts(w2pre  = (conditionw2pre_single + conditionw2pre_multiple)/2 - (conditionw0_single + conditionw0_multiple)/2, 
                       w2post = (conditionw2post_single + conditionw2post_multiple)/2 - (conditionw0_single + conditionw0_multiple)/2 ,
                       w12    =  (conditionw12_single + conditionw12_multiple)/2  - (conditionw0_single + conditionw0_multiple)/2 ,
              
              w0_sm = conditionw0_multiple - conditionw0_single, 
              w2pre_sm = conditionw2pre_multiple - conditionw2pre_single, 
              w2post_sm =conditionw2post_multiple - conditionw2post_single,
              w12_sm = conditionw12_multiple - conditionw12_single,
              
              levels = c("conditionw0_single", 
                         "conditionw0_multiple", 
                         "conditionw2pre_single", 
                         "conditionw2pre_multiple",
                         "conditionw2post_single",
                         "conditionw2post_multiple",
                         "conditionw12_single", 
                         "conditionw12_multiple"))




# Loop over dge lists 
for(i in 1:length(dge_lists)) {
  
  
  v <- voom(dge_lists[[i]], design = dge_lists[[i]]$design)
  
  corfit <- duplicateCorrelation(v, dge_lists[[i]]$design, block = dge_lists[[i]]$samples$subject)  
  
  v <- voom(dge_lists[[i]], design = dge_lists[[i]]$design, block = dge_lists[[i]]$samples$subject, 
            correlation = corfit$consensus)
  
  fit <- lmFit(v, design = dge_lists[[i]]$design, block = dge_lists[[i]]$samples$subject, 
               correlation = corfit$consensus)
  
  contr.fit <-  contrasts.fit(fit, contrasts = contr)
  
  contr.fit2 <- eBayes(contr.fit)
  
  
  # Retrieve estimates 
 estimates <- rbind(data.frame(topTable(contr.fit2, coef = 1, number = Inf, adjust.method = "fdr"),
                               coef = "w2pre", 
                               method = names(dge_lists[i])),
                    
                    data.frame(topTable(contr.fit2, coef = 2, number = Inf, adjust.method = "fdr"),
                               coef = "w2post", 
                               method = names(dge_lists[i])),
                    
                    data.frame(topTable(contr.fit2, coef = 3, number = Inf, adjust.method = "fdr"),
                               coef = "w12", 
                               method = names(dge_lists[i])),
                    
                    data.frame(topTable(contr.fit2, coef = 4, number = Inf, adjust.method = "fdr"),
                               coef = "w0_sm", 
                               method = names(dge_lists[i])),
                    
                    data.frame(topTable(contr.fit2, coef = 5, number = Inf, adjust.method = "fdr"),
                               coef = "w2pre_sm", 
                               method = names(dge_lists[i])),
                    
                    data.frame(topTable(contr.fit2, coef = 6, number = Inf, adjust.method = "fdr"),
                               coef = "w2post_sm", 
                               method = names(dge_lists[i])),
                    
                    data.frame(topTable(contr.fit2, coef = 7, number = Inf, adjust.method = "fdr"),
                               coef = "w12_sm", 
                               method = names(dge_lists[i])))
 
 
  
  # Organize all results in a list 
  
  results <- list(v = v, 
                  corfit = corfit, 
                  fit = fit, 
                  contr.fit = contr.fit2, 
                  estimates = estimates)
  
  
  vl_results[[i]] <- results
  
  names(vl_results)[i] <- names(dge_lists)[i]
  
  message("n loops = ", i)
  
  
}


#### Retrieve results from lmFit-aproach 

lmfit_estimates <- bind_rows(vl_results[[1]]$estimates, 
                             vl_results[[2]]$estimates, 
                             vl_results[[3]]$estimates, 
                             vl_results[[4]]$estimates, 
                             vl_results[[5]]$estimates)


# Save results 
saveRDS(list(lmfit_estimates = lmfit_estimates, 
             lmfit_results = vl_results), 
        file = "./data/derivedData/DE/lm_results.RDS")

lmfit_estimates %>%
  mutate(pthreshold = if_else(adj.P.Val > 0.05, "ns", "s")) %>%
  ggplot(aes(logFC, -log10(P.Value), fill = pthreshold)) + 
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  theme_bw() +
  facet_grid(method ~ coef)

## Clean environment 
rm(vl_results, v, temp, temp2, fit, corfit, contr.fit, contr.fit2, estimates, lmfit_estimates)




##### GLM fit approach #############


# Instead of block in voom/limma/lmFit approach, subjects can be entered to the model 
# as fixed effects (!). See p. 39 in the edge R user guide 
# https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
# To simplify, this is done separatly for each time point

# dispersion estimates must be made with this itkaen into account 



results_glmfit_paired <- list()

for(i in 1:length(dge_lists)) {
  
  
  ## Pairwise comparisons in each time-point
   
  timepoints <- c("w0", "w2pre", "w2post", "w12")
  
  intermediate_results <- list()
  
  
  for(j in 1:length(timepoints)) {
    
    subset  <- dge_lists[[i]][,dge_lists[[i]]$samples$time == timepoints[j]]
      
    
    # set design matrix
    subset$design <- model.matrix( ~ subject + sets, data = subset$samples)
    
    # estimate dispersion 
    subset <- estimateDisp(subset, design = subset$design, robust = TRUE)
    
    # fit model 
    fit <- glmQLFit(subset, design = subset$design)
    
    ftest <- glmQLFTest(fit, coef = "setsmultiple")
    
    intermediate_results[[j]] <- data.frame(topTags(ftest, n = Inf), method = names(dge_lists[i]), coef = paste0(timepoints[j], "_setsmultiple")) %>%
      mutate(gene = rownames(topTags(ftest    , n = Inf)))
    
    message("Time-point ", j, " whitin dge-list ", i)
    
    
  }
  
  results_glmfit_paired[[i]] <- bind_rows(intermediate_results)
  
  
}



glm_estimates <- bind_rows(results_glmfit_paired)



glm_estimates %>%

  mutate(pthreshold = if_else(FDR > 0.05, "ns", "s")) %>%
  ggplot(aes(logFC, -log10(PValue), fill = pthreshold)) + 
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  theme_bw() +
  facet_grid(method ~ coef)

  
  
#### GlM for full model #######
glm_full_results <- list()

for(i in 1:length(dge_lists)) {
  
  
  # set design matrix
  
  dge_lists[[i]]$design <- model.matrix( ~ 0 + subject + time + time:sets, data = dge_lists[[i]]$samples) 
  
  # calculate dispersion 
  
  dge_lists[[i]] <- estimateDisp(dge_lists[[i]], design = dge_lists[[i]]$design)
  
  message("Dispersion calculated for method ", i)
  
  # Differential expression 
  
  # fit model 
  fit <- glmQLFit(dge_lists[[i]], design =  dge_lists[[i]]$design)
  

  
  coefs <- c("timew2pre", "timew2post",  "timew12", "timew0:setsmultiple",
  "timew2pre:setsmultiple",  "timew2post:setsmultiple", "timew12:setsmultiple") 
  
  intermediate_results <- list()
  
  for(j in 1:length(coefs)) { 
  
  
  ftest <- glmQLFTest(fit, coef = coefs[j])
  
  intermediate_results[[j]] <- data.frame(topTags(ftest, n = Inf), 
                                          method = names(dge_lists[i]), 
                                          coef = coefs[j]) %>%
    mutate(gene = rownames(topTags(ftest    , n = Inf)))
  
  
  }
  message("Model fitted/estimates retrieved for dgelist ", i)
  
  glm_full_results[[i]] <- bind_rows(intermediate_results)
  
}





#### In the MSD plot, sex might have been a significant contributor to biological variation. 
# To account for this, sex is added to the design matrix. 



for(i in 1:length(dge_lists)) {
  
  
  dge_lists[[i]]$design <- model.matrix(~ sex + time + time:sets, data = dge_lists[[i]]$samples)

}


########## Dream/variancePartition approach ######################

library('BiocParallel')

param <- SnowParam(2, "SOCK", progressbar=TRUE)
register(param)


dream_results <- list()



for(i in 1:length(dge_lists)){
  
  form <- ~ time + time:sets + (1|subject) 
  
  vobjDream <- voomWithDreamWeights(dge_lists[[i]], form, dge_lists[[i]]$samples)
  
  fitmm <- dream(vobjDream, form, dge_lists[[i]]$samples)
  
  estimates <- rbind(data.frame(topTable(fitmm, coef = 5, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w0_setsmultiple"), 
                     
                     data.frame(topTable(fitmm, coef = 6, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w2pre_setsmultiple"), 
                     
                     data.frame(topTable(fitmm, coef = 7, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w2post_setsmultiple"), 
                     
                     data.frame(topTable(fitmm, coef = 8, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w12_setsmultiple")) 
  
  ## Save results as a list ##
  file.path <- paste0("./data/derivedData/DE/dream_results_", names(dge_lists[i]), ".RDS")
  
  saveRDS(list(vobjDream = vobjDream, 
       fitmm = fitmm, 
       estimates = estimates), file = file.path)
  
  rm(vobjDream, fitmm)
  
  message("n methods estimated = ", i)
  
}


rsem_mm <- readRDS(file = "./data/derivedData/DE/dream_results_rsem.RDS")




estimates %>%
  mutate(pthreshold = if_else(adj.P.Val > 0.05, "ns", "s")) %>%
  ggplot(aes(logFC, -log10(P.Value), fill = pthreshold)) + 
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  theme_bw() +
  facet_grid(method ~ coef)



########## Acute phase estimates ###################

# Subset the dge_lists 

subset_dge <- list()


for(i in 1:length(dge_lists)) {
  
  subset_dge[[i]] <- dge_lists[[i]][,dge_lists[[i]]$samples$time %in% c("w2pre", "w2post")]
  
  names(subset_dge)[i] <- names(dge_lists[i])
  
}



# Estimate dispersion and set design matrix 


for(i in 1:length(subset_dge)) {
  
  # make new factor variables 
  
  subset_dge[[i]]$samples <- subset_dge[[i]]$samples %>%
    mutate(time = factor(time, levels = c("w2pre", "w2post"))) 
  
  # Add design matrix
  subset_dge[[i]]$design <- model.matrix(~ 0 + subject + time * sets, data = subset_dge[[i]]$samples) 

  # Estimate dispersion
  subset_dge[[i]] <- estimateDisp(subset_dge[[i]], design = subset_dge[[i]]$design, robust = TRUE)
  
  message("Dispersion estimated for n = ",i)
  
}



### Explorative plots of dispersion parameters #### 

par(mfrow = c(2, 3))

for(i in 1:length(subset_dge)) {
  
  plotBCV(subset_dge[[i]], main = names(subset_dge[i]))
  
}
par(mfrow = c(1,1))


## Dispersion parameters per method ##

data.frame(method = names(subset_dge), 
           dispersion = c(subset_dge[[1]]$common.dispersion, 
                          subset_dge[[2]]$common.dispersion,
                          subset_dge[[3]]$common.dispersion,
                          subset_dge[[4]]$common.dispersion,
                          subset_dge[[5]]$common.dispersion)) %>%
  mutate(BCV = sqrt(dispersion)) %>%
  print()



glm_acute_results <- list()

for(i in 1:length(subset_dge)) {
  
  # Fit the model
  fit <- glmQLFit(subset_dge[[i]], design = subset_dge[[i]]$design)
  
  
  colnames(fit)
  
  # Retrieve interesting coefficients 
  timew2post          <- glmQLFTest(fit, coef = 25)
  setsmultiple        <- glmQLFTest(fit, coef = 26)
  w2post_setsmultiple <- glmQLFTest(fit, coef = 27)
  
  
  # Get estimates
  coefs <- rbind(data.frame(topTags(timew2post, n = Inf), 
                            method = names(subset_dge[i]), 
                            coef = "timew2post"), 
                 data.frame(topTags(setsmultiple, n = Inf), 
                            method = names(subset_dge[i]), 
                            coef = "setsmultiple"),
                 data.frame(topTags(w2post_setsmultiple, n = Inf), 
                            method = names(subset_dge[i]), 
                            coef = "w2post_setsmultiple"))
  
  # Store a results data frame
  results <- list(fit = fit, 
                  estimates = coefs)
  
  
  # Store results list 
  glm_acute_results[[i]] <- results
  
  # Name with the method
  names(glm_acute_results)[i] <- names(subset_dge[i])
  
  message("Fit estimated for n = ",i)
  
  
}


## Compile all results 

glm_acute_estimates <- bind_rows(glm_acute_results[[1]]$estimates, 
                           glm_acute_results[[2]]$estimates, 
                           glm_acute_results[[3]]$estimates, 
                           glm_acute_results[[4]]$estimates, 
                           glm_acute_results[[5]]$estimates)

## Save results 
saveRDS(list(glm_acute_estimates = glm_acute_estimates, 
             glm_acute_results = glm_acute_results), file = "./data/derivedData/DE/glm_acute_results.RDS")




glm_acute_estimates %>%
  mutate(pthreshold = if_else(FDR > 0.05, "ns", "s")) %>%
  filter(coef == "w2post_setsmultiple") %>%
  
  ggplot(aes(logFC, -log10(PValue), fill = pthreshold)) + 
  geom_point(shape = 21, color = "black", alpha = 0.8) +
  theme_bw() +
  facet_grid(method ~ coef, scales = "free")



########## Dream variancePartition on acute data ###############





### This is where im at ###############

# TODO rewrite to fit acute phase data --- the below is copied from above.
library('BiocParallel')

param <- SnowParam(2, "SOCK", progressbar=TRUE)
register(param)


dream_results <- list()



for(i in 1:length(dge_lists)){
  
  form <- ~ time + time:sets + (1|subject) 
  
  vobjDream <- voomWithDreamWeights(dge_lists[[i]], form, dge_lists[[i]]$samples)
  
  fitmm <- dream(vobjDream, form, dge_lists[[i]]$samples)
  
  estimates <- rbind(data.frame(topTable(fitmm, coef = 5, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w0_setsmultiple"), 
                     
                     data.frame(topTable(fitmm, coef = 6, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w2pre_setsmultiple"), 
                     
                     data.frame(topTable(fitmm, coef = 7, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w2post_setsmultiple"), 
                     
                     data.frame(topTable(fitmm, coef = 8, number = Inf), 
                                method = names(dge_lists[i]), 
                                coef = "w12_setsmultiple")) 
  
  ## Save results as a list ##
  file.path <- paste0("./data/derivedData/DE/dream_results_", names(dge_lists[i]), ".RDS")
  
  saveRDS(list(vobjDream = vobjDream, 
               fitmm = fitmm, 
               estimates = estimates), file = file.path)
  
  rm(vobjDream, fitmm)
  
  message("n methods estimated = ", i)
  
}





# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/ ##






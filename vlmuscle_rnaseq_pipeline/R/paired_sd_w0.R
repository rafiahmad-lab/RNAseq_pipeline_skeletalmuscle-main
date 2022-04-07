###############################################################
#
## Paired analysis of reproducability across pipelines ##


# Load libraries
source("./R/lib_fun.R")
 

### Import count tables and store in a list called "quant"

rsem_rc <- read_delim("./seqDATA/rsem_geneid_expectedcount_matrix.txt", delim = "\t")

kallisto_rc <- read_delim("./seqDATA/kallisto_geneid_readcount_all_samples.tsv", delim = "\t")

salmon_rc <- read_delim("./seqDATA/salmon_geneid_readcount_all_samples.tsv", delim = "\t")

hisat2_rc <- read_delim("./seqDATA/hisat2-htseq_geneid_expectedcount_matrix1.txt", delim = "\t")

star_rc <- read_delim("./seqDATA/star-htseq_geneid_expectedcount_matrix1.txt", delim = "\t")

raw.files <- list(data.frame(rsem_rc),
                  data.frame(kallisto_rc),
                             data.frame(salmon_rc),
                                        data.frame(hisat2_rc),
                                                   data.frame(star_rc))

rm(list = c("rsem_rc", "kallisto_rc", "salmon_rc", "hisat2_rc", "star_rc"))


quant <- list(rsem = as.matrix(raw.files[[1]][,-1]), 
              kallisto = as.matrix(raw.files[[2]][,-1]), 
              salmon = as.matrix(raw.files[[3]][,-1]), 
              hisat = as.matrix(raw.files[[4]][,-1]), 
              star = as.matrix(raw.files[[5]][,-1]))



# Set rownames on quant files
for(i in 1:length(quant)){
  
  rownames(quant[[i]]) <- raw.files[[i]][,1]

  
}

rm(raw.files)



###### Subset count tables ######## 

# Count tables should only include Week 0 samples (these *pairs* will be 
# regarded as biologically equal).

# Extract sample information
# Load leg-condition data
legs_condition <- read.csv("./data/oneThreeSetLeg.csv", sep = ";") %>%
  pivot_longer(names_to = "sets", values_to = "leg", multiple:single) %>%
  print()

# Extract sample information from colnames in the count matrix
samp <- data.frame(id = colnames(quant[[1]])) %>%
  rowwise() %>%
  mutate(id = as.character(id), 
         subject = str_split_fixed(id, "w", 2)[1], 
         leg = substrRight(id, 1), 
         time = gsub(subject, "", id), 
         time = gsub(leg, "", time)) %>%
  inner_join(legs_condition) %>%
  filter(time == "w0") %>%
  print()
  
  
# Subset lists only to include relevant samples
for(i in 1:length(quant)){
  
  quant[[i]] <- quant[[i]][, colnames(quant[[i]]) %in% samp$id]
  
  
}

# Combine meta data and count tables

comp_data <- list(quant = quant, 
                  samp = samp)




# Subset genes for use in analysis

# From each method, remove genes with all zero counts

for(i in 1:length(comp_data$quant)) {
 
  comp_data$quant[[i]] <- comp_data$quant[[i]][rowSums(comp_data$quant[[i]]) > 0,]
  
}


# All operations downstream are to performed on a common set of genes, a vector
# of genes common to all methods is created.

common_genes <- Reduce(intersect, list(rownames(comp_data$quant[[1]]), 
                       rownames(comp_data$quant[[2]]),
                       rownames(comp_data$quant[[3]]),
                       rownames(comp_data$quant[[4]]),
                       rownames(comp_data$quant[[5]])))


# Subsetting all lists to common genes

for(i in 1:length(comp_data$quant)) {
  
  comp_data$quant[[i]] <- comp_data$quant[[i]][row.names(comp_data$quant[[i]]) %in% common_genes,]

}


# Makes a feature table

meta <- data.frame(gene = common_genes, 
                   include = TRUE)



# A list of "housekeeping" genes have been created by Eisenberg & Levanon (2013).
# This list is used to specify what genes to use for calibration in Teng et al (2016).

# Using the list, 

eisenberg <- read_excel("./data/eisenberg2013.xlsx") %>%
  
  print()

housekeeping <- unique(eisenberg$Refseq)

# Using biomaRt to convert to ensemble IDs

mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "www")

hk_refseq.ensembl <- getBM(attributes=c("refseq_mrna", "ensembl_gene_id", "hgnc_symbol"), 
      filters = "refseq_mrna", values = housekeeping, mart = mart)

# Create a logical vector of gene being a house keeping gene (identified in Esienberg 2013)
meta$hk <- meta$gene %in% hk_refseq.ensembl$ensembl_gene_id



# Store meta data set in comp_data-list
comp_data$meta <- meta


## Calculate average count per million 

aCPM <- list()

for(i in 1:length(comp_data$quant)) {
  
  aCPM[[i]] <- data.frame(gene = row.names(comp_data$quant[[i]]),
                     method = names(comp_data$quant)[i], 
                     aCPM = aveLogCPM(comp_data$quant[[i]], normalized.lib.sizes = FALSE))
  
  
}


aCPM <- bind_rows(aCPM) %>%
  pivot_wider(names_from = method, values_from = aCPM) %>%
  print()



comp_data$aCPM <- aCPM

rm(aCPM)





##### Paired deviation metrics expressed per average abundance ######


results_i <- list()

for(j in 1:length(comp_data$quant)) {
  
  
  # set up parallel processing 
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1) # not to overload cpu #### Change this if on server!
  registerDoSNOW(cl)
  
  # Progress bar (should work on linux also)
  iterations <- nrow(comp_data$meta) ### Change this when live-looping!
  pb <- txtProgressBar(max = iterations, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  
  ### foreach loop 
  results <- foreach(i = 1:nrow(comp_data$meta), 
                     .packages = c("tidyverse"), 
                     .options.snow = opts) %dopar% {
   
   
                       results_j <- data.frame(gene = rep(NA, 1), 
                                               method = rep(NA, 1),
                                               SD = rep(NA, 1),
                                               A = rep(NA, 1))  
                       
                       
   
temp <-  data.frame(count = comp_data$quant[[j]][comp_data$meta$gene[i], ]) %>%
    mutate(id = row.names(.)) %>%
    inner_join(comp_data$samp, by = "id") %>%
    dplyr::select(subject, leg, count) %>%
    pivot_wider(names_from = leg, values_from = count) %>%
    mutate(d = sqrt((log2(L + 0.5) - log2(R + 0.5))^2), 
           m = log2(((L + 0.5) + (R + 0.5))/2)) %>%
    group_by() %>%
    summarise(SD = mean(d), 
              A = mean(m)) %>%
    mutate(gene = comp_data$meta[i,1], 
           method = names(comp_data$quant)[j])%>%
   data.frame()
   

results_j[1,1] <- temp[1,3]
results_j[1,2] <- temp[1,4]
results_j[1,3] <- temp[1,1]
results_j[1,4] <- temp[1,2]

results_j

 } 
  
  close(pb)
  stopCluster(cl)
  
 results_i[[j]] <- bind_rows(results)
  
  
}




results_i <- bind_rows(results_i)


paired_sd_data <- results_i %>%
  mutate(Gene = rep(comp_data$meta$gene, 5)) %>%
  inner_join(comp_data$aCPM %>%
               mutate(Gene = gene) %>%
               pivot_longer(names_to = "method", values_to = "aCPM", cols = rsem:star) %>%
               dplyr::select(Gene, method, aCPM)) %>%
  inner_join(comp_data$meta %>%
               dplyr::select(Gene = gene, hk)) %>%
  mutate(method = factor(method, levels = c("rsem", "kallisto", "salmon", "star", "hisat"))) %>%
  print()






saveRDS(paired_sd_data, file = "./data/derivedData/paired_sd_data.RDS")



############ scraps #############

condInfo <- factor(comp_data$samp$leg)

repInfo <- factor(comp_data$samp$subject)

evaluationFeature <- rep(TRUE, nrow(comp_data$meta))
calibrationFeature <- simdata$meta$house & simdata$meta$chr == 'chr1'
unitReference <- 1

?signalCalibrate

sc_obj <- signalCalibrate(quantData = comp_data$quant, 
                          condInfo = condInfo, 
                          repInfo = repInfo, 
                          unitReference = 1, 
                          evaluationFeature = comp_data$meta$include, 
                          calibrationFeature = comp_data$meta$hk)
sc_obj

plotSD(sc_obj, ylim=c(0,3.9))

plotNE(sc_obj)




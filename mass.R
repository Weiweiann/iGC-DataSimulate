require(MASS)
library(data.table)


dataSimulation <- function(n_people=1, n_genes=100, cnv_proportion=0.3, cnv_threhold = 2.5){
  Genes <- vector()
  CNVs <- vector()

  for(i in 1:n_genes){
  # simulate a gene
    
  # setup correlation vector: random from 0~0.3 0.3~0.7 0.7~1
  correlation = c(runif(1,0,0.3), runif(1,0.3,0.7),runif(1,0.7,1))
  sample<-vector()
  
    for(cor in correlation){
      
      # cor_matirx is the matrix that fills with symmetry correlation value.
      cor_matrix = diag(nrow = 2,ncol = 2)
      cor_matrix[upper.tri(cor_matrix)] <- cor
      cor_matrix[lower.tri(cor_matrix)] <- cor
      
      # mu is the mean of simulated cnv and gene
      # the sd of "mvrnorm" is 1, thus we have to adjust the mean to set up the proportion of cnvs
      # increase cnv_default_mean and calculate the threshold until it equals or larger to cnv_threhold
   
      cnv_normal_possibility = 1-cnv_proportion 
      cnv_default_mean = 0
      while(qnorm(cnv_normal_possibility, mean = cnv_default_mean, sd = 1) < cnv_threhold){
        cnv_default_mean = cnv_default_mean + 0.05
      }
      
      mu = c(cnv_threhold, 5) # the second one is the mean value for gene expression, it can be any number since it does not effect the simulation
      
      # simulate n_genes with sample correlation
      sample_per_cor <- mvrnorm(n=n_people, mu = mu, Sigma = cor_matrix , empirical = T)
      colnames(sample_per_cor) = c("CNV", "Gene")
      sample <- rbind(sample,sample_per_cor)
    }
  
  # combine with all the simulate reuslts.
  CNVs <- cbind(CNVs, sample[,1])
  Genes <- cbind(Genes, sample[,2])
  }
  
  Gene_name = paste('Gene',seq(1:(n_genes*length(correlation))))
  Sample_name = paste('Sample', seq(1:n_people))
  
  CNVs <- data.table(Gene_name, CNVs)
  Genes <- data.table(Gene_name, Genes)
  
  colnames(CNVs) <- c('GENE', Sample_name)
  colnames(Genes) <- c('GENE', Sample_name)
  
  simData <- list()
  simData[[1]] <- CNVs
  simData[[2]] <- Genes
  simData[[3]] <- correlation
  return(simData)
}
 
## data[[1]]: CNV, data[[2]]: Gene
data <- dataSimulation(n_people = 1)

# iGC

require(iGC)



# cnv data should be -1, 0, 1 (for iGC)
cnv = data[[1]]
gene = data[[2]]
cnv[,2:30] = 1
cnv[,31:40] = 0
cnv[,41:51] = -1
result <- find_cna_driven_gene(gene_cna = cnv, gene_exp = gene)






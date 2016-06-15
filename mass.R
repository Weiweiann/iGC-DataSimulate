require(MASS)
library(data.table)

dataSimulation <- function(n_people=1, n_genes=100){
  Genes <- vector()
  CNVs <- vector()

  for(i in 1:n_people){
  # simulate a sample 
    
  
  correlation = c(runif(1,0,0.3), runif(1,0.3,0.7),runif(1,0.7,1))
  n_genes = 100
  sample<-vector()
  
  for(cor in correlation){
    cor_matrix = diag(nrow = 2,ncol = 2)
    cor_matrix[upper.tri(cor_matrix)] <- cor
    cor_matrix[lower.tri(cor_matrix)] <- cor
    sample_per_cor <- mvrnorm(n=n_genes, mu = c(5,10), Sigma = cor_matrix , empirical = T) # 只能設定mean, 
    colnames(sample_per_cor) = c("CNV", "Gene")
    sample <- rbind(sample,sample_per_cor)
  }
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
  
  return(simData)
}
 
## data[[1]]: CNV, data[[2]]: Gene
data <- dataSimulation(n_people = 1)

# iGC

require(iGC)



# cnv data should be -1, 0, 1
cnv = data[[1]]
gene = data[[2]]
cnv[,2:30] = 1
cnv[,31:40] = 0
cnv[,41:51] = -1
result <- find_cna_driven_gene(gene_cna = cnv, gene_exp = gene)






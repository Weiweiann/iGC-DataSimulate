require(MASS)
library(data.table)

n_genes = 100
n_people = 100
cnv_proportion = 0.3
cnv_threhold = 2.5
correlation = correlation_simualate(n_genes = n_genes)


dataSimulation <- function(n_people, cnv_proportion, cnv_threhold , correlation){
  Genes <- vector()
  CNVs <- vector()

  n_genes = length(correlation)
  
  for(i in 1:n_genes){
    # simulate a gene
    cor = correlation[i]
   
        
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
  
    # the second one is the mean value for gene expression, 
    # it can be any number since it does not effect the simulation
    mu = c(cnv_threhold, 5) 
    
    # simulate n_people with sample correlation
    one_gene_simulate <- t(mvrnorm(n=n_people, mu = mu, Sigma = cor_matrix , empirical = T))
    row.names(one_gene_simulate) = c("CNV", "GENE")
    
    # combine with all the simulate reuslts.
    # row: sample, col: gene
    CNVs <- rbind(CNVs, one_gene_simulate[1,])
    Genes <- rbind(Genes, one_gene_simulate[2,])
  }
  
  
  Gene_name = paste('GENE',seq(1:n_genes),sep = "")
  Sample_name = paste('Sample', seq(1:n_people), sep="")
  
  CNVs <- data.table(GENE = Gene_name,CNVs)
  colnames(CNVs) <- c('GENE',Sample_name)
  
  Genes <- data.table(GENE = Gene_name,Genes)
  colnames(Genes) <- c('GENE',Sample_name)

  
  simData <- list()
  simData[[1]] <- CNVs
  simData[[2]] <- Genes
  simData[[3]] <- correlation
  
  cat("Number of genes:", n_genes, '\n')
  cat("Number of people", n_people, '\n')
  return(simData)
}
 
## data[[1]]: CNV, data[[2]]: Gene
data <- dataSimulation(n_people = n_people, cnv_proportion = cnv_proportion, 
                       cnv_threhold = cnv_threhold, correlation = correlation)

# iGC

require(iGC)



# cnv data should be -1, 0, 1 (for iGC)
cnv = data[[1]]
gene = data[[2]]

# for iGC
cnv_threhold = 2.5
temp = which(cnv < cnv_threhold, arr.ind = T)
# 隨便亂設0跟1
cnv[temp] = 0
cnv[temp[(dim(temp)[1]/3):dim(temp)[1],]] = -1
cnv[which(cnv >= cnv_threhold, arr.ind = T)] = 1
cnv[,1] =  paste('GENE',seq(1:n_genes),sep="")

igc_result <- find_cna_driven_gene(gene_cna = cnv, gene_exp = gene)


# simulate answer 
answer = data.frame(Gene=  paste('GENE',seq(1:n_genes),sep = ""), Cor =data[[3]])
answer$Cor[which(answer$Cor > 0.3)] = 'true' 
answer$Cor[which(answer$Cor <= 0.3)] = 'false'
colnames(answer) = c('Gene', 'Answer')

# pvalue < 0.5 as True, else as False
igc_T_gene <- igc_result$both$GENE[which(igc_result$both$gain_p_value < 0.05)]
igc_F_gene <- igc_result$both$GENE[which(igc_result$both$gain_p_value >= 0.05)]
igc_answer = rbind( data.frame(Gene = igc_T_gene, Answer= 'true'), data.frame(Gene = igc_F_gene, Answer='false'))
igc_answer = igc_answer[order(igc_answer$Gene),]
answer_sort = answer[order(answer$Gene),] # 因為沒辦法完全按照數字來排列，所以只好統一用奇怪的方式排列了~"~
answer_sort$Answer = factor(answer_sort$Answer)
xtab = table(answer_sort$Answer, igc_answer$Answer)

library(caret)
cat("iGC_sen: ",sensitivity(answer_sort$Answer, igc_answer$Answer), "iGC_spe: ", specificity(answer_sort$Answer, igc_answer$Answer),'\n')


source('./sim_simulate.R')

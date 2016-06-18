## for sim, simulate data can not use data.table, otherwise there will be error.



dataSimulation <- function(n_people=100, n_genes=300, cnv_proportion=0.3, cnv_threhold = 2.5){
  Genes <- vector()
  CNVs <- vector()
  
  # setup correlation vector: random from 0~0.3 0.3~0.7 0.7~1
  if (n_genes %% 3 ==1){
    correlation = c(runif(round(n_genes/3),0,0.3), runif(round(n_genes/3),0.3,0.7),runif(round(n_genes/3)+1,0.7,1))
  }else if(n_genes %% 3 ==2 ){
    correlation = c(runif(round(n_genes/3),0,0.3), runif(round(n_genes/3),0.3,0.7),runif(round(n_genes/3)+2,0.7,1))
  }else{
    correlation = c(runif(round(n_genes/3),0,0.3), runif(round(n_genes/3),0.3,0.7),runif(round(n_genes/3),0.7,1))
  }
  
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
  
  CNVs <- data.frame(GENE = Gene_name,CNVs)
  colnames(CNVs) <- c('GENE',Sample_name)
  
  Genes <- data.frame(GENE = Gene_name,Genes)
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
data <- dataSimulation(n_people = 100, n_genes = 300)


cnv = data[[1]]
gene = data[[2]]

fake_ann <- chrom.table[1:dim(cnv)[1], c('chr','start','end')]
cnv_sim <- data.frame(gene = cnv$GENE,fake_ann, cnv[,-1])
gene_sim <- data.frame(gene = cnv$GENE, fake_ann, gene[,-1])


assemble.data(dep.data = cnv_sim,
              indep.data = gene_sim,
              dep.ann = colnames(cnv_sim)[1:4],
              indep.ann = colnames(gene_sim)[1:4],
              dep.id = "gene",
              dep.chr = "chr",
              dep.pos = "start",
              indep.id = "gene",
              indep.chr = "chr",
              indep.pos = "start",
              overwrite = TRUE,
              run.name = "simulate")

sim_samples = paste('Sample',seq(1:100),sep="")

integrated.analysis(samples = sim_samples,
                    input.regions = c(1,2,3,4),
                    zscores = TRUE,
                    method = "full",
                    run.name = "simulate")






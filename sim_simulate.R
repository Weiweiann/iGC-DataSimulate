## for sim, simulate data can not use data.table, otherwise there will be error.
data("chrom.table")

dataSimulation <- function(n_people,cnv_proportion, cnv_threhold, correlation){
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
  
  CNVs <- data.frame(GENE = Gene_name,CNVs) # data.table  -> data.frame
  colnames(CNVs) <- c('GENE',Sample_name)
  
  Genes <- data.frame(GENE = Gene_name,Genes) # data.table  -> data.frame
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
data <- dataSimulation(n_people = n_people, cnv_proportion = cnv_proportion, cnv_threhold = cnv_threhold, correlation = correlation)

cnv = data[[1]]
gene = data[[2]]

# generate a fake annotation file from data("chrom.table")
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

sim_samples = paste('Sample',seq(1:n_people),sep="")

# fake annotation containes chr 1 to chr 6.
integrated.analysis(samples = sim_samples,
                    zscores = TRUE,
                    method = "full",
                    run.name = "simulate")

# the p-value results are located in ./simulate/full/gpvals.pat.c.*
pvalues <-vector()
for(file in list.files(path = './simulate/full/', pattern = 'gpvals.pat.c*')){
  n =1 
  variable <- paste('pvalue.',n,sep="")
  path = paste('./simulate/full/',file,sep="")
  pvalues <- c(pvalues,dget(path))
}


sim_result = data.frame(Gene = paste('Gene',seq(1:n_genes),sep=''), 
                        pval = pvalues)
sim_result.2 <- data.frame(Gene=sim_result$Gene, Answer='empty',stringsAsFactors = F)
sim_result.2$Answer[which(sim_result$pval < 0.05)] = 'true'
sim_result.2$Answer[which(sim_result$pval >= 0.05)] = 'false'
answer$Answer = factor(answer$Answer)
sim_result.2$Answer = factor(sim_result.2$Answer)
sim_tab <- table(answer$Answer, sim_result.2$Answer)
cat("sim_sen: ",sensitivity(answer$Answer, sim_result.2$Answer), "sim_spe: ", specificity(answer$Answer, sim_result.2$Answer),'\n')





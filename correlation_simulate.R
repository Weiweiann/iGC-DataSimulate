
correlation_simualate <- function(n_genes){
  # setup correlation vector: random from 0~0.3 0.3~0.7 0.7~1
  if (n_genes %% 3 ==1){
    correlation = c(runif(round(n_genes/3),0,0.3), runif(round(n_genes/3),0.3,0.7),runif(round(n_genes/3)+1,0.7,1))
  }else if(n_genes %% 3 ==2 ){
    correlation = c(runif(round(n_genes/3),0,0.3), runif(round(n_genes/3),0.3,0.7),runif(round(n_genes/3)+2,0.7,1))
  }else{
    correlation = c(runif(round(n_genes/3),0,0.3), runif(round(n_genes/3),0.3,0.7),runif(round(n_genes/3),0.7,1))
  }
  return(correlation)
}
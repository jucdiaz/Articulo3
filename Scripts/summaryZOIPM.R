summary.ZOIPM<-function(mod){
  
  estimate <- c(mod$Fixed_Parameters.mu,mod$Fixed_Parameters.sigma
                ,mod$Fixed_Parameters.p0,mod$Fixed_Parameters.p1,mod$Parameters.randoms[,1])
  se       <- sqrt(diag(solve(mod$HM)))
  zvalue   <- estimate / se
  pvalue   <- 2 * pnorm(abs(zvalue), lower.tail=F)
  res      <- cbind(estimate=estimate, se=se, zvalue=zvalue, pvalue=pvalue)
  colnames(res) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
  res      <- as.data.frame(res)
  
  a <- 1:length(mod$Fixed_Parameters.mu)
  b <-
    (length(mod$Fixed_Parameters.mu) + 1):(length(mod$Fixed_Parameters.mu) +
                                             length(mod$Fixed_Parameters.sigma))
  c <-
    (length(mod$Fixed_Parameters.mu) + length(mod$Fixed_Parameters.sigma) +
       1):(
         length(mod$Fixed_Parameters.mu) + length(mod$Fixed_Parameters.sigma) + length(mod$Fixed_Parameters.p0)
       )
  d <-
    (
      length(mod$Fixed_Parameters.mu) + length(mod$Fixed_Parameters.sigma) + length(mod$Fixed_Parameters.p0) +
        1
    ):(
      length(mod$Fixed_Parameters.mu) + length(mod$Fixed_Parameters.sigma) + length(mod$Fixed_Parameters.p0) +
        length(mod$Fixed_Parameters.p1)
    )
  e <-
    (
      length(mod$Fixed_Parameters.mu) + length(mod$Fixed_Parameters.sigma) + length(mod$Fixed_Parameters.p0) +
        length(mod$Fixed_Parameters.p1) + 1
    ):(
      length(mod$Fixed_Parameters.mu) + length(mod$Fixed_Parameters.sigma) + length(mod$Fixed_Parameters.p0) +
        length(mod$Fixed_Parameters.p1) + length(mod$Parameters.randoms[, 1])
    )
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            link[1], "(mu) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res[a,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            link[2], "(sigma) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res[b,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            link[3], "(p0) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res[c,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat(paste("Fixed effects for ",
            link[4], "(p1) \n", sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res[d,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat("---------------------------------------------------------------\n")
  cat(paste("Random effects for mu and sigma \n",sep=''))
  cat("---------------------------------------------------------------\n")
  printCoefmat(res[e,], P.value=TRUE, has.Pvalue=TRUE)
  cat("---------------------------------------------------------------\n")
  cat("---------------------------------------------------------------\n")
  
}

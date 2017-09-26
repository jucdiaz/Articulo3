
print.ZOIPM<-function(mod){
  
  cat("Call:\n")
  print(mod$call)
  cat("\n Results: \n")
  cat("\n Estimated fixed coefficients for h(mu): \n")
  print(mod$Fixed_Parameters.mu)
  cat("\n Estimated fixed coefficients for h(sigma): \n")
  print(mod$Fixed_Parameters.sigma)
  cat("\n Estimated fixed coefficients for h(p0): \n")
  print(mod$Fixed_Parameters.p0)
  cat("\n Estimated fixed coefficients for h(p1): \n")
  print(mod$Fixed_Parameters.p1)
  cat("\n Estimated random coefficients for h(mu) and h(sigma) \n")
  print(mod$Parameters.randoms)
  cat("\n message \n")
  print(mod$message)
  cat("\n time \n")
  print(mod$Time)
  cat("\n iterations \n")
  print(mod$num.iter)
  cat("\n Log-likelihood \n")
  print(mod$logverosimilitud)
  
}

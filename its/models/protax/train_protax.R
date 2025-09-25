protax_dir <- Sys.getenv("PROTAX")
num_taxlevels <- as.integer(Sys.getenv("NUM_TAXLEVELS"))

source(file.path(protax_dir, "amcmc.rcode.txt"))
library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)

num.params=1+4
ind <- 1001:2000

pdf(file = "mcmc_trace.pdf", width = 10, height = 5)

for (level in seq.int(num_taxlevels)[-1]) {

  dat <- read.xdata(sprintf("train%i.xdat", level))
  pp1 <- adaptiveMCMC(dat, num.params, 10000, 2000, 1000, rseed=1, info=1)
  traceplot.all(pp1, ind, num.levels = 1, title = paste0("L", level, "a"))

  k <- which.max(pp1$postli[ind])
  write.postparams(pp1, paste0("mcmc", level, "a"), ind[k])
  
  initstate=initialize.adaptation(pp1$params[2000,])
  pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
  traceplot.all(pp1,ind,num.levels=1, title=paste0("L", level, "b"))
  
  k=which.max(pp1$postli[ind])
  write.postparams(pp1,paste0("mcmc", level),ind[k])
}

dev.off()

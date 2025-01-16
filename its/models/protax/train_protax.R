protax_dir <- Sys.getenv("PROTAX")
num_taxlevels <- as.integer(Sys.getenv("NUM_TAXLEVELS"))

source(file.path(protax_dir, "amcmc.rcode_noweight.txt"))
library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)

num.params=1+4
ind <- 1001:2000

pdf(file = "mcmc_trace.pdf", width = 10, height = 5)

for (level in seq.int(num_taxlevels)) {

  dat <- read.xdata(sprintf("train%i.scxdat", level))
  pp <- adaptiveMCMC(dat, num.params, 10000, 2000, 1000, rseed=1, info=1)
  traceplot.all(pp, ind, num.levels = 1, title = paste0("L", level))

  k <- which.max(pp$postli[ind])
  write.postparams(pp, paste0("mcmc", level), ind[k])
}

dev.off()

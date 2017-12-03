library(MCMCglmm)
library(phytools)
library(rjags)
load.module("glm")

df <- read.csv("noCLSE.csv", header = T)
df.strat <- df[df$Strat.Length != 0,]
phy <- force.ultrametric(collapse.singles(read.tree("Species_tree.txt")))
phy$node.label <- NULL
Ap <- inverseA(phy, nodes = "TIPS")$Ainv

run.mcmc <- function(df, fn) {
  n.samples <- nrow(df)
  n.species <- length(levels(df$Species))
  data.species <- as.integer(df$Species)
  data.pop <- rep(NA, n.samples)
  n.populations <- rep(NA, length(levels(df$Species)))
  for (sp in 1:n.species) {
    i <- df$Species == levels(df$Species)[sp]
    x <- as.factor(as.character(df$Population[i]))
    n.populations[as.integer(sp)] <- length(levels(x))
    data.pop[i] <- as.integer(x)
  }
  data.sl <- as.integer(as.factor(df$Strat.Length))
  data.st <- as.integer(df$Strat.Temp)
  data.it <- as.integer(df$Inc.Temp)
  data.germinated <- df$Total.Germ
  data.viable <- df$Total.Viable

  permute <- sapply(1:n.species, function(i) which(rownames(Ap) == levels(df$Species)[i]))
  Ainv <- matrix(NA, nrow = n.species, ncol = n.species)
  for (i in 1:n.species) {
    for (j in 1:n.species) {
      Ainv[i,j] <- Ap[permute[i], permute[j]]
    }
  }

  jags <- jags.model("glmm.bugs", data = list(n.species = n.species, data.species = data.species, data.pop = data.pop, n.populations = n.populations, data.sl = data.sl, data.st = data.st, data.it = data.it, data.germinated = data.germinated, data.viable = data.viable, Ainv = Ainv), n.chains = 4, n.adapt = 5000)

  update(jags, 5000)
  mcmc <- coda.samples(jags, c("p", names(jags[["state"]]()[[1]])), 50000, thin = 1)
  save(mcmc, file=fn)
}

run.mcmc(df, "glmm.coda")
run.mcmc(df.strat, "glmm.strat.coda")

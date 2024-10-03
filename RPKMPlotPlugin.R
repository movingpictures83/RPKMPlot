
###################################################
### code chunk number 1: cqn.Rnw:7-8
###################################################
options(width=70)


###################################################
### code chunk number 2: load
###################################################
library(cqn)
library(scales)
library(edgeR)

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix = prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }
}

run <- function() {}

output <- function(outputfile) {

montgomery.subset <- read.table(paste(pfix, parameters["dataset", 2], sep="/"), sep=",")
sizeFactors.subset <- as.double(read.csv(paste(pfix, parameters["sizefactors", 2], sep="/")))
uCovar <- read.table(paste(pfix, parameters["covar", 2], sep="/"), sep=",")
cqn.subset <- readRDS(paste(pfix, parameters["cqnsubset", 2], sep="/"))
grp1 <- readLines(paste(pfix, parameters["grp1", 2], sep="/"))
grp2 <- readLines(paste(pfix, parameters["grp2", 2], sep="/"))


### code chunk number 9: normalizedvalues
RPKM.cqn <- cqn.subset$y + cqn.subset$offset

### code chunk number 10: rpkmvalues
RPM <- sweep(log2(montgomery.subset + 1), 2, log2(sizeFactors.subset/10^6))
RPKM.std <- sweep(RPM, 1, log2(uCovar$length / 10^3))

### code chunk number 12: whGenes
whGenes <- which(rowMeans(RPKM.std) >= 2 & uCovar$length >= 100)
M.std <- rowMeans(RPKM.std[whGenes, grp1]) - rowMeans(RPKM.std[whGenes, grp2])
A.std <- rowMeans(RPKM.std[whGenes,])
M.cqn <- rowMeans(RPKM.cqn[whGenes, grp1]) - rowMeans(RPKM.cqn[whGenes, grp2])
A.cqn <- rowMeans(RPKM.cqn[whGenes,])

saveRDS(M.std, paste(outputfile,"M","std","rds",sep="."))
saveRDS(M.cqn, paste(outputfile,"M","cqn","rds",sep="."))
saveRDS(A.std, paste(outputfile,"A","std","rds",sep="."))
saveRDS(A.cqn, paste(outputfile,"A","cqn","rds",sep="."))
saveRDS(whGenes, paste(outputfile,"whGenes","rds",sep="."))

pdf(paste(outputfile,"pdf",sep="."))
### code chunk number 13: maplots
plot(A.std, M.std, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
     main = "Standard RPKM", ylim = c(-4,4), xlim = c(0,12),
     col = alpha("black", 0.25))
plot(A.cqn, M.cqn, cex = 0.5, pch = 16, xlab = "A", ylab = "M",
     main = "CQN normalized RPKM", ylim = c(-4,4), xlim = c(0,12),
     col = alpha("black", 0.25))

}

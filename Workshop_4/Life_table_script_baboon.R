######
# Winter school script baboon
######

#install.packages("MortalityLaws")
#install.packages("popbio")

#Load baboon life table
tab <- read.table("life_table_baboon.csv",sep=',',dec='.',header=TRUE)
head(tab)

#plot the different life table column for each age class
par(mfrow=c(3,1))
plot(tab$qx~tab$Age,ylab="qx (mortality)", xlab="Age (years)",pch=3,cex.lab=1.5)
plot(tab$lx~tab$Age,ylab="lx (survivorship)", xlab="Age (years)",pch=3,cex.lab=1.5)
plot(tab$dx~tab$Age,ylab="dx (age at death distribution)", xlab="Age (years)",pch=3,cex.lab=1.5)


#longevity measures computation

#life expectancies
tab$ex <- numeric(length(tab$Age))
for (i in 1: length(tab$ex)){
  tab$ex[i] <- sum(tab$lx[i:length(tab$lx)])/sum(tab$lx[i]) - 0.5
}
par(mfrow=c(1,1))
plot(tab$ex~tab$Age,ylab="ex (life expectancy", xlab="Age (years)",pch=3,cex.lab=1.5)

#median lifespan
max(tab$Age[(tab$lx-0.5)>0])

#maximum lifespan
max(tab$Age)


#mortality smoothing
#load package
library(MortalityLaws)

#fitting gompertz model to qx data
mod_gompertz <- MortalityLaw(x = tab$Age,
                          qx = tab$qx,law="gompertz")
#check siler model fit (much better)
plot(mod_gompertz)

#fitting siler model to qx data
mod_siler <- MortalityLaw(x = tab$Age,
                          qx = tab$qx,law="siler")
#check siler model fit (much better)
plot(mod_siler)

#compute the life table from the different siler parameters of the fitted model
tab_fit <- LawTable(x = tab$Age, par = mod_siler$coefficients, law = "siler")
tab_fit <- tab_fit[[1]][,c("x","dx","qx","lx","ex")]
tab_fit$lx <- tab_fit$lx/100000
tab_fit$dx <- tab_fit$dx/100000


#plot the different life table column for each age class with predictions
par(mfrow=c(3,1))
plot(tab$qx~tab$Age,ylab="qx (mortality)", xlab="Age (days)",pch=3,cex.lab=1.5)
lines(tab_fit$qx~tab_fit$x)
plot(tab$lx~tab$Age,ylab="lx (survivorship)", xlab="Age (days)",pch=3,cex.lab=1.5)
lines(tab_fit$lx~tab_fit$x)
plot(tab$dx~tab$Age,ylab="dx (age at death distribution)", xlab="Age (days)",pch=3,cex.lab=1.5)
lines(tab_fit$dx~tab_fit$x)

#plotting mx
par(mfrow=c(1,1))
plot(tab$mx~tab$Age,ylab="mx (fertility)", xlab="Age (year)",pch=3,cex.lab=1.5)

#plotting mx lx
par(mfrow=c(1,1))
plot(tab$mx*tab$lx~tab$Age,ylab="lxmx", xlab="Age (year)",pch=3,cex.lab=1.5)

#R0 calculation
R0 <- sum(tab$mx*tab$lx)

#G calculation
G <- sum(tab$Age*tab$mx*tab$lx)/sum(tab$mx*tab$lx)


#load prebreeding matrix fro this example
mat <- read.table('matrix_baboon.csv',sep=',',dec='.')
#matrix reduced last age with zero reproduction removed
matreduced <- mat[1:24,1:24]

#load pobio package
library(popbio)

#compute lambda
lambda(mat)

#compute stage age distribution
age_dist <- stable.stage(mat)

#compute reproductive value (scaled for the first value to be 1)
rep_val <- reproductive.value(matreduced)

#compute sensitivities
sensi <- sensitivity(matreduced)

#compute sensitivities
elas <- elasticity(matreduced)
sum(elas)


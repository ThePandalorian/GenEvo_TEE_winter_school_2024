######
# Winter school script Ageing dolphins
######

#Load age at death data
surv <- read.table("Surv_dolphin.csv",sep=',',dec='.',header=TRUE)
head(surv)



#load maternal age at birth data
repro <- read.table("Repro_dolphin.csv",sep=',',dec='.',header=TRUE)
head(repro)

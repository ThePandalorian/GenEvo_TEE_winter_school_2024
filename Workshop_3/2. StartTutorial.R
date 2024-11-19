rm(list=ls())

options(timeout = 1000) ##Necessary for slow connections#
install.packages("tidyverse")
install.packages("learnr")
rmarkdown::run("Workshop_3/ChangingEnv_tutorial.Rmd")

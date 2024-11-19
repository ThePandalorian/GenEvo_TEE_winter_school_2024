##Script for running the agent based model###
##If you have not done so, install the package tidyverse (install.packages("tidyverse"))
###Set working directory and source script_model#####
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Script_model.R")

####_________Assignment 1: Different response modes_________####
#run each scenario at least twice, using the following seeds
#in case you want to run a scenario more than two times, use a new seed
seed1 = 102
seed2 = 128943
#
###Scenario 1: Fast environmental change (R = 1), highly predictable cues (P = 1)####
scenario1A <- run.model(R = 1, P = 1, seed = seed1, Max_gen = 2000)
plot.phenotypes(scenario1A)
plot.genes(scenario1A)

scenario1B <- run.model(R = 1, P = 1, seed = seed2, Max_gen = 2000)
plot.phenotypes(scenario1B)
plot.genes(scenario1B)


###Scenario 2: Intermediate environmental change (R = 315), highly predictable cues (P = 1)####
scenario2A <- run.model(R = 315, P = 1, seed = seed1, Max_gen = 2000)
plot.phenotypes(scenario2A)
plot.genes(scenario2A)

scenario2B <- run.model(R = 315, P = 1, seed = seed2, Max_gen = 2000)
plot.phenotypes(scenario2B)
plot.genes(scenario2B)

###Scenario 3: Slow environmental change (R = 10.000), highly predictable cues (P = 1)####
scenario3A <- run.model(R = 10000, P = 1, seed = seed1, Max_gen = 20000, plot_gen = 15000)
plot.phenotypes(scenario3A)
plot.genes(data = scenario3A)
scenario3B <- run.model(R = 10000, P = 1, seed = seed2, Max_gen = 20000, plot_gen = 15000)
plot.phenotypes(scenario3B)
plot.genes(data = scenario3B)



###Scenario 4: Fast environmental change (R = 1), somewhat predictable cues (P = 0.4)####
scenario4A <- run.model(R = 1, P = 0.4, seed = seed1, Max_gen = 2000)
plot.phenotypes(scenario4A)
plot.genes(scenario4A)
scenario4B <- run.model(R = 1, P = 0.4, seed = seed2, Max_gen = 2000)
plot.phenotypes(scenario4B)
plot.genes(scenario4B)


###Scenario 5: Intermediate environmental change (R = 315), somewhat predictable cues (P = 0.4)####
scenario5A <- run.model(R = 315, P = 0.4, seed = seed1, Max_gen = 10000)
plot.phenotypes(scenario5A)
plot.genes(scenario5A)
scenario5B <- run.model(R = 315, P = 0.4, seed = seed2, Max_gen = 10000)
plot.phenotypes(scenario5B)
plot.genes(scenario5B)

###Scenario 6: Slow environmental change (R = 10.000), somewhat predictable cues (P = 0.4)####
scenario6A <- run.model(R = 10000, P = 0.4, seed = seed1, Max_gen = 20000)
plot.phenotypes(scenario6A)
plot.genes(scenario6A)
#In case your computer is very slow, you can skip the second calculation for this scenario
scenario6B <- run.model(R = 10000, P = 0.4, seed = seed2, Max_gen = 20000)
plot.phenotypes(scenario6B)
plot.genes(scenario6B)


###Scenario 7: Fast environmental change (R = 1), unpredictable cues (P = 0)####
scenario7A <- run.model(R = 1, P = 0, seed = seed1, Max_gen = 2000)
plot.phenotypes(scenario7A)
plot.genes(scenario7A)
scenario7B <- run.model(R = 1, P = 0, seed = seed2, Max_gen = 2000)
plot.phenotypes(scenario7B)
plot.genes(scenario7B)

###Scenario 8: Intermediate environmental change (R = 315), unpredictable cues (P = 0)####
scenario8A <- run.model(R = 315, P = 0, seed = seed1, Max_gen = 10000)
plot.phenotypes(scenario8A)
plot.genes(scenario8A)

scenario8B <- run.model(R = 315, P = 0, seed = seed2, Max_gen = 10000)
plot.phenotypes(scenario8B)
plot.genes(scenario8B)

###Scenario 9: Slow environmental change (R = 10.000), unpredictable cues (P = 0)####
scenario9A <- run.model(R = 10000, P = 0, seed = seed1, Max_gen = 20000)
plot.phenotypes(scenario9A)
plot.genes(scenario9A)
#In case your computer is very slow, you can skip the second calculation for this scenario
scenario9B <- run.model(R = 10000, P = 0, seed = seed2, Max_gen = 20000)
plot.phenotypes(scenario9B)
plot.genes(scenario9B)



#####_________Assignment 2 (Advanced): Changing environmental conditions_________#####

##First, look at each scenario separately, starting with an initial population with plenty of variation in the 7 traits
R = 10
P3 = 0.3
P2 = 0.2
P1 = 0.1
qval = 2.3

Scenario_P3 <- run.model(P = P3, R = R, seed = seed2, Max_gen = 2000, plot_gen = 2000, 
                         rel_fitness = F, 
                         q = qval)
plot.phenotypes(Scenario_P3)
plot.genes(Scenario_P3, maxtime = 500)
plot.pop(Scenario_P3)

Scenario_P2 <- run.model(P = P2, R = R, seed = seed2, Max_gen = 2000, plot_gen = 2000, 
                         rel_fitness = F,
                         q = qval)
plot.phenotypes(Scenario_P2)
plot.genes(Scenario_P2, maxtime = 500)
plot.pop(Scenario_P2)

Scenario_P1 <- run.model(P = P1, R = R, seed = seed2, Max_gen = 2000, plot_gen = 2000,
                         rel_fitness = F,
                         q = qval)
plot.phenotypes(Scenario_P1)
plot.genes(Scenario_P1, maxtime = 500)
plot.pop(Scenario_P1)

#Now, use the population at the end of the previous run as the initial population of the new run
#thereby simulation how a population will respond to changing conditions

Scenario_P2_prevpop <- run.model(P = P2, R = R, seed = seed2, Max_gen = 2000, plot_gen = 2000, 
                                 rel_fitness = F,
                                 q = qval, inds = Scenario_P3$inds)
plot.phenotypes(Scenario_P2_prevpop)
plot.genes(Scenario_P2_prevpop, maxtime = 500)
plot.pop(Scenario_P2_prevpop)


Scenario_P1_prevpop <- run.model(P = P1, R = R, seed = seed2, Max_gen = 2000, plot_gen = 2000, 
                                 rel_fitness = F, 
                                 q = qval, inds = Scenario_P2_prevpop$inds)

plot.pop(Scenario_P1_prevpop)
plot.genes(Scenario_P1_prevpop)
#plot.phenotypes(Scenario_P1_prevpop)



#####_________Additional assignment (Advanced): Costs of plasticity_________######
#The standard costs of plasticity are 0.02 for being plastic, and 0.01 for each time you switch your phenotype during development
##Program yourself
mxgen = 2000 #Number of generations to run the simulation for. If this is too short, you may not have reached the optimal phenotype. However, if you run for many generations, the simulations take very long
c_plast = 0.02
c_devp = 0.01
R = 1
P = 1
Check_plast <- run.model(R = R, P = P, seed = seed1, Max_gen = mxgen,
                         c_plast = c_plast, 
                         c_devp = c_devp)


#####_________Additional assignment (Advanced): Speed of evolution_________######
#The mutation probability can be changed by changing parameter mu_mut,
#the mutational step can be changed by changing parameter sd_mu
##Program yourself
mxgen = 2000 #Number of generations to run the simulation for. If this is too short, you may not have reached the optimal phenotype. However, if you run for many generations, the simulations take very long

mu_mut = 0.001
sd_mut = 0.05
R = 1
P = 1
Check_mut <- run.model(R = R, P = P, seed = seed1, Max_gen = mxgen,
                       mu_mut = mu_mut, sd_mut = sd_mut)


#####_________Additional assignment (Advanced): Transitions in strategies_________#####
##Program yourself
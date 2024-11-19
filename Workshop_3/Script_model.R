###This script contains all functions of the model and plotting the output####
rm(list=ls())
library("tidyverse")
###Multiplot function####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#####Functions model####
##Environment
env_function_p1 <- function(time, L, A, R) {
  #A and B are scaling coefficients, scaling importance of predictable and unpredictable part
  #A+B = 1
  B = 1 - A
  #L = number of timesteps per gen (lifespan)
  #R = (relative) timescale of env. variation
  env = A * sin(2 * pi * time / (L * R)) + B * runif(1)
  #env = sin(time * 2 * pi)
  return(env)
}
env_function <- Vectorize(env_function_p1)

#Environmental cue#
env_cues_p1 <- function(P, env) {
  mu = P * env
  sigma = (1 - P)/3
  C = rnorm(1, mu, sigma)
  C = max(-1, min(C, 1))
  return(C)
}
env_cues <- Vectorize(env_cues_p1)


##Create individuals##
make.inds <- function(id = NaN, n = NaN, age = NaN, gene_bh = NaN){
  inds <- data.frame(id = id, age = 0, 
                     gene_bh = gene_bh,
                     gene_devp = runif(n, 0, 1),
                     gene_pls = runif(n, 0, 1),
                     gene_I0 = runif(n, -1, 1),
                     gene_I0_alt = runif(n, -1, 1),
                     gene_b = runif(n, -2, 2),
                     gene_b_alt = runif(n, -2, 2),
                     phenotype = NaN,
                     mismatch = 0,
                     fitness = NaN,
                     pheno_update = 0,
                     offspring = 0) # give them starting attribute
  inds
}
##mutate individuals#
mutate.inds <- function(inds, mu_mut, sd_mut) {
  n <- nrow(inds)
  ##Mutate bet hedging trait [=h]
  mutate <- runif(n) < mu_mut
  inds$gene_bh[mutate] <- inds$gene_bh[mutate] + rnorm(sum(mutate), 0, sd_mut)
  inds$gene_bh <- pmin(1, pmax(0, inds$gene_bh))
  
  #mutate plasticity trait [=s]
  mutate <- runif(n) < mu_mut
  inds$gene_pls[mutate] <- inds$gene_pls[mutate] + rnorm(sum(mutate), 0, sd_mut)
  inds$gene_pls <- pmin(1, pmax(0, inds$gene_pls))
  
  #decide plasticty and alternative phenotype of offspring
  pls <<- inds$gene_pls > 0.5
  alt <<- runif(nrow(inds)) < inds$gene_bh #alternative phenotype?
  
  #Mutate intercept traits
  mutate <- runif(n) < mu_mut
  inds$gene_I0[mutate] <- inds$gene_I0[mutate] + rnorm(sum(mutate), 0, sd_mut)
  inds$gene_I0 <- pmin(1, pmax(-1, inds$gene_I0))
  
  mutate <- runif(n) < mu_mut
  inds$gene_I0_alt[mutate] <- inds$gene_I0_alt[mutate] + rnorm(sum(mutate), 0, sd_mut)
  inds$gene_I0_alt <- pmin(1, pmax(-1, inds$gene_I0_alt))
  
  #mutate developmental plasticity trait [trait = a in paper]
  mutate <- runif(n) < mu_mut
  inds$gene_devp[mutate & pls] <- inds$gene_devp[mutate & pls] + rnorm(sum(mutate & pls), 0, sd_mut)
  inds$gene_devp = pmin(1, pmax(0, inds$gene_devp))
  
  #mutate slopes#
  mutate <- runif(n) < mu_mut
  inds$gene_b[mutate & pls] <- inds$gene_b[mutate & pls] + rnorm(sum(mutate & pls), 0, sd_mut)
  inds$gene_b = pmin(2, pmax(-2, inds$gene_b))
  
  mutate <- runif(n) < mu_mut
  inds$gene_b_alt[mutate & pls] <- inds$gene_b_alt[mutate & pls] + rnorm(sum(mutate & pls), 0, sd_mut)
  inds$gene_b_alt = pmin(2, pmax(-2, inds$gene_b_alt))
  
  inds
}
#Increase age individuals#
grow.inds <- function(inds){
  inds$age <- inds$age + 1 # advance age by one
  inds
}
#Change phenotype individuals
pheno_change.inds <- function(inds, Cue) {
  if(inds$age[1] == 0) {
    ###Set plastic traits to zero for non-plastic individuals###
    inds$gene_devp[!pls] <- 0
    inds$gene_b[!pls] <- 0
    inds$gene_b_alt[!pls] <- 0
    ##Phenotype non-plastic individuals#
    inds$phenotype[!pls & !alt] = inds$gene_I0[!pls & !alt]
    inds$phenotype[!pls & alt] = inds$gene_I0_alt[!pls & alt]
    ##phenotype plastic individuals###
    inds$phenotype[pls & !alt] = inds$gene_I0[pls & !alt] + inds$gene_b[pls & !alt] * Cue
    inds$phenotype[pls & alt] = inds$gene_I0_alt[pls & alt] + inds$gene_b_alt[pls & alt] * Cue
  }
  else { ##developmental plasticity
    switch <- pls & runif(nrow(inds)) < inds$gene_devp #reversible switch
    inds$pheno_update[switch] = inds$pheno_update[switch] + 1
    inds$phenotype[switch & !alt] =  inds$gene_I0[switch & !alt] + inds$gene_b[switch & !alt] * Cue
    inds$phenotype[switch & alt] =  inds$gene_I0_alt[switch & alt] + inds$gene_b_alt[switch & alt] * Cue
  }
  inds
}
#calculate mismatch#
mismatch_update.inds <- function(inds, E) {
  mismatch = abs(inds$phenotype - E)
  inds$mismatch = inds$mismatch + mismatch #We're adding them
  inds
}
#calculate individual and pop fitness
fitness.inds <- function(inds, tau, c_plast, c_devp) {
  inds$fitness[!pls] <- exp(-tau * inds$mismatch[!pls])
  inds$fitness[pls] <-exp(-tau * inds$mismatch[pls]) - c_plast - c_devp * inds$pheno_update[pls]
  inds$fitness <- pmax(0, inds$fitness)
  inds
}
#offspring creation#
offspring.inds <- function(inds, rel_fit = T, q = 1, K) {
  if(rel_fit == T) {
    inds$offspring <- rpois(nrow(inds), inds$fitness / avg_fitness)
  } else {
    max_fitness = 1
    inds$offspring <- rpois(nrow(inds), q * inds$fitness / max_fitness)
  }
  ###add individuals if too few (not in absolute fitness scenario)#
  if(rel_fit == T & sum(inds$offspring) < K) {
    empty_spots <- K - sum(inds$offspring)
    add <- sample(1:nrow(inds), empty_spots, replace = F)
    inds$offspring[add] = inds$offspring[add] + 1
  }
  
  if(sum(inds$offspring) > 0) {
    
    #Get the new offspring 
    inds <- inds %>% slice(rep(row_number(), inds$offspring))
    
    #Remove individuals if too many (also in the absolute fitness scenario)
    if (nrow(inds) > K) {
      inds <- inds[sample(nrow(inds), K), ]
    }
    #reset some traits 
    inds$age = 0
    inds$offspring = 0
    inds$fitness = 0
    inds$mismatch = 0
    inds$pheno_update = 0
    inds$phenotype = NaN
    inds
  } else {
    return()
  }
}
#get output data#
data.output <- function(inds, E, time, dataframe) {
  mean_plast <- mean(inds$gene_pls)
  sd_plast <- sd(inds$gene_pls)
  mean_bh <- mean(inds$gene_bh)
  sd_bh <- sd(inds$gene_bh)
  mean_rev <- mean(inds$gene_devp)
  sd_rev <- sd(inds$gene_devp)
  mean_I <- mean(inds$gene_I0)
  sd_I <- sd(inds$gene_I0)
  mean_I2 <- mean(inds$gene_I0_alt)
  sd_I2 <- sd(inds$gene_I0_alt)
  mean_slope = mean(inds$gene_b)
  sd_slope = sd(inds$gene_b)
  mean_slope2 = mean(inds$gene_b_alt)
  sd_slope2 = sd(inds$gene_b_alt)
  pop_size = nrow(inds)
  mean_phenotype1 = mean(inds$phenotype[!alt])
  mean_phenotype2 = mean(inds$phenotype[alt])
  sd_phenotype = sd(inds$phenotype)
  new.row <- data.frame(Time = time, Env = Env, pop_size,
                        mean_phenotype1, 
                        mean_phenotype2,
                        sd_phenotype,
                        mean_plast, sd_plast,
                        mean_bh, sd_bh,
                        mean_rev, sd_rev,
                        mean_I, sd_I,
                        mean_I2, sd_I2,
                        mean_slope, sd_slope,
                        mean_slope2, sd_slope2)
  dataframe <- rbind(dataframe, new.row)
  dataframe
}

###Plotting functions#####
#plot environment#
plot_env <- function(R = 1, A = 1, L = 5, max_time = 100) {
  p1 <- ggplot() + xlim(0, max_time) +
    stat_function(fun = env_function,
                  args = list(R = R, A = A, L = L),
                  n = 100000,  linewidth = 2) +
    theme_bw() +   
    xlab("Time") + ylab("Environment") +
    ylim(c(-1, 1)) +
    NULL
  show(p1)
}
##plot cue#
plot_cue <- function(P) {
  p1 <- ggplot() + xlim(-1, 1) +
    stat_function(fun = env_cues,
                  args = list(P = P),
                  n = 1000,  linewidth = 2) +
    theme_bw() +   
    ylab("Cue") + xlab("Environment") +
    NULL
  show(p1)
}

plot_cue_withenv <- function(R = 1, P = 1, A = 1, L = 5, max_time = 100){
  time <- c(0:max_time)
  env <- env_function(time = time, L = L, A = A, R = R)
  cue <- env_cues(P = P, env = env)
  data <- data.frame(time = time, env = env, cue = cue)
  p1 <- ggplot(data = data, aes(x = time, y = env)) + 
    xlim(0, max_time) +
    geom_path(aes(colour = 'Environment'), linewidth = 2) +
    geom_path(aes(y = cue, colour = 'Cue'), linewidth = 2, alpha = 1) +
    theme_bw() +   
      theme(legend.title = element_blank()) +
    
    xlab("Time") + ylab("Environment/Cue") +
    ylim(c(-1, 1)) +
    NULL
  show(p1)
}



#plot data#
plot.phenotypes <- function(data = NULL, gen_length = 5) {
  inds = data$inds
  mean_plast = mean(inds$gene_pls)
  mean_I = mean(inds$gene_I0)
  mean_I2 = mean(inds$gene_I0_alt)
  mean_slope = mean(inds$gene_b)
  mean_slope2 = mean(inds$gene_b_alt)
  mean_rev = mean(inds$gene_devp)
  mean_bh = mean(inds$gene_bh)
  possibletraits <- mean_I + mean_slope * seq(-1, 1, 0.01)
  possibletraits2 <- mean_I2 + mean_slope2 * seq(-1, 1, 0.01)
  miny = min(possibletraits, possibletraits2, -1)
  maxy = max(possibletraits, possibletraits2, 1)
  
  ###plot data of last generation
  if(mean_plast > 0.5){
    color = colorRamp(c("blue", "red"))(mean_rev)
    col = rgb(color[1], color[2], color[3], maxColorValue = 255)
  } else {
    col = 'black'
  }
  
  if(mean_bh > 0.01 & mean_bh < 0.99) {
    lw1 <- 2 * (1 - mean_bh)   # Width for strategy 1
    lw2 <- 2 * mean_bh         # Width for strategy 2
  } else if (mean_bh <= 0.01) {
    lw1 <- 2 # Width for strategy 1
    lw2 <- 0 # Width for strategy 2
  } else {
    lw1 <- 0 # Width for strategy 1
    lw2 <- 2 # Width for strategy 2
  }
  
  p1 = ggplot() + 
    geom_path(aes(x = seq(-1, 1, 0.01),
                  y = mean_I2 + mean_slope2 * seq(-1, 1, 0.01)), 
              colour = col,
              linewidth = lw2) +
    geom_path(aes(x = seq(-1, 1, 0.01),
                  y = mean_I + mean_slope * seq(-1, 1, 0.01)),
              colour = col, linewidth = lw1) +
    theme_bw() + 
    ylab("Phenotype (Insulation)") + xlab('Environmental cue') +
    ggtitle("Average phenotype(s) at end of simulation") +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    ylim(c(-1, 1)) +
    NULL
  
  if(mean_plast > 0.5) {
    testdata <- data.frame(x = seq(-1, 1, 0.1), y = seq(0, 1, 0.05))
    p1 <- p1 + geom_path(data = testdata, aes(x =x, y = y, colour = y), alpha = 0) + 
      scale_color_gradient2(low = "#0000FF", high = "#FF0000", mid = "#7F007F",
                            name = "Reversibility", midpoint = 0.5,
                            breaks = c(0, 1),
                            labels = c("low", "high")) +
      theme(legend.text = element_text(size = 8), legend.key.width = unit(0.3, 'cm'),
            legend.title = element_text(size = 8)) +
      NULL
  }
  return(p1)
}

#plot individual outcome#
plot.strategy <- function(gene_pls = 1, gene_devp = 0, gene_bh = 0,
                      intercept1 = 0, slope1 = 1,
                      intercept2 = 0, slope2 = 1) {
  if(gene_bh > 1 | gene_bh < 0) {return("gene_bh should have a value between 0 and 1")}
  if(gene_pls > 1 | gene_pls < 0) {return("gene_pls should have a value between 0 and 1")}
  if(gene_devp > 1| gene_devp <0 ) {return("gene_devp should have a value between 0 and 1")}
  
  if(gene_bh > 0.01 & gene_bh < 0.99) {
    l1 <- 2 * (1 - gene_bh)   # Width for strategy 1
    l2 <- 2 * gene_bh         # Width for strategy 2
  } else if (gene_bh <= 0.01) {
    l1 <- 2 # Width for strategy 1
    l2 <- 0 # Width for strategy 2
  } else {
    l1 <- 0 # Width for strategy 1
    l2 <- 2 # Width for strategy 2
  }
  
  if(gene_pls > 0.5){
    color = colorRamp(c("blue", "red"))(gene_devp)
    col = rgb(color[1], color[2], color[3], maxColorValue = 255)
  } else {
    col = 'black'
  }
  
  if(gene_pls > 0.5) {
  p1 = ggplot(data = NULL) + 
    geom_path(aes(x = seq(-1, 1, 0.01),
                  y = intercept1 + slope1 * seq(-1, 1, 0.01)), 
              colour = col,
              linewidth = l1) +
    geom_path(aes(x = seq(-1, 1, 0.01),
                  y = intercept2 + slope2 * seq(-1, 1, 0.01)),
              colour = col, linewidth = l2) +
    theme_bw() + 
    ylab("Phenotype (Insulation)") + xlab('Environmental cue') +
    ggtitle("Phenotype(s)") +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    ylim(c(-1, 1)) +
    NULL
  
    testdata <- data.frame(x = seq(-1, 1, 0.1), y = seq(0, 1, 0.05))
    p1 <- p1 + geom_path(data = testdata, aes(x = x, y = y, colour = y), alpha = 0) + 
      scale_color_gradient2(low = "#0000FF", high = "#FF0000", mid = "#7F007F",
                            name = "Reversibility", midpoint = 0.5,
                            breaks = c(0, 1),
                            labels = c("low", "high")) +
      theme(legend.text = element_text(size = 8), legend.key.width = unit(0.3, 'cm'),
            legend.title = element_text(size = 8)) +
      NULL
  } else {
    p1 = ggplot(data = NULL) + 
      geom_path(aes(x = seq(-1, 1, 0.01),
                    y = intercept1), 
                colour = col,
                linewidth = l1) +
      geom_path(aes(x = seq(-1, 1, 0.01),
                    y = intercept2),
                colour = col, linewidth = l2) +
      theme_bw() + 
      ylab("Phenotype (Insulation)") + xlab('Environmental cue') +
      ggtitle("Phenotype(s)") +
      theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
            axis.text = element_text(size = 8)) +
      ylim(c(-1, 1)) +
      NULL
  }
  show(p1)
}


#number of generations to show -> also in the run scenario!!!
plot.genes <- function(data, maxtime = NULL, gen_length = 5) {
  plotdata <- data$output_data
  
  if(is.null(maxtime)) {
    max_plot <- min(plotdata$Time)
  } else {
    max_plot <- max(min(plotdata$Time), max(plotdata$Time) - maxtime)
  }
  plotdata <- subset(plotdata, Time >= max_plot)
  
  bh <- plotdata[nrow(plotdata), "mean_bh"] 
  l1 <- 2
  l2 <- 2
  
  p1 <- ggplot(data = plotdata, aes(x = Time / gen_length, y = mean_I)) +
    geom_path(aes(y = mean_I2, colour = "2Intercept 2"), linewidth = l2) +
    geom_path(linewidth = l1, aes(colour = "1Intercept 1")) +
    theme_bw() +
    ylim(c(-1,1)) +
    ylab("Intercept1 (black) and Intercept 2 (green)") +
    xlab("Generation") +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    theme(legend.position = "none", #legend.position.inside = c(0.3, 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.background = element_blank()) +
    scale_color_manual(values = c("black", "darkgreen"))
    NULL
    
  minslope <- min(plotdata$mean_slope, plotdata$mean_slope2) 
  maxslope <- max(plotdata$mean_slope, plotdata$mean_slope2)
  ymax <- max(maxslope, 1.5)
  ymin <- min(minslope, -1.5)
  p2 <- ggplot(data = plotdata, aes(x = Time / gen_length, y = mean_slope)) +
    geom_path(aes(y = mean_slope2, colour = "2Slope 2"), linewidth = l2) +
    geom_path(linewidth = l1, aes(colour = "1Slope 1")) +
    theme_bw() +
    ylim(c(ymin, ymax)) +
    ylab("Slope1 (black) and Slope2 (green)") +
    xlab("Generation") +
    scale_color_manual(values = c("black", "darkgreen")) +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
      theme(legend.position = "none", #legend.position.inside = c(0.3, 0.5),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.background = element_blank()) +
    NULL
  p3 <-  ggplot(data = plotdata, aes(x = Time / gen_length, y = mean_rev)) +
    geom_path(aes(y = mean_plast, colour = "2Plasticity"), linewidth = 2) +
    geom_path(linewidth = 2, aes(colour = "1Reversibility")) +
    scale_color_manual(values = c("grey", "black")) +
    geom_hline(aes(yintercept = 0.5), linetype = 'dashed') +
    theme_bw() +
    ylim(c(0, 1)) +
    ylab("Plasticity (black) and reversibility (grey)") +
    xlab("Generation") +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    theme(legend.position = "none", #legend.position.inside = c(0.3, 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.background = element_blank()) +
    NULL
  p4 <- ggplot(data = plotdata, aes(x = Time / gen_length, y = mean_bh)) +
    geom_path( linewidth = 2) +
    theme_bw() +
    ylim(c(0, 1)) +
    ylab("Probability of alternative phenotype") +
    xlab("Generation") +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    NULL
  p_all <- multiplot(p1, p2, p3, p4, cols = 2)
  show(p_all)
}


##plot population size
plot.pop <- function(data, maxtime = NULL, gen_length = 5, K = 5000) {
  plotdata <- data$output_data
  
  if(is.null(maxtime)) {
    max_plot <- min(plotdata$Time)
  } else {
    max_plot <- max(min(plotdata$Time), max(plotdata$Time) - maxtime)
  }
  plotdata <- subset(plotdata, Time >= max_plot)
  

  p1 <- ggplot(data = plotdata, aes(x = Time, y = pop_size)) +
    geom_path(linewidth = 1) +
    theme_bw() +
    ylim(c(0, K)) +
    ylab("Population size") +
    xlab("Generation") +
    theme(plot.title = element_text(size = 8, hjust = 0.5), axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    NULL
  show(p1)
}

##simulation function#####
run.model <- function(P = 1, #reliability of environmental cue
                      R = 1, #period of environmental change (per generation)
                      c_plast = 0.02, #Basic costs of plasticity
                      c_devp = 0.01, #additional costs of switching phenotype during life
                      mu_mut = 0.001, #mutational rate
                      sd_mut = 0.05, #mutational stepsize
                      q = 2.3, #max fecundity (in case of absolute fitness)
                      L = 5, #Lifespan
                      K = 5000, #carrying capacity
                      tau = 0.25, #scaling of fitness
                      seed = 1, #random seed
                      Max_gen = 20000, #generations per simulations
                      plot_gen = NULL, #how many generations to keep track of in output (default is max(30,R))
                      rel_fitness = T, #relative vs absolute fitness
                      h_start = 0, #initial value for bet hedging gene
                      inds = NULL,
                      test = F,
                      progressbar = T) {
  A = 1 #Predictability of enviroment
  name <- paste("Simulation over ", Max_gen, " generations. R = ", R, " P = ", P, " seed = ", seed, "\n" , sep ="")
  n.initial <- K #Initial number of individuals
  perc <- (Max_gen * L) / 20 #5% steps
  max_time <- L * Max_gen #maximum time of simulation
  if(is.null(plot_gen)) {
    plot_gen <- max(30, R * 1)
  }
  #if(is.null(start_output)) {
  start_output <- max(0, max_time - L * plot_gen) #when to generate timeseries data
  #}
  if(R < 5) { #how often to generate timeserie data
    period_output <- 1
  } else {
    period_output <- 5
  }
  set.seed(seed) ##make sure repeatable#
  output_data <- data.frame() #create empty data frame
  if(is.null(inds)) { #create initial individuals
    inds <- make.inds(id=1:n.initial, n = n.initial, gene_bh = h_start) 
  }
  pls <<- inds$gene_pls > 0.5 #determine who is plastic
  alt <<- runif(nrow(inds)) < inds$gene_bh #determine who gets alternative phenotype
  param <- c(R = R, P = P, seed = seed,
             c_plast = c_plast, 
             c_devp = c_devp, 
             mu_mut = mu_mut, 
             sd_mut = sd_mut, 
             q = q,
             L = L,
             K = K, 
             tau = tau) #Save parameters
  time <- 1 #start time
  start.time <- Sys.time() #to calculate runtime
  for(time in 1:max_time) {
    
    #1) update environment + cues
    Env <<- env_function(time = time, R = R, L = L, A = A)
    Cue <- env_cues(P = P, env = Env)
    #2) phenotypic development/update
    inds <- pheno_change.inds(inds, Cue = Cue)
    #3) calculate mismatch and add it
    inds <- mismatch_update.inds(inds, E = Env)
    #update time and age#
    inds <- grow.inds(inds)
    
    #Get population state at certain time points for output#
    if(time %% period_output == 0 & time > start_output) {
      output_data <- data.output(inds, E = Env, time = time, dataframe = output_data)
    }
    
    #at end of each generation:
    if(inds$age[1] == L) {
      #A) Calculate individual and average fitness
      inds <- fitness.inds(inds, tau = tau, c_plast = c_plast, c_devp = c_devp)
      avg_fitness <<- mean(inds$fitness)
      #B) Create offspring
      inds <- offspring.inds(inds, rel_fit = rel_fitness, q = q, K = K)
      
      if(is.null(inds)) {
        cat(paste('\nAt time ', time, " the population is extinct", sep = ""))
        break
      }
      #C) Mutate offspring
      inds <- mutate.inds(inds, mu_mut = mu_mut, sd_mut = sd_mut)
    }
    
    
    #Progress bar#
    if(progressbar) {
    if(!test) {
    if(time %% perc == 0 | time == 1) {
      cat('\014')
      cat(paste0(name, round(time / (Max_gen * L) * 100), '% completed'))
      if (time == Max_gen * L) cat(': Done\n')
    }}
    }
    
  }
  end.time <- Sys.time() #to calculate runtime
  time.taken <<- end.time - start.time #runtime
  data_save <- list(inds = inds, output_data = output_data, param = param)
  
  cat("Duration of run:", time.taken, units(time.taken))
  if(test) {
    return(cat('\014 Model runs succesfully\n'))
    } else {
    return(data_save)
  }
}


###Test if everything works and let the user know####
run.model(P = 1, R = 1, seed = 1, Max_gen = 1, test = T)
cat("Script_model succesfully loaded")

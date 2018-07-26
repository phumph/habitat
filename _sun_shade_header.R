### header file for Sun/Shade manuscript
### Last updated: 08 Feb 2018 by PTH

# load globals:
suppressPackageStartupMessages(require(ggplot2, warn.conflicts = FALSE))
suppressPackageStartupMessages(require(gridExtra, warn.conflicts = FALSE))
suppressPackageStartupMessages(require(ggpubr, warn.conflicts = FALSE))
suppressPackageStartupMessages(require(glmmTMB, warn.conflicts = FALSE))
suppressPackageStartupMessages(require(bbmle, warn.conflicts = FALSE))
suppressPackageStartupMessages(require(dplyr, warn.conflicts = FALSE))

# load plotting theme:
theme_pth1 <- function(base_size = 9, base_family = "sans") {
  theme_bw(
    base_family = base_family,
    base_size = base_size
  ) +
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 0.4),
      #axis.ticks = element_line(size = 0.3),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = rel(0.9)),
      strip.placement = "outside",
      # strip.background = element_rect(fill = "gray95", color = NA),
      panel.spacing = unit(1.5, "lines"),
      legend.position = "right",
      legend.background = element_blank(),
      legend.text = element_text(size = 9),
      legend.text.align = 0,
      legend.key = element_blank()
    )
}

#### define test statistic functions ####
# proportion of zeros (p(0)) returned by simulation results:
prop.zero <- function(sim1, dat1, fact1 = 'site_type', response1 = 'stipples'){
  # define zero prop calculation function
  calc.prop.zero <- function(x){ length(which(x==0))/length(x) }
  
  # join data and sim data:
  dat2 <- cbind(dat1, sim1)
  
  # find observed number of zeros by grouping factor:
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]
  
  # calculate distribution of proportion zeros by grouping factor for all simulations:
  zeros.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, calc.prop.zero)
  zeros.sim.2 <- apply(lvl2[,names(lvl1) %in% names(sim1)], 2, calc.prop.zero)
  
  # calculate proportion of observed zeros by grouping factor:
  lvl1.0 <- calc.prop.zero(lvl1[,paste0(response1)])
  lvl2.0 <- calc.prop.zero(lvl2[,paste0(response1)])
  
  # construct data.frame
  res1 <- data.frame(prop.zero = c(as.vector(zeros.sim.1),as.vector(zeros.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1],length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(prop.zero = c(lvl1.0,lvl2.0),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))
  
  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.0 < zeros.sim.1))/length(zeros.sim.1),3)
  q.2 <- round(length(which(lvl2.0 < zeros.sim.2))/length(zeros.sim.2),3)
  
  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)
  
  # plot:
  plot1 <- ggplot() + 
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = prop.zero), bins = 25, col = "gray40", fill = "gray40") + 
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = prop.zero), col = "red") +
    facet_wrap(~ fact1, scales = "free") + 
    xlab(expression(p(Y=0))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_pth1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list('plot' = plot1, 'data' = res2))  
}


#### prop.zero.plant-level
prop.zero2 <- function(sim1, dat1, fact1 = 'site_type', response1 = 'stipples', plant.level = TRUE, pfact = 'stem_id'){
  # define zero prop calculation function
  calc.prop.zero <- function(x){ length(which(x==0))/length(x) }
  
  # join data and sim data:
  dat2 <- cbind(dat1, sim1)
  
  # if plant.level == TRUE, aggregate by plant id within site to sum up counts:
  if (plant.level == TRUE){
    dat2 <- aggregate(dat2[,names(dat2) %in% c(paste0(response1),grep('sim_', names(sim1), value = T))], by = list(dat2[,paste0(fact1)],dat2[,paste0(pfact)]), FUN="sum")
    names(dat2)[1:2] <- c(paste0(fact1),paste0(pfact))
  } 
  
  # find observed number of zeros by grouping factor:
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]
  
  # calculate distribution of proportion zeros by grouping factor for all simulations:
  zeros.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, calc.prop.zero)
  zeros.sim.2 <- apply(lvl2[,names(lvl1) %in% names(sim1)], 2, calc.prop.zero)
  
  # calculate proportion of observed zeros by grouping factor:
  lvl1.0 <- calc.prop.zero(lvl1[,paste0(response1)])
  lvl2.0 <- calc.prop.zero(lvl2[,paste0(response1)])
  
  # construct data.frame
  res1 <- data.frame(prop.zero = c(as.vector(zeros.sim.1),as.vector(zeros.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1],length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(prop.zero = c(lvl1.0,lvl2.0),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))
  
  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.0 < zeros.sim.1))/length(zeros.sim.1),3)
  q.2 <- round(length(which(lvl2.0 < zeros.sim.2))/length(zeros.sim.2),3)
  
  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)
  
  # plot:
  plot1 <- ggplot() + 
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = prop.zero), bins = 25, col = "gray40", fill = "gray40") + 
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = prop.zero), col = "red") +
    facet_wrap(~ fact1, scales = "free") + 
    xlab(expression(p(Y=0))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_pth1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list('plot' = plot1, 'data' = res2))  
}


# maximum value observed among simulation results:
pp.max <- function(sim1, dat1, fact1 = 'site_type', response1 = 'stipples'){
  
  # join data and sim data:
  dat2 <- cbind(dat1, sim1)
  
  # split data by grouping factor
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]
  
  # calculate maximum value of all simulations:
  max.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, max)
  max.sim.2 <- apply(lvl2[,names(lvl2) %in% names(sim1)], 2, max)
  
  # calculate proportion of observed zeros by grouping factor:
  lvl1.max <- max(lvl1[,paste0(response1)])
  lvl2.max <- max(lvl2[,paste0(response1)])
  
  # construct data.frame
  res1 <- data.frame(max.val = c(as.vector(max.sim.1),as.vector(max.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1], length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(max.val = c(lvl1.max,lvl2.max),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))
  
  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.max < max.sim.1))/length(max.sim.1),3)
  q.2 <- round(length(which(lvl2.max < max.sim.2))/length(max.sim.2),3)
  
  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)
  
  # plot:
  plot1 <- ggplot() + 
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = max.val), bins = 25, col = "gray40", fill = "gray40") + 
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = max.val), col = "red") +
    facet_wrap(~ fact1, scales = "free") + 
    xlab(expression(max(Y))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_pth1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list('plot' = plot1, 'data' = res2))  
}

### calculate means
pp.mean <- function(sim1, dat1, fact1 = 'site_type', response1 = 'stipples'){
  
  # join data and sim data:
  dat2 <- cbind(dat1, sim1)
  
  # split data by grouping factor
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]
  
  # calculate maximum value of all simulations:
  max.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, mean)
  max.sim.2 <- apply(lvl2[,names(lvl2) %in% names(sim1)], 2, mean)
  
  # calculate proportion of observed zeros by grouping factor:
  lvl1.max <- mean(lvl1[,paste0(response1)])
  lvl2.max <- mean(lvl2[,paste0(response1)])
  
  # construct data.frame
  res1 <- data.frame(mean.val = c(as.vector(max.sim.1),as.vector(max.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1], length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(mean.val = c(lvl1.max,lvl2.max),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))
  
  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.max < max.sim.1))/length(max.sim.1),3)
  q.2 <- round(length(which(lvl2.max < max.sim.2))/length(max.sim.2),3)
  
  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)
  
  # plot:
  plot1 <- ggplot() + 
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = mean.val), bins = 25, col = "gray40", fill = "gray40") + 
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = mean.val), col = "red") +
    facet_wrap(~ fact1, scales = "free") + 
    xlab(expression(mean(Y))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_pth1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list('plot' = plot1, 'data' = res2))  
}


pp.mean2 <- function(sim1, dat1, fact1 = 'site_type', response1 = 'stipples', plant.level = FALSE, pfact){
  
  # join data and sim data:
  dat2 <- cbind(dat1, sim1)
  
  if (plant.level == TRUE){
    dat2 <- aggregate(dat2[,names(dat2) %in% c(paste0(response1),grep('sim_', names(sim1), value = T))], by = list(dat2[,paste0(fact1)],dat2[,paste0(pfact)]), FUN="sum")
    names(dat2)[1:2] <- c(paste0(fact1),paste0(pfact))
  } 
  
  # split data by grouping factor
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]
  
  # calculate maximum value of all simulations:
  max.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, mean)
  max.sim.2 <- apply(lvl2[,names(lvl2) %in% names(sim1)], 2, mean)
  
  # calculate proportion of observed zeros by grouping factor:
  lvl1.max <- mean(lvl1[,paste0(response1)])
  lvl2.max <- mean(lvl2[,paste0(response1)])
  
  # construct data.frame
  res1 <- data.frame(mean.val = c(as.vector(max.sim.1),as.vector(max.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1], length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(mean.val = c(lvl1.max,lvl2.max),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))
  
  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.max < max.sim.1))/length(max.sim.1),3)
  q.2 <- round(length(which(lvl2.max < max.sim.2))/length(max.sim.2),3)
  
  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)
  
  # plot:
  plot1 <- ggplot() + 
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = mean.val), bins = 25, col = "gray40", fill = "gray40") + 
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = mean.val), col = "red") +
    facet_wrap(~ fact1, scales = "free") + 
    xlab(expression(mean(Y))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_pth1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list('plot' = plot1, 'data' = res2))  
}


# produce better test statistic: inter-quartile range using R base function IQR()
pp.iqr <- function(sim1, dat1, fact1 = 'site_type', response1 = 'stipples'){
  
  # join data and sim data:
  dat2 <- cbind(dat1, sim1)
  
  # split data by grouping factor
  lvl1 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[1],]
  lvl2 <- dat2[dat2[,paste0(fact1)] == levels(dat2[,paste0(fact1)])[2],]
  
  # calculate maximum value of all simulations:
  iqr.sim.1 <- apply(lvl1[,names(lvl1) %in% names(sim1)], 2, IQR)
  iqr.sim.2 <- apply(lvl2[,names(lvl2) %in% names(sim1)], 2, IQR)
  
  # calculate proportion of observed zeros by grouping factor:
  lvl1.iqr <- IQR(lvl1[,paste0(response1)])
  lvl2.iqr <- IQR(lvl2[,paste0(response1)])
  
  # construct data.frame
  res1 <- data.frame(iqr.val = c(as.vector(iqr.sim.1),as.vector(iqr.sim.2)),
                     fact1 = c(rep(levels(dat2[,paste0(fact1)])[1], length(sim1[1,])),
                               rep(levels(dat2[,paste0(fact1)])[2], length(sim1[1,]))),
                     type = 'sim'
  )
  res2 <- rbind(res1, data.frame(iqr.val = c(lvl1.iqr,lvl2.iqr),
                                 fact1 = c(levels(dat2[,paste0(fact1)])[1],levels(dat2[,paste0(fact1)])[2]),
                                 type = c('obs','obs')
  ))
  
  # calculate p-values for placing them on the plot:
  q.1 <- round(length(which(lvl1.iqr < iqr.sim.1))/length(iqr.sim.1),3)
  q.2 <- round(length(which(lvl2.iqr < iqr.sim.2))/length(iqr.sim.2),3)
  
  # create labels for factor levels passed to ggplot which include p-value:
  thelabs <- c(paste0(levels(dat2[,paste0(fact1)])[1],'\np = ',round(q.1,3)),
               paste0(levels(dat2[,paste0(fact1)])[2],'\np = ',round(q.2,3))
  )
  # add labels to factor level
  res2[,'fact1'] <- factor(res2[,'fact1'], levels = levels(res2[,'fact1']), labels = thelabs)
  
  # plot:
  plot1 <- ggplot() + 
    geom_histogram(data = dplyr::filter(res2, type == 'sim'), aes(x = iqr.val), bins = 25, col = "gray40", fill = "gray40") + 
    geom_vline(data = dplyr::filter(res2, type == 'obs'), aes(xintercept = iqr.val), col = "red") +
    facet_wrap(~ fact1, scales = "free") + 
    xlab(expression(IQR(Y))) +
    #geom_text(data = qv, aes(label = qv), x = xval, y = yval) +
    #ggtitle(data = qv, aes(label = qv)) +
    theme_pth1() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  
  return(list('plot' = plot1, 'data' = res2))  
}


# generate residuals plot:
resid.plot <- function(m1, d1, resids = TRUE, response1 = 'stipples', pearson=TRUE){
  if (pearson==FALSE){
    d2 <- data.frame(d1, resids = resid(m1), fits = predict(m1, type = "response"))
  } else{
    d2 <- data.frame(d1, resids = resid(m1, type = 'pearson'), fits = predict(m1, type = "response"))
  }
  #if (resids = TRUE) {
  # plot residuals vs. fitted
  p1 <- ggplot(d2, aes(x = fits, y = resids, alpha = 0.5)) + 
    geom_point(col = "black") + 
    facet_wrap(~ site_type, scales = "free") +
    xlab("fitted") + ylab("residuals (Pearson)") +
    theme_pth1() + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  
    geom_hline(yintercept = 0, col = "red")
  
  #} else{
  p2 <- ggplot(d2, aes_string(x = 'fits', y = paste0(response1), alpha = 0.5)) + 
    geom_point(col = "black") + 
    facet_wrap(~ site_type, scales = "free") +
    xlab("fitted") + ylab("observed") +
    theme_pth1() + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +  
    geom_abline(intercept = 0, slope = 1, col = "red")
  #  }
  return(list('p1' = p1, 'p2' = p2))
}



#### COMPOSITE PREDICTIVE CHECK FUNCTIONS
pp.plots <- function(thesim, thedat, thefact, theresponse, plant.level = FALSE, plantfact, fourcats = FALSE){
  
  # combine data from sims to data.frame:
  thedat2 <- cbind(thedat,thesim)
  
  # if plant_level == TRUE, aggregate by plantfact:
  if (plant.level == TRUE){
    # modify
    p1 <- prop.zero2(thesim, thedat, fact1 = thefact, response1 = theresponse, plant.level = TRUE, pfact = plantfact)
    p2 <- pp.mean2(thesim, thedat, fact1 = thefact, response1 = theresponse, plant.level = TRUE, pfact = plantfact)
    
    # now aggregate the data for calculations below:
    thedat2b <- aggregate(thedat2[,names(thedat2) %in% c(paste0(theresponse),grep('sim_', names(thesim), value = T))], by = list(thedat2[,paste0(thefact)],thedat2[,paste0(plantfact)]), FUN="sum")
    names(thedat2b)[1:2] <- c(fact1,plantfact)
    thedat3 <- reshape::melt.data.frame(thedat2b, id.vars = c(paste0(thefact),paste0(plantfact)), value.name = paste0(theresponse))
    names(thedat3)[1] <- 'fact'

  } else{
    if (fourcats == FALSE){
    p1 <- prop.zero2(thesim, thedat, fact1 = thefact, response1 = theresponse, plant.level = FALSE)
    p2 <- pp.mean2(thesim, thedat, fact1 = thefact, response1 = theresponse, plant.level = FALSE)
    thedat2b <- thedat2[,names(thedat2) %in% c(paste0(theresponse),paste0(thefact),grep('sim_', names(thesim), value = T))]
    thedat3 <- reshape::melt.data.frame(thedat2b, id.vars = paste0(thefact), value.name = paste0(theresponse))
    names(thedat3)[1] <- 'fact'
    } else if (fourcats == TRUE){
      thedat2b <- thedat2[,names(thedat2) %in% c(paste0(theresponse),paste0(thefact),grep('sim_', names(thesim), value = T))]
      thedat3 <- reshape::melt.data.frame(thedat2b, id.vars = paste0(thefact), value.name = paste0(theresponse))
      names(thedat3)[1] <- 'fact'
    } 
  }
  # generate summary statistics for means and 95% CI on the mean estimate from simulations:
  thedat4 <- dplyr::filter(thedat3, !(variable %in% paste0(theresponse))) %>% group_by(fact, variable) %>% summarise(mu = mean(value))
  thedat5 <- dplyr::group_by(thedat4, fact) %>% summarise(mu2 = mean(mu),
                                                          lower = quantile(mu, probs = 0.025),
                                                          upper = quantile(mu, probs = 0.975))
  
  if (fourcats == FALSE){
  # calculate summary of prop.zero data:
  zerodat <- dplyr::filter(p1$data, type == 'sim') %>% group_by(fact1) %>% summarise(p0 = mean(prop.zero),
                                                                                     lower = quantile(prop.zero, probs = 0.025),
                                                                                     upper = quantile(prop.zero, probs = 0.975))
  } 
  # calculate infection intensity
  # capture only rows that are above zero:
  thedat4b <- dplyr::filter(thedat3, !(variable %in% paste0(theresponse)), value > 0) %>% group_by(fact, variable) %>% summarise(mu = mean(value))
  thedat5b <- dplyr::group_by(thedat4b, fact) %>% summarise(mu2 = mean(mu),
                                                          lower = quantile(mu, probs = 0.025),
                                                          upper = quantile(mu, probs = 0.975))
  
  if (fourcats == FALSE){
  return(list(plots = ggarrange(plotlist = list(p1[[1]],p2[[1]]), ncol = 2, align = 'hv', labels = c('A','B','C','D')),
              means = thedat5,
              p0 = zerodat,
              intensity = thedat5b))
  } else if (fourcats == TRUE){
    return(list(means = thedat5,
                intensity = thedat5b))
  }
}

pred.mean.int <- function(themod, thedat, NSIM = 1000, IDs = 'condition', VAL = 'eggs'){
  require(reshape2)
  require(dplyr)
  sim <- simulate(themod, nsim = NSIM)
  thedat2 <- cbind(thedat[,names(thedat) %in% c(IDs,VAL)],sim)
  thedat3 <- reshape::melt.data.frame(thedat2, id.vars = paste0(IDs), value.name = paste0(VAL))
  thedat4 <- dplyr::group_by(thedat3, condition, variable) %>% summarise(mu = mean(value))
  thedat5 <- dplyr::group_by(thedat4, condition) %>% summarise(mu2 = mean(mu),
                                                               lower = quantile(mu, probs = 0.025),
                                                               upper = quantile(mu, probs = 0.975))
  return(list(each = thedat4, means = thedat5))
}
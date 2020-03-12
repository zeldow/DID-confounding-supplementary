## ----setup, include=FALSE----------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, dpi = 300)
library(tibble)
library(knitr)
library(kableExtra)
library(tidyverse)
library(lme4)
library(multiwayvcov)
library(lmtest)
library(broom)
library(MatchIt)
library(gridExtra)
library(RColorBrewer)

percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "\\%")
}

#nsims <- 5
#n <- 500
nsims <- 400
n <- 800
max.time <- 10 ## do not change (may break other parts)
trt.time <- 6 ## do not change (may break other parts)


## ----functions, include=FALSE------------------------------------------------------------------------------------------
match.on.cov <- function(dat, method = "nearest", m.order = "random", caliper = 0.2, replace=TRUE) {
  tmp <- filter(dat,tp < trt.time) %>% select(id,tp,trt,x)
  tmp2 <- matchit(trt~x,data=tmp,
                  discard="both",method=method,m.order=m.order,caliper=caliper, replace=replace)
  tmp3 <- left_join(match.data(tmp2) %>% select(id), dat, by = "id")
  return(tmp3)
}

match.on.cov2 <- function(dat, method = "nearest", m.order = "random", caliper = 0.2, replace=TRUE) {
  tmp <- filter(dat,tp < trt.time) %>% select(id,tp,trt,x) %>% 
    spread(key=tp,value=x) %>% rename(x1=`1`,x2=`2`,x3=`3`,x4=`4`,x5=`5`)
  tmp2 <- matchit(trt~x1+x2+x3+x4+x5,data=tmp,
                  discard="both",method=method,m.order=m.order,caliper=caliper, replace=replace)
  tmp3 <- left_join(match.data(tmp2) %>% select(id), dat, by = "id")
  return(tmp3)
}

match.on.lev <- function(dat, method = "nearest", m.order = "random", caliper = 0.2, replace=TRUE) {
  tmp <- filter(dat,tp < trt.time) %>% select(id,tp,trt,y)  %>%
    spread(key=tp,value=y) %>% rename(x1=`1`,x2=`2`,x3=`3`,x4=`4`,x5=`5`)
  tmp2 <- matchit(trt ~ x1 + x2 + x3 + x4 + x5, data=tmp,
                  discard="both",method=method,m.order=m.order,caliper=caliper,replace=replace)
  tmp3 <- left_join(match.data(tmp2) %>% select(id), dat, by = "id")
  return(tmp3)
}

match.on.lev2 <- function(dat, method = "nearest", m.order = "random", caliper = 0.2, replace=TRUE) {
  tmp <- filter(dat,tp < trt.time) %>% select(id,tp,trt,y.t)  %>%
    spread(key=tp,value=y.t) %>% rename(x1=`1`,x2=`2`,x3=`3`,x4=`4`,x5=`5`)
  tmp2 <- matchit(trt ~ x1 + x2 + x3 + x4 + x5, data=tmp,
                  discard="both",method=method,m.order=m.order,caliper=caliper,replace=replace)
  tmp3 <- left_join(match.data(tmp2) %>% select(id), dat, by = "id")
  return(tmp3)
}

match.on.tre <- function(dat, method = "nearest", m.order = "random", caliper = 0.2,replace=TRUE) {
  tmp <- filter(dat,tp < trt.time) %>% select(id,tp,trt,y.diff)  %>%
    spread(key=tp,value=y.diff) %>% rename(x1=`1`,x2=`2`,x3=`3`,x4=`4`,x5=`5`) %>% select(-x1)
  tmp2 <- matchit(trt ~ x2 + x3 + x4 + x5, data=tmp,
                  discard="both",method=method,m.order=m.order,caliper=caliper,replace=replace)
  tmp3 <- left_join(match.data(tmp2) %>% select(id), dat, by = "id")
  return(tmp3)
}

match.on.tre2 <- function(dat, method = "nearest", m.order = "random", caliper = 0.2,replace=TRUE) {
  tmp <- filter(dat,tp < trt.time) %>% select(id,tp,trt,y.diff2)  %>%
    spread(key=tp,value=y.diff2) %>% rename(x1=`1`,x2=`2`,x3=`3`,x4=`4`,x5=`5`) %>% select(-x1)
  tmp2 <- matchit(trt ~ x2 + x3 + x4 + x5, data=tmp,
                  discard="both",method=method,m.order=m.order,caliper=caliper,replace=replace)
  tmp3 <- left_join(match.data(tmp2) %>% select(id), dat, by = "id")
  return(tmp3)
}



fit.mods <- function() {
  
  res1 <- lm(y ~ trt*post, dat)
  res1.adj <- coeftest(res1, cluster.vcov(res1, cluster = dat[,"id"]))
  res2 <- lm(y ~ trt*post + x, dat)
  res2.adj <- coeftest(res2, cluster.vcov(res2, cluster = dat[,"id"]))
  res3 <- lm(y ~ trt*post + factor(tp), dat)
  res3.adj <- coeftest(res3, cluster.vcov(res3, cluster = dat[,"id"]))
  res4 <- lm(y ~ trt*post + x + factor(tp), dat)
  res4.adj <- coeftest(res4, cluster.vcov(res4, cluster = dat[,"id"]))
  res5 <- lm(y ~ trt*post + x*factor(tp), dat)
  res5.adj <- coeftest(res5, cluster.vcov(res5, cluster = dat[,"id"]))
  res6 <- lm(y ~ trt*post + factor(tp), matched.lev.dat)
  res6.adj <- coeftest(res6, cluster.vcov(res6, cluster = matched.lev.dat[,"id"]))
  res7 <- lm(y ~ trt*post + factor(tp), matched.tre.dat)
  res7.adj <- coeftest(res7, cluster.vcov(res7, cluster = matched.tre.dat[,"id"]))
  res8 <- lm(y ~ trt*post + factor(tp), matched.cov.dat)
  res8.adj <- coeftest(res8, cluster.vcov(res8, cluster = matched.cov.dat[,"id"]))
  
  
  mean.res <- c(coef(res1)["trt:postTRUE"], coef(res2)["trt:postTRUE"],
                   coef(res3)["trt:postTRUE"], coef(res4)["trt:postTRUE"],
                   coef(res5)["trt:postTRUE"], coef(res6)["trt:postTRUE"],
                   coef(res7)["trt:postTRUE"], coef(res8)["trt:postTRUE"])
  se.res <- c(res1.adj["trt:postTRUE", "Std. Error"],
                   res2.adj["trt:postTRUE", "Std. Error"],
                   res3.adj["trt:postTRUE", "Std. Error"],
                   res4.adj["trt:postTRUE", "Std. Error"],
                   res5.adj["trt:postTRUE", "Std. Error"],
                   res6.adj["trt:postTRUE", "Std. Error"],
                   res7.adj["trt:postTRUE", "Std. Error"],
                   res8.adj["trt:postTRUE", "Std. Error"])
  
  return(list(mean.res, se.res))
}



fit.mods2 <- function() {
  
  res1 <- lm(y.t ~ trt*post, dat)
  res1.adj <- coeftest(res1, cluster.vcov(res1, cluster = dat[,"id"]))
  res2 <- lm(y.t ~ trt*post + x, dat)
  res2.adj <- coeftest(res2, cluster.vcov(res2, cluster = dat[,"id"]))
  res3 <- lm(y.t ~ trt*post + factor(tp), dat)
  res3.adj <- coeftest(res3, cluster.vcov(res3, cluster = dat[,"id"]))
  res4 <- lm(y.t ~ trt*post + x + factor(tp), dat)
  res4.adj <- coeftest(res4, cluster.vcov(res4, cluster = dat[,"id"]))
  res5 <- lm(y.t ~ trt*post + x*factor(tp), dat)
  res5.adj <- coeftest(res5, cluster.vcov(res5, cluster = dat[,"id"]))
  res6 <- lm(y.t ~ trt*post + factor(tp), matched.lev.dat)
  res6.adj <- coeftest(res6, cluster.vcov(res6, cluster = matched.lev.dat[,"id"]))
  res7 <- lm(y.t ~ trt*post + factor(tp), matched.tre.dat)
  res7.adj <- coeftest(res7, cluster.vcov(res7, cluster = matched.tre.dat[,"id"]))
  res8 <- lm(y.t ~ trt*post + factor(tp), matched.cov.dat)
  res8.adj <- coeftest(res8, cluster.vcov(res8, cluster = matched.cov.dat[,"id"]))
  
  
  mean.res <- c(coef(res1)["trt:postTRUE"], coef(res2)["trt:postTRUE"],
                   coef(res3)["trt:postTRUE"], coef(res4)["trt:postTRUE"],
                   coef(res5)["trt:postTRUE"], coef(res6)["trt:postTRUE"],
                   coef(res7)["trt:postTRUE"], coef(res8)["trt:postTRUE"])
  se.res <- c(res1.adj["trt:postTRUE", "Std. Error"],
                   res2.adj["trt:postTRUE", "Std. Error"],
                   res3.adj["trt:postTRUE", "Std. Error"],
                   res4.adj["trt:postTRUE", "Std. Error"],
                   res5.adj["trt:postTRUE", "Std. Error"],
                   res6.adj["trt:postTRUE", "Std. Error"],
                   res7.adj["trt:postTRUE", "Std. Error"],
                   res8.adj["trt:postTRUE", "Std. Error"])
  
  return(list(mean.res, se.res))
}


### plotting functions

makeplot <- function(outcome, dat, ylab, xlab, y.max, y.min, x.axis, title, subtitle, left.title = TRUE, bottom.title = TRUE) {
  
    ggplot(dat, aes(x=x, y=y)) +
      geom_segment( aes(x=x, xend=x, y=0, yend=y), size = 1.2) +
      geom_point( size=2, color="black", fill="black", alpha=0.7, shape=21, stroke=2) +
      theme_minimal() + ylab(ylab) + xlab(xlab) + ggtitle(title, subtitle = subtitle) +
      scale_x_discrete(labels=(c("c" = "Simple", "d" = "CA", "e" = "TVA", 
                                 "f" = "Match (level)", "g" = "Match (trend)", "h" = "Match (cov)" ))) +
      theme(axis.text.x = x.axis,
            axis.text.y = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 10, hjust = 0.5),
            plot.title = element_text(hjust = 0.5),
            text = element_text(face = "bold", size = 14)) + coord_cartesian(ylim = c(y.min, y.max))
}

makegrid <- function(outcomes, scenario = 1:3, 
                     sub = c("Time-Invariant\nCovariate Effect", 
                             "Time-Varying\nCovariate Effect", 
                             "Treatment-Independent\nCovariate"),
                     truth = NULL) {
  len <- length(outcomes)
  if(len%%2 != 0) stop("outcomes variable needs to be even")
  if(is.null(truth)) { truth <- rep(1, times = len)}
  for(i in 1:len) {
    
    if(i %in% c(1, len / 2 + 1)) {
      left.title = TRUE
    } else left.title = FALSE
    if(i == ceiling(len*0.75)) {
      bottom.title = TRUE
    } else bottom.title = FALSE
    
    if(str_sub(names(outcomes)[i], start = 1, end = 2) == "mv") {
      xlab = " "
      if(left.title == TRUE) ylab = "Mean % Bias" else ylab = " "
      dat <- tibble(x = letters[3:8], y = (colMeans(outcomes[[i]][ , 3:8]) - truth[i]) * 100)
      y.max <- 100
      y.min <- -50
      x.axis <- element_blank()
      title <- paste("Scenario", scenario[i])
      subtitle <- sub[i]
    } else {
      if(bottom.title == TRUE) xlab = "Model" else xlab = " "
      if(left.title == TRUE) ylab = "Mean SE" else ylab = " " 
      dat <- tibble(x = letters[3:8], y = colMeans(outcomes[[i]][ , 3:8]))
      y.max <- 0.2
      y.min <- 0
      x.axis <- element_text(angle = 45, hjust = 1, face = "bold", size = 12)
      title <- ""
      subtitle <- ""
    }
    
    assign(paste0("plot", i), makeplot(outcomes[[i]], dat, ylab, xlab, y.max, y.min, x.axis, title, subtitle, left.title, bottom.title))
  }
  
  grid.arrange(ncol = len / 2, 
              plot1, plot2, plot3,
              plot4, plot5, plot6)
}


## ----sim-setup, echo=FALSE---------------------------------------------------------------------------------------------
set.seed(55)


## ----sim-time-invariant-no-confounder, echo=FALSE----------------------------------------------------------------------
mv.no <- matrix(NA, ncol = 8, nrow = nsims)
se.no <- matrix(NA, ncol = 8, nrow = nsims)

for(i in 1:nsims) {

  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=0.5, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           x=rnorm(1, mean = 1.5 - 0.5*trt, sd = 1.5 - 0.5*trt),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(err=rnorm(n*max.time), 
                        y = 1 + x + trt + int + err + treated + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  ## match processing
  matched.cov.dat <- match.on.cov(dat)
  matched.lev.dat <- match.on.lev(dat)
  matched.tre.dat <- match.on.tre(dat)
  
  ## fit regression models
  my.res <- fit.mods()
  
  ## save output
  mv.no[i, ] <- my.res[[1]]
  se.no[i, ] <- my.res[[2]]
}


## ----sim-time-invariant-confounder, echo=FALSE-------------------------------------------------------------------------
mv.yes <- matrix(NA, ncol = 8, nrow = nsims)
se.yes <- matrix(NA, ncol = 8, nrow = nsims)

for(i in 1:nsims) {
  
  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=0.5, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           x=rnorm(1, mean = 1.5 - 0.5*trt, sd = 1.5 - 0.5*trt),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(err=rnorm(n*max.time), 
                        y = 1 + x*tp/10 + trt + int + err + treated + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  
  ## match processing
  matched.cov.dat <- match.on.cov(dat)
  matched.lev.dat <- match.on.lev(dat)
  matched.tre.dat <- match.on.tre(dat)
  
  ## fit regression models
  my.res <- fit.mods()
  
  ## save output
  mv.yes[i, ] <- my.res[[1]]
  se.yes[i, ] <- my.res[[2]]
}


## ----sim-time-invariant-efficiency, echo=FALSE-------------------------------------------------------------------------
mv.eff <- matrix(NA, ncol = 8, nrow = nsims)
se.eff <- matrix(NA, ncol = 8, nrow = nsims)

for(i in 1:nsims) {

  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=0.5, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           x=rnorm(1),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1) # time-varying indicator if treated or not
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(err=rnorm(n*max.time), 
                        y = 1 + x*tp/10 + trt + int + err + treated + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y)) %>% ungroup()
  
  
  ## match processing
  matched.cov.dat <- match.on.cov(dat)
  matched.lev.dat <- match.on.lev(dat)
  matched.tre.dat <- match.on.tre(dat)
  
  ## fit regression models
  my.res <- fit.mods()
  
  ## save output
  mv.eff[i, ] <- my.res[[1]]
  se.eff[i, ] <- my.res[[2]]
}


## ----plot-first, echo=FALSE, fig.height=5, fig.width=7, fig.cap="\\label{fig:time-invariant-results}Simulation Results for Time-Invariant Covariate"----
makegrid(list(mv.no = mv.no, mv.yes = mv.yes, mv.eff = mv.eff, se.no = se.no, se.yes = se.yes, se.eff = se.eff))


## ----sim-time-varying-same, echo=FALSE---------------------------------------------------------------------------------
mv.same <- matrix(NA, ncol = 8, nrow = nsims)
se.same <- matrix(NA, ncol = 8, nrow = nsims)

mv.same2 <- matrix(NA, ncol = 8, nrow = nsims)
se.same2 <- matrix(NA, ncol = 8, nrow = nsims)

for(i in 1:nsims) {

  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=0.5, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           x=rnorm(1, mean = 1.5 - 0.5*trt, sd = 1.5 - 0.5*trt),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1), # time-varying indicator if treated or not
           x=ifelse(tp>=2, lag(x, 1) + (tp-1)/10 * rnorm(1, mean = 1, sd = 0.1), x) 
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(err=rnorm(n*max.time), 
                        y = 1 + x + trt + int + err + treated + ((tp - 2.5)^2)/10,
                        y.t = 1 + x * tp / 10 + trt + int + err + treated + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y), y.diff2 = y.t - lag(y.t)) %>% ungroup() 
  
  ## match processing
  matched.cov.dat <- match.on.cov2(dat)
  matched.lev.dat <- match.on.lev(dat)
  matched.tre.dat <- match.on.tre(dat)
  
  ## fit regression models
  my.res <- fit.mods()
  
  ## save output
  mv.same[i, ] <- my.res[[1]]
  se.same[i, ] <- my.res[[2]]
  
  
  #### processing for matching
  matched.lev.dat <- match.on.lev2(dat)
  matched.tre.dat <- match.on.tre2(dat)
  
  ## fit regression models
  my.res <- fit.mods2()
  
  ## save output
  mv.same2[i, ] <- my.res[[1]]
  se.same2[i, ] <- my.res[[2]]
}


## ----sim-time-varying-cross, echo=FALSE--------------------------------------------------------------------------------
mv.cross <- matrix(NA, ncol = 8, nrow = nsims)
se.cross <- matrix(NA, ncol = 8, nrow = nsims)

mv.cross2 <- matrix(NA, ncol = 8, nrow = nsims)
se.cross2 <- matrix(NA, ncol = 8, nrow = nsims)

for(i in 1:nsims) {

  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=0.5, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           x=rnorm(1, mean = 1.5 - 0.5*trt, sd = 1.5 - 0.5*trt),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1), # time-varying indicator if treated or not
           x=ifelse(tp>=2, lag(x, 1) + (I(trt == 1) - I(trt == 0)) * (tp-1)/10 * rnorm(1, mean = 1, sd = 0.5), x) 
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(err=rnorm(n*max.time), 
                        y = 1 + x + trt + int + err + treated + ((tp - 2.5)^2)/10,
                        y.t = 1 + x * tp / 10 + trt + int + err + treated + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y), y.diff2 = y.t - lag(y.t)) %>% ungroup() 
  
  ## match processing
  matched.cov.dat <- match.on.cov(dat) ## check
  matched.lev.dat <- match.on.lev(dat)
  matched.tre.dat <- match.on.tre(dat)
  
  ## fit regression models
  my.res <- fit.mods()
  
  ## save output
  mv.cross[i, ] <- my.res[[1]]
  se.cross[i, ] <- my.res[[2]]
  
  
  #### processing for matching
  matched.lev.dat <- match.on.lev2(dat)
  matched.tre.dat <- match.on.tre2(dat)
  
  ## fit regression models
  my.res <- fit.mods2()
  
  ## save output
  mv.cross2[i, ] <- my.res[[1]]
  se.cross2[i, ] <- my.res[[2]]
}


## ----sim-time-varying-post, echo=FALSE---------------------------------------------------------------------------------
mv.post <- matrix(NA, ncol = 8, nrow = nsims)
se.post <- matrix(NA, ncol = 8, nrow = nsims)

mv.post2 <- matrix(NA, ncol = 8, nrow = nsims)
se.post2 <- matrix(NA, ncol = 8, nrow = nsims)

for(i in 1:nsims) {

  dat <- expand.grid(id = 1:n, tp = 1:max.time) %>% arrange(id,tp) %>% group_by(id) %>%
    mutate(int=rnorm(1,0,sd=0.25), # random intercept
           p.trt=0.5, # probability of treatment
           trt=rbinom(1, 1, p.trt), # treatment
           x=rnorm(1, mean = 1.5 - 0.5*trt, sd = 1.5 - 0.5*trt),
           post=I(tp >= trt.time), # indicator of post-treatment period
           treated=I(post == 1 & trt == 1), # time-varying indicator if treated or not
           x=ifelse(tp>=2, lag(x, 1) + (tp-1)/10 * rnorm(1, mean = 1, sd = 0.1) - I(trt == 1) * I(tp>6)*(tp)/20, x)
    ) %>% 
    ungroup()
  
  dat <- dat %>% mutate(err=rnorm(n*max.time), 
                        y = 1 + x + trt + int + err + treated + ((tp - 2.5)^2)/10,
                        y.t =  1 + x * tp / 10 + trt + int + err + treated + ((tp - 2.5)^2)/10) %>%
    group_by(id) %>% mutate(y.diff = y - lag(y), y.diff2 = y.t - lag(y.t)) %>% ungroup()
  
  ## match processing
  matched.cov.dat <- match.on.cov2(dat)
  matched.lev.dat <- match.on.lev(dat)
  matched.tre.dat <- match.on.tre(dat)
  
  ## fit regression models
  my.res <- fit.mods()
  
  ## save output
  mv.post[i, ] <- my.res[[1]]
  se.post[i, ] <- my.res[[2]]
  
  
  #### processing for matching
  matched.lev.dat <- match.on.lev2(dat)
  matched.tre.dat <- match.on.tre2(dat)
  
  ## fit regression models
  my.res <- fit.mods2()
  
  ## save output
  mv.post2[i, ] <- my.res[[1]]
  se.post2[i, ] <- my.res[[2]]
}


## ----plot-second, echo=FALSE, fig.height=5, fig.width=7, fig.cap="\\label{fig:time-varying-results}Simulation Results for Time-Varying Covariate"----
makegrid(list(mv.no = mv.same, mv.yes = mv.cross, mv.eff = mv.post, se.no = se.same, se.yes = se.cross, se.eff = se.post), scenario = c("4a", "5a", "6a"), sub = c("Parallel\nEvolution", "Evolution Differs\nBy Group", "Evolution Diverges\nin Post-Period"), truth = c(1, 1, 0.85))


## ----plot-third, echo=FALSE, fig.height=5, fig.width=7, fig.cap="\\label{fig:time-varying2-results}Simulation Results for Time-Varying Covariate"----
makegrid(list(mv.no = mv.same2, mv.yes = mv.cross2, mv.eff = mv.post2, se.no = se.same2, se.yes = se.cross2, se.eff = se.post2), scenario = c("4b", "5b", "6b"), sub = c("Parallel\nEvolution", "Evolution Differs\nBy Group", "Evolution Diverges\nin Post-Period"), truth = c(1, 1, 0.87))


## ----reproducibility, echo=FALSE, eval=FALSE---------------------------------------------------------------------------
## sessionInfo()


## ---- child="appendix-confounding.Rmd", echo=FALSE, eval=FALSE---------------------------------------------------------
## NA


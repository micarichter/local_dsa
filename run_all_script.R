# Code to accompany the manuscript
# Local dynamical survival analysis and epidemic surveillance panels for 
# epidemic tracking and forecasting 
# by Micaela Richter, Eben Kenah, Prescott C. Alexander, and Grzegorz Rempala.
# last modified: July 31, 2025
# This code generates all figures in the manuscript

## Required R packages
require(deSolve)  # version 1.40
require(MASS)     # version 7.3-60.2
require(survival) # version 3.6.4
require(dplyr)    # version 1.1.4
require(tidyr)    # version 1.3.1
require(ggplot2)  # version 3.5.1
require(parallel) # version 4.4.1
require(furrr)    # version 0.3.1
require(cowplot)  # version 1.2.0

set.seed(918) # fixed seed for reproducibility
source("localDSA_functions.R") # back end functions
S_true <- read.csv("s_covs.csv") # read in coverage probability data
I_true <- read.csv("i_covs.csv") # read in coverage probability data


## Generate synthetic data
epidemic <- SEIRepidemic(
  beta = 0.5, delta = 0.3, gamma = 0.2, rhoE = 0.02, rhoI = 0.01, rhoR = 0.02, 
  tmax = 100, tstep = 0.05)

## Read in EpiCast data from Ohio
epidat <- read.csv("ohio_all.csv")

## Figure 3 --------------------------------------------------------------------
# (parameter point estimates (a) and standard errors (b) as the 5-day
# window slides across the epidemic)

res <- DSA_window(epidemic)
# A (parameter point estimates)
res[[1]]
# B (parameter standard errors)
res[[2]]

## Figure 4 --------------------------------------------------------------------
# (predictions based on a 5 day window using simulated data)
fit_synth <- function(n, epidemic, tstart, tstop) {
  
  hsdat <- HSdat(n, epidemic, tstart, tstop)
  hsest <- HSsurv(hsdat)
  DSAest <- DSAmle(hsdat, method = "L-BFGS-B")
  
  # fit and predictions
  mlesamp <- DSApred_mlesamp(DSAest)
  times_pred <- seq(attr(hsdat, "tstart"), 100 , by = 0.1)
  DSApred_mle <- DSApredict(DSAest$point$point, times_pred, mlesamp)
  
  # plot
  DSAplot(DSApred_mle[DSApred_mle$time > attr(hsdat, "tstop"), ], epidemic, xlim = c(0, 100))
  rect(xleft = tstart, xright = tstop, ybottom = -1.2, ytop = 1.2, border = NA,
       col = adjustcolor("gray", alpha = 0.5))
  legend(x = "right", box.lwd = 0.8 ,  
         legend=c("S", "E", "I", "R"),  
         fill = c("#23CE6B","#E3D888", "#79ADDC", "#9D5C63")) 
}

fit_synth(5000, epidemic, 10, 15)

## Table 2 ---------------------------------------------------------------------
# (coverage probabilities for predicted S and I curves)
S_true$type <- "S"
I_true$type <- "I"

cov_tab <- bind_rows(S_true, I_true) %>%
  pivot_wider(names_from = type, values_from = "Coverage Probability")

cov_tab1 <- cov_tab %>% 
  mutate(t_start = tstarts[as.integer(tstart)],
         pred_time = as.numeric(unlist(ptimes)),
         coverage_S = round(S, 3),
         coverage_I = round(I, 3)) %>%
  select(t_start, pred_time, coverage_S, coverage_I) %>%
  distinct(t_start, pred_time, .keep_all = TRUE) %>%
  arrange(t_start, pred_time)

ft <- flextable(cov_tab1)
ft <- set_header_labels(ft, t_start = "tstart", pred_time = "Prediction time",
                        coverage_S = "Coverage probability S", 
                        coverage_I = "Coverage probability I")

ft <- colformat_num(ft, j = c("coverage_S", "coverage_I"), digits = 3)
ft <- merge_v(ft, j = "t_start")
ft <- align(ft, j = c("t_start", "pred_time"), align = "center", part = "all")
ft <- autofit(ft)
ft <- hline(ft, i = c("6", "11", "15"))
ft

## Figure 5 --------------------------------------------------------------------
# (changing parameter estimates in EpiCast data over a sliding window)
res_epicast <- EPI_window(epidat)
# A (parameter point estimates)
res_epicast[[1]]
# B (parameter standard errors)
res_epicast[[2]]

## Figure 6 --------------------------------------------------------------------
# (empirical EpiCast data vs fitted observations in 5 day windows)
width <- 14
tstarts <- seq(40, 199, by = width)
ntstarts <- length(tstarts)
ntimes <- 281
tstart <-  40

# take a sample of 5000
dat <- epidat[sample(nrow(epidat), 5000), ]

DSApredS <- data.frame(matrix(nrow = ntimes, ncol = ntstarts))
DSApredE <- data.frame(matrix(nrow = ntimes, ncol = ntstarts))
DSApredI <- data.frame(matrix(nrow = ntimes, ncol = ntstarts))
DSApredR <- data.frame(matrix(nrow = ntimes, ncol = ntstarts))

fit_S <- data.frame(matrix(nrow = length(seq(tstart, tstart + width, by = 0.1)),
                           ncol = ntstarts))
fit_E <- data.frame(matrix(nrow = length(seq(tstart, tstart + width, by = 0.1)),
                           ncol = ntstarts))
fit_I <- data.frame(matrix(nrow = length(seq(tstart, tstart + width, by = 0.1)),
                           ncol = ntstarts))
fit_R <- data.frame(matrix(nrow = length(seq(tstart, tstart + width, by = 0.1)),
                           ncol = ntstarts))
for (i in 1:ntstarts) {
  try({
    tstart <- tstarts[i]
    
    # data for interval i, fit and obtain MLEs
    hsdat_i <- EPIdat(dat, tstart = tstart, tstop = tstart + width)
    DSAest_i <- DSAmle(hsdat_i, method = "L-BFGS-B")
    pvec <- as.numeric(exp(DSAest_i$point$point))
    beta <- pvec[1]
    delta <- pvec[2]
    gamma <- pvec[3]
    xrhoE <- pvec[4]
    xrhoI <- pvec[5]
    xrhoR <- pvec[6]
    rhoE <- xrhoE / (1 + xrhoE + xrhoI + xrhoR)
    rhoI <- xrhoI / (1 + xrhoE + xrhoI + xrhoR)
    rhoR <- xrhoR / (1 + xrhoE + xrhoI + xrhoR)
    
    # fitted intervals (i)
    fit_i <- SEIRepidemic(beta = beta, delta = delta, 
                          gamma = gamma, rhoE = rhoE, 
                          rhoI = rhoI, rhoR = rhoR, 
                          tmin = tstart, tmax = tstart + width, tstep = 0.1)
    fit_S[, i] <- fit_i$S
    fit_E[, i] <- fit_i$E
    fit_I[, i] <- fit_i$I
    fit_R[, i] <- fit_i$R
    
    # prediction intervals (i+1)
    mlesamp_i <- DSApred_mlesamp(DSAest_i)
    times_pred_i <- seq(tstart, tstart + 2 * width, by = 0.1)
    DSApred_mle_i <- DSApredict(DSAest_i$point$point, times_pred_i, mlesamp_i)
    DSApredS[, i] <- DSApred_mle_i$S
    DSApredE[, i] <- DSApred_mle_i$E
    DSApredI[, i] <- DSApred_mle_i$I
    DSApredR[, i] <- DSApred_mle_i$R
  })
}

emp <- EPIdat(dat, tstart = 0, tstop = 240)
surv <- HSsurv(emp, se.fit = FALSE)


green <- rgb(55/255, 255/255, 20/255, 1)
nblue1 <- rgb(0, 19/255, 222/255, 1)

# clean functions for 4 panel plot
out_merge <- function(time, emp, label) {
  # empirical data
  emp_df <- data.frame(
    time = time,
    empirical = emp,
    source = "empirical",
    compartment = label
  )
}

pred_merge <- function(pred, label) {
  # predicted data
  pred_list <- list()
  for (i in 1:(ntstarts - 1)) {
    times <- seq(tstarts[i] + width, tstarts[i] + 2 * width, by = 0.1)
    preds <- pred[141:281, i]
    pred_list[[i]] <- data.frame(
      time = times,
      prev = preds,
      source = paste0("prediction_", i),
      compartment = label
    )
  }
  pred_df <- bind_rows(pred_list)
}

df_E <- out_merge(surv$Eprev$time, surv$Eprev$prev, "E")
df_I <- out_merge(surv$Iprev$time, surv$Iprev$prev, "I")
df_S <- out_merge(surv$Esurv$time, surv$Esurv$surv, "S")
df_R <- out_merge(surv$Rsurv$time, 1 - surv$Rsurv$surv, "R")

pS <- pred_merge(DSApredS, "S")
pE <- pred_merge(DSApredE, "E")
pI <- pred_merge(DSApredI, "I")
pR <- pred_merge(DSApredR, "R")

emp_combined <- bind_rows(df_S, df_E, df_I, df_R)
pred_combined <- bind_rows(pS, pE, pI, pR)

emp_combined$type <- "Empirical"
pred_combined$type <- "Predicted"

emp_combined$compartment <- factor(emp_combined$compartment, levels = c("S", "E", "I", "R"))
pred_combined$compartment <- factor(pred_combined$compartment, levels = c("S", "E", "I", "R"))

fit_merge <- function(fit, label) {
  # fitted data
  fit_list <- list()
  for (i in 1:(ntstarts)) {
    times <- seq(tstarts[i], tstarts[i] + width, by = 0.1)
    fits <- fit[, i]
    fit_list[[i]] <- data.frame(
      time = times,
      prev = fits,
      source = paste0("fit_", i),
      compartment = label
    )
  }
  fit_df <- bind_rows(fit_list)
}

fS <- fit_merge(fit_S, "S")
fE <- fit_merge(fit_E, "E")
fI <- fit_merge(fit_I, "I")
fR <- fit_merge(fit_R, "R")
fit_combined <- bind_rows(fS, fE, fI, fR)
fit_combined$type <- "Fitted"
fit_combined$compartment <- factor(fit_combined$compartment, levels = c("S", "E", "I", "R"))

ggplot() +
  geom_line(data = emp_combined, aes(x = time, y = empirical, color = type), 
            linewidth = 0.5) +
  geom_line(data = fit_combined, aes(x = time, y = prev, group = source, 
                                     color = type), linewidth = 0.5) +
  facet_wrap(~ compartment, scales = "free_y") +
  geom_vline(xintercept = tstarts, color = "lightgray", linetype = "dotted") +
  xlim(40, 200) +
  scale_color_manual(values = c(
    "Empirical" = nblue1, 
    "Fitted" = green  
  )) +
  labs(x = "Time", y = "Prevalence", color = "") +
  theme_dark() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 14)
  )

## Figure 7 --------------------------------------------------------------------
# (empirical EpiCast data vs predictions in 5 day windows)
ggplot() +
  geom_line(data = emp_combined, aes(x = time, y = empirical, color = type), 
            linewidth = 0.5) +
  geom_line(data = pred_combined, aes(x = time, y = prev, group = source, 
                                      color = type), linewidth = 0.5) +
  facet_wrap(~ compartment, scales = "free_y") +
  geom_vline(xintercept = tstarts, color = "lightgray", linetype = "dotted") +
  xlim(40, 200) +
  scale_color_manual(values = c(
    "Empirical" = nblue1, 
    "Predicted" = green  
  )) +
  labs(x = "Time", y = "Prevalence", color = "") +
  theme_dark() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 14)
  )

## Figure 8 --------------------------------------------------------------------
# (LRT using predictions at 1 week)
lrt_fun <- function() {
  width <- 14
  tstarts <- seq(40, 199, by = width)
  iter <- 100
  LRdat <- data.frame(iter = NA,
                      tstart = NA,
                      H0 = NA,
                      H1 = NA,
                      H2 = NA,
                      ratio1 = NA,
                      ratio2 = NA,
                      LRstat = NA)
  
  hsdat_emp <- EPIdat(epidat, 40, 240)
  
  system.time(for (j in 1:iter) {
    print(j)
    dat_j <- epidat[sample(nrow(epidat), 5000), ]
    
    for (tstart in tstarts) {
      try({
        
        # data for window i
        hsdat_i <- EPIdat(dat_j, tstart, tstart + width)
        
        # fit, window i, obtain MLEs
        DSAest_i <- DSAmle(hsdat_i, method = "L-BFGS-B")
        pvec <- as.numeric(exp(DSAest_i$point$point))
        beta <- pvec[1]
        delta <- pvec[2]
        gamma <- pvec[3]
        xrhoE <- pvec[4]
        xrhoI <- pvec[5]
        xrhoR <- pvec[6]
        rhoE <- xrhoE / (1 + xrhoE + xrhoI + xrhoR)
        rhoI <- xrhoI / (1 + xrhoE + xrhoI + xrhoR)
        rhoR <- xrhoR / (1 + xrhoE + xrhoI + xrhoR)
        
        # calculate E, I, R at end of interval i (use as rhos in interval i+1)
        int_i_pred <- SEIRepidemic(beta = beta, delta = delta, gamma = gamma, 
                                   rhoE = rhoE, rhoI = rhoI, rhoR = rhoR, 
                                   tmin = tstart, tmax = tstart + 2 * width, tstep = 0.05)
        
        rhos <- c(int_i_pred$E[int_i_pred$time == tstart + width], 
                  int_i_pred$I[int_i_pred$time == tstart + width],
                  int_i_pred$R[int_i_pred$time == tstart + width])
        
        xrhos <- rhos / (1 - sum(rhos))
        
        # data for window i+1
        hsdat_iplus1 <- EPIdat(dat_j, tstart + width, tstart + 2 * width)
        
        # calculate loglik null hypothesis
        nullpars_i <- c(DSAest_i$point$point[1:3], log(xrhos))
        logliknull_i <- -nloglikDSA(nullpars_i, hsdat_iplus1, tstep = 0.05)
        
        # fit, obtain MLEs: full model
        DSAest_iplus1 <- DSAmle(hsdat_iplus1, method = "L-BFGS-B", init = nullpars_i)
        
        # calculate loglik alt hypothesis and LRT
        loglikalt_i <- DSAest_iplus1$ll 
        LRstat_ij <- 2 * (loglikalt_i - logliknull_i)
        
        # predict S from end of fitted interval (tstart + 2 * width) to tstart + 4 * width 
        # fitting model from tstart to tstart + 2 * width
        hsdat_full_i <- EPIdat(dat_j, tstart, tstart + 2 * width) 
        DSAest_full_i <- DSAmle(hsdat_full_i, method = "L-BFGS-B")
        
        times_pred_i <- seq(tstart, tstart + 4 * width, by = 0.1)
        DSApred_mle_i <- DSApredict(DSAest_i$point$point, times_pred_i)
        
        # predicted change in H
        S_pred1 <- DSApred_mle_i$S[DSApred_mle_i$time == tstart + 3 * width] # width/2
        S_pred2 <- DSApred_mle_i$S[DSApred_mle_i$time == tstart + 4 * width] # 2
        Hpred0 <- -log(int_i_pred$S[int_i_pred$time == tstart + 2 * width])
        Hpred1 <- -log(S_pred1)
        Hpred2 <- -log(S_pred2)
        change_pred1 <- Hpred1 - Hpred0
        change_pred2 <- Hpred2 - Hpred0
        
        # empirical change in H
        Hemp0 <- -log(1 - ecdf(hsdat_emp$Etime)(tstart + 2 * width))
        Hemp1 <- -log(1 - ecdf(hsdat_emp$Etime)(tstart + 3 * width))
        Hemp2 <- -log(1 - ecdf(hsdat_emp$Etime)(tstart + 4 * width))
        change_emp1 <- Hemp1 - Hemp0
        change_emp2 <- Hemp2 - Hemp0
      })
      
      ratio1_ij <- change_pred1 / change_emp1
      ratio2_ij <- change_pred2 / change_emp2
      
      LRdat <- rbind(LRdat, data.frame(
        iter = j,
        tstart = tstart,
        H0 = Hpred0,
        H1 = Hpred1,
        H2 = Hpred2,
        ratio1 = ratio1_ij,
        ratio2 = ratio2_ij,
        LRstat = LRstat_ij
      ))
    }
  })
}

plot_LR1 <- function() {
  LRdat %>% subset(ratio1 > 0 & LRstat > 1) %>%
    ggplot(aes(x = log(LRstat), y = abs(log(ratio1)))) +
    geom_point(shape = 4, color = "khaki3") +
    geom_smooth(color = "dodgerblue3", show.legend = FALSE) +
    geom_vline(xintercept = log(qchisq(0.975, 6))) +
    coord_cartesian(ylim = c(0, 1.5)) +
    theme_minimal() +
    labs(x = "log(LRT statistic)", 
         y = expression(abs(log(Delta*H[pred] / Delta*H[emp])))) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}

plot_LR2 <- function() {
  LRdat %>% subset(ratio2 > 0 & LRstat > 1) %>%
    ggplot(aes(x = log(LRstat), y = abs(log(ratio2)))) +
    geom_point(shape = 4, color = "khaki3") +
    geom_smooth(color = "dodgerblue3", show.legend = FALSE) +
    geom_vline(xintercept = log(qchisq(0.975, 6))) +
    coord_cartesian(ylim = c(0, 1.5)) +
    theme_minimal() +
    labs(x = "log(LRT statistic)", 
         y = expression(abs(log(Delta*H[pred] / Delta*H[emp])))) +
    theme(
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
}

plot_LR1()

## Figure 9 --------------------------------------------------------------------
# (LRT using predictions at 2 weeks)
plot_LR2()

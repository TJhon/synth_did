# The bootstrap se: Algorithm 2 of Arkhangelsky et al.
bootstrap_se = function(estimate, replications) { sqrt((replications-1)/replications) * sd(bootstrap_sample(estimate, replications)) }
bootstrap_sample = function(estimate, replications) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  weights = attr(estimate, 'weights')
  if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
  theta = function(ind) {
    if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
    else {
      weights.boot = weights
      weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
      do.call(synthdid_estimate, c(list(Y=setup$Y[sort(ind),], N0=sum(ind <= setup$N0), T0=setup$T0, X=setup$X[sort(ind), ,], weights=weights.boot), opts))
    }
  }
  bootstrap.estimates = rep(NA, replications)
  count = 0
  while(count < replications) {
    bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
    if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
  }
  bootstrap.estimates
}



# General case

# https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota_example.csv
librarian::shelf(tidyverse, synthdid, haven)

dir("R", pattern = "utils|solver|vcov|synthdid.R", full.names = T) |> map(source)

quota <- read_dta("https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota.dta?raw=true")

# panel <- quota |> rename(unit = country, time = year, treatment = quota, outcome = womparl)

quota1 <- quota |> filter(!(country %in% c("Algeria", "Kenya", "Samoa", "Swaziland", "Tanzania")))

panel <- quota1 |> rename(unit = country, time = year, treatment = quota, outcome = womparl)

quota_setup <- panel.matrices(panel, "unit", "time", "outcome", "treatment")

att_quota <- staggered_synthdid(quota_setup)

### manually search

units <- panel$unit |> unique()

sample_boot <- sample(units, replace = T)

setup_boot <- panel |> filter(unit %in% sample_boot) |> panel.matrices("unit", "time", "outcome", "treatment")



att_breaks <- att_quota$info_tau$time
att_weight <- att_quota |> pluck("tau_time") |> map(attr, "weight")# |> map(pluck, 'lambda')
names(att_weight) <- att_breaks
# att_omega <- att_quota |> pluck("tau_time") |> map(attr, "weight")# |> map(pluck, 'omega')

br <- names(setup_boot)
setup_eval <- setup_boot[[br[1]]]
att_wt_20 <- att_weight[[br[1]]]
att_wt_20$lambda <- att_wt_20$lambda[1:setup_eval$T0]
att_wt_20$omega <- att_wt_20$omega[1:setup_eval$N0] |> sum_normalize()
synthdid_estimate(setup_eval$Y, setup_eval$N0, setup_eval$T0, weights = att_wt_20)

br <- br[1:(length(br) - 1)]
tau_bt <- rep(0, length(br))
tau_wt <- setup_boot$tau_wt
for(i in 1:length(br)){
  setup_eval <- setup_boot[[br[i]]]
  att_wt_20 <- att_weight[[br[i]]]
  att_wt_20$lambda <- att_wt_20$lambda[1:setup_eval$T0]
  att_wt_20$omega <- att_wt_20$omega[1:setup_eval$N0] |> sum_normalize()
  att <- synthdid_estimate(setup_eval$Y, setup_eval$N0, setup_eval$T0, weights = att_wt_20)
  tau_bt[i] <- as.numeric(att)
}

sum(tau_bt * (tau_wt / sum(tau_wt)))

vcov_stagg <- function(att_estimate, rep_n = 50, method = "bootsrap"){

  sample_rep <- pluck(att_estimate, "tau_time") |> map(attr, "setup") |> map_dbl(pluck, "N0")

  bt_search <- function(){
    panel = pluck(att_estimate, "panel_ref")

    units <- panel$unit |> unique()
    if(method == "bootstrap"){
      sample_boot <- sample(units, replace = T)
    } else if (method == "placebo") {
      sample_boot <- sample(units, sample_rep[1], replace = T)
      print(sample_rep)
    }
    # sample_boot <- sample(units)
    setup_boot <- panel |> filter(unit %in% sample_boot) |> panel.matrices("unit", "time", "outcome", "treatment")

    att_breaks <- att_estimate$info_tau$time
    att_weight <- att_estimate |> pluck("tau_time") |> map(attr, "weight")# |> map(pluck, 'lambda')
    names(att_weight) <- att_breaks
    # att_omega <- att_quota |> pluck("tau_time") |> map(attr, "weight")# |> map(pluck, 'omega')

    br <- names(setup_boot)
    # setup_eval <- setup_boot[[br[1]]]
    # att_wt_20 <- att_weight[[br[1]]]
    # att_wt_20$lambda <- att_wt_20$lambda[1:setup_eval$T0]
    # att_wt_20$omega <- att_wt_20$omega[1:setup_eval$N0] |> sum_normalize()
    # synthdid_estimate(setup_eval$Y, setup_eval$N0, setup_eval$T0, weights = att_wt_20)

    br <- br[1:(length(br) - 1)]
    tau_bt <- rep(0, length(br))
    tau_wt <- setup_boot$tau_wt
    for(i in 1:length(br)){
      setup_eval <- setup_boot[[br[i]]]
      att_wt_20 <- att_weight[[br[i]]]
      att_wt_20$lambda <- att_wt_20$lambda[1:setup_eval$T0]
      att_wt_20$omega <- att_wt_20$omega[1:setup_eval$N0] |> sum_normalize()
      att <- synthdid_estimate(setup_eval$Y, setup_eval$N0, setup_eval$T0, weights = att_wt_20)
      tau_bt[i] <- as.numeric(att)
    }

    att_bt <- sum(tau_bt * (tau_wt / sum(tau_wt)))
    return(att_bt)
  }

  att_b = rep(NA, rep_n)
  for (att in 1:rep_n) {
    att_b[att] = suppressWarnings(try(bt_search(), silent = T))
  }
  att_b <- as.numeric(att_b)

  se_bootstrap <- sqrt(var(att_b, na.rm = T)) * sqrt((rep_n - 1) / rep_n)
  return(list(se_bootstrap = se_bootstrap, att_b))
}

se <- vcov_stagg(att_quota, method = "placebo")
se


att_estimate = att_quota


se[[2]] |> as.numeric() |> sd(na.rm = T)

se[[2]] |> var(na.rm = T) |> std()

bt_search(att_quota)

sim <- 40
bt_att <- rep(0, sim)

se <- map_dbl(1:sim, bt_search)


# stopifnot(is.null(att_lambda[["2000"]]$omega) || length(att_lambda[["2000"]]$omega) == N0)

# |> staggered_synthdid() |> pluck("att")


a <- c(7.303319, 6.078912, 9.45325, 4.854472, 4.94) |> sd()

a * 4 / 5

estimate <- att_quota$tau_time[[1]]

setup = attr(estimate, 'setup')
weights = attr(estimate, 'weights')

attr(att_quota$tau_time[[7]], "setup")$Y |> nrow()

ind <- sample(1:nrow(setup$Y), replace=TRUE)

weights.boot = weights
weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])

synthdid_estimate(setup$Y[sort(ind), ], sum(ind <= setup$N0), setup$T0, weights = weights.boot)




(c(6.719, 14.282) * att_quota$info_tau$tau_wt_time) |> sum()

(c(8.044, 13.945)* att_quota$info_tau$tau_wt_time) |> sum()

sd(c(10.36044, 10.88522))


  sqrt((2-1)/2) * sd(c(10.36044, 10.88522))
#sd(bootstrap_sample(estimate, replications))

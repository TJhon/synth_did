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
quota1 <- quota |> filter(!(country %in% c("Algeria", "Kenya", "Samoa", "Swaziland", "Tanzania")))

panel <- quota1 |> rename(unit = country, time = year, treatment = quota, outcome = womparl)

quota_setup <- panel.matrices(panel, "unit", "time", "outcome", "treatment")

att_quota <- staggered_synthdid(quota_setup)

### manually search

estimate <- att_quota$tau_time[[1]]

setup = attr(estimate, 'setup')
weights = attr(estimate, 'weights')

ind <- sample(1:nrow(setup$Y), replace=TRUE)

weights.boot = weights
weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])

synthdid_estimate(setup$Y[sort(ind), ], sum(ind <= setup$N0), setup$T0, weights = weights.boot)

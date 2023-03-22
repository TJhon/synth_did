# General case

# https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota_example.csv
librarian::shelf(tidyverse, synthdid, haven)

dir("R", pattern = "utils|solver|vcov|synthdid.R", full.names = T) |> map(source)

quota <- read_dta("https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota.dta?raw=true")
quota1 <- quota |> filter(!(country %in% c("Algeria", "Kenya", "Samoa", "Swaziland", "Tanzania")))

panel <- quota1 |> rename(unit = country, time = year, treatment = quota, outcome = womparl)

quota_setup <- panel.matrices(quota1, "country", "year", "womparl", "quota")

att <- staggered_synthdid(quota_setup)


data_0 <-
  # panel |>
  quota_setup[[1]]$data_ref
  # group_by(unit) |>
  # mutate(
  #   treated = max(treatment),
  #   ty = ifelse(treatment == 1, time, NA),
  #   tyear = ifelse(treated == 1, min(ty, na.rm = T), NA)
  # ) |>
  # ungroup() |>
  # arrange(treated, unit, time) |>
  # # mutate(across(c(ty, tyear), replace_na, 0)) |>
  # mutate(across(where(is.numeric), replace_na, 0))
break_points <- data_0 |> pull(tyear) |> sort() |> unique()
break_points <- break_points[-1]

estimate <- att
quota_setup <- quota_setup
se_search <- function(estimate, quota_setup) {
  start_time <- proc.time()
  att_tau <- estimate |> pluck(3)
  data_0 <- quota_setup[[1]]$data_ref

  uniqID <- unique(data_0$unit)
  time_break <- estimate$info_tau$time

  # ind <- 12
  # id <- 2
  N <- length(uniqID)

  jackknife_se <- function(ind, id){
    # ind <- 114
    # year <- 2
    # bt <- break_points[year]
    # weights att
    w_att <- att_tau |> pluck(id) |> attr("weight")
    lambda_aux <- w_att$lambda
    omega_aux <- (w_att$omega)[-ind] |> sum_normalize()
    data_aux <-
      data_0 |> filter(unit != uniqID[ind]) |>
      filter(tyear %in% c(0, time_break[id]))


    # Ynt = el minimo de observaciones por individuo tratado
    # Yng = numero de individuos
    # Yn = numero de fechas
    # yn <- data_aux |> count(time) |> nrow()
    yng <-  data_aux |> count(unit) |> nrow()
    ynt <- data_aux |> count(unit) |> pull(n) |> min()
    yntr <- sum(data_aux$tyear == time_break[id]) / ynt

    npre <- data_aux |> filter(time < time_break[id]) |> count(time) |> nrow()
    npost <- ynt - npre

    Y_aux <-
      data_aux |> select(time, unit, outcome) |> pivot_wider(names_from = time, values_from = outcome) |> select(!unit) |> as.matrix()

    tau_aux <- t(c(-omega_aux, rep(1 / yntr, yntr))) %*% Y_aux %*% c(-lambda_aux, rep(1 / npost, npost))
    tau_wt_aux <- yntr * npost

    tibble(unit = uniqID[ind], time = time_break[id], tau_aux = tau_aux[1], tau_wt_aux)
  }

  # library(progress)
  # try(
  #   pb <- progress::progress_bar$new(format = "[:bar] :percent Elapsed: :elapsed eta: :eta", total = 10)
  # )

  total_iter <- expand.grid(1:length(uniqID), 1:length(time_break))

  result <- tibble()
  for (i in 1:nrow(total_iter)) {
    row_eval <- total_iter[i, ]
    result <- bind_rows(result, jackknife_se(row_eval[[1]], row_eval[[2]]))
    # try(pb$tick())
  }

  # for (i in 1:length(uniqID)) {
  #   for (j in 1:2) {
  #     result <- bind_rows(result, jackknife_se(i, j))
  #   }
  # }
  result_att <- result |> arrange(unit) |> distinct() |> group_by(unit) |>
    mutate(tau_wt = tau_wt_aux / sum(tau_wt_aux), att_aux = tau_aux * tau_wt)

  att_aux <-
    result_att |>
    summarise(att_aux = sum(att_aux))
  se_jackknife <- (((N-1)/N) * (N - 1) * var(att_aux$att_aux)) |> sqrt()
  class(se_jackknife) = "se_jackknife"
  # class(estimate) = 'synthdid_estimate'
  attr(se_jackknife, 'estimator') = "synthdid_estimate"
  attr(se_jackknife, "ATT_info") <- list("att_aux" = att_aux, "att_all_info" = result_att)

  time_diff <- proc.time() - start_time
  print(paste("Elapsed Time", time_diff[3]))
  return(se_jackknife)
}

jk <- se_search(att, quota_setup)

attr(jk, "ATT_info")

uniqID <- unique(data_0$unit)
ind <- 10
year <- 1
bt <- break_points[year]

t2003_wg <- att |> pluck(3, year) |> attr("weight")
lambda_aux <- t2003_wg$lambda
omega_aux <- (t2003_wg$omega)[-ind] |> sum_normalize()

data_aux <-
  data_0 |> filter(unit != uniqID[ind]) |>
    filter(tyear %in% c(0, bt))


# Ynt = el minimo de observaciones por individuo tratado
# Yng = numero de individuos
# Yn = numero de fechas
# yn <- data_aux |> count(time) |> nrow()
yng <-  data_aux |> count(unit) |> nrow()
ynt <- data_aux |> count(unit) |> pull(n) |> min()
yntr <- sum(data_aux$tyear == bt) / ynt

npre <- data_aux |> filter(time < bt) |> count(time) |> nrow()
npost <- ynt - npre

Y_aux <-
  data_aux |> select(time, unit, outcome) |> pivot_wider(names_from = time, values_from = outcome) |> select(!unit) |> as.matrix()

tau_aux <- t(c(-omega_aux, rep(1 / yntr, yntr))) %*% Y_aux %*% c(-lambda_aux, rep(1 / npost, npost))
tau_wt_aux <- yntr * npost
tau_aux

# tibble(unit = uniqID[ind], tau_aux = tau_aux[1], tau_wt_aux)


setup <- panel.matrices(california_prop99)

att_1 <- synthdid_estimate(setup$Y, setup$N0, setup$T0)

placebo <- vcov(att[[3]][[1]], method = "placebo", replications = 10)
placebo2 <- vcov(att[[3]][[2]], method = "placebo", replications = 10)

# sqrt(placebo)*

(sqrt(placebo) + sqrt(placebo2))


(placebo + placebo2) |> sqrt()


# vcov(att_1, method = "placebo", replications = 50)

all_setup <- quota |> panel.matrices("country", "year", "womparl", "quota")

placebo <- all_setup[[1]]

placebo$N0

setup = attr(estimate, 'setup')
# opts = attr(estimate, 'opts')
# weights = attr(estimate, 'weights')

Y <- placebo$Y

N1 = nrow(Y) - placebo$N0

ind <- sample(1:placebo$N0)
N0 = length(ind) - N1

synthdid_estimate(Y[ind, ], N0, placebo$T0)



all_att <- staggered_synthdid(all_setup)

(all_att$tau_time[[1]] |> attr("weight"))$omega



vce <- all_att$tau_time |> map(vcov, method = "bootstrap", replications = 50)

a <- vce |> map_dbl(sqrt)# |> sd(na.rm = T)

(all_att$info_tau$tau_wt_time * a) |> sd(na.rm = T) |> sqrt()

a |> var(na.rm = T) |> sqrt()

sqrt(a * 49 / 50)



placebo_se = function(estimate, replications = 50) {
  setup = attr(estimate, 'setup')
  opts = attr(estimate, 'opts')
  weights = attr(estimate, 'weights')
  N1 = nrow(setup$Y) - setup$N0
  if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
  theta = function(ind) {
    N0 = length(ind)-N1
    weights.boot = weights
    weights.boot$omega = sum_normalize(weights$omega[ind[1:N0]])
    do.call(synthdid_estimate, c(list(Y=setup$Y[ind,], N0=N0,  T0=setup$T0,  X=setup$X[ind, ,], weights=weights.boot), opts))
  }
  # replicate(replications, theta(sample(1:setup$N0)))
  # sqrt( (replications-1)/replications) *
  sd(replicate(replications, theta(sample(1:setup$N0))), na.rm =F)
}

placebo_se(att_1, 50)

all_setup



a <- all_att$tau_time |> map_dbl(placebo_se, replications = 10)

a

(all_att$info_tau$tau_wt_time * a) |> sum()


b <- att$tau_time |> map_dbl(placebo_se)

sim <- function(z){
  print(z)
  c <- quota1 |> pull(country) |> unique() |> sample(112)
  a <- quota1 |>
    filter(country %in% c) |>
    panel.matrices("country", "year", "womparl", "quota") |>
    staggered_synthdid()

  return(a$att)
}
a <- map_dbl(1:50, sim)

(a^2 |> sd()) |> sqrt()






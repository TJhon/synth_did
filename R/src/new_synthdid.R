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
  quota_setup[[1]]$data_ref

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

  total_iter <- expand.grid(1:length(uniqID), 1:length(time_break))

  result <- tibble()
  for (i in 1:nrow(total_iter)) {
    row_eval <- total_iter[i, ]
    result <- bind_rows(result, jackknife_se(row_eval[[1]], row_eval[[2]]))
  }

  result_att <- result |> arrange(unit) |> distinct() |> group_by(unit) |>
    mutate(tau_wt = tau_wt_aux / sum(tau_wt_aux), att_aux = tau_aux * tau_wt)

  att_aux <-
    result_att |>
    summarise(att_aux = sum(att_aux))
  se_jackknife <- (((N-1)/N) * (N - 1) * var(att_aux$att_aux)) |> sqrt()

  return(se_jackknife)
}

jk <- se_search(att, quota_setup)

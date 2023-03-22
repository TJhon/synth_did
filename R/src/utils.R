# General case

# https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota_example.csv
librarian::shelf(tidyverse, synthdid, haven)

dir("R", pattern = "utils|solver|vcov|synthdid.R", full.names = T) |> map(source)

quota <- read_dta("https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota.dta?raw=true")
# quota_sample <- read_csv("https://github.com/d2cml-ai/synthdid/raw/stata_review/vignettes/data/quota_example.csv", show_col_types = F)

# data <- quota_sample |>
#   rename(unit = 1, time = 2, outcome = 3, quota = 4) |>
#   select(1:4)



panel_matrices = function(panel, unit = 1, time = 2, outcome = 3, treatment = 4, treated.last = TRUE) {
  # TODO: add support for covariates X, i.e. could keep all other columns
  keep = c(unit, time, outcome, treatment)
  if (!all(keep %in% 1:ncol(panel) | keep %in% colnames(panel))) {
    stop("Column identifiers should be either integer or column names in `panel`.")
  }

  if(typeof(unit) == "character"){
    panel <- panel |> select(all_of(c(unit, time, outcome, treatment)))
  }else{
    panel <- panel[, c(unit, time, outcome, treatment)]
  }
  panel <- panel |> rename(unit = 1, time = 2, outcome = 3, treatment = 4)
  unit = "unit"; time = "time"; outcome = "outcome"; treatment = "treatment"

  data_0 <-
    panel |>
    group_by(unit) |>
    mutate(
      treated = max(treatment),
      ty = ifelse(treatment == 1, time, NA),
      tyear = ifelse(treated == 1, min(ty, na.rm = T), NA)
    ) |>
    ungroup() |>
    arrange(treated, unit, time) |>
    # mutate(across(c(ty, tyear), replace_na, 0)) |>
    mutate(across(where(is.numeric), replace_na, 0))
  break_points <- data_0 |> pull(tyear) |> sort() |> unique()
  break_points <- break_points[-1]

  get_panel <- function(panel){
    # panel = panel[keep]
    if (!is.data.frame(panel)){
      stop("Unsupported input type `panel.`")
    }

    if (anyNA(panel)) {
      stop("Missing values in `panel`.")
    }

    if (length(unique(panel[, treatment])) == 1) {
      stop("There is no variation in treatment status.")
    }
    if (!all(panel[, treatment] %in% c(0, 1))) {
      stop("The treatment status should be in 0 or 1.")
    }

    val <- as.vector(table(panel[, unit], panel[, time]))
    if (!all(val == 1)) {
      stop("Input `panel` must be a balanced panel: it must have an observation for every unit at every time.")
    }
    panel = panel[order(panel[, unit], panel[, time]), ]
    num.years = length(unique(panel[, time]))
    num.units = length(unique(panel[, unit]))
    Y = matrix(panel[,outcome], num.units, num.years, byrow = TRUE,
               dimnames = list(unique(panel[,unit]), unique(panel[,time])))
    W = matrix(panel[,treatment], num.units, num.years, byrow = TRUE,
               dimnames = list(unique(panel[,unit]), unique(panel[,time])))
    w = apply(W, 1, any)                         # indicator for units that are treated at any time
    T0 = unname(which(apply(W, 2, any))[1]-1)    # last period nobody is treated
    N0 = sum(!w)

    if(! (all(W[!w,] == 0) && all(W[,1:T0] == 0) && all(W[w, (T0+1):ncol(Y)]==1))) {
      stop("The package cannot use this data. Treatment adoption is not simultaneous.")
    }

    unit.order = if(treated.last) { order(W[,T0+1], rownames(Y)) } else { 1:nrow(Y) }
    list(Y = Y[unit.order, ], N0 = N0, T0 = T0, W = W[unit.order, ])

  }
  multiple_breaks <- function(data_ref, b_p){
    panel <- data_ref |> filter(tyear %in% c(b_p, 0)) |>
      select(unit, time, outcome, treatment) |>
      as.data.frame() |>
      get_panel()
    nt <- dim(panel$Y)
    panel$weight_mult <- nt - c(panel$N0,  panel$T0)
    panel$tau_wt <- panel$weight_mult[1] * panel$weight_mult[2]
    return(panel)
  }
  if(length(break_points) < 2){
    return(get_panel(panel))
  }else{
    all_setup <- map(break_points, multiple_breaks, data = data_0)
    names(all_setup) <- break_points
    tau_wt <- map_dbl(all_setup, pluck, "tau_wt")
    all_setup$tau_wt <- tau_wt / sum(tau_wt)
    return(all_setup)
  }
}

# drop if country=="Algeria" | country=="Kenya" | country=="Samoa" |
#   > country=="Swaziland" | country=="Tanzania"


quota1 <- quota |> filter(!(country %in% c("Algeria", "Kenya", "Samoa", "Swaziland", "Tanzania")))

quota_setup <- panel.matrices(quota1, "country", "year", "womparl", "quota") -> s2

att <- staggered_synthdid(quota_setup)

att_weight <- attr(att[[3]][[2]], "weight")
att_omega1 <- att_weight$omega
att_lambda1 <- att_weight$lambda


setup_1 <- quota1 |> filter(country!=country[10]) |> panel.matrices("country", "year", "womparl", "quota")

synthdid_estimate(setup_1[[1]]$Y, setup_1[[1]]$N0, setup_1[[1]]$T0, weights = att_weight)


nice <- function(){
  att_se_2 <- function(variables) {

    s1 <- att$tau_time[[2]]

    setup <- attr(s1, "setup")
    weight <- attr(s1, "weight")
    weight$lambda
    weight$omega
    dim_y <- dim(setup$Y)
    Y <- setup$Y[-dim_y[1], ] # porque estoy evaluando la ultima unidad de tratamiento, recordar que sdid de stata es una de O^4,

    omega_aux <- sum_normalize(weight$omega)
    lambda_aux <- weight$lambda


    N1 <- dim_y - c(length(omega_aux) + 1, length(lambda_aux))

    att_se <-  t(c(-omega_aux, rep(1 / N1[1], N1[1]))) %*% (Y) %*% c(-lambda_aux, rep(1 / N1[2], N1[2]))
  }

  tau_aux <- (c(att$tau_time[[1]] |> as.numeric(), att_se))

  tau_wt_aux <- c(att$info_tau$tau_wt[1], prod(N1))

  tau_wt_aux <- tau_wt_aux / sum(tau_wt_aux)


  tau_aux %*% tau_wt_aux # 6.967746483 1.601097026 ATT_aux

}

country <- quota1$country |> unique()
# quota $country |> unique()

set2 <- quota1 |> filter(country!=country[10]) |> panel.matrices("country", "year", "womparl", "quota")

s1 <- att$tau_time[[1]]

setup <- attr(s1, "setup")
weight <- attr(s1, "weight")
N1 <- setup$Y |> dim() - c(length(omega_aux) + 1, length(lambda_aux))
t(c(-weight$omega, rep(1 / 2, 2))) %*% (setup$Y) %*% c(-weight$lambda, rep(1 / 14, 14))


attr(att$tau_time[[1]], "weight")$lambda

rownames(set2[[2]]$Y)
setup$Y |> rownames()

tau_eval <- function(tau_n){

  att1 <- att$tau_time[[tau_n]]

  setup_tau <- attr(att1, "setup")

  Y_tau <- setup_tau$Y
  N1 <- setup_tau$N0
  T1 <- setup_tau$T0

  unit_drop <- nrow(Y_tau)

  ss <- c()
  for (i in 1:unit_drop) {
    # print(i)
    eval <- synthdid_estimate(Y_tau[-i, ], N1, T1)
    ss[i] <- eval |> as.numeric()
  }
  return(ss)
}

s1 <- tau_eval(1)
s2 <- tau_eval(2)


s2


synthdid_estimate(Y_tau[-1, ], N1, T1)


vcov(att$tau_time[[1]], method = "jackknife")
vcov(att$tau_time[[2]], method = "jackknife")

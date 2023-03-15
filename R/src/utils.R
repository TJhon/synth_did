# General case

# https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota_example.csv
librarian::shelf(tidyverse, synthdid, haven)

quota <- read_dta("https://github.com/d2cml-ai/synthdid/blob/stata_review/vignettes/data/quota.dta?raw=true")
quota_sample <- read_csv("https://github.com/d2cml-ai/synthdid/raw/stata_review/vignettes/data/quota_example.csv", show_col_types = F)

ss <- quota_sample |>
  rename(unit = 1, time = 2, outcome = 3, quota = 4) |>
  select(1:4)

ss |>
  select(1:3) |>
  pivot_wider(names_from = time, values_from = outcome, values_fill = 0) |>
  pivot_longer(!unit, values_to = "outcome", names_to = "time", names_transform = as.numeric) |>
  left_join(ss) |> mutate(across(where(is.numeric), replace_na, 0))


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
quota_setup <- panel_matrices(quota, "country", "year", "womparl", "quota") -> s2


if(is.null(quota_setup$tau_wt)){
  stop()
}

tau_wt <- quota_setup$tau_wt
setup_time <- quota_setup[1:(length(quota_setup) - 1)]

all_att <- function(time, setup){
  att <- synthdid_estimate(setup[[time]]$Y, setup[[time]]$N0, setup[[time]]$T0)
  return(att)
}

break_times <- setup_time |> names()

all_att(break_times[1], setup_time)

# setup_time[[break_times[1]]]$Y
tau <- map(break_times, all_att, setup_time)

tau_dbl <- tau |> map_dbl(as.double)

tau_wt_time <- tau_wt / sum(tau_wt)
att <- tau_dbl %*% tau_wt_time |> as.double()

info_tau <- tibble(time = break_times, tau_wt_time, tau_dbl)

print(att)

return(list(att = att, info_tau = info_tau, tau_time <- tau))


info <- staggered_synthdid(quota_setup)
info |> glimpse()

info$att

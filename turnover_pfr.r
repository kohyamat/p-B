# Load packages -----------------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cowplot)

# Functions ---------------------------
# biomass estimation based on Niiyama et al. (2010)
biomass <- function(dbh, component = "agb") {
  h <- 1 / (1 / (1.61 * dbh) + 1 / 69.) # tree height (m) from dbh (cm)
  w_s <- 0.036 * (dbh * dbh * h)^1.01 # trunk plus branch stem (kg in dry mass)
  w_l <- 1 / (1 / (0.108 * w_s^0.75) + 1 / 105) # leaves (kg in dry mass)
  w_r <- 0.023 * dbh^2.59 # coarse root (kg in dry mss)
  if (component == "agb") {
    return(w_s + w_l)
  } else if (component == "all") {
    return(w_s + w_l + w_r)
  } else if (component == "stem") {
    return(w_s)
  } else if (component == "leaf") {
    return(w_l)
  } else if (component == "root") {
    return(w_r)
  } else {
    stop("Invalid 'component'. Expected one of: 'all', 'agb', 'stem', 'leaf', 'root'")
  }
}


# common turnover rate
turnover <- function(y, z, t) {
  f <- function(rho) {
    sum(y * exp(-rho * t) - z)
  }
  df <- function(rho) {
    sum(-t * y * exp(-rho * t))
  }
  # Newton-Raphson iteration
  rho <- 0.02
  precision <- 1.0e-12 # to stop iteration
  change <- precision + 1.0

  while (change > precision) {
    rho2 <- rho - f(rho) / df(rho)
    change <- abs(rho2 - rho)
    rho <- rho2
  }
  return(rho)
}


# production rate
productivity <- function(dbh1, dbh2, w1, w2, t, dbh_min, mass_min, area) {
  si <- ifelse(dbh1 >= dbh_min & dbh2 >= dbh_min, 1, 0) # survival
  di <- ifelse(dbh1 >= dbh_min & dbh2 < dbh_min, 1, 0) # death
  ri <- ifelse(dbh1 < dbh_min & dbh2 >= dbh_min, 1, 0) # recruitment

  Ns0 <- sum(si, na.rm = TRUE)
  N0 <- Ns0 + sum(di, na.rm = TRUE)
  NT <- Ns0 + sum(ri, na.rm = TRUE)
  Bs0 <- sum(si * w1, na.rm = TRUE)
  BsT <- sum(si * w2, na.rm = TRUE)
  B0 <- Bs0 + sum(di * w1, na.rm = TRUE)
  BT <- BsT + sum(ri * w2, na.rm = TRUE)
  # period-mean biomass and abundance
  Nw <- ifelse(NT != N0, (NT - N0) / log(NT / N0), N0)
  N <- Nw / area # (per ha)
  Bw <- ifelse(BT != B0, (BT - B0) / log(BT / B0), B0)
  B <- Bw / area

  # Standardized maximum tree mass for initial population
  W_max <- as.numeric(quantile(w1[ri != 1], 0.99)) # Mg

  # turnover rates
  r <- turnover(si + ri, si, t)
  m <- turnover(si + di, si, t)
  p <- turnover(w2, si * w1, t)
  l <- turnover(w1, si * w1, t)

  # absolute productivity (Mg per ha per year)
  P <- p * B
  P_simple <- sum(((si + ri) * w2 - si * w1) / t)
  P_simple_Clark <- sum(si * (w2 - w1) / t + ri * (w2 - mass_min) / t)
  P_simple <- P_simple / area
  P_simple_Clark <- P_simple_Clark / area # cf. Clark et al. (2001, Ecol. Appl.)

  return(list(
    "B" = B,
    "N" = N,
    "W_max" = W_max,
    "p" = p,
    "l" = l,
    "r" = r,
    "m" = m,
    "P" = P,
    "P_simple" = P_simple,
    "P_simple_Clark" = P_simple_Clark
  ))
}


# Load data ---------------------------
df1 <- read.csv("data/pfr_observed.csv.gz", header = TRUE)
df2 <- read.csv("data/pfr_identity_free.csv.gz", header = TRUE)

# Boundary size in dbh (cm)
dbh_min <- 1.

data_preparation <- function(d) {
  # remove rows if both dbh1 and dbh2 lower than dbh_min
  d <- filter(d, dbh1 >= dbh_min | dbh2 >= dbh_min)

  # add some columns to the dataframe
  d <- d %>%
    mutate(
      dbh1 = ifelse(dbh1 >= dbh_min, dbh1, 0),
      dbh2 = ifelse(dbh2 >= dbh_min, dbh2, 0),
      w1 = biomass(dbh1) / 1000, # tree biomass in Mg
      w2 = biomass(dbh2) / 1000,
    )

  # Select species with NsT >= 100
  counts <- table(pull(filter(d, dbh1 >= dbh_min & dbh2 >= dbh_min), Cd))
  sp_list <- names(counts)[counts >= 100]

  # rename the remaining species to 'Others'
  d <- d %>%
    mutate(
      Cd = ifelse(Cd %in% sp_list, as.character(Cd), "Others"),
      Cd = factor(Cd, levels = c(sp_list, "Others"))
    )

  # sort Cd to make 'Others' to be the last
  sp_order <- order(d$Cd)
  d <- d[sp_order, ]
  return(d)
}

df1 <- data_preparation(df1)
df2 <- data_preparation(df2)

# Estimate production rate ------------
# for each Cd
dbh_min <- 1.
mass_min <- biomass(2) / 1000 # tree biomass in Mg at dbh = 2 cm
area <- 50 # total area in ha

res1 <- df1 %>%
  group_nest(Cd) %>%
  mutate(y = purrr::map(data, ~ as_tibble(
    with(., productivity(dbh1, dbh2, w1, w2, t, dbh_min, mass_min, area))
  ))) %>%
  select(-data) %>%
  unnest(cols = y)

res2 <- df2 %>%
  group_nest(Cd) %>%
  mutate(y = purrr::map(data, ~ as_tibble(
    with(., productivity(dbh1, dbh2, w1, w2, t, dbh_min, mass_min, area))
  ))) %>%
  select(-data) %>%
  unnest(cols = y)


# write results
outdir <- "output"
dir.create(outdir, showWarnings = FALSE)
write.csv(res1, file.path(outdir, "res_pfr_observed.csv"), row.names = FALSE)
write.csv(res2, file.path(outdir, "res_pfr_identity_free.csv"), row.names = FALSE)


# Plot productivity ~ biomass ---------
# Cd-sum plot-level productivity (Mg per ha per year)
Ps1 <- with(res1, c(sum(p * B), sum(P_simple), sum(P_simple_Clark)))
Ps2 <- with(res2, c(sum(p * B), sum(P_simple), sum(P_simple_Clark)))
names(Ps1) <- c("P", "P_simple", "P_simple_Clark")
names(Ps2) <- c("P", "P_simple", "P_simple_Clark")
print("Observed")
print(Ps1)
print("Identity free")
print(Ps2)

# fit the linear regression model
fit1 <- lm(log(p) ~ log(B), data = filter(res1, Cd != "Others"))
pval1 <- summary(fit1)$coefficients[2, 4]
print("Observed")
print(summary(fit1))

fit2 <- lm(log(p) ~ log(B), data = filter(res2, Cd != "Others"))
pval2 <- summary(fit2)$coefficients[2, 4]
print("Identity-free")
print(summary(fit2))

# draw a plot
xbreaks = c(.01, 1, 100)
p1 <- ggplot(filter(res1, Cd != "Others"), aes(x = B, y = p)) +
  geom_point(alpha = .5) +
  geom_point(data = filter(res1, Cd == "Others"), shape = 18, size = 3,  alpha = .5) +
  annotate(
    "curve",
    x = 20, y = .12, xend = 45.4, yend = .028,
    curvature = -.25,
    arrow = arrow(length = unit(2, "mm"))
  ) +
  annotate("text", x = 10, y = .2, label = "rare species\ncombined", hjust = .5) +
  scale_x_log10(limits = c(.0005, 100), breaks = xbreaks, labels = xbreaks) +
  scale_y_log10(limits = c(.001, .5)) +
  labs(
    title = "(a) observed",
    x = parse(text = "Biomass*','~B~(Mg %.% ha^{-1})"),
    y = parse(text = "Relative~productivity*','~p~(year^{-1})")
  )

p2 <- ggplot(filter(res2, Cd != "Others"), aes(x = B, y = p)) +
  geom_point(alpha = .5) +
  geom_point(data = filter(res2, Cd == "Others"), shape = 18, size = 3,  alpha = .5) +
  scale_x_log10(limits = c(.0005, 100), breaks = xbreaks, labels = xbreaks) +
  scale_y_log10(limits = c(.001, .5)) +
  labs(
    title = "(b) identity-free",
    x = parse(text = "Biomass*','~B~(Mg %.% ha^{-1})"),
    y = parse(text = "Relative~productivity*','~p~(year^{-1})")
  )

if (pval1 < .05) {
  p1 <- p1 + geom_smooth(method = "lm", formula = y ~ x)
}

if (pval2 < .05) {
  p2 <- p2 + geom_smooth(method = "lm", formula = y ~ x)
}

p <- cowplot::plot_grid(p1, p2)

# save a plot
outdir <- "figs"
dir.create(outdir, showWarnings = FALSE)
ggsave(p, file = file.path(outdir, "p-B.png"), width = 8.5, height = 4)

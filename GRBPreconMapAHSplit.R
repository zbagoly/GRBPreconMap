library(mgcv)
library(ggplot2)
library(mapproj)
library(dplyr)
library(patchwork)

# Global theme settings
theme_set(theme_bw(base_size = 14))

# ---------------------------------------------------------
# Config and Helpers
# ---------------------------------------------------------
F_MAIN <- "all.csv"
SCENARIOS <- list(merged = "blue.csv", grb_a = "bluea.csv", grb_h = "blueh.csv")

# Quick loader for GRB data
read_grb <- function(fname) {
  if (!file.exists(fname)) {
    warning(paste("File missing:", fname, "- Generating dummy set."))
    set.seed(123)
    n_sim <- if (grepl("all", fname)) 1000 else 100
    return(data.frame(l = runif(n_sim, -180, 180), b = rnorm(n_sim, 0, 30)))
  }
  
  d <- read.csv(fname, header = FALSE)
  if (ncol(d) < 2) stop("Input must have 2 cols (l, b)")
  
  colnames(d)[1:2] <- c("l", "b")
  # Fix longitudes
  d$l <- ifelse(d$l > 180, d$l - 360, d$l)
  return(d[, c("l", "b")])
}

# Convert Galactic (l, b) to Cartesian unit vectors
to_cartesian <- function(l, b) {
  rl <- l * (pi / 180)
  rb <- b * (pi / 180)
  cbind(cos(rb) * cos(rl), cos(rb) * sin(rl), sin(rb))
}

# ---------------------------------------------------------
# Load Datasets
# ---------------------------------------------------------
dat_all <- read_grb(F_MAIN)
scenarios <- lapply(SCENARIOS, read_grb)
store_res <- list()

# Load background for ROC (Healpix)
if (file.exists("hpixlb5.csv")) {
  bg_raw <- as.matrix(read.csv("hpixlb5.csv", header = FALSE))
  bg_pts <- data.frame(l = bg_raw[, 1], b = bg_raw[, 2])
} else {
  stop("Background file 'hpixlb5.csv' is missing.")
}

# ---------------------------------------------------------
# Core Analysis Loop
# ---------------------------------------------------------
for (s_id in names(SCENARIOS)) {
  cat(sprintf("\nProcessing: %s ...\n", s_id))
  
  dat_sub <- scenarios[[s_id]]
  
  # Prepare vectors
  vec_all <- to_cartesian(dat_all$l, dat_all$b)
  vec_sub <- to_cartesian(dat_sub$l, dat_sub$b)

  # --- 1. Adaptive KDE (Abramson) ---
  
  # Log-Sum-Exp for stability
  lse <- function(x) {
    m <- max(x)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(x - m)))
  }

  # vMF Normalization constant (log)
  # Approx for high kappa, exact for low
  log_C_vmf <- function(k) {
    if (k < 1e-6) return(-log(4 * pi))
    if (k > 50) return(log(k) - log(2 * pi) - k)
    log(k) - log(4 * pi * sinh(k))
  }

  # MLCV Objective (maximize LOO likelihood)
  get_cv_score <- function(k, x) {
    N <- nrow(x)
    S <- tcrossprod(x)
    lC <- log_C_vmf(k)
    
    # Kernel matrix in log space
    lk <- S * k + lC
    diag_val <- k + lC
    
    # Cap to avoid overflow
    lk_safe <- pmin(lk, 700)
    kerns <- exp(lk_safe)
    
    r_sums <- rowSums(kerns)
    self <- exp(min(diag_val, 700))
    
    # LOO calculation
    loo <- r_sums - self
    loo[loo <= 0] <- 1e-300
    
    # Negative Log Likelihood
    nll=-sum(log(loo) - log(N - 1))
    return(nll)
  }

  solve_pilot <- function(x) {
    optimize(get_cv_score, interval = c(5, 700), x = x)$minimum
  }

  # Abramson's variable bandwidths
  get_adaptive_k <- function(x, kp) {
    S <- tcrossprod(x)
    lC <- log_C_vmf(kp)
    expon <- pmin(S * kp + lC, 700)
    
    # Pilot estimates
    p_dens <- rowMeans(exp(expon))
    g <- exp(mean(log(p_dens[p_dens > 0])))
    
    # Scaling law
    loc_k <- kp * (p_dens / g)
    pmin(pmax(loc_k, 1.0), 1000.0)
  }

  # Compute pilots and local kappas
  kp_all <- solve_pilot(vec_all)
  kp_sub <- solve_pilot(vec_sub)
  
  k_vec_all <- get_adaptive_k(vec_all, kp_all)
  k_vec_sub <- get_adaptive_k(vec_sub, kp_sub)

  # Grid projection
  project_map <- function(g_vec, d_vec, k_local) {
    c_mat <- tcrossprod(g_vec, d_vec)
    lC_v <- sapply(k_local, log_C_vmf)
    
    t1 <- sweep(c_mat, 2, k_local, "*")
    l_kern <- sweep(t1, 2, lC_v, "+")
    
    rowMeans(exp(pmin(l_kern, 700)))
  }

  # Setup grid
  g_l <- seq(-180, 180, len = 360)
  g_b <- seq(-89, 89, len = 180)
  pmap <- expand.grid(l = g_l, b = g_b)
  vec_grid <- to_cartesian(pmap$l, pmap$b)

  # Render maps
  pmap$d_bg <- project_map(vec_grid, vec_all, k_vec_all)
  pmap$d_sig <- project_map(vec_grid, vec_sub, k_vec_sub)

  # Normalize for viz (0.5 max target)
  norm_h <- function(v, t=0.5) if(max(v)==0) v else (v/max(v))*t
  pmap$bg_n <- norm_h(pmap$d_bg)
  pmap$sig_n <- norm_h(pmap$d_sig)

  # --- 2. Bayesian Logic (GAM) ---
  
  # Binary labeling
  s_all <- paste(sprintf("%.5f", dat_all$l), sprintf("%.5f", dat_all$b))
  s_sub <- paste(sprintf("%.5f", dat_sub$l), sprintf("%.5f", dat_sub$b))
  dat_all$is_sub <- as.integer(s_all %in% s_sub)

  # Find best 'k' via AIC
  k_opts <- c(10:80, 82, 84, 86, 90, 100, 110, 130, 160, 200)
  aic_res <- data.frame()
  
  cat("   Optimizing splines...\n")
  for (kv in k_opts) {
    try({
      m <- gam(is_sub ~ s(b, l, bs = "sos", k = kv), 
               data = dat_all, family = binomial, method = "REML", gamma = 1.0)
      aic_res <- rbind(aic_res, data.frame(k = kv, AIC = AIC(m)))
    }, silent = TRUE)
  }
  
  # Select best k where AIC drop stabilizes
  best_k <- aic_res$k[max(which(((min(aic_res$AIC) + 0.1) > aic_res$AIC)))]
  
  # Save AIC plot
  pa <- ggplot(aic_res, aes(k, AIC)) + 
    geom_line(col="red") + geom_vline(xintercept=best_k, lty=2) +
    labs(title=paste("AIC Select:", s_id))
  ggsave(paste0("AIC_", s_id, ".pdf"), pa, width=8, height=5)

  # Fit final model
  mod_final <- gam(is_sub ~ s(b, l, bs = "sos", k = best_k), 
                   data = dat_all, family = binomial, method = "REML")
  
  pmap$prob <- predict(mod_final, newdata = pmap, type = "response")

  # --- 3. Precon Calc ---
  pmap$precon <- pmap$prob * pmap$d_bg
  pmap$precon_n <- norm_h(pmap$precon)
  
  store_res[[s_id]] <- pmap

  # --- 4. Visualization ---
  
  plot_moll <- function(d, var, pts, tit, sub) {
    ggplot(d, aes(x = -l, y = b)) + 
      geom_tile(aes(fill = !!sym(var))) + 
      scale_fill_gradientn(colors = c("darkblue","blue","cyan","green","orange","red"), 
                           limits = c(0, 0.6)) +
      geom_point(data = pts$all, col="black", size=0.5, alpha=0.3) +
      geom_point(data = pts$sub, shape=4, col="red", size=2) +
      coord_map("mollweide") + 
      labs(title = tit, subtitle = sub, x=NULL, y=NULL) +
      theme_minimal()
  }

  pt_ov <- list(all = dat_all, sub = dat_sub)
  
  g1 <- plot_moll(pmap, "bg_n", pt_ov, paste("Background:", s_id), "KDE")
  g2 <- plot_moll(pmap, "sig_n", pt_ov, paste("Signal:", s_id), "KDE")
  g3 <- plot_moll(pmap, "prob", pt_ov, paste("Prob:", s_id), "GAM")
  g4 <- plot_moll(pmap, "precon_n", pt_ov, paste("Precon:", s_id), "Combined")

  ggsave(paste0("m1_bg_", s_id, ".pdf"), g1, w=9, h=6)
  ggsave(paste0("m2_sig_", s_id, ".pdf"), g2, w=9, h=6)
  ggsave(paste0("m3_pr_", s_id, ".pdf"), g3, w=9, h=6)
  ggsave(paste0("m4_pre_", s_id, ".pdf"), g4, w=9, h=6)

  # --- 5. Stats & Validation ---
  
  # Point extraction helper
  get_val <- function(grd, pts, v) {
    sapply(1:nrow(pts), function(i) {
      d2 <- (grd$l - pts$l[i])^2 + (grd$b - pts$b[i])^2
      grd[[v]][which.min(d2)]
    })
  }
  
  # Likelihood check
  v_pre <- get_val(pmap, dat_sub, "precon_n")
  v_kde <- get_val(pmap, dat_sub, "sig_n")
  
  ll_rat <- sum(log(v_pre + 1e-10)) - sum(log(v_kde + 1e-10))
  cat(sprintf("   Likelihood Ratio (Precon vs KDE): %.4f\n", ll_rat))

  # ROC Analysis
  n_pre <- get_val(pmap, bg_pts, "precon_n")
  n_kde <- get_val(pmap, bg_pts, "sig_n")
  
  auc <- function(pos, neg) {
    r <- rank(c(pos, neg))
    n_pos <- length(pos); n_neg <- length(neg)
    (sum(r[1:n_pos]) - n_pos*(n_pos+1)/2) / (n_pos*n_neg)
  }
  
  cat(sprintf("   AUC Precon: %.3f | AUC KDE: %.3f\n", 
              auc(v_pre, n_pre), auc(v_kde, n_kde)))
  
  # Manual ROC calc to keep dependencies low
  calc_roc_curve <- function(s, v_pos, v_neg) {
    dat <- data.frame(score = c(v_pos, v_neg), 
                      y = c(rep(1, length(v_pos)), rep(0, length(v_neg))))
    dat <- dat[order(dat$score, decreasing = TRUE), ]
    data.frame(
      fpr = cumsum(1 - dat$y) / sum(1 - dat$y),
      tpr = cumsum(dat$y) / sum(dat$y),
      Model = s
    )
  }

  roc_df <- rbind(
    calc_roc_curve("Precon", v_pre, n_pre),
    calc_roc_curve("KDE", v_kde, n_kde)
  )

  p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, color = Model)) +
    geom_line(size = 1.2) +
    geom_abline(lty = 2, col = "grey") +
    labs(title = paste("ROC:", s_id), x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal()
    
  ggsave(paste0("roc_", s_id, ".pdf"), p_roc, w=6, h=4)
}

# ---------------------------------------------------------
# Split vs Lumped Check
# ---------------------------------------------------------
if (all(c("merged", "grb_a", "grb_h") %in% names(store_res))) {
  cat("\nComparing Split vs Merged Models...\n")
  
  r_m <- store_res[["merged"]]
  r_a <- store_res[["grb_a"]]
  r_h <- store_res[["grb_h"]]
  
  # Reconstruct densities
  d_raw_m <- r_m$d_bg * r_m$prob
  d_raw_s <- r_m$d_bg * (r_a$prob + r_h$prob)
  
  # Extract at validation points (merged set)
  # (Using simple nearest neighbor on the grid)
  v_pts <- SCENARIOS[["merged"]] # filename
  pts_check <- read_grb(SCENARIOS[["merged"]])
  
  # Grid lookup wrapper
  lookup <- function(vals, g_l, g_b, q_l, q_b) {
    sapply(1:length(q_l), function(i) {
      idx <- which.min((g_l - q_l[i])^2 + (g_b - q_b[i])^2)
      vals[idx]
    })
  }
  
  val_m <- lookup(d_raw_m, r_m$l, r_m$b, pts_check$l, pts_check$b)
  val_s <- lookup(d_raw_s, r_m$l, r_m$b, pts_check$l, pts_check$b)
  
  dll <- sum(log(val_s + 1e-10)) - sum(log(val_m + 1e-10))
  cat(sprintf("   Delta Log-Likelihood (Split - Merged): %.4f\n", dll))
  
  # Plot combined split map
  r_m$split_sum <- norm_h(d_raw_s)
  g_split <- plot_moll(r_m, "split_sum", list(all=dat_all, sub=pts_check), 
                       "Combined Split Model", "Sum of GrbA + GrbH")
  ggsave("map_split_combined.pdf", g_split, w=9, h=6)
}

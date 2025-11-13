TV_size = .4
x = sim.tvVAR(burnin = 20, m = n)

compute_freq_diagnostics = function(x, short_M = 5, B = 8) {
  # x: n x p matrix (no burnin)
  n_rows = nrow(x)
  p = ncol(x)
  JJ0 = mvfft(x) / sqrt(n_rows)            # complex DFTs (length n_rows)
  # reorder as in your main code (backward first then 1:(n/2))
  backward = c(((n_rows / 2) + 1):n_rows)
  JJ = rbind(JJ0[backward, ], JJ0[c(1:(n_rows / 2)), ])
  K = nrow(JJ)  # number of freq bins considered (should be n)
  
  # Per-frequency multivariate vector JJ[k, ]  (complex length p)
  # Energy per k = trace( periodogram ) = sum_j |J_j(omega_k)|^2
  Energy = sapply(1:K, function(k) sum(Mod(JJ[k, ])^2))
  
  # Condition-proxy: small-window smoothed spectral matrix eigen-ratio
  cond_proxy = rep(NA, K)
  for (k in 1:K) {
    idx = ((k - short_M):(k + short_M))
    idx = ((idx - 1) %% K) + 1  # wrap-around
    S_hat = matrix(0+0i, nrow = p, ncol = p)
    for (j in idx) {
      v = matrix(JJ[j, ], ncol = 1)
      S_hat = S_hat + (v %*% Conj(t(v)))
    }
    S_hat = S_hat / length(idx)
    S_hat_reg = S_hat + diag(rep(1e-8, p))
    ev = Re(eigen(S_hat_reg, symmetric = FALSE, only.values = TRUE)$values)
    cond_proxy[k] = (max(ev) / pmax(min(ev), 1e-12))
  }
  
  # Time-variation score S_k: split into B blocks and compute blockwise energy,
  # then measure sum squared differences across blocks per-frequency
  block_len = floor(n_rows / B)
  E_block = matrix(0, nrow = K, ncol = B)
  for (b in 1:B) {
    start_i = (b - 1) * block_len + 1
    end_i = if (b < B) (b * block_len) else n_rows
    x_block = x[start_i:end_i, , drop = FALSE]
    JJb = mvfft(x_block) / sqrt(nrow(x_block))
    # reorder for block (approx)
    backward_b = c(((nrow(JJb) / 2) + 1):nrow(JJb))
    JJb2 = rbind(JJb[backward_b, ], JJb[c(1:(nrow(JJb) / 2)), ])
    # map block frequencies to full K by nearest indexing
    idx_map = round(seq(1, nrow(JJb2), length.out = K))
    E_block[, b] = sapply(1:K, function(k) sum(Mod(JJb2[idx_map[k], ])^2))
  }
  TimeVar = apply(E_block, 1, function(v) sum(diff(v)^2))
  
  omega = 2 * pi * (0:(K - 1)) / K
  return(tibble(k = 1:K, omega = omega, Energy = Energy, Cond = cond_proxy, TimeVar = TimeVar))
}

# Select informative ks using composite score
select_informative_ks = function(diagnostics_tbl, top_frac = 0.20,
                                 weight_energy = 1, weight_cond = 1, weight_timevar = 0.5,
                                 cond_penalty_scale = 1) {
  df = diagnostics_tbl %>%
    mutate(E_norm = Energy / max(Energy, na.rm = TRUE),
           Cond_norm = Cond / max(Cond, na.rm = TRUE),
           S_norm = TimeVar / max(TimeVar, na.rm = TRUE),
           Score = weight_energy * E_norm - weight_cond * (Cond_norm^cond_penalty_scale) + weight_timevar * S_norm)
  q = ceiling(nrow(df) * top_frac)
  selected = df %>% arrange(desc(Score)) %>% slice(1:q)
  return(list(selected_ks = selected$k, diagnostics = df, selected_tbl = selected))
}

# Plot diagnostics and highlight selected ks
plot_freq_diagnostics = function(diagnostics_df, selected_ks = NULL, out_prefix = NULL) {
  df = diagnostics_df
  p1 = ggplot(df, aes(x = k, y = Energy)) +
    geom_line() + geom_point(size = 0.7) +
    labs(title = "Frequency index energy", x = "k (frequency index)", y = "Energy") +
    theme_minimal()
  if (!is.null(selected_ks)) {
    p1 = p1 + geom_point(data = df %>% filter(k %in% selected_ks), aes(x=k,y=Energy), colour = "red", size = 1.5)
  }
  print(p1)
  
  p2 = ggplot(df, aes(x = k)) +
    geom_line(aes(y = Cond), colour = "grey40") +
    geom_line(aes(y = TimeVar), colour = "blue") +
    labs(title = "Condition proxy (grey) and Time-variation score (blue)", x = "k", y = "value") +
    theme_minimal()
  print(p2)
  
  if (!is.null(out_prefix)) {
    tryCatch({
      ggsave(sprintf("%s_energy.png", out_prefix), p1, width = 9, height = 3)
      ggsave(sprintf("%s_cond_timevar.png", out_prefix), p2, width = 9, height = 3)
    }, error = function(e) message("Warning: unable to save plots: ", e$message))
  }
}


pilot_n = 2048   # or whichever n you want to pilot on
pilot_M = 60     # pilot M (short smoothing for cond proxy)
#pilot_TV_size = 0.4

# simulate a single series (use same burnin as your pipeline)
x_pilot = sim.tvVAR(burnin = 20, m = pilot_n)

# compute diagnostics
diag_tbl = compute_freq_diagnostics(x_pilot, short_M = 5, B = 8)

# choose top 20% by default; tweak weights if you prefer energy-heavy or timevar-heavy selection
sel = select_informative_ks(diag_tbl, top_frac = 0.001, weight_energy = 1, weight_cond = 0.9, weight_timevar = 0.6)

# visualize
plot_freq_diagnostics(sel$diagnostics, selected_ks = sel$selected_ks, out_prefix = "pilot")

# Save chosen ks to reuse:
chosen_ks = sel$selected_ks



select_freq_band <- function(x, top_prop = 0.1, plot = TRUE) {
  n <- nrow(x)
  p <- ncol(x)
  
  # FFT and normalization
  JJ0 <- mvfft(x) / sqrt(n)
  
  # Backward reordering: π–2π band first
  backward <- c((n / 2 + 1):n)
  JJ <- rbind(JJ0[backward, ], JJ0[1:(n / 2), ])
  
  # Compute total spectral energy (sum over all variables)
  spec_energy <- apply(Mod(JJ)^2, 1, sum)
  
  # Restrict to π–2π band (first n/2 points)
  freq_idx <- 1:(n / 2)
  spec_energy <- spec_energy[freq_idx]
  
  # Sort by energy magnitude
  k_keep <- order(spec_energy, decreasing = TRUE)[1:ceiling(top_prop * length(spec_energy))]
  
  if (plot) {
    plot(freq_idx, rev(spec_energy), type = "l", col = "darkblue",
         xlab = "Frequency Index (Backward Ordered: 2π → π)",
         ylab = "Spectral Energy",
         main = "Spectral Energy over Backward Frequency Indices")
    points(freq_idx[length(freq_idx) + 1 - k_keep], spec_energy[k_keep], col = "red", pch = 19)
  }
  
  return(list(
    JJ = JJ,
    freq_indices = k_keep,
    energy = spec_energy
  ))
}

# Example usage:
set.seed(1)
x <- sim.tvVAR(burnin = 20, m = 256)
res <- select_freq_band(x, top_prop = 0.05)



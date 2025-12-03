smoothed_spectral_density_LOO = function(JJ, k, M, Kernel_func = Kernel_Triangular, leave_out = k) {
  n = nrow(JJ)
  p = ncol(JJ)
  S_k = matrix(0+0i, p, p)
  weight_sum = 0
  
  freq_window = max(1, k - M):min(n, k + M)
  
  for (j in freq_window) {
    if (j == leave_out) next  # leave-one-out step
    weight = Kernel_func((k - j) / M)
    if (weight > 0) {
      I_j = outer(JJ[j,], Conj(JJ[j,]))
      S_k = S_k + weight * I_j
      weight_sum = weight_sum + weight
    }
  }
  if (weight_sum > 0) S_k = S_k / weight_sum
  return(S_k)
}


whittle_loss = function(S_hat, I_k) {
  # Force real symmetric matrices
  S_hat = Re((S_hat + t(S_hat)) / 2)
  I_k   = Re((I_k   + t(I_k)) / 2)
  
  # Regularization for invertibility
  S_hat = S_hat + diag(1e-8, nrow(S_hat))
  
  loss = sum(diag(solve(S_hat) %*% I_k)) + log(det(S_hat))
  
  return( Re(loss) )    # <- MAKE LOSS PURELY REAL
}

adaptive_M_search = function(JJ, 
                             M_start = 2, 
                             M_step = 5, 
                             M_max = floor(0.25 * nrow(JJ)), 
                             Kernel_func = Kernel_Triangular) {
  
  n_freq = floor(nrow(JJ)/2)
  prev_M = NA
  prev_score = Inf
  results = data.frame(M = numeric(0), CV = numeric(0))
  
  M = M_start
  
  repeat {
    cat("Testing M =", M, "\n")
    total_loss = 0
    
    for (k in 1:n_freq) {
      S_hat = smoothed_spectral_density_LOO(JJ, k, M, Kernel_func, leave_out = k)
      I_k   = outer(JJ[k,], Conj(JJ[k,]))
      total_loss = total_loss + whittle_loss(S_hat, I_k)
    }
    
    # store
    results = rbind(results, data.frame(M = M, CV = total_loss))
    
    # stopping rule: CV has increased
    if (!is.na(prev_M) && total_loss > prev_score) {
      cat("\nCV increase detected: stopping.\n")
      break
    }
    
    # move forward
    prev_M = M
    prev_score = total_loss
    M = M + M_step
    
    if (M > M_max) {
      cat("\nReached M_max, stopping.\n")
      break
    }
  }
  
  # choose minimizer so far
  best_row = results[which.min(results$CV), , drop = FALSE]
  cat("\nOptimal M found:", best_row$M, "\n")
  
  return(list(best_M = best_row$M, results = results))
}


M_grid = unique(round(seq(2, .25 * n, length.out = 20)))  
n_freq = floor(nrow(JJ)/2)

CV_scores = numeric(length(M_grid))

for (m_index in seq_along(M_grid)) {
  M = M_grid[m_index]
  total_loss = 0
  
  for (k in 1:n_freq) {
    S_hat = smoothed_spectral_density_LOO(JJ, k, M, Kernel_Triangular, leave_out = k)
    I_k   = outer(JJ[k,], Conj(JJ[k,]))
    
    total_loss = total_loss + whittle_loss(S_hat, I_k)
  }
  
  CV_scores[m_index] = total_loss
  cat("M =", M, "→ CV =", total_loss, "\n")
}


best_M = M_grid[ which.min(CV_scores) ]
cat("\nOptimal bandwidth:", best_M, "\n")


out = adaptive_M_search(JJ, M_start = 5, M_step = 5)
best_M = out$best_M
print(best_M)





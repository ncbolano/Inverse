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
  # Numerical stabilization
  S_hat = Re(S_hat)
  
  # Regularization for invertibility if needed
  S_hat = S_hat + diag(1e-8, nrow(S_hat))
  
  return( sum(diag(solve(S_hat) %*% I_k)) + log(det(S_hat)) )
}


M_grid = unique(round(seq(5, 100, length.out = 15)))  
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


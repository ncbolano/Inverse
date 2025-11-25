estimate_candidate_ks = function(x, M, num_candidates = 5) {
  n = nrow(x)
  p = ncol(x)
  
  # FFT
  JJ = mvfft(x) / sqrt(n)
  
  # Smooth periodogram with your kernel
  freq_grid = 1:floor(n/2)
  scores = numeric(length(freq_grid))
  
  for (k in freq_grid) {
    # Estimate S(k) using your kernel smoother
    S_k = estimate_spectral_density(JJ, k, M, Kernel_function)
    
    # Score this frequency
    trace_S = sum(diag(S_k))
    eigenvals = eigen(S_k, only.values = TRUE)$values
    cond_num = max(eigenvals) / min(eigenvals)
    
    scores[k] = trace_S / cond_num  # or another criterion
  }
  
  # Select top peaks that are spaced apart
  candidate_ks = select_spaced_peaks(scores, freq_grid, min_spacing = M, 
                                     num_candidates)
  
  return(candidate_ks)
}

smoothed_spectral_density = function(JJ, k, M, Kernel_func = Kernel_Triangular) {
  n = nrow(JJ)
  p = ncol(JJ)
  S_k = matrix(0 + 0i, p, p)  # Complex matrix
  weight_sum = 0
  
  # Determine frequency window
  freq_window = max(1, k - M):min(n, k + M)
  
  for (j in freq_window) {
    weight = Kernel_func((k - j) / M)
    if (weight > 0) {
      I_j = outer(JJ[j, ], Conj(JJ[j, ]))
      S_k = S_k + weight * I_j
      weight_sum = weight_sum + weight
    }
  }
  
  # Normalize by sum of weights
  if (weight_sum > 0) {
    S_k = S_k / weight_sum
  }
  
  return(S_k)
}

return_max = function(x, M) {
  max1 = which.max(x)
  
  # Create a copy and set values within 2*M of max1 to -Inf
  x2 = x
  min_dist = 2 * M
  
  # Set forbidden region to -Inf
  forbidden_indices = max(1, max1 - min_dist):min(length(x), max1 + min_dist)
  x2[forbidden_indices] = -Inf
  
  # Find second maximum (at least 2*M away)
  max2 = which.max(x2)
  
  # Validity check of existing other point to prevent breakdown for M > n/4
  if (is.infinite(x2[max2])) {
    warning("No valid second maximum found at distance >= 2*M")
    return(c(max1, NA))
  }
  
  return(c(max1, max2))
}


library(ggplot2)
# Trace example
R <- 500 # number of replications
nu <- 2 # dimension of matrix
p <- 3 # dimension of time series
n <- 2^12
M<-30
TV_size=0.2

set.seed(123)
n = 2048
x = sim.tvVAR(burnin = 500, m = n)

# Compute FFT
JJ = mvfft(x) / sqrt(nrow(x))

# Raw periodogram loading
n_freq = floor(nrow(JJ) / 2)
trace_periodogram = numeric(n_freq)
diag_periodograms = matrix(0, nrow = n_freq, ncol = 3)
largest_eigenvalue = numeric(n_freq)
condition_number = numeric(n_freq)

# Smoothed periodogram loading
trace_smooth = numeric(n_freq)
diag_smooth = matrix(0, nrow = n_freq, ncol = 3)
largest_eig_smooth = numeric(n_freq)
condition_number_smooth = numeric(n_freq)

for (k in 1:n_freq) {
  # Periodogram matrix at frequency k
  I_k = outer(JJ[k, ], Conj(JJ[k, ]))
  
  # Trace (sum of diagonal = total power)
  trace_periodogram[k] = sum(diag(I_k))
  
  # Individual periodograms (diagonal elements)
  diag_periodograms[k, ] = diag(I_k)
  
  # Eigenvalues for other metrics
  eigenvals = eigen(I_k, only.values = TRUE)$values
  eigenvals = Re(eigenvals) 
  
  largest_eigenvalue[k] = max(eigenvals)
  condition_number[k] = max(eigenvals) / (min(eigenvals) + 1e-10)
  
  # Smoothed Periodogram (weighted by kernel)
  S_k = smoothed_spectral_density(JJ, k, M, Kernel_Triangular)
  # Smoothed Versions
  trace_smooth[k] = sum(diag(S_k))
  diag_smooth[k, ] = diag(S_k)
  
  eigenvals_smooth = eigen(S_k, only.values = TRUE)$values
  eigenvals_smooth = Re(eigenvals_smooth)
  largest_eig_smooth[k] = max(eigenvals_smooth)
  condition_number_smooth[k] = max(eigenvals_smooth) / (min(eigenvals_smooth) + 1e-10)
}

frequencies = (1:n_freq) / n  

plot_data = data.frame(
  frequency = frequencies,
  trace = Re(trace_periodogram),
  var1 = Re(diag_periodograms[, 1]),
  var2 = Re(diag_periodograms[, 2]),
  var3 = Re(diag_periodograms[, 3]),
  largest_eig = largest_eigenvalue,
  cond_num = condition_number
)


p1 <- ggplot(plot_data, aes(x = frequency, y = trace)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  labs(
    title = "Trace of Spectral Density Matrix",
    subtitle = "Total power across all 3 variables",
    x = "Frequency (cycles per observation)",
    y = "Trace(I(ω))"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(p1)

composite_value = trace_periodogram / condition_number
which.max(composite_value)

composite_smooth = trace_smooth / condition_number_smooth
chosen_freq = return_max(composite_smooth ,M)

chosen_ks = chosen_freq / n




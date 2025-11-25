smoothed_spectral_density = function(JJ, k, m, Kernel_func = Kernel_Triangular) {
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

return_max = function(x, m) {
  max1 = which.max(x)
  
  # Create a copy and set values within 2*M of max1 to -Inf
  x2 = x
  min_dist =  M
  
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

sim.tvVAR = function(burnin, m, TV_size) {
  A = matrix(c(0.5, 0.2, 0, 0, 0.8, 0, 0, 0.3, 0.6), ncol = 3, byrow = T)
  n = m + burnin
  p = 3
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)
  for (tt in (2:n)) {
    A.t = A
    A.t[1, 1] = st[tt]  # Fixed indexing
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

library(ggplot2)
# Trace example
R <- 500 # number of replications
nu <- 2 # dimension of matrix
p <- 3 # dimension of time series
n <- 2^12
M<-30
TV_size=0.6
m = 8

set.seed(123)
n = 2048
x = sim.tvVAR(burnin = 500, m = n, TV_size)

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
  S_k = smoothed_spectral_density(JJ, k, m, Kernel_Triangular)
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


p1 = ggplot(plot_data, aes(x = frequency, y = trace)) +
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


plot_data_smooth = data.frame(
  frequency = frequencies,
  trace = Re(trace_smooth),
  var1 = Re(diag_smooth[, 1]),
  var2 = Re(diag_smooth[, 2]),
  var3 = Re(diag_smooth[, 3]),
  largest_eig = largest_eig_smooth,
  cond_num = condition_number_smooth
)

p1_smooth = ggplot(plot_data_smooth, aes(x = frequency, y = trace)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  labs(
    title = "Trace of Smoothed Spectral Density Matrix",
    subtitle = "Total power across all 3 variables (kernel-smoothed)",
    x = "Frequency (cycles per observation)",
    y = "Trace(S(ω))"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(p1_smooth)


composite_value = trace_periodogram / condition_number
which.max(composite_value)

composite_smooth = trace_smooth / condition_number_smooth
chosen_freq = return_max(composite_smooth , M)
sort(composite_smooth, decreasing = TRUE )
k_freq = chosen_freq / n


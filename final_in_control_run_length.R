# Load required libraries
library(DRIP)
library(imager)
library(magick)
library(progress)

# Reproducibility
set.seed(0)

# Number of Monte Carlo simulations
total = 1000

# Window size
w = 5

# True noise standard deviation
sig = 0.30

# Precomputed threshold
h = 15.55166

# Image dimensions
n1 = 128
n2 = 128
n = n1 * n2

# Function: Add Gaussian noise to image

add_noise = function(img) {
  m1 = min(img)
  m2 = max(img)
  
  Z = ((2 * img) / (m2 - m1)) - ((m2 + m1) / (m2 - m1))
  Z = Z + matrix(rnorm(n1 * n2, 0, sig), nrow = n1, ncol = n2)
  
  return(Z)
}

# Function: Detect edges using Canny detector

detect_edges = function(Z) {
  D = array(Z, dim = c(n1, n2, 1))
  D = image_read(D)
  D = image_canny(D)
  D = image_data(D)
  D = matrix(as.numeric(D), nrow = n1, ncol = n2)
  return(D)
}

# Load dataset

# Load grayscale image and convert to matrix
img = as.matrix(grayscale(load.image("star_line_appear_128/test1.png")))

# Initial reference sample

# Generate w noisy observations
I = array(dim = c(w, n1, n2))
for (i in 1:w) {
  I[i, , ] = add_noise(img)
}

# Convert each to edge maps
for (i in 1:w) {
  I[i, , ] = detect_edges(I[i, , ])
}

# Monte Carlo estimation of run length

# Progress bar
pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = total,
                      clear = FALSE)

run_length = c()
for (r in 1:total) {
  # Generate initial monitoring window
  
  W = array(dim = c(w, n1, n2))
  for (i in 1:w) {
    W[i, , ] = add_noise(img)
  }
  
  # Convert to edge maps
  for (i in 1:w) {
    W[i, , ] = detect_edges(W[i, , ])
  }
  
  # Sequential monitoring loop
  
  # Continue while statistic < threshold
  flag = TRUE
  
  # Store monitoring statistics
  T = c()
  
  while (flag) {
    # Slide window forward
    
    # Remove oldest observation and append new one
    W[1:(w - 1), , ] = W[2:w, , ]
    W[w, , ] = detect_edges(add_noise(img))
    
    # Combine all detected edges
    temp = matrix(rep(FALSE, n1 * n2), nrow = n1, ncol = n2)
    for (i in 1:w) {
      temp = temp | I[i, , ] | W[i, , ]
    }
    
    # Expand edge support via 3x3 convolution
    M = as.matrix(convolve(as.cimg(temp), as.cimg(matrix(1, 3, 3))))
    M = (M > 0)
    
    # Construct feature matrix using only active pixels
    X = matrix(nrow = 2 * w, ncol = sum(M))
    for (i in 1:w) {
      X[i, ] = as.numeric(I[i, , ][M])
      X[i + w, ] = as.numeric(W[i, , ][M])
    }
    
    # Compute pairwise distance matrix
    L = as.matrix(dist(X))
    
    # Between-group distances
    r12 = sum(L[1:w, (w + 1):(2 * w)])
    
    # Within-group distances
    r11 = sum(L[1:w, 1:w])
    r22 = sum(L[(w + 1):(2 * w), (w + 1):(2 * w)])
    
    # Test statistic
    T = c(T, (2 * r12 - r11 - r22) / (2 * w))
    
    # Continue monitoring until statistic exceeds threshold
    flag = ((2 * r12 - r11 - r22) / (2 * w) < h)
  }
  
  # Store stopping time / run length
  run_length = c(run_length, length(T))
  
  pb$tick()
}

cat(run_length, "\n")
cat(mean(run_length), "\n")

# Save results to CSV
write.table(run_length,
            file = "star_line_grow_128/run_length_30_5.csv",
            row.names = FALSE,
            col.names = FALSE)
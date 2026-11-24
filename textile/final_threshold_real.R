# Load required libraries
library(DRIP)
library(imager)
library(magick)
library(progress)

# Reproducibility
set.seed(0)

# Number of Monte Carlo simulations
total = 100

# Window size
w = 5

# Hyperparameter
lambda = sqrt(0.03685)

# Significance level for threshold
alpha = 0.95

# True noise standard deviation
sig = 0.05

# Image dimensions
n1 = 500
n2 = 500
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

# Load and preprocess base image

# Load grayscale image and convert to matrix
img = as.matrix(grayscale(load.image("test1.png")))

# Generate w noisy versions of the image
I = array(dim = c(w, n1, n2))
for (i in 1:w) {
  I[i, , ] = add_noise(img)
}

# Estimate noise level via median smoothing

# Create cleaned version via repeated median filtering
clean_I = I
for (i in 1:15) {
  for (j in 1:w) {
    clean_I[j, , ] = as.matrix(medianblur(as.cimg(clean_I[j, , ]), n = 3))
  }
}

# Estimate noise by subtracting smoothed images
noise_I = array(dim = c(w, n1, n2))
for (i in 1:w) {
  noise_I[i, , ] = I[i, , ] - clean_I[i, , ]
}

# Empirical noise standard deviation estimate
sig_hat = sqrt(mean(noise_I^2) - mean(noise_I)^2)

# Compute edge maps for cleaned images
temp_I = array(dim = c(w, n1, n2))
for (i in 1:w) {
  temp_I[i, , ] = detect_edges(clean_I[i, , ])
}

# Monte Carlo threshold estimation

# Progress bar
pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = total,
                      clear = FALSE)

thresh_T = c()
for (m in 1:total) {
  # Bootstrap sample from cleaned images
  B = clean_I[sample(1:w, replace = TRUE), , ]
  
  temp_B = array(dim = c(w, n1, n2))
  for (i in 1:w) {
    # Add inflated/deflated Gaussian noise
    B[i, , ] = B[i, , ] + matrix(rnorm(n1 * n2, 0, lambda * sig_hat),
                                 nrow = n1,
                                 ncol = n2)
    
    # Detect edges
    temp_B[i, , ] = detect_edges(B[i, , ])
  }
  
  # Combine all detected edges
  temp = matrix(rep(FALSE, n1 * n2), nrow = n1, ncol = n2)
  for (i in 1:w) {
    temp = temp | temp_I[i, , ] | temp_B[i, , ]
  }
  
  # Expand edge support via 3x3 convolution
  M = as.matrix(convolve(as.cimg(temp), as.cimg(matrix(1, 3, 3))))
  M = (M > 0)
  
  # Construct feature matrix using only active pixels
  X = matrix(nrow = 2 * w, ncol = sum(M))
  for (i in 1:w) {
    X[i, ] = as.numeric(temp_I[i, , ][M])
    X[i + w, ] = as.numeric(temp_B[i, , ][M])
  }
  
  # Compute pairwise distance matrix
  L = as.matrix(dist(X))
  
  # Between-group distances
  r12 = sum(L[1:w, (w + 1):(2 * w)])
  
  # Within-group distances
  r11 = sum(L[1:w, 1:w])
  r22 = sum(L[(w + 1):(2 * w), (w + 1):(2 * w)])
  
  # Test statistic
  thresh_T = c(thresh_T, (2 * r12 - r11 - r22) / (2 * w))
  
  pb$tick()
}

# Final threshold
h = quantile(thresh_T, alpha)
cat("h = ", h, "\n")
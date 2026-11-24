# Load required libraries
library(DRIP)
library(imager)
library(magick)
library(progress)

# Reproducibility
set.seed(0)

# Number of frames from each regime

# Baseline (img1)
f1 = 14
# First change (img2)
f2 = 7
# Second change (img3)
f3 = 7
# Third change (img4)
f4 = 7
# Fourth change (img5)
f5 = 14

# Total number of frames in one simulated sequence
f_star = f1 + f2 + f3 + f4 + f5

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

img1 = as.matrix(grayscale(load.image("star_line_appear_128/test1.png")))
img2 = as.matrix(grayscale(load.image("star_line_appear_128/test2.png")))
img3 = as.matrix(grayscale(load.image("star_line_appear_128/test3.png")))
img4 = as.matrix(grayscale(load.image("star_line_appear_128/test4.png")))
img5 = as.matrix(grayscale(load.image("star_line_appear_128/test5.png")))

# Initial reference sample

# Generate w noisy observations
I = array(dim = c(w, n1, n2))
for (i in 1:w) {
  I[i, , ] = add_noise(img1)
}

# Convert each to edge maps
for (i in 1:w) {
  I[i, , ] = detect_edges(I[i, , ])
}

# Monte Carlo delay estimation

# Progress bar
pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = total,
                      clear = FALSE)

# Early detections
neg_delay = c()
# Valid detection delays
delay = c()

# Repeat until we collect 'total' valid delays
while (length(delay) < total) {
  # Generate full sequence with regime changes
  W = array(dim = c(f_star, n1, n2))
  for (i in 1:f1) {
    W[i, , ] = add_noise(img1)
  }
  for (i in 1:f2) {
    W[i + f1, , ] = add_noise(img2)
  }
  for (i in 1:f3) {
    W[i + f1 + f2, , ] = add_noise(img3)
  }
  for (i in 1:f4) {
    W[i + f1 + f2 + f3, , ] = add_noise(img4)
  }
  for (i in 1:f5) {
    W[i + f1 + f2 + f3 + f4, , ] = add_noise(img5)
  }
  
  # Convert to edge maps
  for (i in 1:f_star) {
    W[i, , ] = detect_edges(W[i, , ])
  }
  
  T = c()
  
  # Slide window of size w across sequence
  for (m in 1:(f_star - w + 1)) {
    # Combine all detected edges
    temp = matrix(rep(FALSE, n1 * n2), nrow = n1, ncol = n2)
    for (i in 1:w) {
      temp = temp | I[i, , ] | W[i + m - 1, , ]
    }
    
    # Expand edge support via 3x3 convolution
    M = as.matrix(convolve(as.cimg(temp), as.cimg(matrix(1, 3, 3))))
    M = (M > 0)
    
    # Construct feature matrix using only active pixels
    X = matrix(nrow = 2 * w, ncol = sum(M))
    for (i in 1:w) {
      X[i, ] = as.numeric(I[i, , ][M])
      X[i + w, ] = as.numeric(W[i + m - 1, , ][M])
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
  }
  
  # Compute detection delay
  if (which(T > h)[1] >= (f1 + 2 - w)) {
    delay = c(delay, which(T > h)[1] - (f1 + 2 - w))
    pb$tick()
  }
  else {
    neg_delay = c(neg_delay, which(T > h)[1] - (f1 + 2 - w))
  }
}

cat(delay, "\n")
cat(mean(delay), "\n")

cat(neg_delay, "\n")
cat(length(neg_delay), "\n")

# Save results
write.table(delay,
            file = "star_line_appear_128/delay_30_5.csv",
            row.names = FALSE,
            col.names = FALSE)
write.table(neg_delay,
            file = "star_line_appear_128/neg_delay_30_5.csv",
            row.names = FALSE,
            col.names = FALSE)
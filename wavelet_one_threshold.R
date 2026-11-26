library(magick)
library(waveslim)
library(imager)
library(progress)

set.seed(0)

total = 1000

w = 5
lambda = sqrt(1.054)
alpha = 0.95

sig = 0.3

n1 = 128
n2 = 128
n = n1 * n2

# Functions

add_noise = function(img) {
  m1 = min(img)
  m2 = max(img)
  
  Z = ((2 * img) / (m2 - m1)) - ((m2 + m1) / (m2 - m1))
  Z = Z + matrix(rnorm(n1 * n2, 0, sig), nrow = n1, ncol = n2)
  
  return(Z)
}

wavelet_one = function(Z) {
  D = dwt(as.numeric(Z), wf = "haar", n.levels = 4)
  
  return(as.numeric(D$s4))
}

# Data

img = as.matrix(grayscale(load.image("blob_pacman_128/test1.png")))

I = array(dim = c(w, n1, n2))

for (i in 1:w) {
  I[i, , ] = add_noise(img)
}

# Threshold

clean_I = I

for (i in 1:15) {
  for (j in 1:w) {
    clean_I[j, , ] = as.matrix(medianblur(as.cimg(clean_I[j, , ]), n = 3))
  }
}

# Estimate noise

noise_I = array(dim = c(w, n1, n2))

for (i in 1:w) {
  noise_I[i, , ] = I[i, , ] - clean_I[i, , ]
}

sig_hat = sqrt(mean(noise_I^2) - mean(noise_I)^2)

# Initial Wavelet

temp_I = matrix(nrow = w,
                ncol = length(wavelet_one(clean_I[1, , ])))

for (i in 1:w) {
  temp_I[i, ] = wavelet_one(clean_I[i, , ])
}

pb = progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = total,
  clear = FALSE
)

thresh_T = c()

for (m in 1:total) {
  
  B = clean_I[sample(1:w, replace = TRUE), , , drop = FALSE]
  
  temp_B = matrix(nrow = w, ncol = ncol(temp_I))
  
  for (i in 1:w) {
    B[i, , ] = B[i, , ] +
      matrix(rnorm(n1 * n2, 0, lambda * sig_hat),
             nrow = n1, ncol = n2)
    
    temp_B[i, ] = wavelet_one(B[i, , ])
  }
  
  # Energy statistic
  
  X = rbind(temp_I, temp_B)
  
  L = as.matrix(dist(X))
  
  r12 = sum(L[1:w, (w + 1):(2 * w)])
  
  r11 = sum(L[1:w, 1:w])
  
  r22 = sum(L[(w + 1):(2 * w), (w + 1):(2 * w)])
  
  T = (2 * r12 - r11 - r22) / (2 * w)
  
  thresh_T = c(thresh_T, T)
  
  pb$tick()
}

# Threshold

h = quantile(thresh_T, alpha)

cat("h = ", h, "\n")

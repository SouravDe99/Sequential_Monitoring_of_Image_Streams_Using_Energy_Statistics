library(DRIP)
library(imager)
library(magick)
library(waveslim)
library(progress)

set.seed(0)

total = 1000

w = 5
sig = 0.3
h = 14.34461         

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

wavelet_two = function(Z) {
  D = dwt.2d(Z, wf = "haar", J = 2)
  return(as.numeric(D$LL2))
}

# Data

img = as.matrix(grayscale(load.image("blob_pacman_128/test1.png")))

I_image = array(dim = c(w, n1, n2))

for (i in 1:w) {
  I_image[i, , ] = add_noise(img)
}

I = matrix(nrow = w,
           ncol = length(wavelet_two(I_image[1, , ])))

for (i in 1:w) {
  I[i, ] = wavelet_two(I_image[i, , ])
}

pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = total,
                      clear = FALSE)

run_length = c()

for (r in 1:total) {
  
  # Final Data
  
  W_image = array(dim = c(w, n1, n2))
  
  for (i in 1:w) {
    W_image[i, , ] = add_noise(img)
  }
  
  W = matrix(nrow = w, ncol = ncol(I))
  
  for (i in 1:w) {
    W[i, ] = wavelet_two(W_image[i, , ])
  }
  
  flag = TRUE
  T = c()
  
  while (flag) {
    
    W[1:(w - 1), ] = W[2:w, ]
    
    new_image = add_noise(img)
    
    W[w, ] = wavelet_two(new_image)
    
    X = rbind(I, W)
    
    L = as.matrix(dist(X))
    
    r12 = sum(L[1:w, (w + 1):(2 * w)])
    
    r11 = sum(L[1:w, 1:w])
    
    r22 = sum(L[(w + 1):(2 * w),
                (w + 1):(2 * w)])
    
    T_current = (2 * r12 - r11 - r22) / (2 * w)
    
    T = c(T, T_current)
    
    flag = (T_current < h)
  }
  
  run_length = c(run_length, length(T))
  
  pb$tick()
}

cat(run_length, "\n")
cat(mean(run_length), "\n")

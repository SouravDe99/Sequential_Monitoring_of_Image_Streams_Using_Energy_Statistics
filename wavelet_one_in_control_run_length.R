library(waveslim)
library(imager)
library(magick)
library(progress)

set.seed(0)

total = 100

w = 5
sig = 0.1
h = 6.785941

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

img = as.matrix( grayscale(load.image("star_line_grow_128/test1.png")))
in_control_images = array(dim = c(w, n1, n2))
for (i in 1:w) {
  in_control_images[i, , ] = add_noise(img)
}
in_control_wavelets = matrix(nrow = w, ncol = length(wavelet_one(in_control_images[1, , ])))
for (i in 1:w) {
  in_control_wavelets[i, ] = wavelet_one(in_control_images[i, , ])
}

pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = total,
                      clear = FALSE)
run_length = c()
for (r in 1:total) {
  # Final Data: generate the initial monitoring window and transform it.
  monitoring_images = array(dim = c(w, n1, n2))
  for (i in 1:w) {
    monitoring_images[i, , ] = add_noise(img)
  }

  monitoring_wavelets = matrix(
    nrow = w,
    ncol = length(wavelet_one(monitoring_images[1, , ]))
  )
  for (i in 1:w) {
    monitoring_wavelets[i, ] = wavelet_one(monitoring_images[i, , ])
  }

  flag = TRUE
  T = c()
  while (flag) {
    monitoring_wavelets[1:(w - 1), ] = monitoring_wavelets[2:w, ]
    new_image = add_noise(img)
    monitoring_wavelets[w, ] = wavelet_one(new_image)

    X = rbind(in_control_wavelets, monitoring_wavelets)
    L = as.matrix(dist(X))
    r12 = sum(L[1:w, (w + 1):(2 * w)])
    r11 = sum(L[1:w, 1:w])
    r22 = sum(L[(w + 1):(2 * w), (w + 1):(2 * w)])
    T_current = (2 * r12 - r11 - r22) / (2 * w)
    T = c(T, T_current)
    flag = (T_current < h)
  }
    run_length = c(run_length,length(T))
    pb$tick()
  }
cat(run_length,"\n")
cat(mean(run_length),"\n")

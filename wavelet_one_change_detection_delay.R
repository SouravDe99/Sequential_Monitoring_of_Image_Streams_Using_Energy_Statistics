library(waveslim)
library(imager)
library(magick)
library(progress)

set.seed(0)

f1 = 14
f2 = 7
f3 = 7
f4 = 7
f5 = 14
f_star = f1 + f2 + f3 + f4 + f5
total = 1000

w = 5
sig = 0.3
h = 14.14299   

n1 = 128
n2 = 128
n = n1 * n2

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

img1 = as.matrix(grayscale(load.image("blob_pacman_128/test1.png")))
img2 = as.matrix(grayscale(load.image("blob_pacman_128/test2.png")))
img3 = as.matrix(grayscale(load.image("blob_pacman_128/test3.png")))
img4 = as.matrix(grayscale(load.image("blob_pacman_128/test4.png")))
img5 = as.matrix(grayscale(load.image("blob_pacman_128/test5.png")))

in_control_images = array(dim = c(w, n1, n2))
for (i in 1:w) {
  in_control_images[i, , ] = add_noise(img1)
}

in_control_wavelets = matrix(nrow = w, ncol = length(wavelet_one(in_control_images[1, , ])))
for (i in 1:w) {
  in_control_wavelets[i, ] = wavelet_one(in_control_images[i, , ])
}

pb = progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = total,
  clear = FALSE
)

neg_delay = c()
delay = c()

while (length(delay) < total) {
  
  monitoring_images = array(dim = c(f_star, n1, n2))
  
  for (i in 1:f1) {
    monitoring_images[i, , ] = add_noise(img1)
  }
  for (i in 1:f2) {
    monitoring_images[i + f1, , ] = add_noise(img2)
  }
  for (i in 1:f3) {
    monitoring_images[i + f1 + f2, , ] = add_noise(img3)
  }
  for (i in 1:f4) {
    monitoring_images[i + f1 + f2 + f3, , ] = add_noise(img4)
  }
  for (i in 1:f5) {
    monitoring_images[i + f1 + f2 + f3 + f4, , ] = add_noise(img5)
  }
  
  monitoring_wavelets = matrix(nrow = f_star, ncol = ncol(in_control_wavelets))
  for (i in 1:f_star) {
    monitoring_wavelets[i, ] = wavelet_one(monitoring_images[i, , ])
  }
  
  T = c()
  
  for (m in 1:(f_star - w + 1)) {
    
    X = rbind(in_control_wavelets, monitoring_wavelets[m:(m + w - 1), ])
    
    L = as.matrix(dist(X))
    
    r12 = sum(L[1:w, (w + 1):(2 * w)])
    r11 = sum(L[1:w, 1:w])
    r22 = sum(L[(w + 1):(2 * w), (w + 1):(2 * w)])
    
    T = c(T, (2 * r12 - r11 - r22) / (2 * w))
  }
  
  detection = which(T > h)
  
  if (length(detection) > 0) {
    
    if (detection[1] >= (f1 + 2 - w)) {
      delay = c(delay, detection[1] - (f1 + 2 - w))
      pb$tick()
    }
    else {
      neg_delay = c(neg_delay, detection[1] - (f1 + 2 - w))
    }
    
  }
}

cat(delay, "\n")
cat(mean(delay), "\n")

cat(neg_delay, "\n")
cat(length(neg_delay), "\n")

# write.table(delay, file = "star_line_grow/delay_wavelet.csv",
#             row.names = FALSE, col.names = FALSE)
# write.table(neg_delay, file = "star_line_grow/neg_delay_wavelet.csv",
#             row.names = FALSE, col.names = FALSE)
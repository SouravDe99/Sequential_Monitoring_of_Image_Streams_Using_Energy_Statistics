library(DRIP)
library(imager)
library(magick)
library(progress)

set.seed(0)

i_star = 15
f1 = 14
f2 = 7
f3 = 7
f4 = 7
f5 = 14
f_star = f1 + f2 + f3 + f4 + f5

w = 5

sig = 0.05

n1 = 500
n2 = 500
n = n1 * n2

# Functions

read_image = function(path) {
  img = image_read(path)
  img = image_convert(img, colorspace = "gray")
  img = image_data(img)
  img = matrix(as.numeric(img), nrow = n1, ncol = n2)
  return(img)
}

add_noise = function(img) {
  m1 = min(img)
  m2 = max(img)
  
  Z = ((2 * img) / (m2 - m1)) - ((m2 + m1) / (m2 - m1))
  Z = Z + matrix(rnorm(n1 * n2, 0, sig), nrow = n1, ncol = n2)
  
  return(Z)
}

detect_edges = function(Z) {
  D = array(Z, dim = c(n1, n2, 1))
  D = image_read(D)
  D = image_canny(D)
  D = image_data(D)
  D = matrix(as.numeric(D), nrow = n1, ncol = n2)
  return(D)
}

# Data

img1 = as.matrix(grayscale(load.image("test1.png")))
img2 = as.matrix(grayscale(load.image("test2.png")))
img3 = as.matrix(grayscale(load.image("test3.png")))
img4 = as.matrix(grayscale(load.image("test4.png")))
img5 = as.matrix(grayscale(load.image("test5.png")))

I = array(dim = c(i_star, n1, n2))
for (i in 1:i_star) {
  I[i, , ] = add_noise(img1)
}

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

# Edge

for (i in 1:i_star) {
  I[i, , ] = detect_edges(I[i, , ])
}

for (i in 1:f_star) {
  W[i, , ] = detect_edges(W[i, , ])
}

# Energy

pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = f_star - w + 1,
                      clear = FALSE)

T = c()
for (m in 1:(f_star - w + 1)) {
  temp = matrix(rep(FALSE, n1 * n2), nrow = n1, ncol = n2)
  for (i in 1:w) {
    temp = temp | I[i, , ] | W[i + m - 1, , ]
  }
  
  M = as.matrix(convolve(as.cimg(temp), as.cimg(matrix(1, 3, 3))))
  M = (M > 0)
  
  X = matrix(nrow = 2 * w, ncol = sum(M))
  for (i in 1:w) {
    X[i, ] = as.numeric(I[i, , ][M])
    X[i + w, ] = as.numeric(W[i + m - 1, , ][M])
  }
  
  L = as.matrix(dist(X))
  
  r12 = sum(L[1:w, (w + 1):(2 * w)])
  r11 = sum(L[1:w, 1:w])
  r22 = sum(L[(w + 1):(2 * w), (w + 1):(2 * w)])
  
  T = c(T, (2 * r12 - r11 - r22) / (2 * w))
  pb$tick()
}

print(T)
plot(T, type = "o")
abline(h = 127.3933, col = "red")
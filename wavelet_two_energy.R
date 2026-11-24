library(magick)
library(waveslim)
library(progress)

i_star = 15
f1 = 14
f2 = 7
f3 = 7
f4 = 7
f5 = 14
f_star = f1 + f2 + f3 + f4 + f5

w = 5

sig = 0.3

n1 = 128
n2 = 128
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

# Data

img1 = read_image("star_line_appear_128/test1.png")
img2 = read_image("star_line_appear_128/test2.png")
img3 = read_image("star_line_appear_128/test3.png")
img4 = read_image("star_line_appear_128/test4.png")
img5 = read_image("star_line_appear_128/test5.png")

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

# Wave

I_W = c()
for (i in 1:i_star) {
  I_W = rbind(I_W, as.numeric(dwt.2d(I[i, , ], wf = "haar", J = 2)$LL2))
}

W_W = c()
for (i in 1:f_star) {
  W_W = rbind(W_W, as.numeric(dwt.2d(W[i, , ], wf = "haar", J = 2)$LL2))
}

# Energy

pb = progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                      total = f_star - w + 1,
                      clear = FALSE)

T = c()
for (m in 1:(f_star - w + 1)) {
  X = c()
  for (i in 1:w) {
    X = rbind(X, I_W[i, ])
  }
  for (i in m:(m + w - 1)) {
    X = rbind(X, W_W[i, ])
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
library(magick)
library(imager)
library(spc4sts)
library(progress)

set.seed(0)

# Parameters

w = 11
sig = 0.3

n1 = 128
n2 = 128
n = n1 * n2

total = 100
fa.rate = 0.01

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
  
  Z = Z + matrix(rnorm(n1 * n2, 0, sig),
                 nrow = n1, ncol = n2)
  
  return(Z)
}

# Data

img1 = read_image("star_line_grow_128/test1.png")

# In-control surface model

model = surfacemodel(
  img = img1,
  nb = 3,
  trim.vars = TRUE,
  cp = 1e-5,
  xval = 5,
  standardize = TRUE,
  subsample = 1,
  verbose = FALSE,
  keep.residuals = TRUE
)

# Phase I Data

I = array(dim = c(n1, n2, total))

pb = progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = total,
  clear = FALSE
)

for (i in 1:total) {
  I[, , i] = add_noise(img1)
  pb$tick()
}

# Threshold estimation

cl = climit(
  imgs = I,
  fa.rate = fa.rate,
  model = model,
  type = 1,
  stat = "ad",
  w = w,
  verbose = TRUE
)

# Threshold

h = as.numeric(cl$localStat$control.limit)

cat("Estimated threshold h =", h, "\n")

# Plot

plot(cl)

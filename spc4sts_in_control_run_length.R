
library(magick)
library(progress)
library(spc4sts)

set.seed(0)

# Parameters

total = 10

w = 11
sig = 0.3

n1 = 128
n2 = 128
n = n1 * n2

# Threshold obtained from the threshold estimation program

h = 452.772

cat("Threshold h =", h, "\n")


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

img = read_image("star_line_grow_128/test1.png")


# In-control surface model

model = surfacemodel(
  img = img,
  nb = 3,
  trim.vars = TRUE,
  cp = 1e-5,
  xval = 5,
  standardize = TRUE,
  subsample = 1,
  verbose = FALSE,
  keep.residuals = TRUE
)


# In-control run length

pb = progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = total,
  clear = FALSE
)

run_length = c()

for (r in 1:total) {
  
  flag = TRUE
  t = 0
  
  while (flag) {
    
    # Generate new in-control image
    
    Z = add_noise(img)
    
    # Anderson-Darling monitoring statistic
    
    result = monitoringStat(
      img = Z,
      model = model,
      type = 1,
      stat = "ad",
      w = w
    )
    
    AD = as.numeric(result$localStat)
    
    # Run length
    
    t = t + 1
    
    flag = (AD <= h)
  }
  
  run_length = c(run_length, t)
  
  pb$tick()
}


# Results

cat("\nRun lengths:\n")
cat(run_length, "\n")

cat("\nMean run length:\n")
cat(mean(run_length), "\n")

cat("\nStandard deviation of run length:\n")
cat(sd(run_length), "\n")


# Plot

plot(
  run_length,
  type = "o",
  pch = 19,
  xlab = "Simulation",
  ylab = "Run length",
  main = "In-control run lengths (spc4sts)"
)

abline(h = 100, lty = 2)


library(magick)
library(progress)
library(spc4sts)

set.seed(0)

# Parameters

total = 10

f1 = 14
f2 = 7
f3 = 7
f4 = 7
f5 = 14

f_star = f1 + f2 + f3 + f4 + f5

w = 11
sig = 0.3

n1 = 128
n2 = 128
n = n1 * n2

# Threshold 

h = 450

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

img1 = read_image("star_line_grow_128/test1.png")
img2 = read_image("star_line_grow_128/test2.png")
img3 = read_image("star_line_grow_128/test3.png")
img4 = read_image("star_line_grow_128/test4.png")
img5 = read_image("star_line_grow_128/test5.png")

I = img1


# In-control surface model

model = surfacemodel(
  img = I,
  nb = 3,
  trim.vars = TRUE,
  cp = 1e-5,
  xval = 5,
  standardize = TRUE,
  subsample = 1,
  verbose = FALSE,
  keep.residuals = TRUE
)


# Change detection delay

delay = c()
neg_delay = c()
missed = 0

pb = progress_bar$new(
  format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
  total = total,
  clear = FALSE
)

for (r in 1:total) {
  
  # Generate final data
  
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
  
  
  # Anderson-Darling monitoring statistic
  
  AD = numeric(f_star)
  
  for (j in 1:f_star) {
    
    result = monitoringStat(
      img = W[j, , ],
      model = model,
      type = 1,
      stat = "ad",
      w = w
    )
    
    AD[j] = as.numeric(result$localStat)
  }
  
  
  # Detection
  
  detection = which(AD > h)
  
  change_point = f1 + 1
  
  if (length(detection) > 0) {
    
    first_detection = detection[1]
    
    if (first_detection >= change_point) {
      
      delay = c(delay, first_detection - change_point)
      
    } else {
      
      neg_delay = c(neg_delay, first_detection - change_point)
      
    }
    
  } else {
    
    missed = missed + 1
    
  }
  
  pb$tick()
}


# Results

cat("\nDetection delays:\n")
cat(delay, "\n")

cat("\nNegative delays (false detections):\n")
cat(neg_delay, "\n")

cat("\nMean detection delay:\n")
cat(mean(delay), "\n")

cat("\nStandard deviation of detection delay:\n")
cat(sd(delay), "\n")

cat("\nNumber of false detections:\n")
cat(length(neg_delay), "\n")

cat("\nFalse detection percentage:\n")
cat((length(neg_delay) / total) * 100, "%\n")

cat("\nNumber of missed detections:\n")
cat(missed, "\n")

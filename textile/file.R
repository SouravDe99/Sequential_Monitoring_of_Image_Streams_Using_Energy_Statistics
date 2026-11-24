library(magick)
library(spc4sts)
library(png)

set.seed(0)

n1 = 500
n2 = 500
n = n1 * n2

read_image = function(path) {
  img = image_read(path)
  img = image_convert(img, colorspace = "gray")
  img = image_data(img)
  img = matrix(as.numeric(img), nrow = n1, ncol = n2)
  return(img)
}

img1 = read_image("trainImg.png")
img2 = imposeDefect(img1*256, loc = c(100, 100), a = 11, b = 11, eps = 0.05, sigma = 0.1)$img/256
img3 = imposeDefect(img1*256, loc = c(100, 100), a = 13, b = 13, eps = 0.05, sigma = 0.1)$img/256
img4 = imposeDefect(img1*256, loc = c(100, 100), a = 16, b = 16, eps = 0.05, sigma = 0.1)$img/256
img5 = imposeDefect(img1*256, loc = c(100, 100), a = 19, b = 19, eps = 0.05, sigma = 0.1)$img/256

writePNG(img1, "test1.png")
writePNG(img2, "test2.png")
writePNG(img3, "test3.png")
writePNG(img4, "test4.png")
writePNG(img5, "test5.png")
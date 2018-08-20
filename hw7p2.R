setwd('C:/Users/Arun/Documents/cs498hw7')
library(matrixStats)
library(jpeg)
library(stringi)

GMM_EM = function(X, cluster.count, converge.threshold = 1e-4) {
  N = dim(X)[1]
  W = matrix(0, N, cluster.count)
  #Take kmean as the initial value for U and Pi
  kmeans_cluster = kmeans(X, cluster.count)
  
  U = kmeans_cluster$centers
  pi = kmeans_cluster$size / N
  
  old_W = W
  log_W = W
  
  for (iter in 1:10000) {
    log_pi = log(pi)
    #E Step
    for (j in 1:cluster.count) {
      D = t(X) - U[j,]
      log_W[,j] = (-0.5 * colSums(D * D)) + log_pi[j]
    }
    log_W = log_W - rowLogSumExps(log_W)
    
    W = exp(log_W)
    
    #M Step
    w_sum = colSums(W)
    U = crossprod(W, X) / w_sum
    pi = w_sum / N
    
    #Check convergence
    #I use mean and small convergence threshold for this model
    d = mean(abs(W - old_W))
    if(interactive()) print(paste("Iteration=", iter, ", D=", d, sep=""))
    if(d < converge.threshold) break()
    
    old_W = W
  }
  
  return (list(W = W, centers = U, pi=pi, iteration = iter))
}

segment_image = function(file.name, cluster.count) {
  img = readJPEG(file.name)
  old_dim = dim(img)
  dim(img) = c(dim(img)[1] * dim(img)[2], 3)
  
  r = GMM_EM(img * 255, cluster.count)
  r$centers = r$centers / 255
  
  for (i in 1:dim(r$W)[1]){
    img[i,] = r$centers[which.max(r$W[i,]),]
  }
  dim(img) = old_dim
  
  return (list(img = img, centers = r$centers, pi=r$pi))
}

plot_image = function(file.name) {
  img = readJPEG(file.name)
  plot(c(0, dim(img)[2]), c(0, dim(img)[1]), type="n", 
       main=paste(name, " (original)"),
       xlab = "X", ylab = "Y")
  rasterImage(img, 0, 0, dim(img)[2], dim(img)[1])
}
image_names = c("RobertMixed03.jpg", "smallstrelitzia.jpg", "smallsunset.jpg")
clusters = c(10, 20, 50)

for (name in image_names) {
  plot_image(name)
  for (c in clusters) {
    r = segment_image(name, c)
    plot(c(0, dim(r$img)[2]), c(0, dim(r$img)[1]), type="n", 
         main=paste(stri_sub(name,1,-4), ", ",c, " Segments", sep=""),
         xlab = "X", ylab = "Y")
    rasterImage(r$img, 0, 0, dim(r$img)[2], dim(r$img)[1]) 
  }
}
test_image = "smallsunset.jpg"
cluster = 20

results = list()

for (i in 1:5) {
  #make different start points for k-means
  set.seed(i * 197208)
  r = segment_image(test_image, cluster)
  plot(c(0, dim(r$img)[2]), c(0, dim(r$img)[1]), type="n", 
       main=paste(stri_sub(name,1,-5), ", Segments=", cluster, sep=""),
       xlab = "X", ylab = "Y")
  rasterImage(r$img, 0, 0, dim(r$img)[2], dim(r$img)[1])
  results[[i]] = list(centers = r$centers, pi = r$pi)
}
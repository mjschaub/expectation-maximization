setwd('D:/CS498/HW7 - EM/')


#problem 1
raw_data_table <-read.csv('docword.nips.txt',header=TRUE,sep=' ')
str(raw_data_table)
dim(raw_data_table)
raw_data_table[745000,]

raw_list <- read.csv('words_used.txt',header=FALSE,stringsAsFactors = FALSE)

num_topics <- 30
#1500 documents, 12419 words, 746316 total words
colnames(raw_data_table) <- c("document", "word", "count")

raw_data_table
document_vec <- matrix(0,1500,12419)
for(x in 1:746316)
{
  document_vec[raw_data_table[x,1], raw_data_table[x,2]] = raw_data_table[x,3]
}
#document_vec = [document #, word #] to find the count


probabilities <- array(0,dim=c(30,12419))
pi <- array(1/30,dim=c(1,30))

for(x in 1:30)
{
  rand_num <- runif(12419)
  prob_cols <- rand_num / sum(rand_num)
  probabilities[x,] <- prob_cols
}


Q_arr <- c()
while(TRUE) #EM steps
{
  #the E step
  log_of_prob <- log(probabilities)
  innermost_summation <- document_vec %*% t(log_of_prob) #[1500,30]
  weights <- array(0,dim=c(1500,30))
  for(x in 1:30)
  {
    weights[,x] <- innermost_summation[,x] + log(pi[x])
  }
  w_ij <- array(0,dim=c(1500,30))
  max_row_weight <- apply(weights,1,max)
  document_weight <- array(0,dim=c(1,1500))
  for(x in 1:1500)
  {
    document_weight[x] <- log(sum(exp(weights[x,]-max_row_weight[x])))#logSumExp(weights[x,] - max_row_weight[x])
  }
  curr_w <- weights - unlist(as.list(max_row_weight - document_weight))
  #unnormalized_w_ij <- exp(curr_w)
  for(x in 1:1500)
  {
    w_ij[x,] <- exp(curr_w)[x,]/sum(exp(curr_w)[x,]) #normalize the weights
  }
  
  curr_Q <- sum(weights*w_ij)
  Q_arr <- c(Q_arr,curr_Q)
  lagrange_mult <- .0001
  #the M step
  for(x in 1:30)
  {
    temp_numerator <- colSums(document_vec * w_ij[,x]) +lagrange_mult
    temp_denominator <- sum(rowSums(document_vec) * w_ij[,x]) + lagrange_mult*12419
    probabilities[x,] <- temp_numerator/temp_denominator
    pi[x] <- sum(w_ij[,x])/1500
  }
  prev_Q <- Q_arr[length(Q_arr)-1]
  if(length(Q_arr) > 1 && (curr_Q - prev_Q) < .001)
    break #stop when the Q doesn't change that much

}

plot(unlist(as.list(pi)), type='l', ylab = "probability", xlab="topic", main = "Probability a Topic is Selected")



#10 most common words
table_of_words <- c()
for(x in 1:30)
{
  tenth_max <- sort(probabilities[x,], decreasing = TRUE)[10]
  new_row <- raw_list[[1]][which(probabilities[x,] >= tenth_max)][1:10]
  table_of_words <- c(table_of_words,new_row)
  
}
table_of_words <- matrix(table_of_words,nrow=30)
table_of_words

###
###
###PROBLEM 2 WAS DONE IN THE OTHER FILE WITH CODE
###
###
#problem 2


# library(jpeg)
# 
# #quad <- readJPEG("UIUC_main_quad.jpg")
# rob_mixed <- readJPEG("RobertMixed03.jpg")
# 
# img255 <- rob_mixed*255
# rob_mixed.r <- matrix(img255[,,1],nrow = 480,ncol=640)
# rob_mixed.g <- matrix(img255[,,2],nrow = 480,ncol=640)
# rob_mixed.b <- matrix(img255[,,3],nrow = 480,ncol=640)
# n <- nrow(rob_mixed.r)*ncol(rob_mixed.r)
# k <- 50
# d <- 3
# 
# rob_mixed.all <- matrix(0,nrow=n,ncol=3)
# for(i in 1:nrow(rob_mixed.r)){
#   for(j in 1:ncol(rob_mixed.r)){
#     rob_mixed.all[640*(i-1)+j,] <- img255[i,j,]
#   }
# }
# 
# 
# 
# strelitzia <- readJPEG("smallstrelitzia.jpg")
# sunset <- readJPEG("smallsunset.jpg")
# 
# 
# EM_picture_method <- function(img, num_clusters)
# {
#   h <- dim(img)[1]
#   w <- dim(img)[2]
#   num_pixels <- h*w
#   
#   img_pixels <- array(0,dim=c(num_pixels,3))
#   for(i in 1:h)
#   {
#     for(j in 1:w)
#     {
#       img_pixels[((i-1)*w)+j,] <- img[i,j,]
#     }
#   }
#   
#   
#   img_pixels <- scale(img_pixels)
#   pi_two <- array(1/num_clusters,dim=c(1,num_clusters))
#   rand_nums <- runif(3*num_clusters)
#   mean_vals <- matrix(rand_nums, nrow=num_clusters)
#   
#   Q_arr_two <- c()
#   while(TRUE)
#   {
#     #E Step
#     innermost_summation_two <- array(0,dim=c(num_pixels, num_clusters))
#     for(x in 1:num_clusters)
#     {
#       distances <- t(t(img_pixels)-mean_vals[x,])
#       innermost_summation_two[,x] <- (-1/2) * rowSums(distances^2)
#     }
# 
#     temp_numerator_two <- exp(innermost_summation_two) %*% diag(pi_two[1:num_clusters])
#     temp_denominator_two <- rowSums(temp_numerator_two)
#     W_ij_two   <- temp_numerator_two/temp_denominator_two
#     curr_Q_two <- sum(innermost_summation_two*W_ij_two)
#     Q_arr_two <- c(Q_arr_two, curr_Q_two)
#     
#     #M step
#     for(x in 1:num_clusters)
#     {
#       temp_numerator_three <- colSums(img_pixels * W_ij_two[,x]) 
#       temp_denominator_three <- sum(W_ij_two[,x]) 
#       mean_vals[x,] <- temp_numerator_three/temp_denominator_three
#       pi_two[x] <- sum(W_ij_two[,x]) / num_pixels
#     }
#     
#     prev_Q_two <- Q_arr_two[length(Q_arr_two)-1]
#     if(length(Q_arr_two) > 1 && (curr_Q_two - prev_Q_two) < .0001)
#       break #stop when the Q doesn't change that much
#   }
#   
#   final_img <- array(0,c(h,w,3))
#   for(i in 1:h)
#   {
#     for(j in 1:w)
#     {
#       max_val_indices <- which(W_ij_two[(i-1)*w + j,] == max(W_ij_two[(i-1)*w + j,]))
#       final_img[i,j,] <- mean_vals[max_val_indices,]*attr(img_pixels, 'scaled:scale') + attr(img_pixels, 'scaled:center')
#     }
#   }
#   return(final_img)
# }
# 
# #writeJPEG(EM_picture_method(quad,10), "quad_segmented10.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,10), "sunset_final_10.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,20), "sunset_final_20.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,50), "sunset_final_50.jpg",quality = 1)
# 
# writeJPEG(EM_picture_method(strelitzia,10), "strelitzia_final_10.jpg",quality = 1)
# writeJPEG(EM_picture_method(strelitzia,20), "strelitzia_final_20.jpg",quality = 1)
# writeJPEG(EM_picture_method(strelitzia,50), "strelitzia_final_50.jpg",quality = 1)
# 
# writeJPEG(EM_picture_method(img255,10), "rob_mixed_final_10_2.jpg",quality = 1)
# writeJPEG(EM_picture_method(rob_mixed,20), "rob_mixed_final_20_2.jpg",quality = 1)
# writeJPEG(EM_picture_method(rob_mixed,50), "rob_mixed_final_50_2.jpg",quality = 1)
# 
# writeJPEG(EM_picture_method(sunset,20), "sunset_segment_1.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,20), "sunset_segment_2.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,20), "sunset_segment_3.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,20), "sunset_segment_4.jpg",quality = 1)
# writeJPEG(EM_picture_method(sunset,20), "sunset_segment_5.jpg",quality = 1)






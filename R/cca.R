# Three functions: 
# First computes the variance v(x) of 2 input matrix and covariance cov(x,y)
# Second computes resulting covariance matrix from list of smaller matrices
# Third computes 

# 1.
##########compute covariances#############
comput_cov <- function(x,y) {
  x <-  scale(x, scale = F)
  y <-  scale(y, scale = F)

  cx <- var(x, na.rm = TRUE, use = "pairwise") 
  cy <- var(y, na.rm = TRUE, use = "pairwise")
  cxy <- cov(x, y, use = "pairwise")

  return(list(cx = cx, cy = cy,cxy = cxy))
}

# 2.
######### Merge cov matrices#########
##########matrices need to have same n of columns and be scaled
merge_cov <- function(listcx = list(), listcy = list(), listcxy = list(list())){
  
  cx <-  Reduce("+" , lapply(listcx,function(x) (dim(x)[1]-1)*cov(x)))
  cov_tot = (Reduce("+", lapply(xl, function(x) dim(x)[1]))-1) *Reduce("+" , lapply(xl,function(x) (dim(x)[1]-1)*cov(x)))
  
  cy <- Reduce("+" , lapply(listcy,function(x) (dim(x)[1]-1)*cov(x)))
  
  cxy = 1/(Reduce("+",lapply(xyl, function(x) dim(x[[1]])[1])) -1 )*Reduce("+",lapply(xyl, function(x) (dim(x[[1]])[1]-1)*cov(x[[1]],x[[2]])))
  
  
  
  return(list(cx = cx, cy = cy, cxy =cxy ))
}

# 3.
#########compute coefficients matrices#########
cov_CCA <- function(cxy,cxx,cyy){
  require(CCA)
  res <- geigen(cxy, cxx, cyy)
  names(res) <- c("cor", "xcoef", "ycoef")
  
  return(list(cor = res$cor, xcoef = res$xcoef, 
              ycoef = res$ycoef))
}


##################### TESTS ######################

# 1) 
#We have a data file, mmreg.dta, with 600 observations on eight variables. The psychological 
#variables are locus_of_control, self_concept and motivation. The academic variables are standardized 
#tests in reading (read), writing (write), math (math) and science (science). Additionally, 
#the variable female is a zero-one indicator variable with the one indicating a female student.

mm <- read.csv("https://stats.idre.ucla.edu/stat/data/mmreg.csv") #600 x 8 matrix

x <- mm[, 1:3]
y <- mm[, 4:8]

res1 = comput_cov(x,y) # computes covariance matrices for both sets

res2 = cov_CCA(res1$cxy, res1$cx, res1$cy) #computes coefficients matrices

scores = comput(x,y, res2) #computes scores and loadings, this function center the data around the mean!

library(CCA)
out = cc(x,y); round(sum(out$scores$xscores - scores$xscores),0); round(sum(out$scores$yscores - scores$yscores),0)
sum(out$xcoef - res2$xcoef) #they provide the same results
#####################################################################################################################

#READ: AT THE END OF ALL OF THIS I TRY THE ACTUAL FUNCTION merge_Cov

#2)
# Test 2: Generate 6 random matrices (x1,x2,x3,y1,y2,y3) to test the merge_Cov function

#define and center x matrices
x1 = matrix(rnorm(200000), ncol = 400); x2 = matrix(rnorm(12000), ncol = 400); x3 = matrix(rnorm(240000), ncol = 400) 
dim(x1); dim(x2);dim(x3)
x1 = scale(x1, scale = F); x2 = scale(x2, scale = F); x3 = scale(x3, scale = F)
xl = list(x1, x2, x3) 
x = rbind(x1,x2,x3); dim(x) #matrix merged by rows


y1 = matrix(rnorm(200000), ncol = 400); y2 = matrix(rnorm(12000), ncol = 400); y3 = matrix(rnorm(24000), ncol = 400) 
dim(y1); dim(y2); dim(y3)
y1 = scale(y1, scale = F); y2 = scale(y2, scale = F); y3 = scale(y3, scale = F)
yl = list(y1, y2, y3) 
y = rbind(y1,y2,y3); dim(x) #matrix merged by rows


cov_tot = 1/(Reduce("+", lapply(xl, function(x) dim(x)[1]))-1) * Reduce("+" , lapply(xl,function(x) (dim(x)[1]-1)*cov(x)))
sum(cov(x) -cov_tot)

           #1/sum rows of each dataset -1                   #sum covariances   #computes covariances multiplied by rows
cov_tot = 1/(Reduce("+", lapply(yl, function(x) dim(x)[1]))-1) * Reduce("+" , lapply(yl,function(x) (dim(x)[1]-1)*cov(x)))

sum(cov(y) -cov_tot) # to prove we reach the same result


########################### Cov(x,y) #############################
#Generate again 6 random matrices to test if cov(x,y) merge function implemented works properly

x1 = matrix(rnorm(150000), ncol = 3000); x2 = matrix(rnorm(90000), ncol = 3000); x3 = matrix(rnorm(120000), ncol = 3000) 
dim(x1); dim(x2);dim(x3)
x1 = scale(x1, scale = F); x2 = scale(x2, scale = F); x3 = scale(x3, scale = F)
xl = list(x1, x2, x3) 
x = rbind(x1,x2,x3); dim(x) #matrix merged by rows

y1 = matrix(rnorm(150000), ncol = 3000); y2 = matrix(rnorm(90000), ncol = 3000); y3 = matrix(rnorm(120000), ncol = 3000) 
dim(y1); dim(y2); dim(y3)
y1 = scale(y1, scale = F); y2 = scale(y2, scale = F); y3 = scale(y3, scale = F)
yl = list(y1, y2, y3) 
y = rbind(y1,y2,y3); dim(x) #matrix merged by rows


xyl = list(list(x1, y1), list(x2,y2), list(x3,y3))

covxy_tot = 1/(Reduce("+",lapply(xyl, function(x) dim(x[[1]])[1])) -1 )*Reduce("+",lapply(xyl, function(x) (dim(x[[1]])[1]-1)*cov(x[[1]],x[[2]])))

sum(cov(x,y)- covxy_tot)


dim(x)
dim(y)

#### I try my functions
covariances = merge_cov(xl, yl, xyl)
cx = covariances$cx
cy = covariances$cy
cxy = covariances$cxy

make_PD = function(C, m){
  D = 0.5 * (C + t(C))
  D =  D + (m - min(eigen(D)$value) * diag(1, dim(D)))
  
  print(sum(C- D))
  return(D)
} # stupid function to make the matrix PD
cx = make_PD(cx, 0.02)
cy = make_PD(cy, 0.02)

res = cov_CCA(cxy, cx, cy) #non positive definite matrices give problems






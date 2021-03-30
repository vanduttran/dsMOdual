
##########compute covariances#############
comput_cov <- function(x,y) {
  x <-  scale(x, scale = F)
  y <-  scale(y, scale = F)

  cx <- var(x, na.rm = TRUE, use = "pairwise") 
  cy <- var(y, na.rm = TRUE, use = "pairwise")
  cxy <- cov(x, y, use = "pairwise")

  return(list(cx = cx, cy = cy,cxy = cxy))
}


#########merge cov matrices#########
##########matrices need to have same n of columns and be scaled
merge_cov <- function(listcx = list(), listcy = list(), listcxy = list(list())){
  
  cx <-  Reduce("+" , lapply(listcx,function(x) (dim(x)[1]-1)*cov(x)))
  
  cy <- Reduce("+" , lapply(listcy,function(x) (dim(x)[1]-1)*cov(x)))
  
  #Function for cxy to be done!
  
  return(list(cx = cx, cy = cy))
}


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


# Test 2: random matrices and merging function required this time

#define and center x matrices
x1 = matrix(rnorm(20), 5); x2 = matrix(rnorm(12), 3)
dim(x1); dim(x2)
x1 = scale(x1, scale = F); x2 = scale(x2, scale = F)
xl = list(x1, x2) 
x = rbind(x1,x2); dim(x) #matrix merged by rows


y1 = matrix(rnorm(20), 5); y2 = matrix(rnorm(12), 3)
dim(x1); dim(x2)
y1 = scale(x1, scale = F); y2 = scale(x2, scale = F)


cov_tot = Reduce("+" , lapply(xl,function(x) (dim(x)[1]-1)*cov(x)))

(dim(x)[1]-1) * cov(x) -cov_tot # to prove we reach the same result

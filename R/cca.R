
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
merge_cov <- function(listcx = list(), listcy = list(), listcxy = list()){
  cx = 0
  for (i in listcx){
    cx = cx + (dim(i)[1]-1)*var(i)
    
  }
  
  cy = 0
  for (i in listcy){
    cy = cy + (dim(i)[1]-1)*var(i)
  }
  
  #Function for cxy to be done!
  
  return(list(cx = cx, cy = cy))
}


#########compute coefficients matrices#########
cov_CCA <- function(cxy,cxx,cyy){
  
  res <- geigen(cxy, cxx, cyy)
  names(res) <- c("cor", "xcoef", "ycoef")
  #scores <- comput(X, Y, res)
  
  return(list(cor = res$cor, xcoef = res$xcoef, 
              ycoef = res$ycoef))
}



my = comput(x,y, res2)
res1 = comput_cov(x,y) #each node

res2 = cov_CCA(res1$cxy, res1$cx, res1$cy) #analyst

scores = comput(x,y, res2) #each node, to compute scores and loadings



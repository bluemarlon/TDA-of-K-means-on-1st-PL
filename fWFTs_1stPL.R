# It gives R code for constructing the fWFTs from a vector of univariate time series with length as a power of 2. 
# And the 1st order persistence landscapes from a vector of discretized univariate continuous function in general. 
################################################ fWFTs construction #################################################
fWFTs <- function(X, num){
  # Thanks to Dr. Stoffer providing the fWFTs code.
  # X: a univariate time series, with length num, 
  # num : the length of X, it should be a power of 2, namely 256, 512, 1024, etc.
  # the output is the fWFTs
  Ipower <- matrix(0, nrow = 1, ncol = 15)
  y <- rep(0, num)
  for(I in 1:num){
    # I <- I+1
    IB <- I-1
    IL <- 1
    while(TRUE){
      IBD <- floor(IB/2) 
      Ipower[IL] <- 1
      if(IB == IBD*2) Ipower[IL] <- 0
      if(IB == 0 | IBD==0) break
      IB <- IBD
      IL <- IL+1
    }
    IP <- 1
    IFAC <- num
    for(t1 in 1:IL){
      IFAC <- floor(IFAC/2)
      IP <- IP+IFAC * Ipower[t1]
    }
    y[IP] <- X[I] 
  }
  X <-  y
  Iter <- floor(log(num, base = 2))
  for(M in 1:Iter){
    if(M==1){nump <- 1}else{nump <- nump*2}
    Mnum <- floor(num/nump)
    Mnum2 <- floor(Mnum/2)
    alph <- 1
    for(MP in 1:nump){
      IB <- (MP-1)*Mnum
      for(MP2 in 1:Mnum2){
        mnum21 <- Mnum2+MP2+IB
        IBA <- IB+MP2
        y[IBA] <- X[IBA]+alph*X[mnum21]
        y[mnum21] <- X[IBA]-alph*X[mnum21]
      }
      alph <- -alph
    }
    r <- 1/sqrt(num)
    for(II in 1:num) X[II] <- y[II]*r
  }
  return(y*num^((log(num, base = 2)-2-1)/2) )
}
################################################ 1stPL construction #################################################
1stPL <- function(functions, PL.L, PL.range = range(functions)){
  # it uses to construct the 1st persistence landscape (PL) from a univariate continuous function
  # functions : a vector of the discretized univariate continuous function values
  # PL.L : the length of 1st PL to be constructed
  # PL.range : the range of 1st PL to be constructed. 
  # When there are multiple time series, the range should be the same for all t.s. in order to be comparable.
  # output will be the 1st PL. 
  dmin = min(x);dmax = max(x)
  ell = 1:PL.L
  land1 = c(apply(cbind(PL.range[1] + (ell-1)*(PL.range[2]-PL.range[1])/(PL.L-1)-dmin,
                dmax-PL.range[1]-(ell-1)*(PL.range[2]-PL.range[1])/(PL.L-1)), 1,
          FUN = function(x){
            min(max(x[1], 0),
                max(x[2], 0) )
          })
  )
  return(land1)
}



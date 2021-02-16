##########################################################################
##########################################################################
##########################################################################
library(LACPD)
library(parallel)
library(trend)
library(strucchange)
library(ecp)
library(bfast)
library(wbs)
library(AdaptiveCpt)

##########################################################################
########################################################################## simulation
##########################################################################
set.seed(1234)
X <- replicate(500,rnorm(200),simplify = FALSE)

##########################################################################
########################################################################## LACPD calculation
##########################################################################

R <- mclapply(X=1:length(X), function(i){

  if(i<length(X)) {
    cat(paste(i),",")
    flush.console()
  } else {
    cat(paste(i),"\n")
    flush.console()
  }

  lacpd_mw(X[[i]],m=100,k=c(2:10),adjust = TRUE,method="BY",history = TRUE)
},mc.cores = detectCores()-4)

##########################################################################
########################################################################## pettitt
##########################################################################

pet <- lapply(X=1:length(X), function(i){
  pettitt.test(X[[i]])
})

##########################################################################
########################################################################## strucchange
##########################################################################

str <- unlist(lapply(X=1:length(X), function(i){
  strucchange::breakpoints(X[[i]]~1,h=0.1)$breakpoints
}))

##########################################################################
########################################################################## Buishand
##########################################################################

br <- lapply(X=1:length(X), function(i){
  br.test(X[[i]])
})

bu <- lapply(X=1:length(X), function(i){
  bu.test(X[[i]])
})

##########################################################################
########################################################################## e.divisive
##########################################################################

e.cp <- mclapply(X=1:length(X), function(i){
    e.divisive(as.matrix(X[[i]],ncol=1))
})

##########################################################################
########################################################################## bfast
##########################################################################

bf <- mclapply(X=1:length(X), function(i){
    bfast(ts(X[[i]]),season = "none",max.iter = 1)$Time
  })

##########################################################################
########################################################################## wbs
##########################################################################

wbs.cp <- mclapply(X=1:length(X), function(i){
  changepoints(wbs(X[[i]]))
})

##########################################################################
########################################################################## AdaptiveCpt
##########################################################################

adapt.cp <- mclapply(X=1:length(X), function(i){
  cpt.mean.test(matrix(X[[i]],ncol = 1),S0=200,
                P=1000,tau0=0.1,B=500,type = 1)
})

######################################################
###################################################### FP
######################################################
alph <- seq(0,1,by=0.01)

LACPD_fp <- c()
for (i in 1:length(alph)) {
  LACPD_fp[i] <-  sum(unlist(lapply(X=1:length(X), function(i)
    R[[i]]$p
  )
  )<=alph[i])/length(X)
}
# [1] 0.006

pet_fp <- c()
for (i in 1:length(alph)) {
  pet_fp[i] <-  sum(unlist(lapply(X=1:length(X), function(i)
    pet[[i]]$p.value
  )
  )<=alph[i])/length(X)
}
# [1] 0.032


br_fp <- c()
for (i in 1:length(alph)) {
  br_fp[i] <-  sum(unlist(lapply(X=1:length(X), function(i)
    br[[i]]$p.value
  )
  )<=alph[i])/length(X)
}

# [1] 0.044

bu_fp <- c()
for (i in 1:length(alph)) {
  bu_fp[i] <-  sum(unlist(lapply(X=1:length(X), function(i)
    bu[[i]]$p.value
  )
  )<=alph[i])/length(X)
}
# [1] 0.052


sum(!is.na(str))/length(X)
# [1] 0.024

sum(!is.na(unlist(bf)))/length(X)
# [1] 0.002


ecp_fp <- c()
for (i in 1:length(alph)) {
  ecp_fp[i] <-  sum(unlist(lapply(X=1:length(X), function(i)
    e.cp[[i]]$p.value[1]
  )
  )<=alph[i])/length(X)
}
# [1] 0.032


wbs.cp.f <- c()
for (i in 1:length(X)) {
  wbs.cp.f[i] <- unlist(wbs.cp[[i]]$cpt.ic$bic.penalty)
}
sum(!is.na(wbs.cp.f))/length(X)
# [1] 0.064

adapt.cp_fp <- c()
for (i in 1:length(alph)) {
  adapt.cp_fp[i] <-  sum(unlist(lapply(X=1:length(X), function(i)
    as.numeric(adapt.cp[[i]]$all.adaptive.pvalue)
  )
  )<=alph[i])/length(X)
}

# [1] 0.054

# plot(alph,LACPD_fp,type = "l",ylab = "FP")
# points(alph,pet_fp,type = "l",ylab = "FP",add=T,col=2)
# points(alph,br_fp,type = "l",ylab = "FP",add=T,col=3)
# points(alph,bu_fp,type = "l",ylab = "FP",add=T,col=4)
# points(alph,ecp_fp,type = "l",ylab = "FP",add=T,col=5)
# points(alph,adapt.cp_fp,type = "l",ylab = "FP",add=T,col=6)

save(alph,LACPD_fp,pet_fp,br_fp
       ,bu_fp,ecp_fp,adapt.cp_fp,file = "FPS.RData")

save.image("outfalseNormal(BY_double_FALSE).RData")

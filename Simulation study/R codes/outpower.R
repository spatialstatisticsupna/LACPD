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

cpvec <- c(40,80,100,120,160)

##########################################################################
########################################################################## LACPD calculation
##########################################################################

R <- list()

for (j in 1:length(cpvec)) {

  Y <- list()
  for (i in 1:length(X)) {
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  R[[j]] <- mclapply(X=1:length(X), function(i){

    if(i<length(X)) {
      cat(paste(i),",")
      flush.console()
    } else {
      cat(paste(i),"\n")
      flush.console()
    }

    lacpd_mw(Y[[i]],m=100,k=c(2:10),adjust = TRUE,history = TRUE,double = FALSE,method="BY")
  },mc.cores = detectCores()-4)
}


##########################################################################
########################################################################## pettitt
##########################################################################

pet <- list()
for (j in 1:length(cpvec)) {

  Y <- list()
  for (i in 1:length(X)) {
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  pet[[j]] <- mclapply(X=1:length(X), function(i){
    pettitt.test(Y[[i]])
  },mc.cores = detectCores()-4)

}

##########################################################################
########################################################################## Buishand
##########################################################################

br <- list()
for (j in 1:length(cpvec)) {

  Y <- list()
  for (i in 1:length(X)) {
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  br[[j]] <- mclapply(X=1:length(X), function(i){
    br.test(Y[[i]])
  },mc.cores = detectCores()-4)

}


bu <- list()
for (j in 1:length(cpvec)) {

  Y <- list()
  for (i in 1:length(X)) {
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  bu[[j]] <- mclapply(X=1:length(X), function(i){
    bu.test(Y[[i]])
  },mc.cores = detectCores()-4)

}

##########################################################################
########################################################################## strucchange
##########################################################################

str <- list()
for (j in 1:length(cpvec)){

  Y <- list()
  for (i in 1:length(X)){
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  str[[j]] <- mclapply(X=1:length(X), function(i){
    strucchange::breakpoints(Y[[i]]~1,breaks=1,h=0.15)$breakpoints
  },mc.cores = detectCores()-4)

}

##########################################################################
########################################################################## e.divisive
##########################################################################

e.cp <- list()
for (j in 1:length(cpvec)){

  Y <- list()
  for (i in 1:length(X)){
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  e.cp[[j]] <- mclapply(X=1:length(X), function(i){
    e.divisive(as.matrix(Y[[i]],ncol=1),k=1)$considered.last
  })

}

##########################################################################
########################################################################## bfast
##########################################################################

bf <- list()
for (j in 1:length(cpvec)){

  Y <- list()
  for (i in 1:length(X)){
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  bf[[j]] <- mclapply(X=1:length(X), function(i){
    bfast(ts(Y[[i]]),season = "none",breaks = 1,max.iter = 1)$Time
  })

}

##########################################################################
########################################################################## wbs
##########################################################################

wbs.cp <- list()
for (j in 1:length(cpvec)){

  Y <- list()
  for (i in 1:length(X)){
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  wbs.cp[[j]] <- mclapply(X=1:length(X), function(i){
    changepoints(wbs(Y[[i]]),Kmax = 1)$cpt.ic$t
  })

}

##########################################################################
########################################################################## AdaptiveCpt
##########################################################################

adapt.cp <- list()
for (j in 1:length(cpvec)){

  Y <- list()
  for (i in 1:length(X)){
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  adapt.cp[[j]] <- mclapply(X=1:length(X), function(i){
    cpt.mean.test(matrix(Y[[i]],ncol = 1),S0=200,
                  P=1000,tau0=0.1,B=500,type = 1)
  })

}


######################################################
###################################################### power
######################################################

alph <- seq(0,1,by=0.01)

lacpd_p <- list()
for (k in 1:length(alph)) {
  lacpd_power <- c()
  for (j in 1:length(cpvec)) {

    lacpd_power[j] <- sum(unlist(lapply(X=1:length(X), function(i)
      R[[j]][[i]]$p
    )
    )<=alph[k])/length(X)
  }
  lacpd_p[[k]] <- lacpd_power
}

pet_p <- list()
for (k in 1:length(alph)) {
  pet_power <- c()
  for (j in 1:length(cpvec)) {

    pet_power[j] <- sum(unlist(lapply(X=1:length(X), function(i)
      pet[[j]][[i]]$p.value
    )
    )<alph[k])/length(X)
  }
  pet_p[[k]] <- pet_power
}

br_p <- list()
for (k in 1:length(alph)) {
  br_power <- c()
  for (j in 1:length(cpvec)) {

    br_power[j] <- sum(unlist(lapply(X=1:length(X), function(i)
      br[[j]][[i]]$p.value
    )
    )<alph[k])/length(X)
  }
  br_p[[k]] <- br_power
}

bu_p <- list()
for (k in 1:length(alph)) {
  bu_power <- c()
  for (j in 1:length(cpvec)) {

    bu_power[j] <- sum(unlist(lapply(X=1:length(X), function(i)
      bu[[j]][[i]]$p.value
    )
    )<alph[k])/length(X)
  }
  bu_p[[k]] <- bu_power
}


str_power <- c()
for (j in 1:length(cpvec)) {
  str_power[j] <- sum(!is.na(unlist(str[[j]])))/length(X)
}


e.cp1 <- list()
for (j in 1:length(cpvec)){

  Y <- list()
  for (i in 1:length(X)){
    Y[[i]] <- c(X[[i]][1:(cpvec[j])],X[[i]][(cpvec[j]+1):200]+1)
  }

  e.cp1[[j]] <- mclapply(X=1:length(X), function(i){
    e.divisive(as.matrix(Y[[i]],ncol=1),k=1)
  })

}

ecp_p <- list()
for (k in 1:length(alph)) {
  ecp_power <- c()
  for (j in 1:length(cpvec)) {
    ecp_power[j] <- sum(unlist(lapply(X=1:length(X), function(i)
      e.cp1[[j]][[i]]$p.value[1]
    )
    )<alph[k])/length(X)
  }
  ecp_p[[k]] <- ecp_power
}

bf_power <- c()
for (j in 1:length(cpvec)) {
  bf_power[j] <- sum(!is.na(as.vector(unlist(bf[[j]]))))/length(X)
}

wbs_power <- c()
for (j in 1:length(cpvec)) {
  wbs_power[j] <- sum(!is.na(as.vector(unlist(wbs.cp[[j]]))))/length(X)
}

adapt.cp_p <- list()
for (k in 1:length(alph)) {
  adapt.cp_power <- c()
  for (j in 1:length(cpvec)) {

    adapt.cp_power[j] <- sum(unlist(lapply(X=1:length(X), function(i)
      as.numeric(adapt.cp[[j]][[i]]$all.adaptive.pvalue)
    )
    )<alph[k])/length(X)
  }
  adapt.cp_p[[k]] <-  adapt.cp_power
}

save(adapt.cp_p,ecp_p ,bu_p,br_p,pet_p,lacpd_p,file = "PS.RData")

######################################################
###################################################### Absolute Bias
######################################################

lacpd_bias <- c()
for (j in 1:length(cpvec)) {

  ins <- unlist(lapply(X=1:length(X), function(i)
    R[[j]][[i]]$cp
  ))
  pv <- unlist(lapply(X=1:length(X), function(i)
    R[[j]][[i]]$p
  ))
  lacpd_bias[j] <- abs((sum(ins[pv<0.05])/length(which(pv<0.05))) - cpvec[j])
}

pet_bias <- c()
for (j in 1:length(cpvec)) {

  ins.pet <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(pet[[j]][[i]]$estimate))
  ))
  pv.pet <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(pet[[j]][[i]]$p.value))
  ))
  pet_bias[j] <- abs((sum(ins.pet[pv.pet<0.05])/length(which(pv.pet<0.05)))- cpvec[j])
}

br_bias <- c()
for (j in 1:length(cpvec)) {

  ins.br <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(br[[j]][[i]]$estimate))
  ))
  pv.br <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(br[[j]][[i]]$p.value))
  ))
 br_bias[j] <- abs((sum(ins.br[pv.br<0.05])/length(which(pv.br<0.05)))- cpvec[j])
}

bu_bias <- c()
for (j in 1:length(cpvec)) {
  ins.bu <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(bu[[j]][[i]]$estimate))
  ))
  pv.bu<- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(bu[[j]][[i]]$p.value))
  ))
  bu_bias[j] <- abs((sum(ins.bu[pv.bu<0.05])/length(which(pv.bu<0.05)))- cpvec[j])
}

str_bias <- c()
for (j in 1:length(cpvec)) {
  str_bias[j] <- abs(mean(unlist(str[[j]]),na.rm = TRUE) - cpvec[j])
}


ecp_bias <- c()
for (j in 1:length(cpvec)) {
  ecp_bias[j] <- abs(mean(unlist(e.cp[[j]]),na.rm = TRUE) - cpvec[j])
}

bf_bias <- c()
for (j in 1:length(cpvec)) {
  bf_bias[j] <- abs(mean(unlist(bf[[j]]),na.rm = TRUE) - cpvec[j])
}

wbs_bias <- c()
for (j in 1:length(cpvec)) {
  wbs_bias[j] <- abs(mean(unlist(wbs.cp[[j]]),na.rm = TRUE) - cpvec[j])
}

adapt.cp_bias <- c()
for (j in 1:length(cpvec)) {

  ins <- unlist(lapply(X=1:length(X), function(i)
    as.numeric(adapt.cp[[j]][[i]]$cpt.est)
  ))
  pv <- unlist(lapply(X=1:length(X), function(i)
    as.numeric(adapt.cp[[j]][[i]]$all.adaptive.pvalue)
  ))
  adapt.cp_bias[j] <- abs((sum(ins[pv<0.05])/length(which(pv<0.05))) - cpvec[j])
}

######################################################
###################################################### MAE
######################################################

lacpd_mae <- c()
for (j in 1:length(cpvec)) {
  pv <- unlist(lapply(X=1:length(X), function(i)
    R[[j]][[i]]$p
  ))
  lacpd_mae[j] <- mean(abs((unlist(lapply(X=1:length(X), function(i)
    abs(R[[j]][[i]]$cp)
  )
  )) - cpvec[j])[pv<0.05])
}


pet_mae <- c()
for (j in 1:length(cpvec)) {
  pv.pet <- unlist(lapply(X=1:length(X), function(i)
    as.numeric(pet[[j]][[i]]$p.value)
  ))
  pet_mae[j] <-mean(abs((unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(pet[[j]][[i]]$estimate))
  )
  )) - cpvec[j])[pv.pet<0.05])
}

br_mae <- c()
for (j in 1:length(cpvec)) {
  pv.br <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(br[[j]][[i]]$p.value))
  ))
  br_mae[j] <-mean(abs((unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(br[[j]][[i]]$estimate))
  )
  )) - cpvec[j])[pv.br<0.05])
}

bu_mae <- c()
for (j in 1:length(cpvec)) {
  pv.bu<- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(bu[[j]][[i]]$p.value))
  ))
  bu_mae[j] <-mean(abs((unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(bu[[j]][[i]]$estimate))
  )
  )) - cpvec[j])[pv.bu<0.05])
}

str_mae <- c()
for (j in 1:length(cpvec)) {
  str_mae[j] <-mean(abs(unlist(str[[j]]) - cpvec[j]),na.rm = TRUE)
}


ecp_mae <- c()
for (j in 1:length(cpvec)) {
  ecp_mae[j] <- mean( abs(unlist(e.cp[[j]])  - cpvec[j]) ,na.rm = TRUE)
}

bf_mae <- c()
for (j in 1:length(cpvec)) {
  bf_mae[j] <- mean( abs(unlist(bf[[j]])  - cpvec[j]) ,na.rm = TRUE)
}

wbs_mae <- c()
for (j in 1:length(cpvec)) {
  wbs_mae[j] <- mean( abs(unlist(wbs.cp[[j]])  - cpvec[j]) ,na.rm = TRUE)
}

adapt.cp_mae <- c()
for (j in 1:length(cpvec)) {
  pv <- unlist(lapply(X=1:length(X), function(i)
    as.numeric(adapt.cp[[j]][[i]]$all.adaptive.pvalue)
  ))
  adapt.cp_mae[j] <- mean(abs((unlist(lapply(X=1:length(X), function(i)
    abs(adapt.cp[[j]][[i]]$cpt.est)
  )
  )) - cpvec[j])[pv<0.05])
}

######################################################
###################################################### variance
######################################################

lacpd_var <- c()
for (j in 1:length(cpvec)) {

  pv <- unlist(lapply(X=1:length(X), function(i)
    R[[j]][[i]]$p
  ))

  av <- sum(unlist(lapply(X=1:length(X), function(i)
    abs(R[[j]][[i]]$cp)
  )
  )[pv<0.05])/length(which(pv<0.05))

  lacpd_var[j] <- mean(

  (unlist(lapply(X=1:length(X), function(i)
    R[[j]][[i]]$cp
  )[pv<0.05]
  )-av)^2)

}


pet_var <- c()
for (j in 1:length(cpvec)) {
  pv.pet <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(pet[[j]][[i]]$p.value))
  ))
  av <- sum(unlist(lapply(X=1:length(X), function(i)
    abs(pet[[j]][[i]]$estimate)
  ))[pv.pet<0.05])/length(which(pv.pet<0.05))

  pet_var[j] <- mean((unlist(lapply(X=1:length(X), function(i)
    pet[[j]][[i]]$estimate
  )
  )-av)^2)

}

br_var <- c()
for (j in 1:length(cpvec)) {
  pv.br <- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(br[[j]][[i]]$p.value))
  ))
  av <- sum(unlist(lapply(X=1:length(X), function(i)
    abs(br[[j]][[i]]$estimate)
  )
  )[pv.br<0.05])/length(which(pv.br<0.05))
  br_var[j] <- mean((unlist(lapply(X=1:length(X), function(i)
    br[[j]][[i]]$estimate
  )
  )-av)^2)

}

bu_var <- c()
for (j in 1:length(cpvec)) {

  pv.bu<- unlist(lapply(X=1:length(X), function(i)
    abs(as.numeric(bu[[j]][[i]]$p.value))
  ))
  av <- sum(unlist(lapply(X=1:length(X), function(i)
    abs(bu[[j]][[i]]$estimate)
  )[pv.bu<0.05])/length(which(pv.bu<0.05)))
  bu_var[j] <- mean((unlist(lapply(X=1:length(X), function(i)
    bu[[j]][[i]]$estimate
  )
  )-av)^2)

}


str_var <- c()
for (j in 1:length(cpvec)) {

  str_var[j] <- mean((unlist(str[[j]])-mean(unlist(str[[j]]),na.rm = TRUE))^2,na.rm = TRUE)

}


ecp_var <- c()
for (j in 1:length(cpvec)) {

  ecp_var[j] <- mean((unlist(e.cp[[j]])-mean(unlist(e.cp[[j]]),na.rm = TRUE))^2,na.rm = TRUE)

}

bf_var <- c()
for (j in 1:length(cpvec)) {

  bf_var[j] <- mean((unlist(bf[[j]])-mean(unlist(bf[[j]]),na.rm = TRUE))^2,na.rm = TRUE)

}

wbs_var <- c()
for (j in 1:length(cpvec)) {

  wbs_var[j] <- mean((unlist(wbs.cp[[j]])-mean(unlist(wbs.cp[[j]]),na.rm = TRUE))^2,na.rm = TRUE)

}

adapt.cp_var <- c()
for (j in 1:length(cpvec)) {

  pv <- unlist(lapply(X=1:length(X), function(i)
    as.numeric(adapt.cp[[j]][[i]]$all.adaptive.pvalue)
  ))

  av <- sum(unlist(lapply(X=1:length(X), function(i)
    as.numeric(adapt.cp[[j]][[i]]$cpt.est)
  )
  )[pv<0.05])/length(which(pv<0.05))

  adapt.cp_var [j] <- mean(

    (unlist(lapply(X=1:length(X), function(i)
      as.numeric(adapt.cp[[j]][[i]]$cpt.est)
    )[pv<0.05]
    )-av)^2)

}

save.image("outpowerNormal(Adjust_TRUE_double_FALSE).RData")

##############################################
############################################## plots
##############################################
############################# data preparation
# load("./outfalseNormal(BY_double_FALSE).RData")

allzs <- list()
for (j in 1:length(cpvec)) {
  allzs[[j]] <- do.call(rbind, lapply(X=1:length(R[[j]]), function(i){
    attr(R[[j]][[i]],"zs")
  }))
}

# ts.plot(colMeans(allzs[[1]]))
# R[[1]][[1]]$s[which.min(colMeans(allzs[[1]]))]

allps <- list()
for (j in 1:length(cpvec)) {
  allps[[j]] <- do.call(rbind, lapply(X=1:length(R[[j]]), function(i){
    attr(R[[j]][[i]],"ps")
  }))
}

# ts.plot(colMeans(allps[[1]]))
# R[[1]][[1]]$s[which.min(colMeans(allps[[1]]))]
# R[[1]][[1]]$s[which.min(colMeans(allps[[2]]))]
# R[[1]][[1]]$s[which.min(colMeans(allps[[3]]))]
# R[[1]][[1]]$s[which.min(colMeans(allps[[4]]))]
# R[[1]][[1]]$s[which.min(colMeans(allps[[5]]))]

allmags <- list()
for (j in 1:length(cpvec)) {
  allmags[[j]] <- do.call(rbind, lapply(X=1:length(R[[j]]), function(i){
    attr(R[[j]][[i]],"mags")
  }))
}

# ts.plot(colMeans(allmags[[1]]))
# R[[1]][[1]]$s[which.max(colMeans(allmags[[1]]))]


S <- R[[1]][[1]]$s

df1.mag <- data.frame(x=S,y=as.numeric(colMeans(allmags[[1]])))
df1.mag <- cbind(df1.mag,ID=rep("40",161))

df2.mag <- data.frame(x=S,y=as.numeric(colMeans(allmags[[2]])))
df2.mag <- cbind(df2.mag,ID=rep("80",161))

df3.mag <- data.frame(x=S,y=as.numeric(colMeans(allmags[[3]])))
df3.mag <- cbind(df3.mag,ID=rep("100",161))

df4.mag <- data.frame(x=S,y=as.numeric(colMeans(allmags[[4]])))
df4.mag <- cbind(df4.mag,ID=rep("120",161))


df5.mag <- data.frame(x=S,y=as.numeric(colMeans(allmags[[5]])))
df5.mag <- cbind(df5.mag,ID=rep("160",161))

df.mag <- rbind(df3.mag,df4.mag,df5.mag,df1.mag,df2.mag)
type <- rep("magnitude",nrow(df.mag))
df.mag <- cbind(df.mag,type)


df1.p <- data.frame(x=S,y=as.numeric(colMeans(allps[[1]])))
df1.p <- cbind(df1.p,ID=rep("40",161))

df2.p <- data.frame(x=S,y=as.numeric(colMeans(allps[[2]])))
df2.p <- cbind(df2.p,ID=rep("80",161))

df3.p <- data.frame(x=S,y=as.numeric(colMeans(allps[[3]])))
df3.p <- cbind(df3.p,ID=rep("100",161))

df4.p <- data.frame(x=S,y=as.numeric(colMeans(allps[[4]])))
df4.p <- cbind(df4.p,ID=rep("120",161))


df5.p <- data.frame(x=S,y=as.numeric(colMeans(allps[[5]])))
df5.p <- cbind(df5.p,ID=rep("160",161))

df.p <- rbind(df3.p,df4.p,df5.p,df1.p,df2.p)
type <- rep("p-value",nrow(df.p))
df.p <- cbind(df.p,type)


df1 <- data.frame(x=S,y=as.numeric(colMeans(allzs[[1]])))
df1 <- cbind(df1,ID=rep("40",161))

df2 <- data.frame(x=S,y=as.numeric(colMeans(allzs[[2]])))
df2 <- cbind(df2,ID=rep("80",161))

df3 <- data.frame(x=S,y=as.numeric(colMeans(allzs[[3]])))
df3 <- cbind(df3,ID=rep("100",161))

df4 <- data.frame(x=S,y=as.numeric(colMeans(allzs[[4]])))
df4 <- cbind(df4,ID=rep("120",161))

df5 <- data.frame(x=S,y=as.numeric(colMeans(allzs[[5]])))
df5 <- cbind(df5,ID=rep("160",161))

df <- rbind(df3,df4,df5,df1,df2)
type <- rep("Z",nrow(df))
df <- cbind(df,type)

S[which.min(colMeans(allps[[1]]))]
S[which.min(colMeans(allps[[2]]))]
S[which.min(colMeans(allps[[3]]))]
S[which.min(colMeans(allps[[4]]))]
S[which.min(colMeans(allps[[5]]))]
S[which(colMeans(allps[[1]])<0.05)]
S[which(colMeans(allps[[2]])<0.05)]
S[which(colMeans(allps[[3]])<0.05)]
S[which(colMeans(allps[[4]])<0.05)]
S[which(colMeans(allps[[5]])<0.05)]

dpzmags <- rbind(df,df.p,df.mag)
dpzmags$fac <- factor(dpzmags$ID,levels=c('40','80','100','120','160'))

library(ggplot2)
library(forcats)
png("normalshift-BY-pzmag.png",width = 2000, height = 2000)
ggplot(dpzmags,aes(x,y))+geom_line(aes(x,y),lwd=2)+
  facet_grid(fct_rev(type) ~ fct_rev(fct_rev(fac)), scales = "free")+
  geom_vline(data = data.frame(xint=41,type="Z",fac="40"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=81,type="Z",fac="80"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=101,type="Z",fac="100"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=121,type="Z",fac="120"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=161,type="Z",fac="160"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=41,type="p-value",fac="40"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=81,type="p-value",fac="80"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=101,type="p-value",fac="100"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=121,type="p-value",fac="120"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=161,type="p-value",fac="160"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=41,type="magnitude",fac="40"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=81,type="magnitude",fac="80"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=101,type="magnitude",fac="100"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=121,type="magnitude",fac="120"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=161,type="magnitude",fac="160"), color=2, aes(xintercept = xint),lwd=1.5)+
  geom_vline(data = data.frame(xint=37,type="p-value",fac="40"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=45,type="p-value",fac="40"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=74,type="p-value",fac="80"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=87,type="p-value",fac="80"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=95,type="p-value",fac="100"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=106,type="p-value",fac="100"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=115,type="p-value",fac="120"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=126,type="p-value",fac="120"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=157,type="p-value",fac="160"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_vline(data = data.frame(xint=164,type="p-value",fac="160"), color="gold", aes(xintercept = xint),lwd=1.5,linetype="dashed")+
  geom_hline(data = data.frame(yint=0.05,type="p-value"), color="brown", aes(yintercept = yint),lwd=1.5)+
  labs(x = "Time",y= "")+
  theme(legend.position="bottom",aspect.ratio=1.5,axis.text=element_text(size=30),
        axis.text.x = element_text(angle = 315),
        axis.title=element_text(size=45),strip.text = element_text(size = 45)
        ,legend.text=element_text(size=rel(3)),legend.title=element_text(size=45))+
  scale_x_continuous(breaks = c(40,80,100,120,160))
dev.off()


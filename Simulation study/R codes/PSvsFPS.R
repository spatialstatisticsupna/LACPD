##############################################
############################################## plots
##############################################
############################# data preparation
load("./PS.RData")
load("./FPS.RData")

lacpd_p <- do.call(rbind,lacpd_p)
pet_p <- do.call(rbind,pet_p )
br_p <- do.call(rbind,br_p)
bu_p <- do.call(rbind,bu_p)
ecp_p <- do.call(rbind,ecp_p)
adapt.cp_p <- do.call(rbind,adapt.cp_p)

######################################################


names_methods <- c("LACPD","Pettitt","Buishand U","e.divisive","AdaptiveCpt")
Fps_all <- c(LACPD_fp,pet_fp,bu_fp,ecp_fp,adapt.cp_fp)

df_40 <- data.frame(ROC=c(lacpd_p[,1],pet_p[,1],bu_p[,1],ecp_p[,1],adapt.cp_p[,1]),time=rep(40,101),method=rep(names_methods,each=101),
                      type=rep("ROC",505))
df_80 <- data.frame(ROC=c(lacpd_p[,2],pet_p[,2],bu_p[,2],ecp_p[,2],adapt.cp_p[,2]),time=rep(80,101),method=rep(names_methods,each=101),
                    type=rep("ROC",505))
df_100 <- data.frame(ROC=c(lacpd_p[,3],pet_p[,3],bu_p[,3],ecp_p[,3],adapt.cp_p[,3]),time=rep(100,101),method=rep(names_methods,each=101),
                     type=rep("ROC",505))
df_120 <- data.frame(ROC=c(lacpd_p[,4],pet_p[,4],bu_p[,4],ecp_p[,4],adapt.cp_p[,4]),time=rep(120,101),method=rep(names_methods,each=101),
                     type=rep("ROC",505))
df_160 <- data.frame(ROC=c(lacpd_p[,5],pet_p[,5],bu_p[,5],ecp_p[,5],adapt.cp_p[,5]),time=rep(160,101),method=rep(names_methods,each=101),
                     type=rep("ROC",505))

df_crit <- rbind(df_40,df_80,df_100,df_120,df_160)
df_crit <- cbind(df_crit,alpha=rep(alph,5),Fp=rep(Fps_all,5))
df_crit <- cbind(df_crit,y=rep(0,nrow(df_crit)))
df_crit <- df_crit[df_crit$Fp<0.05,]


png("FPvsTP.png",width = 3000, height = 1200)
qplot(`Fp`,`ROC`,data = df_crit, size=I(2)
      ,color=method,facets = ~ time, group = method)+geom_line(lwd=2)+
  facet_grid(~time , scales = "free")+ labs(y="Power",x="FP")+
  scale_shape_manual(values=1:nlevels(df_crit$method))+
  theme(legend.position="bottom",legend.key.width = unit(5,"line"),
        aspect.ratio=1.5,axis.text=element_text(size=30),
        axis.text.x = element_text(angle = 315),
        axis.title=element_text(size=45),strip.text = element_text(size = 45),
        legend.text=element_text(size=rel(5)),legend.title=element_text(size=55))+
  # scale_x_continuous(breaks = c(40,80,100,120,160))+
  guides(shape = guide_legend(override.aes = list(size = 10,alpha = 5)))
dev.off()



##################################### AUC

round(sintegral(LACPD_fp[LACPD_fp<0.05],lacpd_p[,1][LACPD_fp<0.05],201)$value,3)
round(sintegral(pet_fp[pet_fp<0.05],pet_p[,1][pet_fp<0.05],201)$value,3)
round(sintegral(adapt.cp_fp[adapt.cp_fp<0.05],adapt.cp_p[,1][adapt.cp_fp<0.05],201)$value,3)
round(sintegral(br_fp[br_fp<0.05],br_p[,1][br_fp<0.05],201)$value,3)
round(sintegral(bu_fp[bu_fp<0.05],bu_p[,1][bu_fp<0.05],201)$value,3)
round(sintegral(ecp_fp[ecp_fp<0.05],ecp_p[,1][ecp_fp<0.05],201)$value,3)


round(sintegral(LACPD_fp[LACPD_fp<0.05],lacpd_p[,2][LACPD_fp<0.05],201)$value,3)
round(sintegral(pet_fp[pet_fp<0.05],pet_p[,2][pet_fp<0.05],201)$value,3)
round(sintegral(adapt.cp_fp[adapt.cp_fp<0.05],adapt.cp_p[,2][adapt.cp_fp<0.05],201)$value,3)
round(sintegral(br_fp[br_fp<0.05],br_p[,2][br_fp<0.05],201)$value,3)
round(sintegral(bu_fp[bu_fp<0.05],bu_p[,2][bu_fp<0.05],201)$value,3)
round(sintegral(ecp_fp[ecp_fp<0.05],ecp_p[,2][ecp_fp<0.05],201)$value,3)



round(sintegral(LACPD_fp[LACPD_fp<0.05],lacpd_p[,3][LACPD_fp<0.05],201)$value,3)
round(sintegral(pet_fp[pet_fp<0.05],pet_p[,3][pet_fp<0.05],201)$value,3)
round(sintegral(adapt.cp_fp[adapt.cp_fp<0.05],adapt.cp_p[,3][adapt.cp_fp<0.05],201)$value,3)
round(sintegral(br_fp[br_fp<0.05],br_p[,3][br_fp<0.05],201)$value,3)
round(sintegral(bu_fp[bu_fp<0.05],bu_p[,3][bu_fp<0.05],201)$value,3)
round(sintegral(ecp_fp[ecp_fp<0.05],ecp_p[,3][ecp_fp<0.05],201)$value,3)



round(sintegral(LACPD_fp[LACPD_fp<0.05],lacpd_p[,4][LACPD_fp<0.05],201)$value,3)
round(sintegral(pet_fp[pet_fp<0.05],pet_p[,4][pet_fp<0.05],201)$value,3)
round(sintegral(adapt.cp_fp[adapt.cp_fp<0.05],adapt.cp_p[,4][adapt.cp_fp<0.05],201)$value,3)
round(sintegral(br_fp[br_fp<0.05],br_p[,4][br_fp<0.05],201)$value,3)
round(sintegral(bu_fp[bu_fp<0.05],bu_p[,4][bu_fp<0.05],201)$value,3)
round(sintegral(ecp_fp[ecp_fp<0.05],ecp_p[,4][ecp_fp<0.05],201)$value,3)


round(sintegral(LACPD_fp[LACPD_fp<0.05],lacpd_p[,5][LACPD_fp<0.05],201)$value,3)
round(sintegral(pet_fp[pet_fp<0.05],pet_p[,5][pet_fp<0.05],201)$value,3)
round(sintegral(adapt.cp_fp[adapt.cp_fp<0.05],adapt.cp_p[,5][adapt.cp_fp<0.05],201)$value,3)
round(sintegral(br_fp[br_fp<0.05],br_p[,5][br_fp<0.05],201)$value,3)
round(sintegral(bu_fp[bu_fp<0.05],bu_p[,5][bu_fp<0.05],201)$value,3)
round(sintegral(ecp_fp[ecp_fp<0.05],ecp_p[,5][ecp_fp<0.05],201)$value,3)

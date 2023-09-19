# Functions_Cell_ADPRS_ver7_3pc_09022023_1.R

#####
library(fmsb)

## ROSMAP
amyloid_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("amyloid_sqrt~",PRSs[i],"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("amyloid_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

tangles_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("tangles_sqrt~",PRSs[i],"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("tangles_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

cogdec_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("cogng_demog_slope~",PRSs[i],"+apoe4n+apoe2n+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("cogng_demog_slope~apoe4n+apoe2n+gbatch+EV01+EV02+EV03")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

ADdem_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("factor(ADdem)~",PRSs[i],"+apoe4n+apoe2n+age_death+msex+educ+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("factor(ADdem)~apoe4n+apoe2n+age_death+msex+educ+gbatch+EV01+EV02+EV03")
    m_base<-glm(as.formula(fmla2),data=subset(data),family="binomial")
    tryCatch({
      m<-glm(as.formula(fmla1),data=subset(data),family="binomial")
      m1<-summary(m)
      a<-exp(cbind(OR=coef(m),confint(m)))
      A[i,2]<-(a[2,1])
      A[i,3]<-paste0((a[2,2])," to ",(a[2,3]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-NagelkerkeR2(m)$R2
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-NagelkerkeR2(m)$R2 - NagelkerkeR2(m_base)$R2
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

dp_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("plaq_d_sqrt~",PRSs[i],"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("plaq_d_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

np_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("plaq_n_sqrt~",PRSs[i],"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("plaq_n_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

nft_PRScs_ROSMAP<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("nft_sqrt~",PRSs[i],"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    fmla2<-paste0("nft_sqrt~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

boot_R2<-function(data,trait,predictors){
  A<-matrix(NA,length(predictors),2)
  for (i in c(1:length(predictors))){
    var=trait
    prs=predictors[i]
    R2.fun <- function(data, idx){
      df <- data[idx, ]
      fmla0<-paste0(var,"~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
      fmla1<-paste0(var,"~",prs,"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
      m0<-summary(lm(as.formula(fmla0),data=subset(df)))
      m1<-summary(lm(as.formula(fmla1),data=subset(df)))
      m1$adj.r.squared - m0$adj.r.squared
    }
    A[i,1]<-gsub("_.*","",prs)
    set.seed(2023)
    A[i,2]<-sd(boot(data,R2.fun,R=1000)$t)
  }
  A<-data.frame(A)
  names(A)<-c("PRS","R2_sd")
  return(A)
}

boot_R2.cog<-function(data,predictors){
  A<-matrix(NA,length(predictors),2)
  for (i in c(1:length(predictors))){
    prs=predictors[i]
    R2.fun <- function(data, idx){
      df <- data[idx, ]
      fmla0<-paste0("cogng_demog_slope~apoe4n+apoe2n+gbatch+EV01+EV02+EV03")
      fmla1<-paste0("cogng_demog_slope~",prs,"+apoe4n+apoe2n+gbatch+EV01+EV02+EV03")
      m0<-summary(lm(as.formula(fmla0),data=subset(df)))
      m1<-summary(lm(as.formula(fmla1),data=subset(df)))
      m1$adj.r.squared - m0$adj.r.squared
    }
    A[i,1]<-gsub("_.*","",prs)
    set.seed(2023)
    A[i,2]<-sd(boot(data,R2.fun,R=1000)$t)
  }
  A<-data.frame(A)
  names(A)<-c("PRS","R2_sd")
  return(A)
}

boot_R2.AD<-function(data,predictors){
  A<-matrix(NA,length(predictors),2)
  for (i in c(1:length(predictors))){
    prs=predictors[i]
    R2.fun.AD <- function(data, idx){
      df <- data[idx, ]
      fmla0<-paste0("factor(ADdem)~apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
      fmla1<-paste0("factor(ADdem)~",prs,"+apoe4n+apoe2n+age_death+msex+gbatch+EV01+EV02+EV03")
      m0<-glm(as.formula(fmla0),data=subset(df),family="binomial")
      m1<-glm(as.formula(fmla1),data=subset(df),family="binomial")
      NagelkerkeR2(m1)$R2 - NagelkerkeR2(m0)$R2
    }
    A[i,1]<-gsub("_.*","",prs)
    set.seed(2023)
    A[i,2]<-sd(boot(data,R2.fun.AD,R=1000)$t)
  }
  A<-data.frame(A)
  names(A)<-c("PRS","R2_sd")
  return(A)
}

## A4
Ab_PRScs<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("Ab_Composite_Summary~",PRSs[i],"+apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3")
    fmla2<-paste0("Ab_Composite_Summary~apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

FTP_composite_PRScs<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("FTP_Composite_Unweighted~",PRSs[i],"+apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3")
    fmla2<-paste0("FTP_Composite_Unweighted~apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(data))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

PACC_PRScs<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("PACC_RN~",PRSs[i],"+apoe4n+apoe2n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3")
    fmla2<-paste0("PACC_RN~apoe4n+apoe2n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(a_d1))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

HV_PRScs<-function(PRS,data){
  PRSs=PRS
  A<-matrix(NA,length(PRSs),9)
  for (i in 1:length(PRSs)){
    A[i,1]<-PRSs[i]
    fmla1<-paste0("HippVol~",PRSs[i],"+apoe4n+apoe2n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3")
    fmla2<-paste0("HippVol~apoe4n+apoe2n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3")
    m0<-summary(lm(as.formula(fmla2),data=subset(data)))
    tryCatch({
      m<-lm(as.formula(fmla1),data=subset(a_d1))
      m1<-summary(m)
      A[i,2]<-(coef(m1)[2,1])
      A[i,3]<-paste0((confint(m)[2,1])," to ",(confint(m)[2,2]))
      A[i,4]<-(coef(m1)[2,3])
      A[i,5]<-(coef(m1)[2,4])
      A[i,6]<-(m1$adj.r.squared)
      A[i,7]<-(-log10(as.numeric(A[i,5])))
      A[i,8]<-df.residual(m)+length(coef(m))
      A[i,9]<-(m1$adj.r.squared-m0$adj.r.squared)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  A<-as.data.frame(A)
  names(A)<-c("PRS","Beta","95% CI","t","p","adj.R2","log10P","n","delta.adj.R2")
  A$PRS<-as.character(A$PRS)
  A$PRS<-gsub("_.*","",A$PRS)
  A$PRS<-factor(A$PRS,levels=A$PRS)
  A$log10P<-as.numeric(as.character(A$log10P))
  A$Beta<-as.numeric(as.character(A$Beta))
  A$t<-as.numeric(as.character(A$t))
  return(A)
}

boot_R2.A4<-function(data,trait,predictors){
  A<-matrix(NA,length(predictors),2)
  for (i in c(1:length(predictors))){
    var=trait
    prs=predictors[i]
    R2.fun <- function(data, idx){
      df <- data[idx, ]
      fmla0<-paste0(var,"~apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3")
      fmla1<-paste0(var,"~",prs,"+apoe4n+apoe2n+PTAGE+PTGENDER+pc1+pc2+pc3")
      m0<-summary(lm(as.formula(fmla0),data=subset(df)))
      m1<-summary(lm(as.formula(fmla1),data=subset(df)))
      m1$adj.r.squared - m0$adj.r.squared
    }
    A[i,1]<-gsub("_.*","",prs)
    set.seed(2023)
    A[i,2]<-sd(boot(data,R2.fun,R=1000)$t)
  }
  A<-data.frame(A)
  names(A)<-c("PRS","R2_sd")
  return(A)
}

boot_R2.A4.hv<-function(data,predictors){
  A<-matrix(NA,length(predictors),2)
  for (i in c(1:length(predictors))){
    prs=predictors[i]
    R2.fun <- function(data, idx){
      df <- data[idx, ]
      fmla0<-paste0("HippVol~apoe4n+apoe2n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3")
      fmla1<-paste0("HippVol~",prs,"+apoe4n+apoe2n+PTAGE+PTGENDER+ICV+pc1+pc2+pc3")
      m0<-summary(lm(as.formula(fmla0),data=subset(df)))
      m1<-summary(lm(as.formula(fmla1),data=subset(df)))
      m1$adj.r.squared - m0$adj.r.squared
    }
    A[i,1]<-gsub("_.*","",prs)
    set.seed(2023)
    A[i,2]<-sd(boot(data,R2.fun,R=1000)$t)
  }
  A<-data.frame(A)
  names(A)<-c("PRS","R2_sd")
  return(A)
}

boot_R2.A4.pacc<-function(data,predictors){
  A<-matrix(NA,length(predictors),2)
  for (i in c(1:length(predictors))){
    prs=predictors[i]
    R2.fun <- function(data, idx){
      df <- data[idx, ]
      fmla0<-paste0("PACC_RN~apoe4n+apoe2n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3")
      fmla1<-paste0("PACC_RN~",prs,"+apoe4n+apoe2n+PTAGE+PTGENDER+PTEDUCAT+pc1+pc2+pc3")
      m0<-summary(lm(as.formula(fmla0),data=subset(df)))
      m1<-summary(lm(as.formula(fmla1),data=subset(df)))
      m1$adj.r.squared - m0$adj.r.squared
    }
    A[i,1]<-gsub("_.*","",prs)
    set.seed(2023)
    A[i,2]<-sd(boot(data,R2.fun,R=1000)$t)
  }
  A<-data.frame(A)
  names(A)<-c("PRS","R2_sd")
  return(A)
}

#########

pfromt<-function(t,d){
  a<-2*pt(-abs(t),df=d)
  print(a)
}
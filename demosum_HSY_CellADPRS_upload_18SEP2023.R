## Demographics Summary for A4/ROSMAP data subsets
demosum_A4<-function(data){

  #Baseline Age (NP_bl)
  age_bl_m<-mean(data$PTAGE, na.rm=TRUE)
  age_bl_sd<-sd(data$PTAGE, na.rm=TRUE)
  age_bl_n<-sum(!is.na(data$PTAGE))
  
  #Sex: female
  femalenum<-sum(data$PTGENDER==2, na.rm=TRUE)
  femalepercent<-sum(data$PTGENDER==2, na.rm=TRUE)*100/sum(data$PTGENDER==1 | data$PTGENDER==2, na.rm=TRUE)
  female_n<-sum(!is.na(data$PTGENDER))
  
  #Education
  educ_m<-mean(data$PTEDUCAT, na.rm=TRUE)
  educ_sd<-sd(data$PTEDUCAT, na.rm=TRUE)
  educ_n<-sum(!is.na(data$PTEDUCAT))
  
  # E4 carrier rate
  E4num<-sum(data$apoe4n>0, na.rm=TRUE)
  E4percent<-sum(data$apoe4n>0, na.rm=TRUE)*100/sum(!is.na(data$apoe4n))
  E4_n<-sum(!is.na(data$apoe4n))
  
  #Baseline Ab
  Ab_m<-mean(data$Ab_Composite_Summary, na.rm=TRUE)
  Ab_sd<-sd(data$Ab_Composite_Summary, na.rm=TRUE)
  Ab_n<-sum(!is.na(data$Ab_Composite_Summary))
  
  #Baseline MMSE
  mmse_m<-median(data$MMSE_RN_RW, na.rm=TRUE)
  mmse_sd<-IQR(data$MMSE_RN_RW)
  mmse_n<-sum(!is.na(data$MMSE_RN_RW))
  
  # age_bl_m, femalenum, whitenum, educ_m, numnpobs_m, pib_m, hv_m, wmh_m, pacc_m
  
  #output
  x1<-c("Baseline Age","Female (%)","Education(yrs)",
       "E4 carrier (%)","Baseline Ab","Baseline MMSE")
  x2<-c(age_bl_m, femalenum, educ_m, E4num, Ab_m, mmse_m)
  x3<-c(age_bl_sd, femalepercent, educ_sd, E4percent, Ab_sd, mmse_sd)
  x4<-c(age_bl_n, female_n, educ_n, E4_n, Ab_n, mmse_n)
  y<-data.frame("Variables"=x1,"Mean/Number"=x2,"SD/Percent"=x3,"n"=x4)
  return(y)
  
}

demosum_ROSMAP<-function(data){
  
  #Age at death
  age_bl_m<-mean(data$age_death, na.rm=TRUE)
  age_bl_sd<-sd(data$age_death, na.rm=TRUE)
  age_bl_n<-sum(!is.na(data$age_death))
  
  #Sex: female
  femalenum<-sum(data$msex==0, na.rm=TRUE)
  femalepercent<-sum(data$msex==0, na.rm=TRUE)*100/sum(!is.na(data$msex))
  female_n<-sum(!is.na(data$msex))
  
  #Education
  educ_m<-mean(data$educ, na.rm=TRUE)
  educ_sd<-sd(data$educ, na.rm=TRUE)
  educ_n<-sum(!is.na(data$educ))
  
  # E4 carrier rate
  E4num<-sum(data$apoe4n>0, na.rm=TRUE)
  E4percent<-sum(data$apoe4n>0, na.rm=TRUE)*100/sum(!is.na(data$apoe4n))
  E4_n<-sum(!is.na(data$apoe4n))
  
  # pathoAD
  pathoADnum<-sum(data$pathoAD==1,na.rm=TRUE)
  pathoADpercent<-sum(data$pathoAD==1,na.rm=TRUE)/sum(!is.na(data$pathoAD))
  pathoAD_n<-sum(!is.na(data$pathoAD))
  
  #Last MMSE
  mmse_m<-median(data$cts_mmse30_lv, na.rm=TRUE)
  mmse_sd<-IQR(data$cts_mmse30_lv, na.rm=TRUE)
  mmse_n<-sum(!is.na(data$cts_mmse30_lv))

  #output
  x1<-c("Age at Death","Female (%)","Education(yrs)",
        "E4 carrier (%)","pathoAD","Last MMSE")
  x2<-c(age_bl_m, femalenum, educ_m, E4num, pathoADnum, mmse_m)
  x3<-c(age_bl_sd, femalepercent, educ_sd, E4percent, pathoADpercent, mmse_sd)
  x4<-c(age_bl_n, female_n, educ_n, E4_n, pathoAD_n, mmse_n)
  y<-data.frame("Variables"=x1,"Mean/Number"=x2,"SD/Percent"=x3,"n"=x4)
  return(y)
  
}



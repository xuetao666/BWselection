
rm(list=ls(all=TRUE))  #same to clear all in stata
cat("\014")
##Import the NHANES data for testing:
library(NHANES)
library(dplyr)
library(lspline)
library(stringi)
data("NHANES")

options(warn=0)

##Stepwise selection:
# data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
# spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",s

backward=function(data,qvarlist,lvarlist,spvarlist,spclist,catvarlist,
                  outcome,type="lm",sig=0.05,complete_case=FALSE){
  
  ##Macros for test, delete later:   
  # data=NHANES
  # 
  # qvarlist<-c("BPSysAve","SleepHrsNight","TotChol")
  # lvarlist<-c("Age","BPDiaAve","Weight","Height")
  # spvarlist<-c("Pulse","DirectChol")
  # spclist<-c(70,1.2)
  # catvarlist<-c("Depressed","Marijuana","Gender")
  # 
  # type="lm"
  # outcome="BMI"
  # sinkfile="test1"
  # complete_case=TRUE
  # sig=0.05
  
  data=data[,c(qvarlist,lvarlist,spvarlist,catvarlist,outcome)]
  
  ##UPdate quadratic terms:
  for(var in qvarlist){
    data[,paste0(var,"_quar2")]=data[,var]^2
  }
  
  ##Test sig, delete later
  ##summary(lm(BMI~Age+BPDiaAve+DirectChol+BPSysAve_2,data=data))
  
  
  varlist=c()
  qvarlist2<-c()
  for(var in qvarlist){
    varlist=c(varlist,var,paste0(var,"_quar2"))
    qvarlist2=c(qvarlist2,var,paste0(var,"_quar2"))
  }
  
  
  for(var in c(lvarlist)){
    varlist=c(varlist,var)
  }
  
  
  catvarlist2=c()
  for(var in c(catvarlist)){
    varlist=c(varlist,paste0("as.factor(",var,")"))
    catvarlist2=c(catvarlist2,paste0("as.factor(",var,")"))
  }
  
  
  spvarlist2=c()
  for(var in spvarlist){
    p=which(spvarlist %in% var)
    np=spclist[p]
    varlist=c(varlist,paste0("lspline(",var,",",np,")"))
    spvarlist2=c(spvarlist2,paste0("lspline(",var,",",np,")"))
  }
  
  
  #Remove missings in the full model
  if(complete_case==TRUE){
    data=data[complete.cases(data),]
  } 
  
  
  #Unweighted Selection:#
  ################################################################################
  ################################################################################
  
  
  #sink(file=sinkfile,append = FALSE)
  
  max_p=1
  itest=1
  warn=FALSE
  full_model=NULL
  aic_pre=Inf
  removelist=NULL
  
  while(max_p>sig){
    
    if(length(full_model)!=0){ ##Get the previous model
      model_pre=full_model 
      fit_pre=fit
      aic_pre=aic_fit
    }
    
    full_model<-paste(varlist,"+",collapse=" ")
    full_model<-substr(full_model,1,nchar(full_model)-1)
    if(type=="lm"){ #For linear regression:
      full_model<-paste0(outcome,"~",full_model)
      fit<-lm(as.formula(full_model),data=data)
      aic_fit=extractAIC(fit)[2]
    }
    if(itest==1){
      fit_or=fit
      full_model_or=full_model
    }
    if(warn==FALSE){
      if(itest!=1){
        fit_row_n=dim(fit[["model"]])[1]
        if(fit_row_n!=fit_row){
          warning("Different rows in different model, result might not reliable, could use(complete_case=TRUE)")
          warn=TRUE
        }
      }
      fit_row=dim(fit[["model"]])[1]
    }
    
    ##Compare the AICs: If later model have higher AICs, we use the previous model and stop the selection process
    if(aic_fit>aic_pre){  
      fit=fit_pre
      full_model=model_pre
      break
    }
    
    # } else if(type=="surv"){  ##Fix later:
    #   full_model<-paste0("Surv(",timevar,", ",deathvar,") ~ ",full_model)
    #   fit<-coxph(as.formula(full_model),data=data)
    # }
    if(length(removelist)!=0){
      print(paste0("Remove Var:",removelist,collapse =","))
    }
    
    print("----------------------------------------------------------------------")
    print(paste0("Test:",itest))
    print(full_model)
    print(paste0("AIC: ",extractAIC(fit)[2]))
    
    ##Overall fit:
    summary(fit)
    
    #get pvalues
    print(coef(summary(fit)))
    if(type=="lm"){
      coef_p=coef(summary(fit))[,"Pr(>|t|)"]
    }
    #coef_p=coef(summary(fit))[,"Pr(>|z|)"]
    
    remove=FALSE
    while(remove==FALSE){
      max_p=max(coef_p)
      
      names_p=names(coef_p)
      names_p=str_remove(names_p," ")
      names(coef_p)=names_p
      max_p_name=names_p[coef_p==max(coef_p)] ##For detecting factors in the formula names
      
      ##remove the numbers and the category indicator of the varnames
      cat_index=gregexpr(")",max_p_name)[[1]][1]
      if(cat_index!=-1){
        max_p_name_c=substr(max_p_name,1,cat_index)
      } else {
        max_p_name_c=max_p_name
      }
      
      
      ##remove the factor from the regression list:
      
      if(max_p>sig){
        ##For category variables, didn't remove if any category is significant:
        if(sum(c(catvarlist2,qvarlist2,spvarlist2) %in% max_p_name_c)>0 & !stri_detect_fixed(max_p_name_c,"_quar2") ){ ##Variable in the category list, but not qvarlist second term
          #Get p values for all category:
          values_p=coef_p[stri_detect_fixed(names_p,max_p_name_c)]
          min_list_p=min(values_p)
          
          if(min_list_p<=sig){ #If the values are less or equals to sig, keep this one and found the next:
            coef_p=coef_p[!stri_detect_fixed(names_p,max_p_name)]
            
            remove=FALSE
            
          } else if(min_list_p>sig){ #If the minimum value > sig, remove all the categories
            max_p_name_c=str_remove(max_p_name_c,"_quar2")
            ##remove from the varlist:
            removelist=varlist[stri_detect_fixed(varlist,max_p_name_c)]
            varlist=varlist[!stri_detect_fixed(varlist,max_p_name_c)]
            
            remove=TRUE
          } ## end of min_lis_p
        } else { ## Not in any chained removed list, remove directly
          removelist=varlist[stri_detect_fixed(varlist,max_p_name_c)]
          varlist=varlist[!stri_detect_fixed(varlist,max_p_name_c)]
          
          remove=TRUE
        }
        
      } else if(max_p<=sig){ ## end close for max_p>sig
        ##Break the loop
        break
      }
      
    }
    print("----------------------------------------------------------------------")
    itest=itest+1
  }
  print(summary(fit))
  re_list=list()
  re_list[["fit"]]=fit
  re_list[["Orignal_formula"]]=full_model_or
  re_list[["Final_formula"]]=full_model
  return(re_list)
}


data=NHANES
# 
qvarlist<-c("BPSysAve","SleepHrsNight","TotChol")
lvarlist<-c("Age","BPDiaAve","Weight","Height")
spvarlist<-c("Pulse","DirectChol")
spclist<-c(70,1.2)
catvarlist<-c("Depressed","Marijuana","Gender")

fit_test1=backward(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
          spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",sig=0.05,complete_case=TRUE)
summary(fit_test1$fit)


test_data=NHANES[,c(qvarlist,lvarlist,spvarlist,catvarlist,"BMI")]
test_data=test_data[complete.cases(test_data),]

or_fit=lm(BMI~BPSysAve + I(BPSysAve^2) + SleepHrsNight + 
            I(SleepHrsNight^2) + TotChol + I(TotChol^2) + 
            Age + BPDiaAve + Weight + Height + as.factor(Depressed) + 
            as.factor(Marijuana) + as.factor(Gender) + lspline(Pulse,70) + 
            lspline(DirectChol,1.2),data=test_data)

fit_2=step(or_fit, direction = "backward", trace=TRUE, test = "F" ) 
sink()

##Test with stepwise selection:

fit_test=fit_or

fit_2=step(fit_test, direction = "backward", trace=TRUE, test = "F" ) 
summary(fit_2)
summary(fit)
summary(fit_or)


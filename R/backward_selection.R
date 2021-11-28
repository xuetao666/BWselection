#'BWselection
#'
#'Perform Backward selection on regression models
#'
#'@param data input Data used for the regression
#'@param qvarlist Variable list in vector format, including all the variables needed for the quadratic form
#'@param lvarlist Variable list in vector format, including all the variables in the continous linear form
#'@param spvarlist variable list in vector format, including all the variables for the spline form
#'@param spclist if spvarlist is not empty, specify the cutpoints for the spline term
#'@param catvarlist variable list in the vector format, including all the categorial variables
#'@param outcome Outcome variable for the regression model
#'@param type type of regression model, default is lm
#'@param sig significant level used for exit the model, default is 0.05
#'@param complete_case whether to use the complete case for all the variables in the model, defualt is FALSE
#'
#'
#'
#'@return a list of original formula, final fit model and final fit formula
#'
#'@examples
#'#Load NHANES DATA
#'library(NHANES)
#'data("NHANES")
#'data=NHANES
#'
#'#Created variable lists for the model selection
#'qvarlist<-c("BPSysAve","SleepHrsNight","TotChol")
#'lvarlist<-c("Age","BPDiaAve","Weight","Height")
#'spvarlist<-c("Pulse","DirectChol")
#'spclist<-c(70,1.2)
#'catvarlist<-c("Depressed","Marijuana","Gender")
#'
#'#Run the model selection, this will return a list
#'fit_test1=BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
#'                      spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",sig=0.05,complete_case=TRUE)
#'
#'
#'#In order to compare the result, we use the step function build in R.
#'#However, the function didn't work with the Incomplete case data, Thus, we need to clean it first and fit the model by hand
#'
#'#Clean the data:
#'test_data=NHANES[,c(qvarlist,lvarlist,spvarlist,catvarlist,"BMI")]
#'test_data=test_data[complete.cases(test_data),]
#'
#'#Fit the same original model by hand:(Note that here, we create quandartic forms as I(BPSysAve^2), same as result of BPSysAve_quar2)
#'or_fit=lm(BMI~BPSysAve + I(BPSysAve^2) + SleepHrsNight +
#'            I(SleepHrsNight^2) + TotChol + I(TotChol^2) +
#'            Age + BPDiaAve + Weight + Height + as.factor(Depressed) +
#'            as.factor(Marijuana) + as.factor(Gender) + lspline(Pulse,70) +
#'            lspline(DirectChol,1.2),data=test_data)
#'
#'#Final model fit
#'fit_2=step(or_fit, direction = "backward", trace=TRUE, test = "F" )
#'
#'##Compare the two results:
#'summary(fit_2)
#'summary(fit_test1$fit)
#'
#'@export
#'

BWselection<-function(data,qvarlist=NULL,lvarlist=NULL,spvarlist=NULL,spclist=NULL,catvarlist=NULL,
                  outcome,type="lm",sig=0.05,complete_case=FALSE){
  library(dplyr)
  library(lspline)
  library(stringi)
  library(stringr)

  # ##Macros for test, delete later:
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

  ##Warning1: If no covariates selected:
  if(length(qvarlist)+length(lvarlist)+length(spvarlist)+length(catvarlist)==0) stop("No covariates inputed.")
  ##Warning2: If no spclist input but have spvarlist input:
  if(length(spvarlist)!=0 & length(spclist)==0) stop("No split cutpoint input")

  ##Generate overall varlist:
  varlist=c()
  #Quandratic terms:
  qvarlist2<-c()
  for(var in qvarlist){
    varlist=c(varlist,var,paste0("I(",var,"^2)"))
    qvarlist2=c(qvarlist2,var,paste0("I(",var,"^2)"))
  }
  #Continous terms:
  for(var in c(lvarlist)){
    varlist=c(varlist,var)
  }
  #Categorical terms:
  catvarlist2=c()
  for(var in c(catvarlist)){
    varlist=c(varlist,paste0("as.factor(",var,")"))
    catvarlist2=c(catvarlist2,paste0("as.factor(",var,")"))
  }
  #Split terms:
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


  #Selection begain#
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

    full_model=paste(varlist,"+",collapse=" ")
    full_model=substr(full_model,1,nchar(full_model)-1)
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
        if(sum(c(catvarlist2,qvarlist2,spvarlist2) %in% max_p_name_c)>0 & !stri_detect_fixed(max_p_name_c,"^2") ){ ##Variable in the category list, but not qvarlist second term
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
#sink()
#

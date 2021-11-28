#'BackWard Selection
#'
#'Perform Backward selection on regression models
#'
#'@param data input Data used for the regression
#'@param qvarlist Variable list in vector format, including all the variables needed for the quadratic form
#'@param lvarlist Variable list in vector format, including all the variables in the continuous linear form
#'@param spvarlist variable list in vector format, including all the variables for the spline form
#'@param spclist if spvarlist is not empty, specify the outpoints for the spline term
#'@param catvarlist variable list in the vector format, including all the categorical variables
#'@param outcome Outcome variable for the regression model
#'@param type type of regression model, default is lm, logistic regression use type="logistic"
#'@param sig significant level used for exit the model, default is 0.05
#'@param complete_case whether to use the complete case for all the variables in the model, default is TRUE
#'@param track Whether print the selection process, could select from "Simple","All", "None"
#'
#'
#'@return a list of original formula, final fit model and final fit formula
#'
#'@examples
#'
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
#'@export
#'

BWselection<-function(data,qvarlist=NULL,lvarlist=NULL,spvarlist=NULL,spclist=NULL,catvarlist=NULL,
                  outcome,type="lm",sig=0.05,complete_case=TRUE,track="Simple"){
  library(dplyr)
  library(lspline)
  library(stringi)
  library(stringr)

  data=data[,c(qvarlist,lvarlist,spvarlist,catvarlist,outcome)] #Remove missings in the full model
  if(complete_case==TRUE){
    data=data[complete.cases(data),]
  }
  #------Error Messages--------------------:
  #-----------------------------------------------------------------------------
  if(length(qvarlist)+length(lvarlist)+length(spvarlist)+length(catvarlist)==0) stop("No covariates inputed.") #Error1: If no covariates selected:
  if(length(spvarlist)!=0 & length(spclist)==0) stop("No split cutpoint input")   #Error2: If no spclist input but have spvarlist input:
  if(track!="Simple" & track!="None" & track!="All") stop("track type not recognized")   #Error3: Type of track not recognized:
  if(type!="lm" & type!="logistic") stop("type of regression not recognized")   #Error4: Type of regression not recognized:
  #------Generate overall varlist----------:
  #-----------------------------------------------------------------------------
  varlist=c()
  qvarlist2<-c()   #Quandratic terms:
  for(var in qvarlist){
    varlist=c(varlist,var,paste0("I(",var,"^2)"))
    qvarlist2=c(qvarlist2,var,paste0("I(",var,"^2)"))
  }
  for(var in c(lvarlist)){   #Continous terms:
    varlist=c(varlist,var)
  }
  catvarlist2=c()   #Categorical terms:
  for(var in c(catvarlist)){
    varlist=c(varlist,paste0("as.factor(",var,")"))
    catvarlist2=c(catvarlist2,paste0("as.factor(",var,")"))
  }
  spvarlist2=c() #Split terms:
  for(var in spvarlist){
    p=which(spvarlist %in% var)
    np=spclist[p]
    varlist=c(varlist,paste0("lspline(",var,",",np,")"))
    spvarlist2=c(spvarlist2,paste0("lspline(",var,",",np,")"))
  }
  #------Selection began---------------:
  #-----------------------------------------------------------------------------
  max_p=1
  itest=1
  warn=FALSE
  full_model=NULL
  aic_pre=Inf
  removelist=NULL
  while(max_p>sig){ #If there is non-significant terms in the model:
    if(length(full_model)!=0){ ##Get the previous model
      model_pre=full_model
      fit_pre=fit
      aic_pre=aic_fit
    }
    full_model=paste(varlist,"+",collapse=" ") #Define full model formula
    full_model=substr(full_model,1,nchar(full_model)-1)
    if(type=="lm"){ #For linear regression:
      full_model<-paste0(outcome,"~",full_model)
      fit<-lm(as.formula(full_model),data=data)
    } else if(type=="logistic"){
      full_model<-paste0(outcome,"~",full_model)
      fit<-glm(as.formula(full_model),data=data,family = "binomial")
    }
    aic_fit=extractAIC(fit)[2]
    if(itest==1){ #Save the original model in advance
      fit_or=fit
      full_model_or=full_model
    }
    if(warn==FALSE){ #Check if there the observations in the model is different
      if(itest!=1){
        fit_row_n=dim(fit[["model"]])[1]
        if(fit_row_n!=fit_row){
          warning("Different rows in different model, result might not reliable, could use(complete_case=TRUE)")
          warn=TRUE
        }
      }
      fit_row=dim(fit[["model"]])[1]
    }
    if(aic_fit>aic_pre){ #Compare the AICs: If later model have higher AICs, we use the previous model and stop the selection process
      fit=fit_pre
      full_model=model_pre
      break
    }
    if(track!="None"){ #Print selection process
      if(length(removelist)!=0){ #Output remove list and other information
        print(paste0("Remove Var:",removelist,collapse =","))
      }
      print("----------------------------------------------------------------------")
      print(paste0("Test:",itest))
      print(full_model)
      print(paste0("AIC: ",extractAIC(fit)[2]))
      if(track=="All"){
        print(summary(fit))  #Overall fit:
      }
      print("----------------------------------------------------------------------")
    }
    coef_p=coef(summary(fit))[,ncol(coef(summary(fit)))] #Get coefficients
    remove=FALSE
    while(remove==FALSE){ #Loop until at least one term is removed from the selection process
      max_p=max(coef_p) #Get the terms of max p value
      names_p=names(coef_p)
      names_p=str_remove(names_p," ")
      names(coef_p)=names_p
      max_p_name=names_p[coef_p==max(coef_p)] #For detecting factors in the formula name
      cat_index=gregexpr(")",max_p_name)[[1]][1] #remove the numbers and the category indicator of the varnames
      if(cat_index!=-1){
        max_p_name_c=substr(max_p_name,1,cat_index)
      } else {
        max_p_name_c=max_p_name
      }
      if(max_p>sig){ #remove the factor from the regression list:
        #For category variables, didn't remove if any category is significant:
        if(sum(c(catvarlist2,qvarlist2,spvarlist2) %in% max_p_name_c)>0 & !stri_detect_fixed(max_p_name_c,"^2") ){ #Variable in the category list, but not qvarlist second term
          values_p=coef_p[stri_detect_fixed(names_p,max_p_name_c)] #Get p values for all category:
          min_list_p=min(values_p)
          if(min_list_p<=sig){ #If the values are less or equals to sig, keep this one and found the next:
            coef_p=coef_p[!stri_detect_fixed(names_p,max_p_name)]
            remove=FALSE
          } else if(min_list_p>sig){ #If the minimum value > sig, remove all the categories
            max_p_name_c=str_remove(max_p_name_c,"_quar2")
            removelist=varlist[stri_detect_fixed(varlist,max_p_name_c)] #remove from the varlist:
            varlist=varlist[!stri_detect_fixed(varlist,max_p_name_c)]
            remove=TRUE
          } # end of min_lis_p
        } else { # Not in any chained removed list, remove directly
          removelist=varlist[stri_detect_fixed(varlist,max_p_name_c)]
          varlist=varlist[!stri_detect_fixed(varlist,max_p_name_c)]
          remove=TRUE
        }
      } else if(max_p<=sig){ # end close for max_p>sig
        break #Break the loop
      } #End for max_p>sig
    } #End for loop term of remove selection
    itest=itest+1
  }
  print("Final Model")
  print(summary(fit))
  re_list=list()
  re_list[["fit"]]=fit
  re_list[["Orignal_formula"]]=full_model_or
  re_list[["Final_formula"]]=full_model
  return(re_list)
}
#sink()
#

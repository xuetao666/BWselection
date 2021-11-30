test_that("BWselection works", {
  library(NHANES)
  data("NHANES")

  #Created variable lists for the model selection
  qvarlist<-c("BPSysAve","SleepHrsNight","TotChol")
  lvarlist<-c("Age","BPDiaAve","Weight","Height")
  spvarlist<-c("Pulse","DirectChol")
  spclist<-c(70,1.2)
  catvarlist<-c("Depressed","Marijuana","Gender")

  #Start testing:
  #1.Test if the final selection correct
  expect_equal(BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                           spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",
                           sig=0.05,complete_case=TRUE)$Final_formula,
               "BMI~BPSysAve + I(BPSysAve^2) + SleepHrsNight + I(SleepHrsNight^2) + TotChol + I(TotChol^2) + Age + BPDiaAve + Weight + Height + lspline(DirectChol,1.2) ")
  expect_equal(BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                           spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",
                           sig=0.05,complete_case=TRUE,track="All")$Final_formula,
               "BMI~BPSysAve + I(BPSysAve^2) + SleepHrsNight + I(SleepHrsNight^2) + TotChol + I(TotChol^2) + Age + BPDiaAve + Weight + Height + lspline(DirectChol,1.2) ")
  expect_equal(BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                           spclist=spclist,outcome="Depressed",type="logistic",
                           sig=0.05,complete_case=TRUE)$Final_formula,
               "Depressed~BPSysAve + I(BPSysAve^2) + SleepHrsNight + I(SleepHrsNight^2) + TotChol + I(TotChol^2) + Age + BPDiaAve + Weight + Height + lspline(Pulse,70) + lspline(DirectChol,1.2) " )
  #2.Test if error message correctlly outputed:
  expect_error(BWselection(data=NHANES,outcome="BMI",type="lm",
                             sig=0.05,complete_case=TRUE),"No covariates inputed.")
  expect_error(BWselection(data=NHANES,outcome="BMI",type="lm",spvarlist=spvarlist,lvarlist=lvarlist,
                           sig=0.05,complete_case=TRUE),"No split cutpoint input")
  expect_error(BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                           spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",
                           sig=0.05,track="Simplel",complete_case=TRUE),"track type not recognized")
  expect_error(BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                          spclist=spclist,catvarlist=catvarlist,outcome="BMI",
                          sig=0.05,type="Survival",complete_case=TRUE),"type of regression not recognized")
  #3. Warnings:
  expect_warning(BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,spvarlist=spvarlist,
                             spclist=spclist,catvarlist=catvarlist,outcome="BMI",type="lm",
                             sig=0.05,complete_case=FALSE))

  #4. Test other branches:

  expect_equal(BWselection(data=NHANES,lvarlist=c("Age","BPDiaAve","Weight","Height","BPSysAve"),outcome="BMI",type="lm",
                           sig=0.05,complete_case=TRUE)$Final_formula,"BMI~Age + BPDiaAve + Weight + Height ")

  expect_equal(BWselection(data=NHANES,catvarlist = c("Depressed","Marijuana","Gender","AgeDecade","HomeRooms"),
                           outcome="BMI",type="lm",track="All")$Final_formula,
               "BMI~as.factor(Depressed) + as.factor(Marijuana) + as.factor(AgeDecade) + as.factor(HomeRooms) ")

  qvarlist<-c("BPSysAve","SleepHrsNight","TotChol")
  lvarlist<-c("Age","Weight","Height")
  BWselection(data=NHANES,qvarlist=qvarlist,lvarlist=lvarlist,outcome="BMI",type="lm",
              sig=0.05,complete_case=FALSE)

})

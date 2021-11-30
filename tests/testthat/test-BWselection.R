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
  #2.Test if error message correctlly outputed:
  expect_error(BWselection(data=NHANES,outcome="BMI",type="lm",
                             sig=0.05,complete_case=TRUE),"No covariates inputed.")

})

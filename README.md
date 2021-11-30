  <!-- badges: start -->
  [![R-CMD-check](https://github.com/xuetao666/BWselection/workflows/R-CMD-check/badge.svg)](https://github.com/xuetao666/BWselection/actions)
  <!-- badges: end -->

# BWslection
## Overview:
BWselection is a 


This is a backward selection function based on AIC and F-test P-values. Current version could handle linear regression and logistic regression. 
Compared to the "step" function in base R, this version of backward selection have the following improvements:
1. In the long-term list of variables needed to be selected from, this version only need to input the variable list instead of the formalized formula. For example, in the step version, one need to fit the model first, then do the selection, in the model fit step, formulas like "A~B+C+D+E+..." needed to be writened down carefully, including the split terms. If one wants to change the orginal model, it is hard to modifiy and could higher chances to make mistakes. However, in our model, you only need to specify the variable type and list in the function, the function will generate the formula for you.
2. This version handles quandratic terms better. The "step" function might remove the linear term of quandratic form and leave the quandratic term in which is not correct

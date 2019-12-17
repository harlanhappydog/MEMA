######

# devtools::install_github("quentingronau/bridgesampling@master")
# install.packages("rstan")

library(bridgesampling)
library(rstan)

library(RCurl)

script <- getURL("https://raw.githubusercontent.com/harlanhappydog/MEMA/master/MEMA_functions.R", ssl.verifypeer = FALSE)

eval(parse(text = script))


library(httr)

schooldata <- read.csv("https://github.com/harlanhappydog/MEMA/raw/master/13schools.csv")

head(schooldata)


## conduct meta-analysis of simple linear regressions:

simple_lm <- simple_linear_regression_MA(y = schooldata$sci, X = cbind(schooldata$rdg), studyID = schooldata$sch, mcmc_iter=3000)

simple_lm

## conduct meta-analysis of simple linear regressions (with multivariable framework):

schooldata_simple <- linear_regression_MA(y = schooldata$sci, X = cbind(schooldata$rdg), studyID = schooldata$sch, mcmc_iter=3000)  

schooldata_simple            

## conduct meta-analysis of multivariable linear regressions:

schooldata_multi <- linear_regression_MA(y = schooldata$sci, X = cbind(schooldata$rdg, schooldata$math), studyID = schooldata$sch, mcmc_iter=3000)

schooldata_multi           





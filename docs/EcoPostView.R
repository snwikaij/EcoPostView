## ----include = FALSE----------------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
#install.packages("devtools")
library(devtools)
#install.github("snwikaij/EcoPostView")
library(EcoPostView)

## ----data---------------------------------------------------------------------
data(example1)
head(example1)


## ----meta1--------------------------------------------------------------------
mod1 <- meta(estimate=example1$est,         #Model estimate
             stderr=example1$se,            #Standard error of the model estimate
             parameter=example1$parameter,  #Model parameter (b0 or b1)
             predictor=example1$predictor,  #Predictor variable  (independent variable)
             link_function=example1$link,   #Link function
             grouping=example1$group,       #Group
             Nsamp=example1$n,              #Sample size
             method=2)                      #Adjustment method (0=none, 1=Egger's (1/se), 2=Peters (1/n))

#remove model to keep the environment clean
rm(mod1)


## ----meta1 warning------------------------------------------------------------
#mod1$model$JAGS_model

## ----meta2--------------------------------------------------------------------
mod2 <-meta(estimate=example1$est,        
            stderr=example1$se,            
            parameter=example1$parameter,  
            predictor=example1$predictor,  
            link_function=example1$link,   
            grouping=example1$group,       
            Nsamp=example1$n,              
            method=2,
            n_chain = 6)                   #Increasing the number of chains from 2 to 6   

#remove model to keep the environment clean
rm(mod2)


## ----meta2 example priors-----------------------------------------------------
only_priors <- meta(estimate=example1$est,        
            stderr=example1$se,            
            parameter=example1$parameter,  
            predictor=example1$predictor,  
            link_function=example1$link,   
            grouping=example1$group,       
            Nsamp=example1$n,            
            method=2,
            n_chain = 6,
            get_prior_only=TRUE) #Only show the structure of the priors

print(only_priors)

#remove data frame of priors to keep environment clean
rm(only_priors)



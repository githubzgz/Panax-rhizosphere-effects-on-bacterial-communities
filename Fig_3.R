#Analysis of beta-diversity

#1.Extraction of within-pair community dissimilarity
diff_extract <- function(dis){
  #dis -- Bray-Curtis dissimilarity (object:dist)
  require(stringr)
  require(tidyr)
  require(dplyr)
  dis_data <- dis %>% as.matrix %>% data.frame
  result <- data.frame(sample=str_sub(names(dis_data)[1:81],end = -4),diff=NA)
  dis_data <- dis_data[1:81,82:162]
  for(i in 1:nrow(result)){
    result[i,2] = dis_data[rownames(dis_data)%in%str_subset(rownames(dis_data),result[i,1]),
                           names(dis_data)%in%str_subset(names(dis_data),result[i,1])]
  }
  result <- result %>% mutate(.,name = .$sample) %>% separate(.,col = 'name',into = c('plant','year','replicate'),sep = '_')
  result$year_num <- NA
  for(i in 1:nrow(result)){
    if(result[i,'year']=='two') result[i,'year_num']=2
    if(result[i,'year']=='three') result[i,'year_num']=3
    if(result[i,'year']=='four') result[i,'year_num']=4
  }
  return(result)
}

#2.Multiple variates regression and variance decomposition
library(relaimpo)
library(leaps)
#Using the PG community model as an example
gin_bray_lm <- regsubsets(diff~.,data=gin_differ_data,nvmax = 10)
plot(gin_bray_lm) #Examine the R2 distribution
gin_bray_lm_summary <- gin_bray_lm %>% summary

#Select the optimal model on the basis of adjusted R2 and BIC:
gin_bray_lm_summary$adjr2 
gin_bray_lm_summary$bic

#Fitting the model using the selected variables:
gin_bray_lm <- lm(diff~AK+NIN+Mg,data=gin_differ_data)

#Variance decomposition
calc.relimp(gin_bray_lm)
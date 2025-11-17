#### Greetings, 
###This script was made to be open-access
### Our data is in open-access too (available online --> OSF)
### For any information concerning the data / Data analysis, please contact the authors
### Any feedback is always appreciated.


library("specificity")
library("car")
library("psych")
library("stats")


#Computing a function to detect observations associated to rsd > 3 (pre-registered criteria of outliers exclusion) [Many thanks to François Ric for providing this function]
outliers <- function (formule,DF,order="sdr") {
  monModel <- formule
  fit <- lm(monModel,DF)
  
  DF<-model.frame(monModel,DF)
  
  DF$sdr<-round(rstudent(fit),2)
  DF$fstar<-round(DF$sdr^2,2)
  DF$cookd<-round(cooks.distance(fit),2)
  DF$leverage<-round(hatvalues(fit),2)
  
  if(order=="sdr") {
    DF[order(-DF$fstar),]
  } else {
    if(order=="cookd") {
      DF[order(-DF$cookd),]
    } else {
      DF[order(-DF$leverage),]
      
    }}
}

#Opening data (depends on your personnal directory)
MirrorDataMCS2019 <- read.csv("~/Mirror effect Octobre 2018/MirrorDataMCS2019.csv")
#Renaming data for an easier use
df<-MirrorDataMCS2019

#Pre-registered data exclusion:
#One non-native french speaker :
df<-df[-20,]

#Coding factors as factors
df$Mirror<-as.factor(df$Mirror)
df$discrepancies<-as.factor(df$discrepancies)
df$genre<-as.factor(df$genre)




###########################PRE REGISTERED ANALYSES ################################




#testing residuals normality:
shapiro.test(resid(lm(suicidewords~Mirror, data = df))) #Not normal

# With Log transformation
shapiro.test(resid(lm(log(suicidewords)~Mirror, data = df))) #Normal

#Variances:
leveneTest(log(df$suicidewords), df$Mirror)

###############################Regression with NEUTRAL RT
#Log-tranforming data
df$suicidewords<-log(df$suicidewords)

#Model called "modneu" for model neutral:
modneu<-lm(df$suicidewords~df$Mirror*df$discrepancies+df$neutralwords)

#Cleaning the DF
outliers(modneu, df) #Checking studentized residuals
dfclean<-df[-69,] #Removing line 69 associated to studentized residual >3 in a new data set called "dfclean" for clean dataframe
modneu<-lm(dfclean$suicidewords~dfclean$Mirror*dfclean$discrepancies+dfclean$neutralwords)
shapiro.test(resid(modneu)) #Normal distribution

#Obtaining results
summary(modneu)

####### WARNING: The floowing lines will provide confidence intervals for the log transformed DV. To obtain the reported untrosformed DV, run again the script from the beginning except line 67 (df$suicidewords<-log(df$suicidewords))
#To obtain 95% one-sided confidence interval, we take the 90% confidence interval and only keep the lower bound while replacing the upperbound with infinity 
confint(modneu, level = 0.90)
#To obtain the non-direction effects 95% CI:
confint(modneu)
#Obtaining partial eta squared:
eta2(modsui)


##################Regression controlling for NEGATIVE RT 
#Same as before
modneg<-lm(df$suicidewords~df$Mirror*df$discrepancies+df$negativewords)
outliers(modneg, df)
#On peut enlever un participant
dfclean<-df[-69,]
modneg<-lm(suicidewords~Mirror*discrepancies+negativewords, data= dfclean)
summary(modneg)
confint(modneg, level = 0.90)
confint(modneg)
eta2(modneg)


###########Exploratory analyses ############
####It is necessary to re-start from the beginning because 
#Renaming data for an easier use
df<-MirrorDataMCS2019

#Pre-registered data exclusion:
#One non-native french speaker :
df<-df[-20,]


#Coding factors as factors
df$Mirror<-as.factor(df$Mirror)
df$discrepancies<-as.factor(df$discrepancies)
df$genre<-as.factor(df$genre)


#Computing datasets excluding outliers using the MAD2 method

###################Suicide words RT
MADsuicide<-2*mad(df$suicidewords) #calculating 2 MAD
MADmaxsuicide<-median(df$suicidewords)+MADsuicide #Computing the upper bound of 2 MAD from the median
MADminsuicide<-median(df$suicidewords)-MADsuicide #Lower bound


###################Neutral words RT
MADneutre<-2*mad(df$neutralwords) 
MADmaxneu<-median(df$neutralwords)+MADneutre
MADminneu<-median(df$neutralwords)-MADneutre

#Computing data set using MAD2 method: Making a dataframe called "madneutral" corresponding to the dataframe "df" but only with the valus of suicide RT and neutral RT between the upper and lower bound
madneutral<-df[-(which(df$suicidewords < MADminsuicide | df$suicidewords > MADmaxsuicide |
                           df$neutralwords < MADminneu | df$neutralwords > MADmaxneu)),]

#testing hypotheses:
#H1
madmodneu<-lm(madneutral$suicidewords~madneutral$Mirror*madneutral$discrepancies+madneutral$neutralwords)
shapiro.test(resid(madmodneu))
summary(madmodneu)
#Failing to replicate 1st hypothesis
confint(madmodneu, level = .9)
confint(madmodneu)
eta2(madmodneu)



###########################SET SUICIDE NEGATIF
#Same steps
#Compute 2 MAD value
MADneg<-2*mad(df$negativewords) 
# Compute the upper limit for 2 MAD frome the median
MADmaxneg<-median(df$negativewords)+MADneg
#Lower limit
MADminneg<-median(df$negativewords)-MADneg
#Make the dataset
madnegative<-df[-(which(df$suicidewords < MADminsuicide | df$suicidewords > MADmaxsuicide |
                        df$negativewords < MADminneg | df$negativewords > MADmaxneg)),]


#testing hypotheses:
#H2
madmodneg <- lm(suicidewords ~ Mirror*discrepancies  + negativewords, data=madnegative)
shapiro.test(resid(madmodneg))
summary(madmodneg)
confint(madmodneg, level = .9)
confint(madmodneg)
eta2(madmodneg)
#Successfully replicated H2 when using MAD2 method


#################  To obtain the untransformed 95% CI :

MirrorDataMCS2019 <- read.csv("~/Mirror effect Octobre 2018/MirrorDataMCS2019.csv")
df<-MirrorDataMCS2019
df<-df[-20,]
df$Mirror<-as.factor(df$Mirror)
df$discrepancies<-as.factor(df$discrepancies)
df$genre<-as.factor(df$genre)

#H1 Confirmatory :
dfclean<-df[-69,] #Removing line 69 associated to studentized residual >3 in a new data set called "dfclean" for clean dataframe
modneu<-lm(dfclean$suicidewords~dfclean$Mirror*dfclean$discrepancies+dfclean$neutralwords)
#To obtain 95% one-sided confidence interval, we take the 90% confidence interval and only keep the lower bound while replacing the upperbound with infinity 
confint(modneu, level = 0.90)
#To obtain the non-direction effects 95% CI:
confint(modneu)


#H2 Confirmatory
dfclean<-df[-69,]
modneg<-lm(suicidewords~Mirror*discrepancies+negativewords, data= dfclean)
confint(modneg, level = 0.90)
confint(modneg)




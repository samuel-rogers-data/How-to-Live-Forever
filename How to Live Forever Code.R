# Sam Rogers 
# How to Live Forever 
# Language: R
####################### Initial Import Data (Below)  ######################## ######################################################################
# load required packages using a bulk check/install function from professor
packages <- c("cowplot", "googleway", "ggplot2", "ggrepel",  "ggspatial", "lwgeom", "sf", "rnaturalearth",
              "rnaturalearthdata","maps","dplyr", "sqldf","readr","reshape","neuralnet","dplyr","e1071","kernlab",
              "googleway","rgeos","rgdal","maptools","scales","shiny","Rcmdr","tidyverse","RcmdrMisc")
## function to check for libraries and install if necessary.
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

#This script imports data from the chsi dataset and combines several sheets to
#make it easier to work with. You need to run it from within the folder that
#contains all the individual csv files

# use sqldf for joins
library(sqldf)
library(readr)

#!!change to your own directory where the CSV files are located!!
setwd("C:/Users/sarogers.TX.AP/Syracuse University/Group-Data Science MS Courses - Documents/Additionals/PastMSCourses/IST 687 Applied Data Science/Project/chsi_dataset")


#this is the one that explains what columns mean
key = read_csv("DATAELEMENTDESCRIPTION.csv")

#Healthypeople is a reference sheet; we shouldn't join it. Instead, we'll use it
#later for comparison to standards
healthypeople = read_csv("HEALTHYPEOPLE2010.csv")

#These 5 will be combined
demographics = read_csv("DEMOGRAPHICS.csv")
vulnpopsandenvhealth = read_csv("VUNERABLEPOPSANDENVHEALTH.csv")
riskfactors = read_csv("RISKFACTORSANDACCESSTOCARE.csv")
leadingcausesofdeath = read_csv("LEADINGCAUSESOFDEATH.csv")
summarymeasures = read_csv("SUMMARYMEASURESOFHEALTH.csv")

combined <- sqldf('
select 
  d.county_fips_code  ,d.state_fips_code  ,d.chsi_county_name  ,d.chsi_state_name  ,d.chsi_state_abbr  ,d.strata_id_number  ,d.strata_determining_factors  ,d.number_counties  ,d.population_size  ,d.population_density  ,d.poverty  ,d.Age_19_Under  ,d.Age_19_64 ,d.Age_65_84, d.Age_85_and_Over, d.White, d.Black, d.Native_American, d.Asian, d.Hispanic, v.No_HS_Diploma, v.Unemployed, v.Sev_Work_Disabled, v.Major_Depression, v.Recent_Drug_Use, v.Ecol_Rpt, v.Ecol_Rpt_Ind, v.Ecol_Exp, v.Salm_Rpt, v.Salm_Rpt_Ind , v.Salm_Exp, v.Shig_Rpt, v.Shig_Rpt_Ind, v.Shig_Exp, v.Toxic_Chem, v.Carbon_Monoxide_Ind, v.Nitrogen_Dioxide_Ind, v.Sulfur_Dioxide_Ind, v.Ozone_Ind, v.Particulate_Matter_Ind,  v.Lead_Ind, v.EH_Time_Span, r.No_Exercise, r.Few_Fruit_Veg, r.Obesity, r.High_Blood_Pres, r.Smoker, r.Diabetes, r.Uninsured, r.Elderly_Medicare, r.Disabled_Medicare, r.Prim_Care_Phys_Rate, r.Dentist_Rate, r.Community_Health_Center_Ind, r.HPSA_Ind, l.A_Wh_Comp, l.A_Bl_Comp, l.A_Ot_Comp, l.A_Hi_Comp, l.A_Wh_BirthDef, l.A_Bl_BirthDef, l.A_Ot_BirthDef, l.A_Hi_BirthDef, l.B_Wh_Injury, l.B_Bl_Injury, l.B_Ot_Injury, l.B_Hi_Injury, l.B_Wh_Cancer, l.B_Bl_Cancer, l.B_Ot_Cancer, l.B_Hi_Cancer, l.B_Wh_Homicide, l.B_Bl_Homicide, l.B_Ot_Homicide, l.B_Hi_Homicide, l.C_Wh_Injury, l.C_Bl_Injury, l.C_Ot_Injury, l.C_Hi_Injury, l.C_Wh_Homicide, l.C_Bl_Homicide, l.C_Ot_homicide, l.C_Hi_Homicide, l.C_Wh_Suicide, l.C_Bl_Suicide, l.C_Ot_Suicide, l.C_Hi_Suicide, l.C_Wh_Cancer, l.C_Bl_Cancer, l.C_Ot_Cancer, l.C_Hi_Cancer, l.D_Wh_Injury, l.D_Bl_Injury, l.D_Ot_Injury, l.D_Hi_Injury, l.D_Wh_Cancer, l.D_Bl_Cancer, l.D_Ot_Cancer, l.D_Hi_Cancer, l.D_Wh_HeartDis, l.D_Bl_HeartDis, l.D_Ot_HeartDis, l.D_Hi_HeartDis, l.D_Wh_Suicide, l.D_Bl_Suicide, l.D_Ot_Suicide, l.D_Hi_Suicide, l.D_Wh_HIV, l.D_Bl_HIV, l.D_Ot_HIV, l.D_Hi_HIV, l.D_Wh_Homicide, l.D_Bl_Homicide, l.D_Ot_Homicide, l.D_Hi_Homicide, l.E_Wh_Cancer, l.E_Bl_Cancer, l.E_Ot_Cancer, l.E_Hi_Cancer, l.E_Wh_HeartDis, l.E_Bl_HeartDis, l.E_Ot_HeartDis, l.E_Hi_HeartDis, l.F_Wh_HeartDis, l.F_Bl_HeartDis, l.F_Ot_HeartDis, l.F_Hi_HeartDis, l.F_Wh_Cancer, l.F_Bl_Cancer, l.F_Ot_Cancer, l.F_Hi_Cancer, l.LCD_Time_Span, s.ALE, s.US_ALE, s.All_Death  ,s.US_All_Death, s.Health_Status, s.US_Health_Status, s.Unhealthy_Days
,s.US_Unhealthy_Days 
from demographics d 
left join vulnpopsandenvhealth v on v.county_fips_code=d.county_fips_code and v.state_fips_code=d.state_fips_code
left join riskfactors r on r.county_fips_code=d.county_fips_code and r.state_fips_code=d.state_fips_code
left join leadingcausesofdeath l on l.county_fips_code=d.county_fips_code and l.state_fips_code=d.state_fips_code
left join summarymeasures s on s.county_fips_code=d.county_fips_code and s.state_fips_code=d.state_fips_code
' )

# We know that negative numbers are really NAs; let's replace them
combined[,-1:-7] <-  data.frame(lapply(combined[,-1:-7], function(x){
  as.numeric(gsub("-1*|-2*|-9*",NA,x))
}))

str(combined)
write_csv(combined,'combined.csv')

#replace NAs with mean of column auto
na_mean_swap <- function(x) {
  replace(x, is.na(x),mean(as.numeric(x),na.rm=TRUE))
}

mean_clean <- cbind(combined[,1:7],replace(combined[,-1:-7],TRUE, lapply(combined[,-1:-7], na_mean_swap)))
str(mean_clean)
################### Initial Import Data ^ #############################
################################################################

############### Plot Average Life Expectancy by State ####################
# Sam Rogers
#assistance from https://www.datanovia.com/en/blog/ggplot-themes-gallery/
Life_Exp_States <- ggplot(data = mean_clean, aes(x=CHSI_State_Abbr, y = ALE, color = ALE, group=CHSI_State_Abbr)) + 
  geom_point() +
  scale_color_continuous(name = "Average Age") +
  theme_classic() +
  xlab(label = "States") +
  ylab(label = "Average Life Expectancy") + 
  scale_y_reverse() +
  ggtitle(label = "State Average Life Expectancy") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Life_Exp_States
################### plot ALE Under 1st Quartile #######################
#################################################################
# Sam Rogers
Life_Exp_75 <- ggplot(data = mean_clean[mean_clean$ALE<75,], aes(x=CHSI_State_Abbr, y = ALE, color = ALE, group=CHSI_State_Abbr)) + 
  geom_point() +
  scale_color_continuous(name = "Average Age") +
  theme_classic() +
  xlab(label = "States") +
  ylab(label = "Average Life Expectancy") + 
  scale_y_reverse() +
  ggtitle(label = "State ALE Below Nation's 1st Quartile") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Life_Exp_75

################ Plot of Smokers relating to ALE and Exercise ################ #####################################################################
# Sam Rogers
# Shows a plot of no exercise to ALE. Less smoking in green, heavy smoking in blue. 
#determine the range of smokers 
summary(combined$Smoker)
#plot 
#assistance from https://www.datanovia.com/en/blog/ggplot-themes-gallery/
gtt <- ggplot(data = combined, aes(x=No_Exercise, y = ALE,   color = Smoker)) + 
  theme_classic() +  
  geom_point() +  
  xlab(label = "Percent not Exercising") +
  scale_x_reverse() + #x_axis reversed to show positive correlation
  ylab(label = "Average Life Expectancy") + 
  scale_colour_gradient(limits=c(3.6,46.2), low="green", high="red") + #Smoker Range
  ggtitle(label = "Negative Affects of Smoking and Not Working Out on ALE")
gtt

#Create formula to later take statistics on the plot
#These are all class or personal knowledge 
ggtt <- lm(formula= ALE ~ No_Exercise, data = combined)    #ALE(dep)
ggts <- lm(formula= ALE ~ Smoker, data = combined)         #ALE(dep)
# shows adj R square of ALE to smoking and No exercise plots 
summary(ggtt)
summary(ggts)

############# Barchart Summary of Variables Most Affecting Death #############
#####################################################################
#Sam Rogers
#this is from the section which takes population and percent of population into #consideration. It's a true average of the population as values are originally presented #as percentage of county.
#This calculates the count then takes percent of the US
DF_Bar <- mean_clean[,c('Poverty','Unemployed','Sev_Work_Disabled','Recent_Drug_Use',' No_Exercise','Obesity','High_Blood_Pres','Smoker','Diabetes','Uninsured','Elderly_Medicare', 'Disabled_Medicare')]
PLot_Bar <- replace(DF_Bar, TRUE, lapply(DF_Bar, na_mean_swap))
#calculate mean of each column
PLot_Bar <- colMeans(PLot_Bar)
#assistance from Quick-R by DataCamp
#used melt from library(reshape)
#this makes each row a unique id-variable combination
PLot_Bar <- melt(PLot_Bar)
PLot_Bar$Variables <- rownames(PLot_Bar)
#used Cookbook for R, Colors(ggplot2) for assistance
#toned the colors down to be a bit darker and match the colors of our other graphs
BarPlotFin<- ggplot(data=PLot_Bar) + geom_bar(aes(x = Variables,y=value,fill=Variables),stat="identity") + coord_flip() +scale_fill_hue(l=45)

###############################################################
######################Suicide to Toxic Chemicals###########
#Randall
scatterplot(Suicide~Toxic_Chem, regLine=TRUE, smooth=FALSE,
            boxplots=FALSE, xlab="Toxic Chemicals", ylab="Suicide",
            data=sContribs.df)
cor(sContribs.df[,c("Suicide","Toxic_Chem")], use="complete")

################Suicide to Poverty#######################
scatterplot(Suicide~Poverty,
            regLine=TRUE, smooth=FALSE,
            boxplots=FALSE,
            xlab="Poverty",
            ylab="Suicide", data=sContribs.df)



####Correlations Matrix runs, and Visualizations
#First run removed environmental factors that had no correlative value

sContribs.corr3 = cor(ALEdimin.df[,c( "No_HS_Diploma", "Unemployed","Sev_Work_Disabled", "Major_Depression", "Recent_Drug_Use","No_Exercise", "Few_Fruit_Veg", "Obesity", "Smoker", "Diabetes", "Uninsured", "Elderly_Medicare", "Disabled_Medicare", "High_Blood_Pres","Toxic_Chem", "Carbon_Monoxide_Ind", "Ozone_Ind", "Particulate_Matter_Ind", "Lead_Ind", "Disabled_Medicare", "Prim_Care_Phys_Rate", "Dentist_Rate", "Community_Health_Center_Ind","ALE","Suicide", "Premature", "Under_18", "Over_40", "Unmarried", "Late_Care", "Infant_Mortality", "IM_Hisp", "IM_Neonatal", "IM_Postneonatal", "Brst_Cancer", "Col_Cancer", "CHD", "Homicide", "Lung_Cancer", "MVA", "Stroke", "Injury" ,  "Total_Deaths", "All_Death", "Unhealthy_Days")], use="complete")
#### Second run; removed categorical descriptions: US ALE, US UNHEALTHY DAYS
sContribs.corr3 = cor(ALEdimin.df[,c( "No_HS_Diploma", "Unemployed","Sev_Work_Disabled", "Major_Depression", "Recent_Drug_Use","No_Exercise", "Few_Fruit_Veg", "Obesity", "Smoker", "Diabetes", "Uninsured", "Elderly_Medicare", "Disabled_Medicare", "High_Blood_Pres","Toxic_Chem", "Carbon_Monoxide_Ind", "Ozone_Ind", "Particulate_Matter_Ind", "Lead_Ind", "Disabled_Medicare", "Prim_Care_Phys_Rate", "Dentist_Rate", "Community_Health_Center_Ind","ALE","Suicide", "Premature", "Under_18", "Over_40", "Unmarried", "Late_Care", "Infant_Mortality", "IM_Hisp", "IM_Neonatal", "IM_Postneonatal", "Brst_Cancer", "Col_Cancer", "CHD", "Homicide", "Lung_Cancer", "MVA", "Stroke", "Injury" ,  "Total_Deaths", "All_Death", "Unhealthy_Days")], use="complete")

#####look at a single correlation between ALE and Suicide 
sContribs.corr3.1 = cor(ALEdimin.df[c("ALE", "Suicide")], use='complete')
View(sContribs.corr3.1)

##########Correlations matrix visualization ALE correalting risk factors 
sContribs.corrPlt3 = corrplot::corrplot(sContribs.corr3, method = ('pie'), type = ('upper'), add = FALSE, col = NULL, bg = "white", title = "ALE correlating risk factors", is.corr=TRUE, diag = TRUE, outline = FALSE, addgrid.col = TRUE,addCoef.col = NULL,addCoefasPercent = TRUE, sig.level = 0.05)
plot(sContribs.corrPlt3)
##########Correlations matrix visualization ALE correalting risk factors next iteration
sContribs.corrPlt3.1 = corrplot::corrplot(sContribs.corr3.1, method = ('pie'), type = ('upper'), add = FALSE, col = NULL, bg = "white", title = "ALE correlating risk factors", is.corr=TRUE, diag = TRUE, outline = FALSE, addgrid.col = TRUE,addCoef.col = NULL,addCoefasPercent = TRUE, sig.level = 0.05)
######Visualization for interactive maps##### 
##Required Libraries are as follows:
packages <- c("readr","sqldf","cowplot", "googleway", "ggplot2", "ggrepel",  "ggspatial", "lwgeom", "sf", "rnaturalearth",
              "rnaturalearthdata","maps","dplyr","raster","stringr","spData","leaflet","sp","tidyr","devtools","tmap","tmaptools","corrplot","rgeos","rgdal","maptools","scales","shiny","Rcmdr")

ggplot(data = world) +
  geom_sf() +
  geom_sf(data=noHSDIPLOMA1.map,aes(fill=No_HS_Diploma),color=gray(0.1)) +
  geom_sf(data=states,fill=NA,size=0.8,color=gray(0.1)) #+
#coord_sf(xlim=xlimit,ylim=ylimit, expand=FALSE)
plot(noHSDIPLOMA1.map)

nohsD.map1 = tm_shape(noHSDIPLOMA1.map) + tm_fill("No_HS_Diploma", legend.show = TRUE, id = "county", palette = "Greens" ) +  tm_polygons("county", "orange", style = "cont",  n = 2, alpha = 1 , stretch.palette = TRUE, popup.vars = NA, convert2density = TRUE, midpoint = FALSE)
tmap_mode("view")
nohsD.map1
#####Bar Graphs for Demographics######
#Transform demographics for valid columns
# subset data to only include some columns
demographics <- subset(demographics, select = 
                         c(CHSI_County_Name:CHSI_State_Abbr, Population_Size, 
                           Population_Density, Poverty, Age_19_Under, Age_19_64, 
                           Age_65_84, Age_85_and_Over, White, Black, Native_American, Asian, Hispanic))
# clean names of the subset
nms_demo_dat <- c("county.name","state.name","state.abbr","pop.size","pop.density",
                  "poverty","age.19_under","age.19_64","age.65_84","age.85_over",
                  "white","black","nat.amer","asian","hispanic")
# change col names
names(demographics)<-nms_demo_dat
#This data is a representation at the county level
#we have to do some work to find it at state level
mean(demographics$pop.size)
max(demographics$pop.size)
min(demographics$pop.size)
sd(demographics$pop.size)
#made with help from R Studio Community
library(tidyverse) #load library
# build vectors  
county <- c(demographics$county.name)
state <- c(demographics$state.name)
pop <- c(demographics$pop.size)
# now we assemble into data frame from vectors
df <- data.frame(county, state, pop)
# assembled pipe 
stateandmeanpop <- df %>% group_by(state) %>% summarize(mean_pop = mean(pop))
#Now we have a dataframe we can make different charts from
#with each states mean population
#This is the true mean of the population in our dataset
mean(stateandmeanpop$mean_pop)
#This is the range of the population in our dataset
max(stateandmeanpop$mean_pop)
min(stateandmeanpop$mean_pop)
sd(stateandmeanpop$mean_pop)
#B.	Show mean age, and ranges (a distribution would be a good visual for this as well)
#We can do the same with age. 
demographics$age.19_underraw = (demographics$pop.size*demographics$age.19_under/100)
demographics$age.19_64raw = (demographics$pop.size*demographics$age.19_64/100)
demographics$age.65_84raw = (demographics$pop.size*demographics$age.65_84/100)
demographics$age.85_overraw = (demographics$pop.size*demographics$age.85_over/100)
xviiii_U <- mean(demographics$age.19_under)
max(demographics$age.19_under)
min(demographics$age.19_under)

xvii_lxiv <- mean(demographics$age.19_64)
max(demographics$age.19_64)
min(demographics$age.19_64)

lxv_lxxxiv <- mean(demographics$age.65_84)
max(demographics$age.65_84)
min(demographics$age.65_84)

lxxxv_O <- mean(demographics$age.85_over)
max(demographics$age.85_over)
min(demographics$age.85_over)
meanofage <- c(xviiii_U, xvii_lxiv, lxv_lxxxiv, lxxxv_O)
barplot (meanofage,
         main = "National Age %",
         xlab = "Mean % Age",
         ylab = "Population %",
         names.arg = c("Under19", "19-64", "65-84", "Over85"),
         col = "blue",
         horiz = FALSE)
#C.	Show ethnicity nationally (e.g. 87% white, etc.) - bar chart for visual
#How to create raw population numbers for age and ethnicity
#First we have to add new columns in the data.frame with the raw data
demographics$whiteraw = (demographics$pop.size*demographics$white/100)
demographics$blackraw = (demographics$pop.size*demographics$black/100)
demographics$hispanicraw = (demographics$pop.size*demographics$hispanic/100)
demographics$asianraw = (demographics$pop.size*demographics$asian/100)
demographics$nat.amerraw = (demographics$pop.size*demographics$nat.amer/100)

#Create mean vectors for each ethnicity
w <- mean(demographics$whiteraw)
b <- mean(demographics$blackraw)
h <- mean(demographics$hispanicraw)
a <- mean(demographics$asianraw)
n.a <- mean(demographics$nat.amerraw)
#Create dataframe
meanofrace <- c(w,b,h,a,n.a)
#Create bar chart 
barplot (meanofrace,
         main = "National Mean Ethnicity",
         xlab = "Mean Ethnicity",
         ylab = "Population",
         names.arg = c("White", "Black", "Hispanic", "Asian", "NativeAmerican"),
         col = "darkred",
         horiz = FALSE)

##################
# linear model & svm
##################
# Linear model
# After you've imported CHSI
library(neuralnet)
library(dplyr)
library(e1071)
library(kernlab)

ggplot(data=mean_clean,aes(x='Life expectancy',y=ALE)) + geom_boxplot()  

vuln <- sqldf('
          	select
            	county_fips_code
            	,state_fips_code
            	,chsi_county_name
            	,chsi_state_name
            	,population_size
            	,poverty
            	,No_HS_Diploma
            	,Unemployed
            	,Sev_Work_Disabled
            	,Major_Depression
            	,Recent_Drug_Use
            	,No_Exercise
            	,Few_Fruit_Veg
            	,Obesity
            	,High_Blood_Pres
            	,Smoker
            	,Diabetes
            	,Uninsured
            	,Elderly_Medicare
            	,Disabled_Medicare
            	,Prim_Care_Phys_Rate
            	,Dentist_Rate
            	,ALE
            	,All_Death
            	,US_All_Death
          	from combined')

str(vuln)
#Convert raw count data into percentages so we can compare it across counties with different populations
vuln$No_HS_Diploma <- vuln$No_HS_Diploma/vuln$Population_Size*100
vuln$Unemployed <- vuln$Unemployed/vuln$Population_Size*100
vuln$Sev_Work_Disabled <-vuln$Sev_Work_Disabled/vuln$Population_Size*100
vuln$Major_Depression <- vuln$Major_Depression/vuln$Population_Size*100
vuln$Recent_Drug_Use <- vuln$Recent_Drug_Use/vuln$Population_Size*100
vuln$Uninsured <- vuln$Uninsured/vuln$Population_Size*100
vuln$Elderly_Medicare <- vuln$Elderly_Medicare/vuln$Population_Size*100
vuln$Disabled_Medicare <- vuln$Disabled_Medicare/vuln$Population_Size*100
vuln$DeathRate <- vuln$All_Death/vuln$Population_Size*100

cor(vuln[,c("ALE","Dentist_Rate","Diabetes","Disabled_Medicare",
            "Elderly_Medicare","Few_Fruit_Veg","High_Blood_Pres","Major_Depression",
            "No_Exercise","No_HS_Diploma","Obesity","Population_Size","Poverty",
            "Prim_Care_Phys_Rate","Recent_Drug_Use","Sev_Work_Disabled","Smoker",
            "Unemployed","Uninsured")], use="complete")

LinearModel.1 <- lm(ALE ~ Dentist_Rate + Diabetes + Disabled_Medicare +
                      Elderly_Medicare + Few_Fruit_Veg + High_Blood_Pres + Major_Depression +
                      No_Exercise + No_HS_Diploma + Obesity + Poverty + Prim_Care_Phys_Rate +
                      Recent_Drug_Use + Sev_Work_Disabled + Smoker + Unemployed + Uninsured,
                    data=vuln)
summary(LinearModel.1)

v_lm <- lm(ALE ~ Poverty +No_HS_Diploma +Unemployed +Sev_Work_Disabled +Major_Depression +Recent_Drug_Use +No_Exercise +Few_Fruit_Veg +Obesity +High_Blood_Pres +Smoker +Diabetes +Uninsured +Elderly_Medicare +Disabled_Medicare +Prim_Care_Phys_Rate +Dentist_Rate
           ,data=vuln)
summary(v_lm)

v_lm1 <- lm(ALE ~ Poverty + Unemployed + Sev_Work_Disabled + Recent_Drug_Use +No_Exercise + Obesity + High_Blood_Pres + Smoker + Diabetes + Uninsured + Elderly_Medicare + Disabled_Medicare
            ,data=vuln)
summary(v_lm1)

######################
# svm
######################
#add ID number so we can easily split into train/test
vuln <- vuln %>% mutate(id=row_number())
vuln <- na.omit(vuln)

#put 70 percent of the data into a training dataset
rtrain <- vuln %>% sample_frac(0.7)

#put the rest into a testing dataset
rtest <- anti_join(vuln,rtrain, by='id')

svmOutput <- ksvm(ALE ~ Poverty + Unemployed + Sev_Work_Disabled + Recent_Drug_Use +No_Exercise + Obesity + High_Blood_Pres + Smoker + Diabetes + Uninsured + Elderly_Medicare + Disabled_Medicare
                  , data = rtrain,
                  kernel = "rbfdot", # kernel function that projects the low dimensional problem into higher dimensional space
                  kpar = "automatic",# kpar refer to parameters that can be used to control the radial function kernel(rbfdot)
                  C = 10, # C refers to "Cost of Constrains"
                  cross = 10, # use 10 fold cross validation in this model
                  prob.model = TRUE # use probability model in this model
)
# check the model
svmOutput

# 2) Test the model with the testData data set
svmPred <- predict(svmOutput, # use the built model "svmOutput" to predict
                   rtest, # use testData to generate predictions
                   type = "votes" # request "votes" from the prediction process
)


# create a comparison dataframe that contains the exact "Ozone" value and the predicted "Ozone" value
# use for RMSE calc

compTable <- data.frame(rtest[,c('ALE')], svmPred[,1])
colnames(compTable) <- c("actual","Pred")


compTable$diff <- abs(compTable$actual-compTable$Pred)
compTable$err <- (1 -compTable$diff/compTable$actual)*100
compTable
# compute the Root Mean Squared Error
ksvm_rms <- sqrt(mean((compTable$actual-compTable$Pred)^2)) #A smaller value indicates better model performance.

cat('Actual average:  	',mean(compTable$actual),' years',
    '\nActual deviation:	',sd(compTable$actual),' years',
    '\nPredicted average:   ',mean(compTable$Pred),' years',
    '\nPredicted deviation: ',sd(compTable$Pred),' years',
    '\nAverage difference:  ',mean(compTable$diff),' years',
    '\nPercent accuracy:	',mean(compTable$err),'%',
    '\nRMS is           	',ksvm_rms)


ggplot(data=compTable,aes(x=1:nrow(compTable),y=actual,color='Actual')) + geom_line() + geom_line(aes(y=Pred,color="Predicted"),alpha=0.5) + xlab('rownum') + scale_color_manual(values=c('darkblue','red')) +labs(color='Type')

ggplot(data=compTable) + geom_boxplot(aes(x='Actual',y=actual,fill='Actual')) + geom_boxplot(aes(x='Predicted',y=Pred,fill='Predicted')) + labs(fill='Data type',y='ALE',x='')

sd(rtest[,c('ALE')])
summary(rtest[,c('ALE')])


#######################
# Add neuralnet
#######################

#sampling example from https://medium.com/@HollyEmblem/training-and-test-dataset-creation-with-dplyr-41d9aa7eab31
#Assign row number to an id column within dataframe; this allows us to easily
#split the data
str(rtrain)
#create neural net model using training data
vuln_net <- neuralnet(ALE ~ Poverty + Unemployed + Sev_Work_Disabled + Recent_Drug_Use +No_Exercise + Obesity + High_Blood_Pres + Smoker + Diabetes + Uninsured + Elderly_Medicare + Disabled_Medicare
                      ,data=rtrain, hidden=2,lifesign='minimal',linear.output=FALSE,threshold = 0.5)

#here are coeeficients
vuln_net$result.matrix

#here's what the net looks like
plot(vuln_net,rep='best')

#predict results of the testing dataset using the model
v.results <- predict(vuln_net,rtest[,c('Poverty','Unemployed','Sev_Work_Disabled','Recent_Drug_Use','No_Exercise','Obesity','High_Blood_Pres','Smoker','Diabetes','Uninsured','Elderly_Medicare','Disabled_Medicare' )])
#create comparison dataframe with the real data and the predicted data
v.compare <- data.frame(actual = rtest[,c('ALE')],predicted=v.results)
typeof(v.compare)
# add diff column
v.compare$diff <- abs(v.compare$actual - v.compare$predicted )
# add accuracy column. Since we're comparing percentages, I *think* we just need
# to take 1 minus the diff. Is that right?
v.compare$accuracy <- 100 - abs(v.compare$actual - v.compare$predicted )
mean(v.compare$accuracy)
sd(v.compare$accuracy)
v.compare

#bar plot for Sam

dgraph <- vuln[,c('Poverty','Unemployed','Sev_Work_Disabled','Recent_Drug_Use','No_Exercise','Obesity','High_Blood_Pres','Smoker','Diabetes','Uninsured','Elderly_Medicare','Disabled_Medicare')]

dgraph <- colMeans(dgraph)

dgraph <- melt(dgraph)
dgraph$variable <- rownames(dgraph)

ggplot(data=dgraph) + geom_bar(aes(x = variable,y=value,fill=variable),stat="identity") + coord_flip()
#Code found on Github by Todd and provided by Luke Miller
#Code for Poverty maps 
#a copy of the function
lower.df = function(v) 
{
  if(is.character(v)) return(tolower(v)) 
  else return(v)
}
# subset data to only include some columns
demographics <- subset(demographics, select = 
                         c(CHSI_County_Name:CHSI_State_Abbr, Population_Size, 
                           Population_Density, Poverty, Age_19_Under, Age_19_64, 
                           Age_65_84, Age_85_and_Over, White, Black, Native_American, Asian, Hispanic))
# clean names of the subset
nms_demo_dat <- c("county.name","state.name","state.abbr","pop.size","pop.density",
                  "poverty","age.19_under","age.19_64","age.65_84","age.85_over",
                  "white","black","nat.amer","asian","hispanic")
# change col names
names(demographics)<-nms_demo_dat
head(demographics)
demographics$state.name <- lower.df(demographics$state.name)
demographics$county.name <- lower.df(demographics$county.name)
demographics$state.abbr <- lower.df(demographics$state.abbr)
# get state map
us_state <- map_data("state")
# change us_state_map names
names(us_state)<- c("long","lat","group","order","state.name","subregion")
# merge us_state_map data frame with demographics_dat by state.name
# only include matching records
Po_data<-merge(us_state, demographics, by ='state.name')
# preserve order
Po_data<-Po_data[order(Po_data$order),]
# remove subregion column
Po_data$subregion<-NULL
# split %'s into 6 cuts 
Po_data$poverty <- cut_interval(Po_data$poverty, 6)
# state data
state_df <- map_data("state")
# create dataframe with county information from maps 
# Longitude and Latitude information here
county_df<-map_data("county")
# change names of county_df
names(county_df)<- c("long","lat","group","order","state.name","county.name")
# check out state.abb and state.name Datasets 
# will add state.abbr to county_df based on match
head(state.abb)
head(state.name)
# add a column with state abbreviations based on matching between 
# county_df$state_name and lowercase state.name dataset
county_df$state.abbr<- state.abb[match(x = county_df$state.name, tolower(state.name))]
# remove state_name column since have abbreviations
county_df$state.name<-NULL

# make all names lowercase to match
county_df <- data.frame(lapply(county_df, lower.df))

#--------------------------------------------------------------
# will use this to zoom in on % poverty at County level
#--------------------------------------------------------------
# merge county_df and demographics_dat by county.name and state.abbr
Pop_map <- merge(x = county_df, y = demographics, by=c("county.name","state.abbr"))
# retain order
Pop_map <- Pop_map[order(Pop_map$order), ]
# add breaks for ranges in new map risk_map
Pop_map$poverty <- cut_interval(Pop_map$poverty ,6)
#-----------------------------------------------------------------------------------------------
# plot Poverty across US states - not at county level
#-----------------------------------------------------------------------------------------------
# create dataframe to add state map abbreviations on map
# state.center - x=long, y=lat. state.abb is a list containing state abbreviations
state.info <- data.frame(state.center, state.abb)
# lower names
state.info$state.abb <- tolower(state.info$state.abb)
# add group info
state.info$group <- Po_data$group[match(x = state.info$state.abb, Po_data$state.abbr)]
# remove ak and hi (no group)
state.info <- state.info[!is.na(state.info$group),]
# map of poverty at state level - org palette
# doesnt include state names
Pop_map_all_org <- ggplot(Po_data, aes(long, lat, group = group)) +
  geom_polygon(aes(fill = poverty)) +
  geom_polygon(data = state_df, colour = "black", fill = NA, size =0.2) + 
  geom_polygon(data = county_df, colour = "snow", fill = NA, size =0.1) + 
  geom_text(data = state.info, aes(x=x, y=y, label = state.abb, group = group), colour ='black') +
  theme_classic() +  
  theme(legend.position="right") + 
  xlab(label = "") +
  ylab(label = "") +
  scale_x_continuous(expand = c(0,0)) + # expand size of map along x axis 
  scale_y_continuous(expand = c(0,0)) + # expand size of map along y axis
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + 
  ggtitle(label = "US Map showing Poverty Intervals at the State Level") +
  scale_fill_brewer(type='div', palette = 'RdBu', name = "% Poverty")
Pop_map_all_org

#Poverty Calculated at the county level

poverty_Cmap_purd <- ggplot(Pop_map, aes(long, lat, group = group)) +
  theme(panel.background = element_rect(fill = "snow")) +
  geom_polygon(aes(fill = poverty), colour = alpha("snow", 1/2), size = 0.2) + 
  geom_polygon(data = state_df, colour = "black", fill = NA,size =0.2) + 
  geom_polygon(data = county_df, colour = "grey", fill = NA, size =0.1) + 
  theme_classic() +  
  theme(legend.position="bottom") + 
  xlab(label = "") +
  ylab(label = "") +
  scale_x_continuous(expand = c(0,0)) + # expand size of map along x axis 
  scale_y_continuous(expand = c(0,0)) + # expand size of map along y axis
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank()) + 
  ggtitle(label = "US Map showing poverty \nIncludes County information") +
  scale_fill_brewer(type='seq',palette = 'PuRd', name ="% Poverty") 
poverty_Cmap_purd 

#######################
# SF map plotting groundwork
#######################
# Luke Miller
# 2 September 2019
# following map drawing tutorial
# from https://www.r-spatial.org/r/2018/10/25/ggplot2-sf-2.html
# Modifying to visualize project data

#load required packages using a bulk check/install function from professor
packages <- c("cowplot", "googleway", "ggplot2", "ggrepel",  "ggspatial", "lwgeom", "sf", "rnaturalearth", "rnaturalearthdata","maps",'dplyr','reshape')

package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

theme_set(theme_bw())

#################
# Plot world map
#################

world <- ne_countries(scale='medium',returnclass='sf')
#class(world)

ggplot(data=world) +
  geom_sf() +
  coord_sf()

##################
# Plot state map
##################

#turn map state into a shapefile
states <- st_as_sf(map("state",plot=FALSE,fill=TRUE))
#add coordinates of 'centroid' so we can later plot names
states <- cbind(states, st_coordinates(st_centroid(states)))
#label state names as uppercase
states$name <- toupper(states$ID)

#add us so we can use the limits
us <- map_data("county")[,c('long','lat')]
xlimit <- c(max(us$long)+2,min(us$long)-2)
ylimit <- c(max(us$lat)+2,min(us$lat)-2)


ggplot(data=world) +
  geom_sf() +
  geom_sf(data=states,fill=NA) +
  geom_label(data=states,aes(X,Y, label=name), size=3) +
  coord_sf(xlim=xlimit,ylim=ylimit, expand=FALSE)

#############################
# Zoom in on chosen location
#############################

zoomAmount <- 5

#knoxville x=-83.99,y=35.82 Pick your favorite city
centerx <- -83.99
centery <- 35.82

ylimit <- c(centery-zoomAmount, centery+zoomAmount)
xlimit <- c(centerx-zoomAmount*2, centerx+zoomAmount*2)

ggplot() +
  geom_sf(data=states) +
  geom_label(data=states,aes(X,Y, label=name), size=3) +
  coord_sf(xlim=xlimit,ylim=ylimit, expand=FALSE)

#################
# Counties
#################

counties <- st_as_sf(map("county", plot=FALSE, fill=TRUE))

ggplot(data = world) +
  geom_sf() +
  geom_sf(data=states,size=0.8) +
  geom_sf(data=counties,fill=NA,color=gray(0.5)) +
  coord_sf(xlim=xlimit,ylim=ylimit, expand=FALSE)

############################
# Counties with smoking data
############################

#must have loaded chsi already
smokers <- combined[,c('CHSI_State_Name','CHSI_County_Name','Smoker','Poverty')]
smokers$Smoker <- as.numeric(smokers$Smoker)/100
smokers$state <- tolower(smokers$CHSI_State_Name)
smokers$county<- tolower(smokers$CHSI_County_Name)
smokers$ID <- paste(smokers$state,smokers$county,sep=',')

counties <- st_as_sf(map("county", plot=FALSE, fill=TRUE))

#create new object that combines the counties sf with the smoker dataset
s.map <- left_join(counties,smokers)

zoomAmount <- 3
#knoxville x=-83.99,y=35.82 Pick your favorite city
centerx <- -83.99
centery <- 35.82

ylimit <- c(centery-zoomAmount, centery+zoomAmount)
xlimit <- c(centerx-zoomAmount*2, centerx+zoomAmount*2)

#inherits xlimit and ylimit from above; change as desired
ggplot(data = world) +
  geom_sf() +
  geom_sf(data=s.map,aes(fill=Smoker),color=gray(0.1)) +
  geom_sf(data=states,fill=NA,size=0.8,color=gray(0.1)) +
  coord_sf(xlim=xlimit,ylim=ylimit, expand=FALSE) + labs(title='US Smoking Percentage by County')

##################################
# ALE
##################################

ALE <- combined[,c('CHSI_State_Name','CHSI_County_Name','ALE')]
ALE$state <- tolower(ALE$CHSI_State_Name)
ALE$county<- tolower(ALE$CHSI_County_Name)
ALE$ID <- paste(ALE$state,ALE$county,sep=',')

counties <- st_as_sf(map("county", plot=FALSE, fill=TRUE))

#create new object that combines the counties sf with the smoker dataset
s.map <- left_join(counties,ALE)

#zoomAmount <- 3
##knoxville x=-83.99,y=35.82 Pick your favorite city
#centerx <- -83.99
#centery <- 35.82
#
#ylimit <- c(centery-zoomAmount, centery+zoomAmount)
#xlimit <- c(centerx-zoomAmount*2, centerx+zoomAmount*2)

#inherits xlimit and ylimit from above; change as desired
ggplot(data = world) +
  geom_sf() +
  geom_sf(data=s.map,aes(fill=ALE),color=gray(0.1)) +
  geom_sf(data=states,fill=NA,size=0.8,color=gray(0.1)) +
  coord_sf(xlim=xlimit,ylim=ylimit, expand=FALSE) + labs(title='US life expectancy by county')


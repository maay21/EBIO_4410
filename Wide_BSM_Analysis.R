##################     Wide Blood Glucose Data Analysis       #########################################
##################            Marina Ayala                    #########################################
##################  Independent Project, The Safran-Flaxman Lab   #####################################


############   Data Upload   ######################################################################################################

# uploading the data file
wide_bsm <- read.csv("Bsm_wide.csv", stringsAsFactors = TRUE)

############## Install Packages ####################################################################################################
# downloading the dplyr package in order to use 'filter', 'mutate', and various other functions
install.packages("dplyr")

# downloading ggplot to create seamless graphs
install.packages("ggplot2")

# downloading tidyverse for the time variable
install.packages("tidyverse")

# downloading lubridate for the time variable
install.packages("lubridate")

# downloadding emmeans 
install.packages("emmeans")

# downloading car
install.packages("car")

# downloading circular
install.packages('circular')

# downloading lmerTest
install.packages('lmerTest')

# downloading broom
install.packages('broom')

# downloading lme4
install.packages('lme4')

# downloading dotwhisker
install.packages('dotwhisker')



# calling on the packages from the R library 
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(emmeans)
library(car)
library(circular)
library(lmerTest)
library(broom)
library(lme4)
library(dotwhisker)



################### Filtering Data #################################################################################################
# filtering out the birds that cannot be used for analysis because they do not have a paired glucose measurement

# this bird needs to be dropped from all analysis,because this bird was attending multiple nests 
# renaming wide_bsm to wide_data
wide_data <- wide_bsm %>% filter_all(all_vars(Band != "2850-57510"))
# this bird was used in 2019 spatial analysis, not for the experiment
wide_data <- wide_data %>% filter_all(all_vars(Band != "2640-97570"))
# these birds cannot be used for change data since they were not recaptured 
wide_data <- wide_data %>% filter_all(all_vars(Band != "2350-21097"))
wide_data <- wide_data %>% filter_all(all_vars(Band != "2640-97545"))
wide_data <- wide_data %>% filter_all(all_vars(Band != "2640-97594"))
wide_data <- wide_data %>% filter_all(all_vars(Band != "2640-97598"))

# check the data to make sure that there is is correct number of rows using the length function
length(wide_data$Band)

# filtering out any and all NAs and converting them to strings
wide_data <- wide_data %>% filter(is.na(wide_data)|wide_data!="str")

################ Creating the main dataframe- (band.df) ########################################################################################
# creating a dataframe for the band numbers so I am able keep track of the band numbers
band.df <- as.data.frame(wide_data$Band)
# renaming the column so that I know it denotes band number
band.df <- band.df %>% rename(Band_Number = `wide_data$Band` )

# adding blood glucose 1 and 2 to the band dataframe
band.df <- band.df %>% mutate(wide_data$BloodGluc_1)
band.df <- band.df %>% rename(BloodGluc_1 = 'wide_data$BloodGluc_1')

band.df <- band.df %>% mutate(wide_data$BloodGluc_2)
band.df <- band.df %>% rename(BloodGluc_2 = 'wide_data$BloodGluc_2')


# calculating the change in blood glucose equation
# BloodGluc2 - BloodGluc1 = Delta Blood Gluc
# using that equation to add a new column in band.df that is denoted change_in_blood_gluc
band.df <- band.df %>% mutate(wide_data$BloodGluc_2 - wide_data$BloodGluc_1)
band.df <- band.df %>% rename(change_in_blood_gluc = 'wide_data$BloodGluc_2 - wide_data$BloodGluc_1')

# adding tag / no tag to the band.df and renaming the column to tag/no tag
band.df <- band.df %>% mutate(wide_data$Tagged)
band.df <- band.df %>% rename (tag_no_tag = 'wide_data$Tagged')

# adding mass 1 and mass 2 to the dataframe
band.df <- band.df %>% mutate(wide_data$Mass_1)
band.df <- band.df %>% rename(Mass_1 = 'wide_data$Mass_1')

band.df <- band.df %>% mutate(wide_data$Mass_2)
band.df <- band.df %>% rename(Mass_2 = 'wide_data$Mass_2')

# finding change in mass and adding that variable to band.df
band.df <- band.df %>% mutate(wide_data$Mass_2 - wide_data$Mass_1)
band.df <- band.df %>% rename(change_in_mass = 'wide_data$Mass_2 - wide_data$Mass_1')

# fitting the equation in so I can add Brood Treatment to the dataframe band.df
# ΔBlood Glucose ~ tag + broodTrt + tag*broodTrt + mass1 + year + time +(1|site)
# adding brood size treatment to band.df
band.df <- band.df %>% mutate(wide_data$BroodTrt)
# renaming the column
band.df <- band.df %>% rename(BroodTrt = 'wide_data$BroodTrt')

# adding time_two, date_two, and time_date_two to band.df
band.df <- band.df %>% mutate(wide_data$Time_2)
band.df <- band.df %>% mutate(wide_data$Date_2)


#################### Univariate Statistics #########################################################################################
######## Blood Glucose #############################################################################################################
# first variable to create a table and plot for is blood glucose 1
BloodGluc1_stat <- band.df %>%
  summarise(BG = sum(!is.na(BloodGluc_1)),
            avg.BG = round (mean(BloodGluc_1, 
                                          na.rm = BG),2),
            stdev.BG = round (sd(BloodGluc_1, 
                                          na.rm = BG), 2),
            med.BG = round(median(BloodGluc_1,
                                           na.rm = BG), 2),
            min.BG = round(min(BloodGluc_1,
                                        na.rm = BG), 2),
            max.BG = round(max(BloodGluc_1,
                                        na.rm = BG), 2))
# plotting blood glucose 1
hist(band.df$BloodGluc_1, breaks = 15, xlim = c(140,280), ylim =c(0,8), main = "Frequency of Blood Glucose Levels (μg/dL) Before Treatment", xlab = "Blood Glucose Levels Before Treatment (μg/dL)", col = "turquoise")

# blood glucose 2
BloodGluc2_stat <- band.df %>%
  summarise(BG = sum(!is.na(BloodGluc_2)),
            avg.BG = round (mean(BloodGluc_2, 
                                 na.rm = BG),2),
            stdev.BG = round (sd(BloodGluc_2, 
                                 na.rm = BG), 2),
            med.BG = round(median(BloodGluc_2,
                                  na.rm = BG), 2),
            min.BG = round(min(BloodGluc_2,
                               na.rm = BG), 2),
            max.BG = round(max(BloodGluc_2,
                               na.rm = BG), 2))
# plotting blood glucose 2
hist(band.df$BloodGluc_2, breaks = 15, xlim = c(150,350), ylim = c(0,6), main = "Frequency of Blood Glucose Levels (μg/dL) After Treatment", xlab= "Blood Glucose Levels After Treatment (μg/dL)", col = "orange")

# change in blood glucose
BloodGlucChange_stat <- band.df %>%
  summarise(BG = sum(!is.na(change_in_blood_gluc)),
            avg.BG = round (mean(change_in_blood_gluc, 
                                 na.rm = BG),2),
            stdev.BG = round (sd(change_in_blood_gluc, 
                                 na.rm = BG), 2),
            med.BG = round(median(change_in_blood_gluc,
                                  na.rm = BG), 2),
            min.BG = round(min(change_in_blood_gluc,
                               na.rm = BG), 2),
            max.BG = round(max(change_in_blood_gluc,
                               na.rm = BG), 2))
# plotting change in blood glucose intitally to understand the distribution
hist(band.df$change_in_blood_gluc, breaks = 8, xlim= c(-100,150), ylim = c(0,10), main= "Frequency of the Change in Blood Glucose Levels (μg/dL) Throughout Treatment", xlab= "Change in Blood Glucose Levels (μg/dL)", col= "lavender")

######### Mass ##################################################################################################################################################################################################################################
# creating a a stat table for mass 1
Mass_1_stat <- band.df %>%
  summarise(mass = sum(!is.na(Mass_1)),
            avg.mass = round (mean(Mass_1, 
                                 na.rm = mass),2),
            stdev.mass = round (sd(Mass_1, 
                                 na.rm = mass), 2),
            med.mass = round(median(Mass_1,
                                  na.rm = mass), 2),
            min.mass = round(min(Mass_1,
                               na.rm = mass), 2),
            max.mass = round(max(Mass_1,
                               na.rm = mass), 2))
# plotting mass 1
hist(band.df$Mass_1, breaks = 12, xlim=c(15,21), ylim=c(0,12), main="Frequency of Mass (g) Before Treatment", xlab = "Mass (g)", col = "turquoise")

# creating a a stat table for mass 2
Mass_2_stat <- band.df %>%
  summarise(mass = sum(!is.na(Mass_2)),
            avg.mass = round (mean(Mass_2, 
                                   na.rm = mass),2),
            stdev.mass = round (sd(Mass_2, 
                                   na.rm = mass), 2),
            med.mass = round(median(Mass_2,
                                    na.rm = mass), 2),
            min.mass = round(min(Mass_2,
                                 na.rm = mass), 2),
            max.mass = round(max(Mass_2,
                                 na.rm = mass), 2))
# plotting mass 2
hist(band.df$Mass_2, breaks = 8, xlim=c(14,21), ylim=c(0,10), main= "Frequency of Mass (g) After Treatment", xlab= " Mass (g)", col= "orange")

# change in mass stat table
Change_mass_stat <- band.df %>%
  summarise(change_mass = sum(!is.na(change_in_mass)),
            avg.change_mass = round (mean(change_in_mass, 
                                   na.rm = change_mass),2),
            stdev.change_mass = round (sd(change_in_mass, 
                                   na.rm = change_mass), 2),
            med.change_mass = round(median(change_in_mass,
                                    na.rm = change_mass), 2),
            min.change_mass = round(min(change_in_mass,
                                 na.rm = change_mass), 2),
            max.change_mass = round(max(change_in_mass,
                                 na.rm = change_mass), 2))
# plotting the change in mass
hist(band.df$change_in_mass, breaks= 4, xlim=c(-4,4), ylim= c(0,15), main="Frequency of Change in Mass Throughout Treatment", xlab ="Change in Mass (g)", col = "lavender")

####### Tag / No Tag ##############################################################################################
# creating a graph of tag: no need for stat table because the variable is factored
plot(band.df$tag_no_tag, main = "Frequency of Individuals with a Tag (Y) or No Tag (N) Throughout Treatment", ylab = "Frequency", xlab= "No Tag (N) or Tag (Y)", col= "lavender")


####### Brood Treatment #############################################################################################################
plot(band.df$BroodTrt, main = "Frequency of Individuals with Brood Treatment Enlarged (E) or Reduced (R) Throughout Treatment", ylab = "Frequency", xlab= "Brood Treatment Enlarged (E) or Reduced (R)", col= "lavender")


############################## Bivariate Statistics ##################################################################################
############################# Associations Between Variables #########################################################################
## a1) Unadjusted model: Tag
# Blood glucose levels before treatment and tag 
T.gluc1.tag <- lm(band.df$BloodGluc_1 ~ band.df$tag_no_tag, 
                  data = subset(band.df,
                                tag_no_tag == 'Y' & 
                                  !is.na(x = band.df)))

summary(T.gluc1.tag) # print model summary (ln scale)
confint(T.gluc1.tag) # 95% CIs (ln scale)
plot(T.gluc1.tag) # view fitted vs residuals

Anova(T.gluc1.tag, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.gluc1.tag.mmean <- emmeans(T.gluc1.tag, 
                             'tag_no_tag')
summary(T.gluc1.tag.mmean)

## a2) Unadjusted model: Tag
# Blood glucose levels after treatment and tag 
T.gluc2.tag <- lm(band.df$BloodGluc_2 ~ band.df$tag_no_tag, 
                  data = subset(band.df,
                                tag_no_tag == 'Y' & 
                                  !is.na(x = band.df)))

summary(T.gluc2.tag) # print model summary (ln scale)
confint(T.gluc2.tag) # 95% CIs (ln scale)
plot(T.gluc2.tag) # view fitted vs residuals

Anova(T.gluc2.tag, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.gluc2.tag.mmean <- emmeans(T.gluc2.tag, 
                             'tag_no_tag')
summary(T.gluc2.tag.mmean)

## a3) Unadjusted model: Tag
# Change in blood glucose levels throughout treatment and tag 
T.chang.gluc.tag <- lm(band.df$change_in_blood_gluc ~ band.df$tag_no_tag, 
                       data = subset(band.df,
                                     tag_no_tag == 'Y' & 
                                       !is.na(x = band.df)))

summary(T.chang.gluc.tag) # print model summary (ln scale)
confint(T.chang.gluc.tag) # 95% CIs (ln scale)
plot(T.chang.gluc.tag) # view fitted vs residuals

Anova(T.chang.gluc.tag, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.chang.gluc.tag.mmean <- emmeans(T.chang.gluc.tag, 
                                  'tag_no_tag')
summary(T.chang.gluc.tag.mmean)


#b1) Unadjusted model: Mass
# Change in blood glucose levels throughout treatment and mass before treatment
T.chang.gluc.mass1 <- lm(band.df$change_in_blood_gluc ~ band.df$Mass_1, 
                 data = subset(band.df,
                               Mass_1 == 'Y' & 
                                 !is.na(x = band.df)))

summary(T.chang.gluc.mass1) # print model summary (ln scale)
confint(T.chang.gluc.mass1) # 95% CIs (ln scale)
plot(T.chang.gluc.mass1) # view fitted vs residuals

Anova(T.chang.gluc.mass1, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.chang.gluc.mass1.mmean <- emmeans(T.chang.gluc.mass1, 
                            'Mass_1')
summary(T.chang.gluc.mass1.mmean)

#b1) Unadjusted model: Mass
# Change in blood glucose levels throughout treatment and mass after treatment
T.chang.gluc.mass2 <- lm(band.df$change_in_blood_gluc ~ band.df$Mass_2, 
                         data = subset(band.df,
                                       Mass_2 == 'Y' & 
                                         !is.na(x = band.df)))

summary(T.chang.gluc.mass2) # print model summary (ln scale)
confint(T.chang.gluc.mass2) # 95% CIs (ln scale)
plot(T.chang.gluc.mass2) # view fitted vs residuals

Anova(T.chang.gluc.mass2, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.chang.gluc.mass2.mmean <- emmeans(T.chang.gluc.mass2, 
                                    'Mass_2')
summary(T.chang.gluc.mass2.mmean)

#b1) Unadjusted model: Mass
# Change in blood glucose levels throughout treatment and change in mass throughout treatment
T.chang.gluc.chang.mass <- lm(band.df$change_in_blood_gluc ~ band.df$change_in_mass, 
                         data = subset(band.df,
                                       change_in_mass == 'Y' & 
                                         !is.na(x = band.df)))

summary(T.chang.gluc.chang.mass) # print model summary (ln scale)
confint(T.chang.gluc.chang.mass) # 95% CIs (ln scale)
plot(T.chang.gluc.chang.mass) # view fitted vs residuals

Anova(T.chang.gluc.chang.mass, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.chang.gluc.chang.mass.mmean <- emmeans(T.chang.gluc.chang.mass, 
                                    'change_in_mass')
summary(T.chang.gluc.chang.mass.mmean)


#c1) Unadjusted model: Brood Treatment
# Blood glucose levels before treatment and brood treatment
T.gluc1.brdtrt <- lm(band.df$BloodGluc_1 ~ band.df$BroodTrt, 
                         data = subset(band.df,
                                       BroodTrt == 'Y' & 
                                         !is.na(x = band.df)))

summary(T.gluc1.brdtrt) # print model summary (ln scale)
confint(T.gluc1.brdtrt) # 95% CIs (ln scale)
plot(T.gluc1.brdtrt) # view fitted vs residuals

Anova(T.gluc1.brdtrt, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.gluc1.brdtrt.mmean <- emmeans(T.gluc1.brdtrt, 
                                    'BroodTrt')
summary(T.gluc1.brdtrt.mmean)

#c2) Unadjusted model: Brood Treatment
# Blood glucose levels after treatment and brood treatment
T.gluc2.brdtrt <- lm(band.df$BloodGluc_2 ~ band.df$BroodTrt, 
                     data = subset(band.df,
                                   BroodTrt == 'Y' & 
                                     !is.na(x = band.df)))

summary(T.gluc2.brdtrt) # print model summary (ln scale)
confint(T.gluc2.brdtrt) # 95% CIs (ln scale)
plot(T.gluc2.brdtrt) # view fitted vs residuals

Anova(T.gluc2.brdtrt, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.gluc2.brdtrt.mmean <- emmeans(T.gluc2.brdtrt, 
                                'BroodTrt')
summary(T.gluc2.brdtrt.mmean)

#c3) Unadjusted model: Brood Treatment
# Change in blood glucose levels throughout treatment and brood treatment
T.chang.gluc.brdtrt <- lm(band.df$change_in_blood_gluc ~ band.df$BroodTrt, 
                     data = subset(band.df,
                                   BroodTrt == 'Y' & 
                                     !is.na(x = band.df)))

summary(T.chang.gluc.brdtrt) # print model summary (ln scale)
confint(T.chang.gluc.brdtrt) # 95% CIs (ln scale)
plot(T.chang.gluc.brdtrt) # view fitted vs residuals

Anova(T.chang.gluc.brdtrt, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.chang.gluc.brdtrt.mmean <- emmeans(T.chang.gluc.brdtrt, 
                                'BroodTrt')
summary(T.chang.gluc.brdtrt.mmean)

######################## Bivariate Plots ############################################################################################
# scatter plot using ggplot to compare change in mass (x) vs. change in blood glucose (y)
# plot without an influence from other factors using base R to construct a scatterplot

# using lm to find the intercept in order to plot regression line 
lm(band.df$change_in_mass ~ band.df$change_in_blood_gluc)
plot(band.df$change_in_mass, band.df$change_in_blood_gluc, main = " Change in Blood Glucose Levels (μg/dL) vs. Change in Mass (g) Throughout Treatment", xlab = "Change in Mass (g)", ylab = "Change in Blood Glucose (μg/dL)",
     abline( -0.894375,0.006594))


# plot with tag/no tag as the factored variable
cm_v_cbg_plot1 <- ggplot(band.df, aes(x= change_in_mass, y=change_in_blood_gluc, color= tag_no_tag))+ 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4,) +
 geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth = 5)
print(cm_v_cbg_plot1)

# no outlier here, added dots to see data distribution #

# plot with brood treatment as the factored variable
cm_v_cbg_plot2 <- ggplot(band.df, aes(x=change_in_mass, y=change_in_blood_gluc, color=BroodTrt)) +
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, bins = 20)
print(cm_v_cbg_plot2)

# there is an outlier here, added dots to see data distribution #

####################### Tables #######################################################################################################
##### Glucose ##############################################################################################################
# calculating the mean and standard deviation for the change in blood glucose
mean_chng_glucose <- mean(band.df$change_in_blood_gluc)
mean_chng_glucose

sd_chng_glucose <- sd(band.df$change_in_blood_gluc)
sd_chng_glucose

# calculating the mean and standard deviation for blood glucose 1 (initial, before treatment)
mean_glucose_1 <- mean(band.df$BloodGluc_1)
mean_glucose_1

sd_glucose_1 <- sd(band.df$BloodGluc_1)
sd_glucose_1
  
# calculating the mean and standard deviation for blood glucose 2 (final, after treatment)
mean_glucose_2 <- mean(band.df$BloodGluc_2)
mean_glucose_2

sd_glucose_2 <- sd(band.df$BloodGluc_2)
sd_glucose_2

#### Mass ###################################################################################################################
# calculating the mean and standard deviation for change in mass
mean_chng_mass <- mean(band.df$change_in_mass)
mean_chng_mass

sd_chng_mass <- sd(band.df$change_in_mass)
sd_chng_mass

# calculating the mean and standard deviation for mass 1 (initial mass, before treatment)
mean_mass_1 <- mean(band.df$Mass_1)
mean_mass_1

sd_mass_1 <- sd(band.df$Mass_1)
sd_mass_1

# calculating the mean and standard deviation for mass 2 (final mass, after treatment)
mean_mass_2 <- mean(band.df$Mass_2)
mean_mass_2

sd_mass_2 <- sd(band.df$Mass_2)
sd_mass_2
#### Brood Treatment ######################################################################################################
E <- band.df$BroodTrt == 'E'
sum(E)
# percentage calculation
20/37 
#### Tag ##################################################################################################################
TA <- band.df$tag_no_tag == 'Y'
sum(TA)
# percentage calculation
20/37
####################### Data Cleaning for Time #########################################################################################################
# accounting for time
# reformatting the time variable #
## Format Date_One date in wide_data    
wide_data <- wide_data %>%
  mutate(Date_One = as.Date(wide_data$Date_1,
                            format = '%m/%d/%y'))

## Format the times as 24hr 
wide_data$Time_One<-format(strptime(wide_data$Time_1, "%I:%M %p"), format="%H:%M:%S")

## create date time variable and store as a datetime class
wide_data$date_time_1 <- 
  as.POSIXct(paste(wide_data$Date_One,
                   wide_data$Time_One), 
             format = '%Y-%m-%d %H:%M:%S')

# formatting Date_two in wide_data
wide_data <- wide_data %>%
  mutate(Date_Two = as.Date(wide_data$Date_2,
                            format = '%m/%d/%y'))

# format time 2 as 24 hr
wide_data$Time_Two<- format(strptime(wide_data$Time_2, "%I:%M %p"), format = "%H:%M%S")

################### Formatting Hours and Univariate Stats for Time ##############################################################
# create a date time variable and store it as a datetime class
wide_data$date_time_2 <- 
  as.POSIXct(paste(wide_data$Date_Two,
                   wide_data$Time_Two),
             format = '%Y-%m-%d %H:%M%S')

# extract those from the wide_data frame and add them to the band.df
band.df <- band.df %>% mutate(wide_data$date_time_1)
band.df <- band.df %>% rename(date_time_1 = 'wide_data$date_time_1')

band.df <- band.df %>% mutate(wide_data$date_time_2)
band.df <- band.df %>% rename(date_time_2 = 'wide_data$date_time_2')

# using lubridate to use to convert time
## create date time variable and store as a datetime class (wide_data for before treatment)
wide_data$date_time_1 <- 
  as.POSIXct(paste(wide_data$Date_One,
                   wide_data$Time_One), 
             format = '%Y-%m-%d %H:%M:%S')

wide_data$hour_1 <- hour(wide_data$date_time_1)

wide_data <- wide_data %>%
  mutate(day_night_1 = ifelse(hour_1 >= 5.45 & hour_1 <= 8.22, 'day', 'night'))

## creating a datetime variable for the second set of times and dates and store as a datetime class (band.df for after treatment)
band.df$date_time_2 <-
  as.POSIXct(paste(band.df$Date_Two,
                   band.df$Time_Two),
             format = '%Y-%m-%d %H:%M:%S')

band.df$hour_2 <- hour(band.df$date_time_2)

band.df <- band.df %>%
  mutate(day_night_2 = ifelse(hour_2 >= 5.45 & hour_2 <= 8.22, 'day', 'night'))

# taking the sin of the hour data and storing that as a variable
sin_hour_1 <- sin(wide_data$hour_1)
sin_hour_2 <- sin(band.df$hour_2)


#################### Circular plots for time variable (Univariate) #################################################
# plotting the hour data on a timeclock using the R package circular
# coercing the data into a circular dataframe

# plotting the circular dataframe onto a graph for before treatment (accounting for sin)
plot_sin_hour_1 <- as.data.frame.circular(sin_hour_1)
plot.circular(plot_sin_hour_1, col= "dodgerblue1", main = "Time of Day That Data Samples Were Collected Before Treatment")

# plotting the circular dataframe onto a graph for after treatment (accounting for sin)
plot_sin_hour_2 <- as.data.frame.circular(sin_hour_2)
plot.circular(plot_sin_hour_2, col= "orange", main = "Time of Day That Data Samples Were Collected After Treatment")

# plotting the circular dataframe without accounting for sin for hour_1 (before trt)
plot_hour_1 <- as.data.frame.circular(wide_data$hour_1)
plot.circular(plot_hour_1, col= "dodgerblue1", main = "Time of Day That Data Samples Were Collected Before Treatment")

# plotting the circular dataframe without accounting for sin for hour_2 (after trt)
plot_hour_2 <- as.data.frame.circular(band.df$hour_2)
plot.circular(plot_hour_1, col= "orange", main = "Time of Day That Data Samples Were Collected After Treatment")


##################### Bivariate Plots for Time  ##################################################################################################
# plotting inital blood glucose (before trt) vs.  initial time of day (before trt)
# saving hour_1 to band.df 
band.df <- band.df %>% mutate(wide_data$hour_1)
band.df <- band.df %>% rename(hour_1 = 'wide_data$hour_1')

# converting day_night_1 to a factor
as.factor(wide_data$day_night_1)
band.df <- band.df %>% mutate(wide_data$day_night_1)
band.df <- band.df %>% rename(day_night_1 = 'wide_data$day_night_1')

# converting day_night_2 to a factor as well
as.factor(band.df$day_night_2)

# plotting blood glucose (before treatment) vs. time (before treatment) bivariate
bg_bt_tod_plot <- ggplot(band.df, aes(x= BloodGluc_1, color= day_night_1))+ 
  geom_boxplot()
print(bg_bt_tod_plot)

# plotting blood glucose (after treatment) vs. time (after treatment) bivariate
bg_at_tod_plot <- ggplot(band.df, aes(x= BloodGluc_2, color= day_night_2))+ 
  geom_boxplot()
print(bg_at_tod_plot) 

# plotting blood glucose 1 (before trt) vs. mass 1 (before trt) vs. time of day 1 before experiment
bg_m_tod_plot1 <- ggplot(band.df, aes(x= Mass_1, y=BloodGluc_1, color= day_night_1))+ 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth = 3)
print(bg_m_tod_plot1)

# plotting blood glucose 2 (after trt) vs mass 2 (after trt) vs time of day 2 after expeirment
bg_m_tod_plot2 <- ggplot(band.df, aes(x = Mass_2 , y = BloodGluc_2, color = day_night_2)) +
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1, binwidth = 4)
print(bg_m_tod_plot2)

# saving sin hour 1 and 2 to the band.df
band.df <- band.df %>% mutate(sin_hour_1)
band.df <- band.df %>% mutate(sin_hour_2)

# graphing the radian hours and blood glucose (rough plot) (before trt)
plot(band.df$BloodGluc_1, sin_hour_1)
# graphing the radian hours and blood glucose (rough plot) (after trt)
plot(band.df$BloodGluc_2, sin_hour_2)


# detailed plot of radian hours for Blood_Gluc_1 (before trt)
ggplot(data = band.df, aes(x = BloodGluc_1, y = sin_hour_1)) +
  geom_point(col= 'dodgerblue') +
  ggtitle("Initial Blood Glucose (Before Treatment) vs. Time of Day Before Treatment") +
  xlab("Initial Blood Glucose (Before Treatment") +
  ylab("Time of Day (Radian Hours)") +
  theme(plot.title = element_text(hjust = 0.5)) 

# detailed plot of radian hours for Blood_Gluc_2 (after trt)
ggplot(data = band.df, aes(x = BloodGluc_2, y = sin_hour_2)) +
  geom_point(col= 'orange') +
  ggtitle("Final Blood Glucose (After Treatment) vs. Time of Day After Treatment") +
  xlab("Final Blood Glucose (After Treatment)") +
  ylab("Time of Day (Radian Hours)") +
  theme(plot.title = element_text(hjust = 0.5)) 

################### Univariate Stats for Time of Day as Day / Night ###################################################################################################################################
Sin_Hour_1_stat <- band.df %>%
  summarise(BG = sum(!is.na(sin_hour_1)),
            avg.BG = round (mean(sin_hour_1, 
                                 na.rm = BG),2),
            stdev.BG = round (sd(sin_hour_1, 
                                 na.rm = BG), 2),
            med.BG = round(median(sin_hour_1,
                                  na.rm = BG), 2),
            min.BG = round(min(sin_hour_1,
                               na.rm = BG), 2),
            max.BG = round(max(sin_hour_1,
                               na.rm = BG), 2))
# plotting sin_hour_1
hist(band.df$sin_hour_1, breaks = 5, xlim = c(-1.5, 1.5), ylim =c(0,25), main = "Frequency of Radian Hours Before Treatment", xlab = "Radian Hours Before Treatment", col = "turquoise")

Sin_Hour_2_stat <- band.df %>%
summarise(BG = sum(!is.na(sin_hour_2)),
          avg.BG = round (mean(sin_hour_2, 
                               na.rm = BG),2),
          stdev.BG = round (sd(sin_hour_2, 
                               na.rm = BG), 2),
          med.BG = round(median(sin_hour_2,
                                na.rm = BG), 2),
          min.BG = round(min(sin_hour_2,
                             na.rm = BG), 2),
          max.BG = round(max(sin_hour_2,
                             na.rm = BG), 2))
# plotting sin_hour_2
hist(band.df$sin_hour_2, breaks = 5, xlim = c(-1.5, 1.5), ylim =c(0,25), main = "Frequency of Radian Hours After Treatment", xlab = "Radian Hours After Treatment", col = "orange")

###################### Assocaitions Between the Time Variable and Other Variables ################################################################################################################
# t1) Unadjusted model: Time
# Blood glucose levels before treatment and time
T.gluc1.time <- lm(band.df$BloodGluc_1 ~ band.df$day_night_1, 
                  data = subset(band.df,
                                day_night_1 == 'day' & 
                                  !is.na(x = band.df)))

summary(T.gluc1.time) # print model summary (ln scale)
confint(T.gluc1.time) # 95% CIs (ln scale)
plot(T.gluc1.time) # view fitted vs residuals

Anova(T.gluc1.time, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.gluc1.time.mmean <- emmeans(T.gluc1.time, 
                             'day_night_1')
summary(T.gluc1.time.mmean)

# t2) Unadjusted model: Time
# Blood glucose levels after treatment and time
T.gluc2.time2 <- lm(band.df$BloodGluc_2 ~ band.df$day_night_2, 
                   data = subset(band.df,
                                 day_night_2 == 'day' & 
                                   !is.na(x = band.df)))

summary(T.gluc2.time2) # print model summary (ln scale)
confint(T.gluc2.time2) # 95% CIs (ln scale)
plot(T.gluc2.time2) # view fitted vs residuals

Anova(T.gluc2.time2, type = 'II') # type II SS from Car package

# Use emmeans to estimate marginal means
T.gluc2.time2.mmean <- emmeans(T.gluc2.time2, 
                              'day_night_2')
summary(T.gluc2.time2.mmean)

# t3) Change in Blood Glucose and Time of Day 


################### Table Time #########################################################################################################################################################################################
# sin hour 1 mean and sd
mean(band.df$sin_hour_1)
sd(band.df$sin_hour_2)

# sin hour 2 mean and sd
mean(band.df$sin_hour_2)
sd(band.df$sin_hour_2)
######################## Fitting the Model #############################################################################################################################################################################
#### ΔBlood Glucose ~ tag + broodSize + tag*broodSize + mass1 + year + time +(1|site)
# adding site to band.df
band.df <- band.df %>% mutate(wide_data$Site)
band.df <- band.df %>% rename(Site = 'wide_data$Site')

# adding brood size 1 and 2 to band.df
band.df <- band.df %>% mutate(wide_data$BroodSize_1)
band.df <- band.df %>% rename(BroodSize_1 = 'wide_data$BroodSize_1')

band.df <- band.df %>% mutate(wide_data$BroodSize_2)
band.df <- band.df %>% rename(BroodSize_2 = 'wide_data$BroodSize_2')

# creating a change in brood size column
change_BroodSize <-band.df$BroodSize_2 - band.df$BroodSize_1
band.df <- band.df %>% mutate(change_BroodSize)

# creating a change in time (radian hours) column
change_time_rad <- band.df$sin_hour_2 - band.df$sin_hour_1
band.df <- band.df %>% mutate(change_time_rad)

# adding year to band.df by extracting the date then the year
wide_data$Date <- 
  as.POSIXct(paste(wide_data$Date_One,
                   wide_data$Date_One), 
             format = '%Y-%m-%d')

wide_data$Date <- Date(wide_data$Date)

wide_data$year <- 
  as.POSIXct(paste(wide_data$Date), 
             format = '%Y')

wide_data$year <- year(wide_data$year)

# adding year to band.df
band.df <- band.df %>% mutate(wide_data$year)
band.df <- band.df %>% rename(year = 'wide_data$year')



########################### Mixed Model #######################################################################################
# model 1- brood size is inital before BSM manipulation, time is radians for inital measurement
Model <- lmer(change_in_blood_gluc ~ tag_no_tag + BroodSize_1 + tag_no_tag*BroodSize_1 + Mass_1 + year + sin_hour_1 +
(1|Site), data= band.df)
print(Model)

# model 2- brood size is final after BSM manipulation, time is radians for inital measurement
Model_2 <- lmer(change_in_blood_gluc ~ tag_no_tag + BroodSize_2 + tag_no_tag*BroodSize_2 + Mass_1 + year + sin_hour_1 + 
(1|Site), data= band.df)
print(Model_2)

# model 3- brood size is change in brood size throughout experiment, time is radians for inital measurement
Model_3 <- lmer(change_in_blood_gluc ~ tag_no_tag + change_BroodSize + tag_no_tag*change_BroodSize + Mass_1 + year + sin_hour_1 + 
                  (1|Site), data= band.df)
print(Model_3)

# model 4- brood size is change in brood size throughout experiment, time is in radians for change in time throughout experiment
Model_4 <- lmer(change_in_blood_gluc ~ tag_no_tag + change_BroodSize + tag_no_tag*chang_BroodSize + Mass_1 + year + change_time_rad + 
                  (1|Site), data= band.df)
print(Model_4)

# model 5- brood treatment (E,R), time is in radians for initial measurement
Model_5 <-lmer(change_in_blood_gluc ~ tag_no_tag + BroodTrt + tag_no_tag*BroodTrt + Mass_1 + year + sin_hour_1 +
                 (1|Site), data = band.df)
print(Model_5)

################################ Residuals ######################################################################################################################3
# creating a residual model of blood_gluc_1 vs. time in order to extrapolate the residuals
# testing this with blood glucose 1 (before trt) and sin_hour_1 (before trt)
resid_model <- lm(BloodGluc_1 ~ sin_hour_1, data= band.df)
Time_1_blood_gluc <- residuals(resid_model)


# testing this with blood glucose 2 (after trt) and sin_hour_2 (after trt)
resid_model_2 <- lm(BloodGluc_2 ~ sin_hour_2, data = band.df)
Time_2_blood_gluc <- residuals(resid_model_2)

# creating a dataframe for the residuals only with band number
resid.df <- data.frame(band.df$Band_Number)
resid.df <- resid.df %>% rename(Band_Number = `band.df.Band_Number` )

# adding Time_1_blood_gluc and Time_2_blood_gluc to that dataframe
resid.df <- resid.df %>% mutate(Time_1_blood_gluc)
resid.df <- resid.df %>% mutate(Time_2_blood_gluc)

# calculating change in glucose (with time accounted for already)
change_time_blood_gluc <- Time_1_blood_gluc - Time_2_blood_gluc

# adding change_time_blood_gluc to resid.df and band.df
resid.df <- resid.df %>% mutate(change_time_blood_gluc)
band.df <- band.df %>% mutate(change_time_blood_gluc)

# residual mixed model
factor_year<- as.factor(band.df$year)
band.df <- band.df %>% mutate(factor_year)

# checking to make sure year is a factor
class(band.df$factor_year)

# fitting the residual model
residual_model <- lmer(change_time_blood_gluc ~ tag_no_tag + BroodTrt + Mass_1 + factor_year +
                                   (1|Site), data = band.df)
print(residual_model)
summary(residual_model)

# calculating the confidence intervals 
confint(residual_model)

# deleting the year column and just keeping the factor_year column
band.df <- select(band.df, -year)

######################################### Final Plots ##############################################################################################
# creating the final plots for the poster

# Use emmeans to estimate marginal means for change in blood gluc for 2019-2020
residual_model_emmean <- emmeans(residual_model, 
                                      'change_in_blood_gluc')
summary(residual_model)


# Graph of change_in_blood_gluc vs. year   
## a) create a new dataframe of means and SD for graphing
mean.sd.year <- data.frame(term = c('2019',
                                             '2020'),
                                    estimate = c(32.1604, 13.7279),
                                    # Note forcing dwplot to graph SD by specificy the upper and
                                    # lower bounds instead of automatically calculating actual
                                    # 95% CI from SE
                                    conf.low = c((6.426533 + 0)),
                                    conf.high = c(( 57.5961 + 0)))


## b) Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
dwplot(mean.sd.year,
       vline = geom_vline(xintercept = 0, colour = 'gray20', 
                          linetype = 2), # line at zero behind coefs
       dot_args = list(size = 3),
       whisker_args = list(size = 1),
       dodge_size = 1) +
  labs(title = 'Difference in Means for Change in Blood Glucose ',
       subtitle = ('')) +
  theme(plot.title = element_text(hjust = 0.5)) + # center title
  theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
  # bold and size title and axes labels
  theme(text = element_text(size=12, face =  'bold')) +
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'lightgrey', 
                             size = 1, linetype = 'solid')) + 
  theme(axis.text.x = element_text(face='bold', color='black', 
                                   size=18, angle=0),
        axis.text.y = element_text(face='bold', color='black', 
                                   size=18, angle=0))  +
  scale_color_grey (start = 0, end = 0) + # make color estimates black
  # color the dot and whiskers
  scale_color_manual(values=c('red3')) + 
  xlab('Change in Blood Glucose') +
  ylab('Year') +
  xlim(c(0,100))


mean.sd.tag <- data.frame(term = c('Y',
                                      'N'),
estimate = c(-7.5111,13.4994),
# Note forcing dwplot to graph SD by specificy the upper and
# lower bounds instead of automatically calculating actual
# 95% CI from SE
conf.low = c((-32.807256 + 0)),
conf.high = c(( 17.48130 + 0)))


## b) Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
dwplot(mean.sd.tag,
       vline = geom_vline(xintercept = 0, colour = 'gray20', 
                          linetype = 2), # line at zero behind coefs
       dot_args = list(size = 3),
       whisker_args = list(size = 1),
       dodge_size = 1) +
  labs(title = 'Difference in Means for Change in Blood Glucose with or without GPS tag',
       subtitle = ('')) +
  theme(plot.title = element_text(hjust = 0.5)) + # center title
  theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
  # bold and size title and axes labels
  theme(text = element_text(size=12, face =  'bold')) +
  # remove background color
  theme(panel.background = element_rect(fill = 'white')) +
  # add major axes
  theme(axis.line = element_line(colour = 'lightgrey', 
                                 size = 1, linetype = 'solid')) + 
  theme(axis.text.x = element_text(face='bold', color='black', 
                                   size=18, angle=0),
        axis.text.y = element_text(face='bold', color='black', 
                                   size=18, angle=0))  +
  scale_color_grey (start = 0, end = 0) + # make color estimates black
  # color the dot and whiskers
  scale_color_manual(values=c('red3')) + 
  xlab('Change in Blood Glucose') +
  ylab('Tag') +
  xlim(-40,40)




# tag vs. blood glucose
tag_final <- ggplot(band.df, aes(tag_no_tag, change_time_blood_gluc, fill= tag_no_tag)) + 
geom_boxplot(filled.contour( x= tag_no_tag, y= change_time_blood_gluc, color.palette = 'Set2' ))

tag_final +  theme(legend.position = "none")+
  theme(axis.line = element_line(colour = 'black', size = 0.5, linetype = 'solid')) + 
  theme(axis.text.x = element_text(color='black', size= 12, angle=0)) + 
  theme(axis.text.y = element_text(color='black', size=12, angle=0)) +
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x=element_blank()) 



# BSM vs. blood glucose
brood_final <- ggplot(band.df, aes(BroodTrt, change_time_blood_gluc, fill= BroodTrt)) + 
geom_boxplot()
  
brood_final +  theme(legend.position = "none")+
  theme(axis.line = element_line(colour = 'black', size = 0.5, linetype = 'solid')) + 
  theme(axis.text.x = element_text(color='black', size= 12, angle=0)) + 
  theme(axis.text.y = element_text(color='black', size=12, angle=0)) +
  theme(axis.title.y = element_blank())+ 
  theme(axis.title.x=element_blank()) 



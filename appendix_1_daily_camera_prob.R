##########################################################
#
# Code to determine the daily probability that a camera trap is active
# from the urban biodiversity monitoring project
#
# Written by Mason Fidino 
#
#

### Note: Although camera trap seasons are 28 days long, cameras can be put up
###       slightly before that window

# load packages, download them if you do not have them

package_load<-function(packages = c("dplyr", "magrittr", "curl",
                                    "data.table"), quiet=TRUE, 
                       verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

package_load()

# We need to download 2 tables (all_possible_dates and check_log)
# that contains daily information on when a camera trap
# was active and if it malfunctioned. Reading the data from 
# the github repository for this project. Therefore,
# this script requires you to be online.

##############
# Note: This is getting saved in your working directory, slightly different methods
#       for mac and linux users. I did not test this on other systems, so I am not
#       100% certain of the code for them for downloading files.
##############

if(Sys.info()[['sysname']] == "Windows"){
  # download all_possible_dates
download.file("https://raw.github.com/mfidino/Fidino_Simonis_Magle_Bayesian_metacom/master/all_possible_dates.txt",
                      destfile = "all_possible_dates.txt", method = "wininet")
  # download check_log
download.file("https://raw.github.com/mfidino/Fidino_Simonis_Magle_Bayesian_metacom/master/check_log.txt",
                destfile = "check_log.txt", method = "wininet")
print(paste("files saved in", getwd(), sep = " "))
}else{
  download.file("https://raw.github.com/mfidino/Fidino_Simonis_Magle_Bayesian_metacom/master/all_possible_dates.txt",
                destfile = "all_possible_dates.txt", method = "curl")
  download.file("https://raw.github.com/mfidino/Fidino_Simonis_Magle_Bayesian_metacom/master/check_log.txt",
                destfile = "check_log.txt", method = "curl")
  print(paste("files saved in", getwd(), sep = " "))
} 

# read in the downloaded data
all_possible_dates <- read.table("all_possible_dates.txt", header = TRUE, sep = "\t")
check_log <- read.table("check_log.txt", header = TRUE, sep = "\t")

# a vector of ones to aid with summation
all_possible_dates$ones <- 1

# 1 = camera active, NA = inactive
the_na <- which(is.na(all_possible_dates$active)) 

# make NA moments 0 on columns 'ones'
all_possible_dates$ones[the_na] <- 0 

# basically the opposite of vector 'ones'
all_possible_dates$dead <-  1 - all_possible_dates$ones 


# Using dplyr and pipes from magrittr here to determine which cameras were inactive
# the entire 28 day season. There are also sites that have either been retired
# or had not been started yet (28 inactive days). 
# Since we are interested in probability
# that an active camera trap is active per day we need to exclude these.

off_all_season <- all_possible_dates %>%
  select(SurveyID, ones, dead) %>% # select these three columns
  group_by(SurveyID) %>% #group by surveyid
  summarise(days_off = sum(dead)) %>% #sum dead for each surveyid
  filter(days_off==28) # get sites that were not active the whole season
  

# remove the sites not active the whole season
working_cam <- all_possible_dates[-which(all_possible_dates$SurveyID %in% off_all_season$SurveyID),] %>%
  select(SurveyID, ones, dead, Date, Season)


# from check_log, get camera surveyids that went inoperable at some point(2,4,and 5)
# in the CameraConditon column

in_operable <- check_log %>%
  select(SurveyID, CameraCondition, CheckDate) %>%
  filter(CameraCondition==2 | CameraCondition==4 | CameraCondition==5)

# make an inoperable category and fill it with 1s if a site was inoperable
# during that season
working_cam$inop <- 0
working_cam$inop[which(working_cam$SurveyID %in% in_operable$SurveyID)] <- 1


# determine the number of active and inactive days per site
# per season

# 
work_sum <- data.table(working_cam)[,rle(ones), by = SurveyID]


# change the rle of ones from work_sum into a survey_id by number of sample days
# matrix. Rows = site * season, columns = days
big_matrix <- matrix(unlist(mapply(rep, x = work_sum$values, each = work_sum$lengths )),
                     ncol = 28, byrow = TRUE)

# take average of each column

daily_prob <- apply(big_matrix, 2, mean)


####################################################
### step 2. assess fit of 28 vs 30 draws
####################################################

# get number of days each camera trap was active and summarize it
n_days <- table(apply(big_matrix, 1, sum))

# pull 803 random draws at the daily probabilities (get_count), do this 1000 times


### 28 draws
draw28 <- matrix(rbinom(803*1000, 28, daily_prob), ncol = 1000, nrow = 803)

# make a contingency table for each simulation
draw28_contingency <- apply(draw28, 2, table)

# make vector to store rmse for 28 draws
rmse28 <- numeric(ncol(draw28))

# blank contingency table with 29 possibilities (0 to 28 days) 
observed <- matrix(0, ncol = 29, nrow = 1)

# put the observed data in blank table
observed[,as.numeric(names(n_days))+1] <- as.numeric(n_days)

# calculate rmse for each simulation
for(i in 1:length(draw28_contingency)){

# blank contingency table for predicted data
predicted <- matrix(0, ncol = 29, nrow = 1)

# add simulated data to predicted
predicted[,as.numeric(names(draw28_contingency[[i]]))+1] <- as.numeric(draw28_contingency[[i]])

# calculate rmse and store it in rmse28
rmse28[i] <- sqrt(mean((predicted - observed)^2))
}


### do the same thing for 30 draws
# each column is 803 pulls
draw30 <- matrix(rbinom(803*1000, 30, c(daily_prob, 0.85, 0.85)), ncol = 1000, nrow = 803)
draw30[draw30>28] <- 28

# create contingency table for each simulation
draw30_contingency<- apply(draw30, 2, table)

# make vector to store rmse for 28 draws
rmse30 <- numeric(ncol(draw30))


# calculate rmse for each simulation
for(i in 1:length(draw30_contingency)){
  
  # blank contingency table for predicted data
  predicted <- matrix(0, ncol = 29, nrow = 1)
  
  # add simulated data to predicted
  predicted[,as.numeric(names(draw30_contingency[[i]]))+1] <- as.numeric(draw30_contingency[[i]])
  
  # calculate rmse and store it in rmse28
  rmse30[i] <- sqrt(mean((predicted - observed)^2))
}

# compare the 2 methods

quantile(rmse28, probs = c(0.025, 0.5, 0.975))

quantile(rmse30, probs = c(0.025, 0.5, 0.975))

# rmse 30 fits much better

###################################
# step 3. simulate j matrix
###################################


# function to simulate jmat with daily probabilities

sim_jmat <- function(nrep = 28, nsite = 80, nyear = 6, nspecies = 3){
    days_array <- array(0, dim = c(nrep, nsite, nyear))
    # the hard coded probability that a camera was active that day
    # calculated from our own data. 
    # note here that we actually sample for 30 days to increase
    # the number of camera traps that worked for 28 days, but then
    # constrain to 28, which is similar to our own camera trap stuff.
    camera_prob <- c(0.45, 0.72, 0.81, 0.84, 0.83, 0.81, 0.82, 0.81, 0.86, 0.88, 
                     0.88, 0.87, 0.87, 0.89, 0.92, 0.94, 0.93, 0.92, 0.90, 0.89,
                     0.89, 0.89, 0.89, 0.88, 0.88, 0.87, 0.86, 0.84, 0.85, 0.85)
    
    # figure out the number of days each camera trap was active
    n_samp <- rbinom(nsite * nyear, length(camera_prob), camera_prob)
    n_samp[n_samp>28] <- 28 # > 28 change to 28
    # randomly sample which days were active for a particular camera.
    # Right now this is done with equal probability.
    n_list <- lapply(n_samp, function(x) sort(c(1, sample(2:28, x-1))))
    
    to_start <- 0
    # for every time a carmera trap was active
    for(i in 1:length(n_list)){
      # put a 1 on the days active on days_array
      days_array[to_start+n_list[[i]]] <- 1
      to_start <- to_start+ nrep # keep trucking through the array
    }
    # reorganize array to line up with our other arrays
    days_array <- aperm(days_array, c(2,3,1))
    # rep days_array by number of species, for jags 
    jmat_expanded <- array(rep(days_array, each = nspec),
                           dim = c(nspec, nsite, nyear, nrep))
    
    # make the jmat 
    jmat <- array(0, dim = c(nspec, nsite, nyear))
    # sum through to give # of days active
    jmat <- apply(jmat_expanded, c(2,3), rowSums)
    
  return(jmat)
    
}






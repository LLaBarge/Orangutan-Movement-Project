# Lc project script for 2012 long calls
# script for creating AKDEs and then polygons for Long call project
library(move)
#devtools::install_github("ctmm-initiative/ctmm")
library(ctmm)
library(readr)
#setwd("O:/working/processed/Long Call Project - Script and data")
setwd("C:/Users/lrlab/OneDrive/Desktop/LC_Project")
#edit file that will go onto movebank to only have data from months
# with enough location pts/individual to include in
# occurrence distributions >50 locations/month
Tuanan_GPS_2012 <- read_csv("Tuanan GPS 2012.csv")
as.data.frame(Tuanan_GPS_2012)
Tuanan_GPS_2012$ID_Month<-paste(Tuanan_GPS_2012$id, Tuanan_GPS_2012$Month)
Tuanan_filtered<-subset(Tuanan_GPS_2012, with(Tuanan_GPS_2012, ID_Month %in% names(which(table(ID_Month)>=50))))
str(Tuanan_filtered)
nrow(Tuanan_GPS_2012)


#write csv and add to movebank manually for ctmm
write.csv(Tuanan_filtered, file="Tuanan GPS 2012_Corrected.txt")


# log into movebank for file
# movebank contains all the GPS data and a separate file
# named "Bornean Orangutan 2012" was batched edited to contain only that year
login <- movebankLogin(username="llabarge", password="Samango1992!")

orang_study <- getMovebankData(study="Bornean orangutan, 2011-2013", login=login)

# create as.telemetry object for ctmm models
orangs<-as.telemetry(orang_study, projection = "+init=epsg:32750 +proj=utm +zone=50 +units=m +south")

#run continuous-time movement models, select the best for each individual
# then use the best model to run autocorrelated KDEs with an error of 25m assigned
SVF <- list()
for(i in 1:length(orangs)){
  print(i)
  SVF[[i]] <- variogram(orangs[[i]])}
names(SVF) <- names(orangs)


## sampling duration/sample size per length of time studied too small for
## Bobo, Tina, Kentung, Fatih, Helium, Tuko, Manggo, Jordan, Pinky, Dayak, Preman
## fit all models together, but only use Chili, Niko, Teju, Tomi, Wodan for creating AKDEs

## fit models to all animals
# AKDEs for largest area to create available points
#occurrance distributions for risk of encounter
FIT <- list()
for(i in 1:length(orangs)){
  print(i)
  GUESS <- ctmm.guess(orangs[[i]],CTMM=ctmm(error = 15), interactive=FALSE)
  FIT[[i]] <- ctmm.select(orangs[[i]],GUESS, verbose=TRUE,trace=2)
}

# extent based on Tuanan map

Tuanan_r<-raster()
extent(Tuanan_r)<-extent(212812, 219348, 9758957, 9773670)
res(Tuanan_r)<-25
projection(Tuanan_r) <- "+init=epsg:32750 +proj=utm +zone=50 +units=m +south"


AKDE <- list()
for(i in 1:length(orangs)){
  print(i)
  AKDE[[i]] <- akde(orangs[[i]],FIT[[i]][[1]], grid=Tuanan_r)
}
names(AKDE) <- names(orangs)



# some individuals appear so infrequently that ODs are likely not representative
# of competition risk, these are:
#Cinta, Tina, Kentung, Fatih, Tuko, Manggo, Jordan, Pinky Luca, Dolay, Preman
# others were not present in the study area for more than 2 months
#in 2012: Bobo, Fayesh, Flumnmu, Fugit, Talia, Tuko

OD_2012 <- getMovebankData(study="Tuanan 2012 Monthly GPS, Bornean Orangutan, Central Kalimatan, Indonesia,", login=login)
as.data.frame(Tuanan_GPS_2012_Corrected)
str(Tuanan_GPS_2012_Corrected)
strptime(Tuanan_GPS_2012_Corrected$timestamp, format = "%Y-%m-%d %H:%M:%S")

Tuanan_GPS_2012_Corrected$timestamp<-as.POSIXlt.character(Tuanan_GPS_2012_Corrected$timestamp, format="%Y-%m-%d %H:%M")

# create as.telemetry object for ctmm models
OD_tel<-as.telemetry(Tuanan_GPS_2012_Corrected, projection = "+init=epsg:32750 +proj=utm +zone=50 +units=m +south")

#fit movement models to each OD telemetry

OD.FIT <- list()
for(i in 1:length(OD_tel)){
  print(i)
  OD.GUESS <- ctmm.guess(OD_tel[[i]], CTMM=ctmm(error = 15), interactive=FALSE)
  OD.FIT[[i]] <- ctmm.select(OD_tel[[i]],OD.GUESS, verbose=TRUE,trace=2)
}

names(OD.FIT) <- names(OD_tel)

OD <- list()
for(i in 1:length(OD_tel)){
  print(i)
  OD[[i]] <- occurrence(OD_tel[[i]],OD.FIT[[i]][[1]])
}
names(OD) <- names(OD_tel)


save(AKDE, file="akdes.Rdata")
save(FIT, file="modFIT.Rdata")
load(file="modFIT.Rdata")
load(file="akdes.Rdata")
load(file="OD.FIT.Rdata")
save(OD.FIT, file="OD.FIT.Rdata")
save(OD, file="OD.Rdata")
load(file="OD.Rdata")
View(OD)
# create a dataframe of GPS points for "available" locations - to use later

Avail_pts<- do.call(rbind, lapply(OD_tel, function(x) as.data.frame(x)))
Avail_pts$ID <- rownames(Avail_pts)
library(tm)
removeNumbers(Avail_pts$ID)
Avail_pts$ID<- gsub("\\..*","",Avail_pts$ID)
# select males of interest (those that were range resident in 2012)

Avail_pts<-Avail_pts[Avail_pts$ID %in% c("Chili", "Niko",  "Tomi", "Teju" , "Wodan"), ]
str(Avail_pts)
PT<-rep("Available",times=2510)
Avail_pts$PT<-PT


# export rasters with probability distribution function as a proxy for competition risk
# create separate per individual/year to keep track
library(raster)

#create individual raster files for competition risk / prob of
# encounter with potential receptive females
# use nearest OD distribution to the date of data collection
# individuals chosen by looking at Tuanan presence records

library(spatialEco)
library(climateStability)
ZeroRast<-Tuanan_r
vals <- NA
ZeroRast <- setValues(ZeroRast, vals)

Chili<-raster(AKDE[["Chili"]] , DF="PDF")
Chili_Own<-raster.transformation(Chili, trans = "norm", smin = 0, smax = 1)
sp::plot(Chili_Own)

Chili_rast7<-raster(OD[["Chili08"]] , DF="PDF")
Chili_rast7<-rescale0to1(Chili_rast7)
Chili_rast8<-raster(OD[["Chili08"]] , DF="PDF")
Chili_rast8<-rescale0to1(Chili_rast8)
Chili_rast9<-raster(OD[["Chili09"]] , DF="PDF")
Chili_rast9<-rescale0to1(Chili_rast9)
Chili_rast10<-raster(OD[["Chili10"]] , DF="PDF")
Chili_rast10<-rescale0to1(Chili_rast10)
Chili_rast11<-raster(OD[["Chili11"]] , DF="PDF")
Chili_rast11<-rescale0to1(Chili_rast11)
Chili_rast12<-raster(OD[["Chili11"]] , DF="PDF")
Chili_rast12<-rescale0to1(Chili_rast12)


Dayak_rast5<-raster(OD[["Dayak05"]], DF="PDF")
Dayak_rast5<-rescale0to1(Dayak_rast5)
Dayak_rast9<-raster(OD[["Dayak05"]], DF="PDF")
Dayak_rast9<-rescale0to1(Dayak_rast9)
Dayak_rast10<-raster(OD[["Dayak05"]], DF="PDF")
Dayak_rast10<-rescale0to1(Dayak_rast10)
Dayak_rast11<-raster(OD[["Dayak10"]], DF="PDF")
Dayak_rast11<-rescale0to1(Dayak_rast11)

Desy_rast<-raster(AKDE[["Desy"]] , DF="PDF")
Desy_rast<-raster.transformation(Desy_rast, trans = "norm", smin = 0, smax = 1)

Henk_rast3<-raster(OD[["Henk03"]], DF="PDF")
Henk_rast3<-rescale0to1(Henk_rast3)
Henk_rast10<-raster(OD[["Henk10"]], DF="PDF")
Henk_rast10<-rescale0to1(Henk_rast10)
Henk_rast12<-raster(OD[["Henk12"]], DF="PDF")
Henk_rast12<-rescale0to1(Henk_rast12)

Jinak_rast<-raster(AKDE[["Jinak"]] , DF="PDF")
Jinak_rast<-raster.transformation(Jinak_rast, trans = "norm", smin = 0, smax = 1)

Juni_rast<-raster(AKDE[["Juni"]] , DF="PDF")
Juni_rast<-raster.transformation(Juni_rast, trans = "norm", smin = 0, smax = 1)

Inul_rast<-raster(AKDE[["Inul"]] , DF="PDF")
Inul_rast<-raster.transformation(Inul_rast, trans = "norm", smin = 0, smax = 1)

Katmandun_rast4<-raster(OD[["Katmandun04"]] , DF="PDF")
Katmandun_rast4<-rescale0to1(Katmandun_rast4)
Katmandun_rast5<-raster(OD[["Katmandun04"]] , DF="PDF")
Katmandun_rast5<-rescale0to1(Katmandun_rast5)
Katmandun_rast6<-raster(OD[["Katmandun06"]] , DF="PDF")
Katmandun_rast6<-rescale0to1(Katmandun_rast6)

Kerry_rast<-raster(AKDE[["Kerry"]] , DF="PDF")
Kerry_rast<-raster.transformation(Kerry_rast, trans = "norm", smin = 0, smax = 1)

Kondor_rast<-raster(AKDE[["Kondor"]] , DF="PDF")
Kondor_rast<-raster.transformation(Kondor_rast, trans = "norm", smin = 0, smax = 1)

Mindy_rast<-raster(AKDE[["Mindy"]] , DF="PDF")
Mindy_rast<-raster.transformation(Mindy_rast, trans = "norm", smin = 0, smax = 1)

Niko_rast3<-raster(OD[["Niko.3"]] , DF="PDF")
Niko_rast3<-rescale0to1(Niko_rast3)
Niko_rast5<-raster(OD[["Niko.5"]] , DF="PDF")
Niko_rast5<-rescale0to1(Niko_rast5)
Niko_rast7<-raster(OD[["Niko.5"]] , DF="PDF")
Niko_rast7<-rescale0to1(Niko_rast7)
Niko_rast9<-raster(OD[["Niko.9"]] , DF="PDF")
Niko_rast9<-rescale0to1(Niko_rast9)

Otto_rast3<-raster(OD[["Otto.3"]] , DF="PDF")
Otto_rast3<-rescale0to1(Otto_rast3)
Otto_rast8<-raster(OD[["Otto.8"]] , DF="PDF")
Otto_rast8<-rescale0to1(Otto_rast8)

Teju_rast11<-raster(OD[["Teju.11"]] , DF="PDF")
Teju_rast11<-rescale0to1(Teju_rast11)
Teju_rast12<-raster(OD[["Teju.12"]] , DF="PDF")
Teju_rast12<-rescale0to1(Teju_rast12)

Tomi_rast3<-raster(OD[["Tomi.4"]] , DF="PDF")
Tomi_rast3<-rescale0to1(Tomi_rast3)
Tomi_rast4<-raster(OD[["Tomi.4"]] , DF="PDF")
Tomi_rast4<-rescale0to1(Tomi_rast4)
Tomi_rast5<-raster(OD[["Tomi.5"]] , DF="PDF")
Tomi_rast5<-rescale0to1(Tomi_rast5)
Tomi_rast6<-raster(OD[["Tomi.6"]] , DF="PDF")
Tomi_rast6<-rescale0to1(Tomi_rast6)
Tomi_rast7<-raster(OD[["Tomi.7"]] , DF="PDF")
Tomi_rast7<-rescale0to1(Tomi_rast7)
Tomi_rast8<-raster(OD[["Tomi.8"]] , DF="PDF")
Tomi_rast8<-rescale0to1(Tomi_rast8)
Tomi_rast9<-raster(OD[["Tomi.9"]] , DF="PDF")
Tomi_rast9<-rescale0to1(Tomi_rast9)
Tomi_rast10<-raster(OD[["Tomi.10"]] , DF="PDF")
Tomi_rast10<-rescale0to1(Tomi_rast10)
Tomi_rast11<-raster(OD[["Tomi.11"]] , DF="PDF")
Tomi_rast11<-rescale0to1(Tomi_rast11)
Tomi_rast12<-raster(OD[["Tomi.12"]] , DF="PDF")
Tomi_rast12<-rescale0to1(Tomi_rast12)


Wodan_rast3<-raster(OD[["Wodan.5"]] , DF="PDF")
Wodan_rast3<-rescale0to1(Wodan_rast3)
Wodan_rast4<-raster(OD[["Wodan.5"]] , DF="PDF")
Wodan_rast4<-rescale0to1(Wodan_rast4)
Wodan_rast5<-raster(OD[["Wodan.5"]] , DF="PDF")
Wodan_rast5<-rescale0to1(Wodan_rast5)
Wodan_rast6<-raster(OD[["Wodan.6"]] , DF="PDF")
Wodan_rast6<-rescale0to1(Wodan_rast6)
Wodan_rast7<-raster(OD[["Wodan.7"]] , DF="PDF")
Wodan_rast7<-rescale0to1(Wodan_rast7)
Wodan_rast8<-raster(OD[["Wodan.8"]] , DF="PDF")
Wodan_rast8<-rescale0to1(Wodan_rast8)
Wodan_rast9<-raster(OD[["Wodan.9"]] , DF="PDF")
Wodan_rast9<-rescale0to1(Wodan_rast9)
Wodan_rast10<-raster(OD[["Wodan.11"]] , DF="PDF")
Wodan_rast10<-rescale0to1(Wodan_rast10)
Wodan_rast11<-raster(OD[["Wodan.11"]] , DF="PDF")
Wodan_rast11<-rescale0to1(Wodan_rast11)
Wodan_rast12<-raster(OD[["Wodan.12"]] , DF="PDF")
Wodan_rast12<-rescale0to1(Wodan_rast12)


# get long call data on movebank
login <- movebankLogin(username="llabarge", password="Samango1992!")
##getMovebankStudies(login)

##searchMovebankStudies(x="Bornean", login=login)

LongCalls<-getMovebankData(study="Spontaneous Long Calls, Bornean Orangutan Flanged males 2012, Tuanan", login=login)

# format long call data so can be merged as "cases" along with the
# "control" available points
## of these "spontaneous" long calls, 1011 came from the PAM system
LC_tel<- as.telemetry(LongCalls, projection = "+init=epsg:32750 +proj=utm +zone=50 +south")


LC_pts<- do.call(rbind, lapply(LC_tel, function(x) as.data.frame(x)))
str(LC_pts)
PT<-rep("LongCall",times=1234)
LC_pts$PT<-PT

library(stringr)
LC_pts$ID <- rownames(LC_pts)
unique(LC_pts$ID)
library(tm)
removeNumbers(LC_pts$ID)
LC_pts$ID<- gsub("\\..*","",LC_pts$ID)
LC_pts<-LC_pts[!grepl("Helium", LC_pts$ID),]
LC_pts<-LC_pts[!grepl("Preman", LC_pts$ID),]
LC_pts<-LC_pts[!grepl("Dayak", LC_pts$ID),]
LC_points<-LC_pts[,c("x", "y", "PT", "timestamp", "ID")]
rownames(LC_points) <- NULL
as.data.frame(LC_points)


library(tidyverse)
library(lubridate)
# need to filter dataframe to remove duplicate pts or those within 25 minutes of one another
elapsed <- function(x)
{
  y <- abs(as.duration(x[2:length(x)] %--% x[1:(length(x)-1)]))
  y >= 25*60
}

LC_points %>%
  group_split(ID) %>%  
  map_dfr(~ .[c(T, if (nrow(.) > 1) elapsed(.$timestamp)),])
str(LC_points)
# removes around 200 points that were part of the same event leaving 1070
## this was because PAM system often recorded the same long calls as observation


# add a month column (for later indexing for raster extract)
# bind the long calls (cases) and available pts (controls) together

Avail_pts$Month <- str_sub(Avail_pts$timestamp, start = 6L, end = 7L)
LC_com<-LC_points%>% dplyr::select('ID', 'PT', 'Month', 'x', 'y', 'timestamp')

Av_com<-Avail_pts%>% dplyr::select('ID', 'PT', 'Month','x', 'y', 'timestamp')
LC_Control<-rbind(LC_com, Av_com)
as.data.frame(LC_Control)



# denote 1 for long call, 0 for available
library(plyr)
LC_Control$PT <- revalue(LC_Control$PT, c("LongCall"=1))
LC_Control$PT <- revalue(LC_Control$PT, c("Available"=0))
Chili_subset8<-LC_Control[LC_Control$ID == "Chili" & LC_Control$Month == "08", ]
Chili_subset9<-LC_Control[LC_Control$ID == "Chili" & LC_Control$Month == "09", ]
Chili_subset10<-LC_Control[LC_Control$ID == "Chili" & LC_Control$Month == "10", ]
Chili_subset11<-LC_Control[LC_Control$ID == "Chili" & LC_Control$Month == "11", ]
Niko_subset3<-LC_Control[LC_Control$ID == "Niko" & LC_Control$Month == "03", ]
Niko_subset5<-LC_Control[LC_Control$ID == "Niko" & LC_Control$Month == "05", ]
Niko_subset9<-LC_Control[LC_Control$ID == "Niko" & LC_Control$Month == "09", ]
Teju_subset11<-LC_Control[LC_Control$ID == "Teju" & LC_Control$Month == "11", ]
Teju_subset12<-LC_Control[LC_Control$ID == "Teju" & LC_Control$Month == "12", ]
Tomi_subset4<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "04", ]
Tomi_subset5<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "05", ]
Tomi_subset6<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "06", ]
Tomi_subset7<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "07", ]
Tomi_subset8<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "08", ]
Tomi_subset9<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "09", ]
Tomi_subset10<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "10", ]
Tomi_subset11<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "11", ]
Tomi_subset12<-LC_Control[LC_Control$ID == "Tomi" & LC_Control$Month == "12", ]
Wodan_subset5<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "05", ]
Wodan_subset6<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "06", ]
Wodan_subset7<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "07", ]
Wodan_subset8<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "08", ]
Wodan_subset9<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "09", ]
Wodan_subset11<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "11", ]
Wodan_subset12<-LC_Control[LC_Control$ID == "Wodan" & LC_Control$Month == "12", ]

# create monthly long call occurrence distributions by using long call data and creating an ID_date column
# download previously month-separated data from movebank and fit movement models and ODs

library(move)
library(ctmm)
LC_to_df<-getMovebankData(study="Long Calls, Bornean Orangutan, Pongo pygmaeus, Tuanan, Central Kalimatan, Indonesia", login=login)
LC_to_df<-as.telemetry(LC_to_df, projection = "+init=epsg:32750 +proj=utm +zone=50 +south")
LC_df<- do.call(rbind, lapply(LC_to_df, function(x) as.data.frame(x)))
head(LC_df)

library(stringr)
LC_df$ID <- rownames(LC_df)
library(tm)
removeNumbers(LC_df$ID)
LC_df$ID<- gsub("\\..*","",LC_df$ID)
LC_df$ID<- gsub("LC","",LC_df$ID)
LC_df$month <- format(LC_df$timestamp, "%m")
LC_df$year <- format(LC_df$timestamp, "%y")
rownames(LC_df)<-NULL

#column for ID_date to make monthly occurrence distributions
LC_df$ID_date<-gsub(" ", "", paste(LC_df$ID,"_", LC_df$month))
LC_<-subset(LC_df, year == "12")

LC_to_<-LC_[LC_$longitude >= "113.0" & LC_$longitude <= "115.0",]
# need to filter dataframe to remove duplicate pts or those within 25 minutes of one another
library(tidyverse)
library(lubridate)

elapsed <- function(x)
{
  y <- abs(as.duration(x[2:length(x)] %--% x[1:(length(x)-1)]))
  y >= 10*60
}

LC_to_movebank<-LC_to_ %>%
  group_split(ID_date) %>%  
  map_dfr(~ .[c(T, if (nrow(.) > 1) elapsed(.$timestamp)),])

# remove monthly data with less than 10 observations
LC_to_movebank<-LC_to_movebank[with(LC_to_movebank, ID_date %in% names(which(table(ID_date)>=10))), ]
str(LC_to_movebank)
#select only variables needed for Movebank
library(dplyr)
LC_to_movebank$Timestamp<-LC_to_movebank$timestamp
LC_to_movebank$id<-LC_to_movebank$ID_date

LC_movebank<-LC_to_movebank%>% dplyr::select('id', 'latitude', 'longitude', 'Timestamp')
rownames(LC_movebank) <- NULL
LongCalls_to_movebank<-as.data.frame(LC_movebank)
head(LongCalls_to_movebank)
# write txt to upload to movebank
write.csv(LongCalls_to_movebank, file="LongCalls_to_movebank.txt")

LongCalls_Monthly<-as.telemetry(object=LongCalls_to_movebank, projection = "+init=epsg:32750 +proj=utm +zone=50 +south")

# monthly long call occurrence distributions - represents competitive risk
# may be more accurrate than AKDEs/ODs of competitor risk
LC.FIT <- list()
for(i in 1:length(LongCalls_Monthly)){
  print(i)
  LC.GUESS <- ctmm.guess(LongCalls_Monthly[[i]], interactive=FALSE)
  LC.FIT[[i]] <- ctmm.select(LongCalls_Monthly[[i]],LC.GUESS, verbose=TRUE,trace=2)
}

names(LC.FIT) <- names(LongCalls_Monthly)

LC_Monthly <- list()
for(i in 1:length(LongCalls_Monthly)){
  print(i)
  LC_Monthly[[i]] <- occurrence(LongCalls_Monthly[[i]],LC.FIT[[i]][[1]])}
names(LC_Monthly) <- names(LongCalls_Monthly)

save(LC_Monthly, file="LC_Monthly.Rdata")



# create rasters of each of these
library(spatialEco)
View(LC_Monthly)

#March
Henk_03_LC<-raster(LC_Monthly[["Henk_03"]] , DF="PDF")
Henk_03_LC<-raster.transformation(Henk_03_LC, trans = "norm", smin = 0, smax = 1)
Niko_03_LC<-raster(LC_Monthly[["Niko_03"]] , DF="PDF")
Niko_03_LC<-raster.transformation(Niko_03_LC, trans = "norm", smin = 0, smax = 1)
Otto_03_LC<-raster(LC_Monthly[["Otto_03"]] , DF="PDF")
Otto_03_LC<-raster.transformation(Otto_03_LC, trans = "norm", smin = 0, smax = 1)
Tomi_03_LC<-raster(LC_Monthly[["Tomi_03"]] , DF="PDF")
Tomi_03_LC<-raster.transformation(Tomi_03_LC, trans = "norm", smin = 0, smax = 1)
Wodan_03_LC<-raster(LC_Monthly[["Wodan_03"]] , DF="PDF")
Wodan_03_LC<-raster.transformation(Wodan_03_LC, trans = "norm", smin = 0, smax = 1)

#April
Katmandun_04_LC<-raster(LC_Monthly[["Katmandun_04"]], DF="PDF")
Katmandun_04_LC<-raster.transformation(Katmandun_04_LC, trans = "norm", smin = 0, smax = 1)
Tomi_04_LC<-raster(LC_Monthly[["Tomi_04"]], DF="PDF")
Tomi_04_LC<-raster.transformation(Tomi_04_LC, trans = "norm", smin = 0, smax = 1)
Wodan_04_LC<-raster(LC_Monthly[["Wodan_04"]] , DF="PDF")
Wodan_04_LC<-raster.transformation(Wodan_04_LC, trans = "norm", smin = 0, smax = 1)

#May
Niko_05_LC<-raster(LC_Monthly[["Niko_05"]], DF="PDF")
Niko_05_LC<-raster.transformation(Niko_05_LC, trans = "norm", smin = 0, smax = 1)
Tomi_05_LC<-raster(LC_Monthly[["Tomi_05"]], DF="PDF")
Tomi_05_LC<-raster.transformation(Tomi_05_LC, trans = "norm", smin = 0, smax = 1)
Wodan_05_LC<-raster(LC_Monthly[["Wodan_05"]], DF="PDF")
Wodan_05_LC<-raster.transformation(Wodan_05_LC, trans = "norm", smin = 0, smax = 1)

#June
Dayak_06_LC<-raster(LC_Monthly[["Dayak_06"]], DF="PDF")
Dayak_06_LC<-raster.transformation(Dayak_06_LC, trans = "norm", smin = 0, smax = 1)
Henk_06_LC<-raster(LC_Monthly[["Henk_06"]], DF="PDF")
Henk_06_LC<-raster.transformation(Henk_06_LC, trans = "norm", smin = 0, smax = 1)
Katmandun_06_LC<-raster(LC_Monthly[["Katmandun_06"]], DF="PDF")
Katmandun_06_LC<-raster.transformation(Katmandun_06_LC, trans = "norm", smin = 0, smax = 1)
Wodan_06_LC<-raster(LC_Monthly[["Wodan_06"]], DF="PDF")
Wodan_06_LC<-raster.transformation(Wodan_06_LC, trans = "norm", smin = 0, smax = 1)


#July
Henk_07_LC<-raster(LC_Monthly[["Henk_07"]], DF="PDF")
Henk_07_LC<-raster.transformation(Henk_07_LC, trans = "norm", smin = 0, smax = 1)
Tomi_07_LC<-raster(LC_Monthly[["Tomi_07"]], DF="PDF")
Tomi_07_LC<-raster.transformation(Tomi_07_LC, trans = "norm", smin = 0, smax = 1)
Wodan_07_LC<-raster(LC_Monthly[["Wodan_07"]], DF="PDF")
Wodan_07_LC<-raster.transformation(Wodan_07_LC, trans = "norm", smin = 0, smax = 1)


#August
Chili_08_LC<-raster(LC_Monthly[["Chili_08"]] , DF="PDF")
Chili_08_LC<-raster.transformation(Chili_08_LC, trans = "norm", smin = 0, smax = 1)
Henk_08_LC<-raster(LC_Monthly[["Henk_08"]] , DF="PDF")
Henk_08_LC<-raster.transformation(Henk_08_LC, trans = "norm", smin = 0, smax = 1)
Niko_08_LC<-raster(LC_Monthly[["Niko_08"]] , DF="PDF")
Niko_08_LC<-raster.transformation(Niko_08_LC, trans = "norm", smin = 0, smax = 1)
Otto_08_LC<-raster(LC_Monthly[["Otto_08"]] , DF="PDF")
Otto_08_LC<-raster.transformation(Otto_08_LC, trans = "norm", smin = 0, smax = 1)
Teju_08_LC<-raster(LC_Monthly[["Teju_08"]], DF="PDF")
Teju_08_LC<-raster.transformation(Teju_08_LC, trans = "norm", smin = 0, smax = 1)
Wodan_08_LC<-raster(LC_Monthly[["Wodan_08"]], DF="PDF")
Wodan_08_LC<-raster.transformation(Wodan_08_LC, trans = "norm", smin = 0, smax = 1)

#september
Chili_09_LC<-raster(LC_Monthly[["Chili_09"]] , DF="PDF")
Chili_09_LC<-raster.transformation(Chili_09_LC, trans = "norm", smin = 0, smax = 1)
Dayak_09_LC<-raster(LC_Monthly[["Dayak_09"]], DF="PDF")
Dayak_09_LC<-raster.transformation(Dayak_09_LC, trans = "norm", smin = 0, smax = 1)
Helium_09_LC<-raster(LC_Monthly[["Helium_09"]], DF="PDF")
Helium_09_LC<-raster.transformation(Helium_09_LC, trans = "norm", smin = 0, smax = 1)
Henk_09_LC<-raster(LC_Monthly[["Henk_09"]], DF="PDF")
Henk_09_LC<-raster.transformation(Henk_09_LC, trans = "norm", smin = 0, smax = 1)
Niko_09_LC<-raster(LC_Monthly[["Niko_09"]], DF="PDF")
Niko_09_LC<-raster.transformation(Niko_09_LC, trans = "norm", smin = 0, smax = 1)
Tomi_09_LC<-raster(LC_Monthly[["Tomi_09"]], DF="PDF")
Tomi_09_LC<-raster.transformation(Tomi_09_LC, trans = "norm", smin = 0, smax = 1)
Wodan_09_LC<-raster(LC_Monthly[["Wodan_09"]], DF="PDF")
Wodan_09_LC<-raster.transformation(Wodan_09_LC, trans = "norm", smin = 0, smax = 1)

#October
Chili_10_LC<-raster(LC_Monthly[["Chili_10"]]  , DF="PDF")
Chili_10_LC<-raster.transformation(Chili_10_LC, trans = "norm", smin = 0, smax = 1)
Dayak_10_LC<-raster(LC_Monthly[["Dayak_10"]], DF="PDF")
Dayak_10_LC<-raster.transformation(Dayak_10_LC, trans = "norm", smin = 0, smax = 1)
Helium_10_LC<-raster(LC_Monthly[["Helium_10"]], DF="PDF")
Helium_10_LC<-raster.transformation(Helium_10_LC, trans = "norm", smin = 0, smax = 1)
Henk_10_LC<-raster(LC_Monthly[["Henk_10"]], DF="PDF")
Henk_10_LC<-raster.transformation(Henk_10_LC, trans = "norm", smin = 0, smax = 1)
Teju_10_LC<-raster(LC_Monthly[["Teju_10"]], DF="PDF")
Teju_10_LC<-raster.transformation(Teju_10_LC, trans = "norm", smin = 0, smax = 1)
Tomi_10_LC<-raster(LC_Monthly[["Tomi_10"]], DF="PDF")
Tomi_10_LC<-raster.transformation(Tomi_10_LC, trans = "norm", smin = 0, smax = 1)
Wodan_10_LC<-raster(LC_Monthly[["Wodan_10"]], DF="PDF")
Wodan_10_LC<-raster.transformation(Wodan_10_LC, trans = "norm", smin = 0, smax = 1)

#November
Chili_11_LC<-raster(LC_Monthly[["Chili_11"]]  , DF="PDF")
Chili_11_LC<-raster.transformation(Chili_11_LC, trans = "norm", smin = 0, smax = 1)
Dayak_11_LC<-raster(LC_Monthly[["Dayak_11"]]  , DF="PDF")
Dayak_11_LC<-raster.transformation(Dayak_11_LC, trans = "norm", smin = 0, smax = 1)
Helium_11_LC<-raster(LC_Monthly[["Helium_11"]]  , DF="PDF")
Helium_11_LC<-raster.transformation(Helium_11_LC, trans = "norm", smin = 0, smax = 1)
Niko_11_LC<-raster(LC_Monthly[["Niko_11"]]  , DF="PDF")
Niko_11_LC<-raster.transformation(Niko_11_LC, trans = "norm", smin = 0, smax = 1)
Preman_11_LC<-raster(LC_Monthly[["Preman_11"]] , DF="PDF")
Preman_11_LC<-raster.transformation(Preman_11_LC, trans = "norm", smin = 0, smax = 1)
Teju_11_LC<-raster(LC_Monthly[["Teju_11"]], DF="PDF")
Teju_11_LC<-raster.transformation(Teju_11_LC, trans = "norm", smin = 0, smax = 1)
Tomi_11_LC<-raster(LC_Monthly[["Tomi_11"]], DF="PDF")
Tomi_11_LC<-raster.transformation(Tomi_11_LC, trans = "norm", smin = 0, smax = 1)
Wodan_11_LC<-raster(LC_Monthly[["Wodan_11"]], DF="PDF")
Wodan_11_LC<-raster.transformation(Wodan_11_LC, trans = "norm", smin = 0, smax = 1)

#December
Dayak_12_LC<-raster(LC_Monthly[["Dayak_12"]]  , DF="PDF")
Dayak_12_LC<-raster.transformation(Dayak_12_LC, trans = "norm", smin = 0, smax = 1)
Henk_12_LC<-raster(LC_Monthly[["Henk_12"]], DF="PDF")
Henk_12_LC<-raster.transformation(Henk_12_LC, trans = "norm", smin = 0, smax = 1)
Niko_12_LC<-raster(LC_Monthly[["Niko_12"]]  , DF="PDF")
Niko_12_LC<-raster.transformation(Niko_12_LC, trans = "norm", smin = 0, smax = 1)
Preman_12_LC<-raster(LC_Monthly[["Preman_12"]] , DF="PDF")
Preman_12_LC<-raster.transformation(Preman_12_LC, trans = "norm", smin = 0, smax = 1)
Teju_12_LC<-raster(LC_Monthly[["Teju_12"]], DF="PDF")
Teju_12_LC<-raster.transformation(Teju_12_LC, trans = "norm", smin = 0, smax = 1)
Tomi_12_LC<-raster(LC_Monthly[["Tomi_12"]], DF="PDF")
Tomi_12_LC<-raster.transformation(Tomi_12_LC, trans = "norm", smin = 0, smax = 1)
Wodan_12_LC<-raster(LC_Monthly[["Wodan_12"]], DF="PDF")
Wodan_12_LC<-raster.transformation(Wodan_12_LC, trans = "norm", smin = 0, smax = 1)


# create spatial points data frames to extract from specific rasters


CRS<-CRS("+init=epsg:32750 +proj=utm +zone=50 +units=m +south")
library(sp)
library(raster)
#subset for Chili in August 2012
xy <- Chili_subset8[,c("x","y")]
Chili_subset8_sp <- SpatialPointsDataFrame(coords = xy, data = Chili_subset8, proj4string = CRS)
Chili_Own <- raster::extract(Chili_Own, Chili_subset8_sp[,c("x", "y")])
Chili_Dayak<-raster::extract(ZeroRast, Chili_subset8_sp[,c("x", "y")])
Chili_Henk<-raster::extract(ZeroRast, Chili_subset8_sp[,c("x", "y")])
Chili_Niko<-raster::extract(ZeroRast, Chili_subset8_sp[,c("x", "y")])
Chili_Otto8<-raster::extract(Otto_rast8, Chili_subset8_sp[,c("x", "y")])
Chili_Katmandun<-raster::extract(ZeroRast, Chili_subset8_sp[,c("x", "y")])
Chili_Teju<-raster::extract(ZeroRast, Chili_subset8_sp[,c("x", "y")])
Chili_Tomi<-raster::extract(Tomi_rast8, Chili_subset8_sp[,c("x", "y")])
Chili_Wodan<-raster::extract(Wodan_rast8, Chili_subset8_sp[,c("x", "y")])
Chili_Juni<-raster::extract(Juni_rast, Chili_subset8_sp[,c("x", "y")])
Chili_Inul<-raster::extract(Inul_rast, Chili_subset8_sp[,c("x", "y")])
Chili_Kerry <-raster::extract(Kerry_rast, Chili_subset8_sp[,c("x", "y")])
Chili_Mindy <-raster::extract(Mindy_rast, Chili_subset8_sp[,c("x", "y")])
Niko_08_LC<-raster::extract(Niko_08_LC, Chili_subset8_sp[,c("x", "y")])
Teju_08_LC<-raster::extract(Teju_08_LC, Chili_subset8_sp[,c("x", "y")])
Chili_8<-cbind(Chili_subset8_sp, Chili_Own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto8, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy, Niko_08_LC, Teju_08_LC)
as.data.frame(Chili_8)
names(Chili_8)[7] <- "Own"
names(Chili_8)[8] <- "Dayak"
names(Chili_8)[9] <- "Henk"
names(Chili_8)[10] <- "Niko"
names(Chili_8)[11] <- "Otto"
names(Chili_8)[12] <- "Katmandun"
names(Chili_8)[13] <- "Teju"
names(Chili_8)[14] <- "Tomi"
names(Chili_8)[15] <- "Wodan"
names(Chili_8)[16] <- "Juni"
names(Chili_8)[17] <- "Inul"
names(Chili_8)[18] <- "Kerry"
names(Chili_8)[19] <- "Mindy"
names(Chili_8)[20]<-"Niko_LC"
names(Chili_8)[21]<-"Teju_LC"
Chili_8<-as.data.frame(Chili_8@data)

xy <- Chili_subset9[,c("x","y")]
Chili_subset9_sp <- SpatialPointsDataFrame(coords = xy, data = Chili_subset9, proj4string = CRS)
Chili_Own <- raster::extract(Chili_Own, Chili_subset9_sp[,c("x", "y")])
Chili_Dayak<-raster::extract(ZeroRast, Chili_subset9_sp[,c("x", "y")])
Chili_Henk<-raster::extract(ZeroRast, Chili_subset9_sp[,c("x", "y")])
Chili_Niko<-raster::extract(Niko_rast9, Chili_subset9_sp[,c("x", "y")])
Chili_Otto9<-raster::extract(ZeroRast, Chili_subset9_sp[,c("x", "y")])
Chili_Katmandun<-raster::extract(ZeroRast, Chili_subset9_sp[,c("x", "y")])
Chili_Teju<-raster::extract(ZeroRast, Chili_subset9_sp[,c("x", "y")])
Chili_Tomi<-raster::extract(Tomi_rast9, Chili_subset9_sp[,c("x", "y")])
Chili_Wodan<-raster::extract(Wodan_rast9, Chili_subset9_sp[,c("x", "y")])
Chili_Juni<-raster::extract(Juni_rast, Chili_subset9_sp[,c("x", "y")])
Chili_Inul<-raster::extract(Inul_rast, Chili_subset9_sp[,c("x", "y")])
Chili_Kerry <-raster::extract(Kerry_rast, Chili_subset9_sp[,c("x", "y")])
Chili_Mindy <-raster::extract(Mindy_rast, Chili_subset9_sp[,c("x", "y")])


Chili_9<-cbind(Chili_subset9_sp, Chili_Own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto9, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
as.data.frame(Chili_9)
names(Chili_9)[7] <- "Own"
names(Chili_9)[8] <- "Dayak"
names(Chili_9)[9] <- "Henk"
names(Chili_9)[10] <- "Niko"
names(Chili_9)[11] <- "Otto"
names(Chili_9)[12] <- "Katmandun"
names(Chili_9)[13] <- "Teju"
names(Chili_9)[14] <- "Tomi"
names(Chili_9)[15] <- "Wodan"
names(Chili_9)[16] <- "Juni"
names(Chili_9)[17] <- "Inul"
names(Chili_9)[18] <- "Kerry"
names(Chili_9)[19] <- "Mindy"
Chili_9<-as.data.frame(Chili_9@data)



xy <- Chili_subset10[,c("x","y")]
Chili_subset10_sp <- SpatialPointsDataFrame(coords = xy, data = Chili_subset10, proj4string = CRS)
Chili_Own <- raster::extract(Chili_Own, Chili_subset10_sp[,c("x", "y")])
Chili_Dayak<-raster::extract(ZeroRast, Chili_subset10_sp[,c("x", "y")])
Chili_Henk<-raster::extract(Henk_rast10, Chili_subset10_sp[,c("x", "y")])
Chili_Niko<-raster::extract(ZeroRast, Chili_subset10_sp[,c("x", "y")])
Chili_Otto<-raster::extract(ZeroRast, Chili_subset10_sp[,c("x", "y")])
Chili_Katmandun<-raster::extract(ZeroRast, Chili_subset10_sp[,c("x", "y")])
Chili_Teju<-raster::extract(ZeroRast, Chili_subset10_sp[,c("x", "y")])
Chili_Tomi<-raster::extract(Tomi_rast10, Chili_subset10_sp[,c("x", "y")])
Chili_Wodan<-raster::extract(Wodan_rast9, Chili_subset10_sp[,c("x", "y")])
Chili_Juni<-raster::extract(Juni_rast, Chili_subset10_sp[,c("x", "y")])
Chili_Inul<-raster::extract(Inul_rast, Chili_subset10_sp[,c("x", "y")])
Chili_Kerry <-raster::extract(Kerry_rast, Chili_subset10_sp[,c("x", "y")])
Chili_Mindy <-raster::extract(Mindy_rast, Chili_subset10_sp[,c("x", "y")])
Chili_10<-cbind(Chili_subset10_sp, Chili_Own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
as.data.frame(Chili_10)
names(Chili_10)[7] <- "Own"
names(Chili_10)[8] <- "Dayak"
names(Chili_10)[9] <- "Henk"
names(Chili_10)[10] <- "Niko"
names(Chili_10)[11] <- "Otto"
names(Chili_10)[12] <- "Katmandun"
names(Chili_10)[13] <- "Teju"
names(Chili_10)[14] <- "Tomi"
names(Chili_10)[15] <- "Wodan"
names(Chili_10)[16] <- "Juni"
names(Chili_10)[17] <- "Inul"
names(Chili_10)[18] <- "Kerry"
names(Chili_10)[19] <- "Mindy"
Chili_10<-as.data.frame(Chili_10@data)

xy <- Chili_subset11[,c("x","y")]
Chili_subset11_sp <- SpatialPointsDataFrame(coords = xy, data = Chili_subset11, proj4string = CRS)
Chili_Own <- raster::extract(Chili_Own, Chili_subset11_sp[,c("x", "y")])
Chili_Dayak<-raster::extract(ZeroRast, Chili_subset11_sp[,c("x", "y")])
Chili_Henk<-raster::extract(ZeroRast, Chili_subset11_sp[,c("x", "y")])
Chili_Niko<-raster::extract(ZeroRast, Chili_subset11_sp[,c("x", "y")])
Chili_Otto<-raster::extract(ZeroRast, Chili_subset11_sp[,c("x", "y")])
Chili_Katmandun<-raster::extract(ZeroRast, Chili_subset11_sp[,c("x", "y")])
Chili_Teju<-raster::extract(Teju_rast11, Chili_subset11_sp[,c("x", "y")])
Chili_Tomi<-raster::extract(Tomi_rast11, Chili_subset11_sp[,c("x", "y")])
Chili_Wodan<-raster::extract(Wodan_rast11, Chili_subset11_sp[,c("x", "y")])
Chili_Juni<-raster::extract(Juni_rast, Chili_subset11_sp[,c("x", "y")])
Chili_Inul<-raster::extract(Inul_rast, Chili_subset11_sp[,c("x", "y")])
Chili_Kerry <-raster::extract(Kerry_rast, Chili_subset11_sp[,c("x", "y")])
Chili_Mindy <-raster::extract(Mindy_rast, Chili_subset11_sp[,c("x", "y")])
Chili_11<-cbind(Chili_subset11_sp, Chili_Own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
as.data.frame(Chili_11)
names(Chili_11)[7] <- "Own"
names(Chili_11)[8] <- "Dayak"
names(Chili_11)[9] <- "Henk"
names(Chili_11)[10] <- "Niko"
names(Chili_11)[11] <- "Otto"
names(Chili_11)[12] <- "Katmandun"
names(Chili_11)[13] <- "Teju"
names(Chili_11)[14] <- "Tomi"
names(Chili_11)[15] <- "Wodan"
names(Chili_11)[16] <- "Juni"
names(Chili_11)[17] <- "Inul"
names(Chili_11)[18] <- "Kerry"
names(Chili_11)[19] <- "Mindy"
Chili_11<-as.data.frame(Chili_11@data)


Chili_df<- rbind(Chili_8, Chili_9, Chili_10, Chili_11)
#any NA points within UDs would be 0
Chili_df[is.na(Chili_df)] <- 0
write.csv(Chili_df, file="Chili_LC_DF.csv")

library(brms)
library(rstan)
library(bayesplot)
library(DHARMa)
library(tidybayes)
library(performance)

Chili_df$PT<-as.numeric(Chili_df$PT)
table(is.na(Chili_df))

# test spatial autocorrelation
library(DHARMa)
library(stringr)

Chili_df$coords <- paste(Chili_df$x,", ",Chili_df$y)
coords <- c(unique(Chili_df$coords))
x_unique <- c(str_extract(coords, "^.+(?=,)"))
y_unique <- c(str_extract(coords, "(?<=, ).+$"))
fittedModel <- glm(PT ~ 1, data = Chili_df, family="binomial")
sims<-simulateResiduals(fittedModel)
simsrecalc<-recalculateResiduals(sims,group = Chili_df$coords)
testSpatialAutocorrelation(simsrecalc, x = x_unique, y = y_unique)

# no evidence of spatial autocorrelation
# test temporal autocorrelation
res = recalculateResiduals(sims, group = Chili_df$timestamp)
testTemporalAutocorrelation(res, time = unique(Chili_df$timestamp))

# evidence that there is temporal autocorrelation
# need to include within models

library(brms)
Chili_mod1 <- brms::brm(PT ~ Own+Dayak+ Henk +Niko +Katmandun +Otto +Teju+Tomi+Wodan + Juni +Inul+Kerry+Mindy, data = Chili_df, family = bernoulli(link="logit"),iter = 1000, warmup = 500, chains = 1, cores = 4, init=0, control = list(adapt_delta = 0.99, max_treedepth=15), save_pars = save_pars(all = TRUE))









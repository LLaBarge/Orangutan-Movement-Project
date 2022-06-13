
# Lc project script for 2012 long calls
# script for creating AKDEs and then polygons for Long call project
library(move)
#devtools::install_github("ctmm-initiative/ctmm")
library(ctmm)
library(readr)
setwd("C:/Users/lrlab/OneDrive/Desktop/LC_Project")
#edit file that will go onto movebank to only have data from months
# with enough location pts/individual to include in 
# occurrence distributions >50 locations/month
Tuanan_GPS_2012 <- read_csv("Tuanan GPS 2012.csv")
as.data.frame(Tuanan_GPS_2012)
Tuanan_GPS_2012$ID_Month<-paste(Tuanan_GPS_2012$id, Tuanan_GPS_2012$Month)
Tuanan_filtered<-subset(Tuanan_GPS_2012, with(Tuanan_GPS_2012, ID_Month %in% names(which(table(ID_Month)>=50))))
nrow(Tuanan_filtered)
nrow(Tuanan_GPS_2012)

#write csv and add to movebank manually for ctmm
#write.csv(Tuanan_filtered, file="Tuanan GPS 2012_Corrected.csv")


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

# Don't run models - load saved file (will take 1+hrs)
#save(AKDE, file="akdes.Rdata")
#save(FIT, file="modFIT.Rdata")
load(file="modFIT.Rdata")
load(file="akdes.Rdata")

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
  AKDE[[i]] <- akde(orangs[[i]],FIT[[i]][[1]], grid=list(extent=extent, dr=25))
}
names(AKDE) <- names(orangs)



# some individuals appear so infrequently that ODs are likely not representative
# of competition risk, these are:
#Cinta, Tina, Kentung, Fatih, Tuko, Manggo, Jordan, Pinky Luca, Dolay, Preman
# others were not present in the study area for more than 2 months
#in 2012: Bobo, Fayesh, Flumnmu, Fugit, Talia, Tuko

OD_2012 <- getMovebankData(study="Bornean Orangutan, Monthly data for occurrence distributions", login=login)

# create as.telemetry object for ctmm models
OD_tel<-as.telemetry(OD_2012, projection = "+init=epsg:32750 +proj=utm +zone=50 +units=m +south")

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



save(OD.FIT, file="OD.FIT.Rdata")
save(OD, file="OD.Rdata")
load(file="OD.Rdata")

# create a dataframe of GPS points for "available" locations - to use later

Avail_pts<- do.call(rbind, lapply(OD_tel, function(x) as.data.frame(x)))
head(Avail_pts)
Avail_pts$ID <- rownames(Avail_pts)
library(tm)
removeNumbers(Avail_pts$ID)
Avail_pts$ID<- gsub("\\..*","",Avail_pts$ID)
unique(Avail_pts$ID)
# select males of interest (those that were range resident in 2012)

Avail_pts<-Avail_pts[Avail_pts$ID %in% c("Chili", "Niko",  "Tomi", "Teju" , "Wodan"), ]
str(Avail_pts)
PT<-rep("Available",times=2510)
Avail_pts$PT<-PT
head(Avail_pts)

# export rasters with probability distribution function as a proxy for competition risk
# create separate per individual/year to keep track
library(raster)

#create individual raster files for competition risk / prob of 
# encounter with potential receptive females
# use nearest OD distribution to the date of data collection
# individuals chosen by looking at Tuanan presence records


library(climateStability)
View(OD)
ZeroRast<-Tuanan_r
vals <- NA
ZeroRast <- setValues(ZeroRast, vals)
View(AKDE)
Chili_<-raster(AKDE[["Chili"]] , DF="PDF")
Chili_Own<-rescale0to1(Chili_)
sp::plot(Chili_Own)
Chili_rast7<-raster(OD[["Chili.8"]] , DF="PDF")
Chili_rast7<-rescale0to1(Chili_rast7)
Chili_rast7<-raster(OD[["Chili.8"]] , DF="PDF")
Chili_rast7<-rescale0to1(Chili_rast7)
Chili_rast8<-raster(OD[["Chili.8"]] , DF="PDF")
Chili_rast8<-rescale0to1(Chili_rast8)
Chili_rast9<-raster(OD[["Chili.9"]] , DF="PDF")
Chili_rast9<-rescale0to1(Chili_rast9)
Chili_rast10<-raster(OD[["Chili.10"]] , DF="PDF")
Chili_rast10<-rescale0to1(Chili_rast10)
Chili_rast11<-raster(OD[["Chili.11"]] , DF="PDF")
Chili_rast11<-rescale0to1(Chili_rast11)
Chili_rast12<-raster(OD[["Chili.11"]] , DF="PDF")
Chili_rast12<-rescale0to1(Chili_rast12)


Dayak_rast5<-raster(OD[["Dayak.5"]], DF="PDF")
Dayak_rast5<-rescale0to1(Dayak_rast5)
Dayak_rast9<-raster(OD[["Dayak.5"]], DF="PDF")
Dayak_rast9<-rescale0to1(Dayak_rast9)
Dayak_rast10<-raster(OD[["Dayak.5"]], DF="PDF")
Dayak_rast10<-rescale0to1(Dayak_rast10)
Dayak_rast11<-raster(OD[["Dayak.5"]], DF="PDF")
Dayak_rast11<-rescale0to1(Dayak_rast11)

Desy_rast<-raster(AKDE[["Desy"]] , DF="PDF")
Desy_rast<-rescale0to1(Desy_rast)

Henk_rast3<-raster(OD[["Henk.3"]], DF="PDF")
Henk_rast3<-rescale0to1(Henk_rast3)
Henk_rast10<-raster(OD[["Henk.10"]], DF="PDF")
Henk_rast10<-rescale0to1(Henk_rast10)
Henk_rast12<-raster(OD[["Henk.12"]], DF="PDF")
Henk_rast12<-rescale0to1(Henk_rast12)

Jinak_rast<-raster(AKDE[["Jinak"]] , DF="PDF")
Jinak_rast<-rescale0to1(Jinak_rast)

Juni_rast<-raster(AKDE[["Juni"]] , DF="PDF")
Juni_rast<-rescale0to1(Juni_rast)

Inul_rast<-raster(AKDE[["Inul"]] , DF="PDF")
Inul_rast<-rescale0to1(Inul_rast)

Katmandun_rast4<-raster(OD[["Katmandun.4"]] , DF="PDF")
Katmandun_rast4<-rescale0to1(Katmandun_rast4)
Katmandun_rast5<-raster(OD[["Katmandun.4"]] , DF="PDF")
Katmandun_rast5<-rescale0to1(Katmandun_rast5)
Katmandun_rast6<-raster(OD[["Katmandun.6"]] , DF="PDF")
Katmandun_rast6<-rescale0to1(Katmandun_rast6)

Kerry_rast<-raster(AKDE[["Kerry"]] , DF="PDF")
Kerry_rast<-rescale0to1(Kerry_rast)

Kondor_rast<-raster(AKDE[["Kondor"]] , DF="PDF")
Kondor_rast<-rescale0to1(Kondor_rast)

Mindy_rast<-raster(AKDE[["Mindy"]] , DF="PDF")
Mindy_rast<-rescale0to1(Mindy_rast)

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
library(move)
login <- movebankLogin(username="llabarge", password="Samango1992!")
getMovebankStudies(login)

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
head(LC_pts)

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
# need to filter dataframe to remove duplicate pts or those within 20 minutes of one another
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
LC_points$Month <- str_sub(LC_points$timestamp, start = 6L, end = 7L)

Avail_pts$Month <- str_sub(Avail_pts$timestamp, start = 6L, end = 7L)
LC_com<-LC_points%>% dplyr::select('ID', 'PT', 'Month', 'x', 'y', 'timestamp')
Av_com<-Avail_pts%>% dplyr::select('ID', 'PT', 'Month','x', 'y', 'timestamp')
LC_Control<-rbind(LC_com, Av_com)
as.data.frame(LC_Control)
head(LC_Control)


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

# create spatial points data frames to extract from specific rasters

CRS<-CRS("+init=epsg:32750 +proj=utm +zone=50 +units=m +south")
library(sp)
library(raster)
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
Chili_8<-cbind(Chili_subset8_sp, Chili_Own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto8, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
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
Chili_9<-cbind(Chili_subset9_sp, Chili_own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto9, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
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
Chili_own <- raster::extract(Chili_own, Chili_subset10_sp[,c("x", "y")])
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
Chili_10<-cbind(Chili_subset10_sp, Chili_own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
as.data.frame(Chili_10)
names(Chili_10)[7] <- "own"
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
Chili_own <- raster::extract(Chili_own, Chili_subset11_sp[,c("x", "y")])
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
Chili_11<-cbind(Chili_subset11_sp, Chili_own, Chili_Dayak, Chili_Henk, Chili_Niko, Chili_Otto, Chili_Katmandun, Chili_Teju, Chili_Tomi, Chili_Wodan, Chili_Juni, Chili_Inul, Chili_Kerry, Chili_Mindy)
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







# take the raw Passive acoustic monitoring data (provided by B. Spillmann)
# and subset for the known individuals
setwd("C:/Users/lrlab/OneDrive/Desktop/LC_Project/GPS_and_LongCallsGiven")
library(readr)

LC_PAM_data <- read_csv("LongCalls_AAM_Spillman.csv")
head(LC_PAM)
unique(LC_PAM_data$caller_ID.def)
LC_PAM<-LC_PAM_data[LC_PAM_data$caller_ID.def %in% c("Helium", "Henk" ,  "Preman", "Katmandun" , "MAX" ,"Max" ,"Chili" ,"Luca","Dayak","Teju","Flunmu", "Wodan", "Niko", "Otto", "Tomi"), ]
unique(LC_PAM$caller_ID.def)
LC_PAM$caller_ID.def <- gsub("MAX", "Max", LC_PAM$caller_ID.def)
LC_PAM$ID<-LC_PAM$caller_ID.def
LC_PAM$x<-LC_PAM$caller_Pos.x
LC_PAM$y<-LC_PAM$caller_Pos.y
LC_PAM$Date <- as.Date(as.character(LC_PAM$Date), format='%Y%m%d')
LC_PAM$timestamp<- as.POSIXct(paste(LC_PAM$Date, LC_PAM$caller_Time), format="%Y-%m-%d %H:%M:%S")



#select only variables needed for Movebank
library(dplyr)

LC_PAM_to_Movebank<-LC_PAM%>% dplyr::select('ID', 'x', 'y', 'timestamp')


# extract month to a new column and
# create an ID_date column to create monthly occurrence distributions for long calls


LC_PAM_to_Movebank$month <- format(LC_PAM_to_Movebank$timestamp, "%m")  
head(LC_PAM_to_Movebank)
LC_PAM_to_Movebank$ID_date <- paste(LC_PAM_to_Movebank$ID,"_", LC_PAM_to_Movebank$month)
str(LC_PAM_to_Movebank)
# remove any 0s from locations

LC_PAM_to_Movebank <- LC_PAM_to_Movebank[!LC_PAM_to_Movebank$x == "0" | !LC_PAM_to_Movebank$y == "0" , ]


# remove any NA values
LC_PAM_to_Movebank<-na.omit(LC_PAM_to_Movebank)


# transform to geographic coords to fit ctmm/movebank format
library(sp)
coordinates(LC_PAM_to_Movebank) <- c("x", "y")
proj4string(LC_PAM_to_Movebank) <- CRS("+init=epsg:32750 +proj=utm +zone=50 +units=m +south")
LC_PAM_to_Movebank <- spTransform(LC_PAM_to_Movebank, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
LC_PAM_to_Movebank<-as.data.frame(LC_PAM_to_Movebank)
# join with observed long calls

PTS_2012_thru_2015 <- read_csv("C:/Users/lrlab/OneDrive/Desktop/GPS201718/CSVs/PTS_2012_thru_2015.csv")
head(LC_PAM_to_Movebank)
head(PTS_2012_thru_2015)

PTS_2012_thru_2015$x<-PTS_2012_thru_2015$long
PTS_2012_thru_2015$y<-PTS_2012_thru_2015$lat
PTS_2012_thru_2015$ID<-PTS_2012_thru_2015$id
PTS_2012_thru_2015$timestamp<-PTS_2012_thru_2015$date

# add a month column to create another ID_date column

PTS_2012_thru_2015$month <- format(PTS_2012_thru_2015$timestamp, "%m")  
PTS_2012_thru_2015$ID_date <- paste(PTS_2012_thru_2015$ID,"_", PTS_2012_thru_2015$month)

# filter to date range of when PAM system was operating
PTS_2012<-PTS_2012_thru_2015[PTS_2012_thru_2015$timestamp >= "2012-03-01" & PTS_2012_thru_2015$timestamp <= "2012-12-31",]

# select same variables as data from PAM system
OBS_LC<-PTS_2012%>% dplyr::select('ID', 'x', 'y', 'timestamp', 'ID_date')


PAM_LC<-LC_PAM_to_Movebank%>% dplyr::select('ID', 'x', 'y', 'timestamp', 'ID_date')
head(OBS_LC)
head(PAM_LC)
LC_MOVE<- rbind(OBS_LC, PAM_LC)


# check for 0s or NAs in geometry before manual upload to movebank
table(is.na(LC_MOVE))
View(LC_MOVE)

# did not filter for points near in time because this is done later
# save to csv to upload manually to movebank
write.csv(LC_MOVE, file="LC_MOVE.csv")


# goes under the name 	"	Monthly Long Call Data, 2012, Bornean Orangutans"


# take the raw Passive acoustic monitoring data (provided by B. Spillmann)
# and subset for the known individuals
library(readr)
LC_PAM_data <- read_csv("C:/Users/lrlab/OneDrive/Desktop/LC_Project/LC_PAMsystem_forediting_R.csv")
head(LC_PAM)
unique(LC_PAM_data$caller_ID.def)
LC_PAM<-LC_PAM_data[LC_PAM_data$caller_ID.def %in% c("Helium", "Henk" ,  "Preman", "Katmandun" , "MAX" ,"Max" ,"Chili" ,"Luca","Dayak","Teju","Flunmu", "Wodan", "Niko", "Otto", "Tomi"), ]
unique(LC_PAM$caller_ID.def)
LC_PAM$caller_ID.def <- gsub("MAX", "Max", LC_PAM$caller_ID.def)
LC_PAM$ID<-LC_PAM$caller_ID.def
LC_PAM$x<-LC_PAM$caller_Pos.x
LC_PAM$y<-LC_PAM$caller_Pos.y
LC_PAM$timestamp<-LC_PAM$Timestamp.1

#select only variables needed for Movebank
library(dplyr)

LC_PAM_to_Movebank<-LC_PAM%>% dplyr::select('ID', 'Month...5', 'x', 'y', 'timestamp')

# create an ID_date column to create monthly occurrence distributions for long calls
LC_PAM_to_Movebank$ID_date <- paste(LC_PAM_to_Movebank$ID,LC_PAM$Month...5)
str(LC_PAM_to_Movebank)
# save to csv to upload manually to movebank
write.csv(LC_PAM_to_Movebank, file="LC_PAM_to_Movebank.csv")

# goes under the name 	"Monthly LC, Bornean Orangutan, Flanged Males, Tuanan, Central Kalimatan, Indonesia"

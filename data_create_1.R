
require(plyr)
require(dplyr)
require(rvest)
require(stringr)
require(haven)
require(survey)

set.seed(2)
#builds the data

#STATES TO INCLUDE
STATES_INCLUDE=c("Maryland","Virginia","Oregon","Washington")
#if you want all states, should be c(".")

#How Many Cases?
HOW_MANY_CASES<-20000

#if you want to keep all cases, set this to false. How many cases wont be used. However, make sure you remember to weight cases still. 
SAMPLECASES<-TRUE


fips<-read_html("http://www.epa.gov/enviro/html/codes/state.html")%>% html_table()
STATES_INCLUDE_FIPS=fips[[1]]$"FIPS Code"[which(str_detect(fips[[1]]$"State Name", ignore.case(paste(paste("^",STATES_INCLUDE,sep=""), collapse="|"))))]



#reads in population sas files
#using sas files helps to preserve variable types
sasa_p<-read_sas("psam_pusa.sas7bdat")
sasb_p<-read_sas("psam_pusb.sas7bdat")
sasab_p<-rbind(sasa_p,sasb_p)
rm(sasa_p)
rm(sasb_p)
sasa_h<-read_sas("psam_husa.sas7bdat")
sasb_h<-read_sas("psam_husb.sas7bdat")
sasab_h<-rbind(sasa_h,sasb_h)


sasfull<-join(sasab_p,sasab_h, by=c("SERIALNO"))
nrow(sasfull)
rm(sasab_h)
rm(sasab_p)
library(RCurl)
match<-read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/11MUycR2TYT9r8P9bfTEYJgB6t-f_5UFFlzr4RDFCVYw/pub?gid=357603510&single=true&output=csv")))

rm(sasa_h)
rm(sasb_h)

colnames(sasfull)<-make.unique(colnames(sasfull))

#Choose states to include. This includes Washington, Idaho, Virginia, and Maryland
require(rgdal)
statelist<-as.character(STATES_INCLUDE_FIPS)
acsselect<-filter(sasfull, ST%in%statelist)

geos<-lapply(statelist, function(X){readOGR("pumamap",paste("tl_2013_",X,"_puma10",sep=""))})


sprbind.unique<-function(..., makeUniqueIDs=TRUE){
  dots = list(...)
  names(dots) <- NULL
  lst = lapply(dots, function(x) as(x, "SpatialPolygons"))
  lst$makeUniqueIDs = makeUniqueIDs
  pl = do.call(rbind.SpatialPolygons, lst)
  df = do.call(rbind, lapply(dots, function(x) x@data))
  SpatialPolygonsDataFrame(pl, df)
}

geoshape<-do.call(sprbind.unique,geos)
rm(geos)

#plots PUMAS of included states as a first test

plot(geoshape)


#unload full dataset, keeping only desired states
rm(sasfull)

 
#Here, we begin to recode variables. 
library(RCurl)
cbook<-read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/1mfZy2fap2eqLIH9_vrbFs35hQMfGnQdKp39XEjjPpVw/pub?gid=114600137&single=true&output=csv")))

                
acsrecode<-function(X){
  filton<-filter(cbook, Variable==X)
  data.frame("level"=str_trim(as.character(unlist(unlist(strsplit(as.character(filton$Code), ".",fixed=T))))) %>% .[seq(1,length(.),2)],"label"=str_trim(unlist(unlist(strsplit(as.character(filton$Code), ".",fixed=T)))) %>% .[seq(2,length(.),2)],'variable'=X)
}



acsrecode.df<-function(X,INframe){
  INframe[X]<-join(data.frame("level"=as.factor(INframe[,X])),acsrecode(X), by='level',type="left",match="first")$label
  print(summary(INframe[X]))
  INframe
}



options( "survey.replicates.mse" = TRUE ) 

acsrepw<-svrepdesign(repweights = acsselect[204:283], weights =acsselect$PWGTP, combined.weights = TRUE, type="JK1",scale=4/80,rscales=rep(1,80),data=acsselect)

#test_res<-sample_n(acsselect,70000,weight=acsselect$PWGTP)

#mean(as.numeric(test_res$SEX))-1
if(SAMPLECASES==TRUE){
wmem<-weights(acsrepw)
acsrepsamp=list()
wmem[which(wmem<0)]<-0
validweights<-do.call("rbind",lapply(1:80, function(i){if(TRUE%in%c(wmem[,i]<0)==FALSE){i}}))
for(i in validweights){
  sampset<-which(wmem[,i]>0)
acsrepsamp[[i]]<-sample_n(acsselect[sampset,],round(HOW_MANY_CASES/length(validweights)),weight=wmem[sampset,i],replace=FALSE)
}
acsrepsamp<-do.call("rbind",acsrepsamp)
} else {

  acsrepsamp<-acsselect

cat("WARNING: IF YOU SET SAMPLING TO FALSE, YOU NEED TO THINK ABOUT SURVEY WEIGHTS. \n The object acsrepw is a svy object that uses replicate weights.\n I recommend using it.\n However, the recoding of variables operates on a normal dataset.\n Once the recoding sequence is done, simply remake the object using the code in line 103,\n but replacing data with acsrepsamp.")
readline("Press <return to continue.")

}
head(acsrepsamp)

#now, we select only those variables that will be useful to students. all of this speeds up recoding
ncol(acsrepsamp)
acsrepsamp2<-acsrepsamp[,-c(6,128:202,204:283,284,286,288,290,289,290,286,288,434:513,382:433)]
rm(acsrepsamp)
colnames(acsrepsamp)
colnames(acsrepsamp2)
colnames(acsrepsamp2)
summary(acsrepsamp2)
acsrepsamp2.strip<-acsrepsamp2
levels(as.factor(acsrepsamp2$MIGPUMA))
for(VARNAME in colnames(acsrepsamp2)){
  print(VARNAME)
  if(VARNAME%in%c("MIGPUMA","POWPUMA","PUMA","SERIALNO")==FALSE){
    if(class(acsrepsamp2[,VARNAME])=="character"){
      acsrepsamp2=acsrecode.df(VARNAME,acsrepsamp2)
    }
  }
}

head(acsrepsamp2)
library(haven)

fm<-read.csv("farmersmarket.csv")
fm<-filter(fm, State%in%STATES_INCLUDE)
fm<-fm[which(is.na(fm$x)==FALSE),]
coordinates(fm)<-c("x","y")
geoshape<-spTransform(geoshape, CRS("+proj=longlat +datum=WGS84"))
proj4string(fm)<-proj4string(geoshape)
require(rgeos)
acsrepsamp2$FARMMARK<-NA
proj4string(fm)
#distancetoafarmersmmarket


for(i in 1:nrow(acsrepsamp2)){
  try(acsrepsamp2$FARMMARK[i]<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==acsrepsamp2$PUMA[i]),n=1,type="random"),fm, byid=TRUE)))
}
for(i in which(is.na(acsrepsamp2$FARMMARK))
){
  try(acsrepsamp2$FARMMARK[i]<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==acsrepsamp2$PUMA[i]),n=1,type="random"),fm, byid=TRUE)))
}
colnames(geoshape@data)
test<-ddply(acsrepsamp2, .(PUMA), summarize,avgdist=mean(FARMMARK))
head(test)
colnames(test)[1]<-"PUMACE10"
require(ggplot2)
geoshape$id<-row.names(geoshape)
gfort<-fortify(geoshape)
gfort<-join(gfort,geoshape@data, by="id")

gfort<-join(gfort, test, by="PUMACE10")

ggplot(gfort)+geom_polygon(aes(x=long,y=lat, group=group, fill=avgdist))+theme_bw()

anair<-read.csv("annual_all_2014.csv")

anair<-filter(anair, State.Name%in%c("Washington","Oregon","Virginia","Maryland"), Parameter.Name=="PM2.5 - Local Conditions")
head(anair)

coordinates(anair)<-c("Longitude","Latitude")
proj4string(anair)<-CRS("+proj=longlat +datum=NAD27")
anair<-spTransform(anair,proj4string(geoshape))
head(anair@data)
acsrepsamp2$PM25ARM<-NA
for(i in 1:nrow(acsrepsamp2)){
  try(temp<-gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==acsrepsamp2$PUMA[i]),n=1,type="random"),anair, byid=TRUE))
  
  try(acsrepsamp2$PM25ARM[i]<-anair@data$Arithmetic.Mean[which(temp==min(temp)) %>% sample(.,size=1)])
}

for(i in which(is.na(acsrepsamp2$PM25ARM))){
  try(temp<-gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==acsrepsamp2$PUMA[i]),n=1,type="random"),anair, byid=TRUE))
  
  try(acsrepsamp2$PM25ARM[i]<-anair@data$Arithmetic.Mean[which(temp==min(temp)) %>% sample(.,size=1)])
}

test2<-ddply(acsrepsamp2, .(PUMA), summarize,PM25=mean(PM25ARM))
head(test2)
colnames(test2)[1]<-"PUMACE10"
gfort<-join(gfort, test2, by="PUMACE10")
head(gfort)
ggplot(gfort)+geom_polygon(aes(x=long,y=lat,group=group, fill=PM25))+theme_bw()
anair@data$Arithmetic.Mean

summary(acsrepsamp2)


#TRI DATA currently only for WA/OR/ID
tri<-read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/1yJpc6qPw4JmloJQEqcIQwoC7ZO5KeWbBJRSluH8eho0/pub?output=csv")))
head(tri)
tri$Lat<-as.numeric(as.character(tri$Lat))
tri$Long<-as.numeric(as.character(tri$Long))
tri<-tri[which(is.na(tri$Long)==FALSE),]
tri<-filter(tri, Lat<=50)
point 
plot(tri)
library(kriging)
head(tri@data)
tri$Volume<-as.numeric(gsub(",","",as.character(tri$Volume)))
tri<-tri[which(is.na(tri$Volume)==FALSE),]
tri<-tri[which(duplicated(data.frame(tri$Lat,tri$Long))==FALSE),]

tri.sp<-tri
coordinates(tri.sp)<-c("Long","Lat")
proj4string(tri.sp)<-proj4string(geoshape)
temp<-over(tri.sp,geoshape)$PUMACE10

require(reshape2)
temp<-melt(table(temp))
colnames(temp)<-c("PUMA","TRIFACILITES")
acsrepsamp2<-join(acsrepsamp2, temp, by="PUMA")
#Raw count of the number of Toxic Release Inventory facilities within the PUMA

for(i in 1:nrow(acsrepsamp2)){
  try(temp<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==acsrepsamp2$PUMA[i]),n=1,type="random"),tri.sp, byid=TRUE)))
  try(acsrepsamp2$TRIFACILITES[i]<-temp)
}

for(i in which(is.na(acsrepsamp2$TRIFACILITES))){
  try(temp<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==acsrepsamp2$PUMA[i]),n=1,type="random"),tri.sp, byid=TRUE)))
  try(acsrepsamp2$TRIFACILITES[i]<-temp)
}

temp<-ddply(acsrepsamp2, .(PUMA), summarize,TRIFAC=mean(TRIFACILITES))

colnames(temp)[1]<-"PUMACE10"
gfort<-join(gfort, temp, by="PUMACE10")
head(gfort)
ggplot(gfort)+geom_polygon(aes(x=long,y=lat,group=group, fill=TRIFAC))+theme_bw()

##Distance to TRI Facility.

head(acsrepsamp2)

acsrepsamp2$FARMMARK.external<-acsrepsamp2$FARMMARK
acsrepsamp2$PM25ARM.external<-acsrepsamp2$PM25ARM
acsrepsamp2$TRIFACILITIES.external<-acsrepsamp2$TRIFACILITES

acsrepsamp2)
ncol(acsrepsamp2)
acsrepsamp2<-acsrepsamp2[,1:223]
colnames(acsrepsamp2)[which(colnames(acsrepsamp2)%in%c("FARMMARK","PM25ARM","TRIFACILITES"))]<-c("FARMMARK.external","PM25ARM.external","TRIFACILITIES.external")


summary(acsrepsamp2$ST)
summary(filter(acsrepsamp2, ST=="Virginia/VA")$PM25ARM.external)
levels(acsrepsamp2$ST)

eqfile<-read.csv(textConnection(getURL("http://www2.census.gov/geo/docs/reference/puma/2010_PUMA_Names.txt")),colClasses="character")
head(eqfile)
colnames(eqfile)<-c("ST","PUMA","PUMANAMES")
eqfile$ST<-as.numeric(eqfile$ST)
eqfile<-acsrecode.df("ST",eqfile)
acsrepsamp2<-join(acsrepsamp2, eqfile, by=c("ST","PUMA"))

require(sp)
counties<-readOGR("counties","cb_2014_us_county_5m")
counties<-subset(counties, counties$STATEFP%in%STATES_INCLUDE_FIPS)
load("35019-0001-Data.rda")
counties$COUNTYFP<-as.numeric(counties$COUNTYFP)
crimes<-subset(da35019.0001, da35019.0001$FIPS_ST%in%STATES_INCLUDE_FIPS)
head(crimes)
crimes$STATEFP<-crimes$FIPS_ST
crimes$COUNTYFP<-crimes$FIPS_CTY
filter(crimes, STATEFP==53)
crimecount<-join(counties@data, crimes, by=c("STATEFP","COUNTYFP"))
counties@data<-crimecount[,c(colnames(counties@data),"GRNDTOT","DRUGTOT","DRGPOSS","MJPOSS","WEAPONS")]

counties<-spTransform(counties, proj4string(geoshape))

countpops<-read.csv("PEP_2014_PEPANNRES.csv")
countpops$GEOID<-countpops$GEO.id2

counties@data<-plyr::join(counties@data,countpops, by="GEOID")
diabetes<-read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/15wBfoHlhun60Y21BxSic3Vpf3tig6NmR14GlF9T5_JM/pub?gid=1458307865&single=true&output=csv")))
diabetes<-diabetes[,c("FIPS.Codes","percent")]
colnames(diabetes)<-c("GEOID",":DIABETES")
counties@data<-join(counties@data, diabetes, by="GEOID")

acsrepsamp2$COUNTY.external<-NA
acsrepsamp2$GRNDTOT.external<-NA
acsrepsamp2$DRUGTOT.external<-NA
acsrepsamp2$DRGPOSS.external<-NA
acsrepsamp2$MJPOSS.external<-NA
acsrepsamp2$WEAPONS.external<-NA
acsrepsamp2$DIABETES.external<-NA
plot(geoshape)
plot(counties)
acsrepsamp2$st2<-acsrepsamp2.strip$ST
for(i in 1:nrow(acsrepsamp2)){
  temp<-over(spsample(subset(geoshape, paste(geoshape$STATEFP10,geoshape$PUMACE10)==paste(acsrepsamp2$st2[i],acsrepsamp2$PUMA[i])),n=100,type="random"),counties)[,c("NAME","GRNDTOT","DRUGTOT","DRGPOSS","MJPOSS","WEAPONS",":DIABETES")]
  temp<-sample_n(temp,1, weight=temp$respop72013) 
  acsrepsamp2$COUNTY.external[i]<-temp$NAME
  acsrepsamp2$GRNDTOT.external[i]<-temp$GRNDTOT
  acsrepsamp2$DRUGTOT.external[i]<-temp$DRUGTOT
  acsrepsamp2$DRGPOSS.external[i]<-temp$DRGPOSS
  acsrepsamp2$MJPOSS.external[i]<-temp$MJPOSS
  acsrepsamp2$WEAPONS.external[i]<-temp$WEAPONS
  acsrepsamp2$DIABETES.external[i]<-as.numeric(as.character(temp$":DIABETES"))
}
head(acsrepsamp2)
acsrepsamp2.strip[,c("FARMMARK.external","PM25ARM.external","TRIFACILITIES.external","COUNTY.external","PUMANAMES","GRNDTOT.external","DRUGTOT.external","DRGPOSS.external","MJPOSS.external","WEAPONS.external","DIABETES.external")]<-acsrepsamp2[,c("FARMMARK.external","PM25ARM.external","TRIFACILITIES.external","COUNTY.external","PUMANAMES","GRNDTOT.external","DRUGTOT.external","DRGPOSS.external","MJPOSS.external","WEAPONS.external","DIABETES.external")]
write_sav(acsrepsamp2, "PBAF_DRAFT4.sav")
write_sav(acsrepsamp2.strip, "PBAF_DRAFT4strip.sav")

summary(acsrepsamp2$DIABETES.external)
acsrepsamp2temp<-acsrepsamp2
acsrepsamp2<-droplevels(acsrepsamp2)
ddply(.data=acsrepsamp2,.(ST), mean=mean(DIABETES.external,na.rm=T),summarize)


?ddply
temp
tail(acsrepsamp2$DIABETES.external)
acsrepsamp2$DIABETES.external
rm(acsrepsamp)
rm(acsselect)
rm(cbook)
rm(gfort)
rm(match)
rm(temp)
rm(est)
rm(test2)
rm(tri)
rm(validweights)
rm(anair)

rm(fm)
rm(wmem)
rm(sampset)




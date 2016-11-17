library(RCurl)
cbook<-read.csv(textConnection(getURL("https://docs.google.com/spreadsheets/d/1mfZy2fap2eqLIH9_vrbFs35hQMfGnQdKp39XEjjPpVw/pub?gid=114600137&single=true&output=csv")))
head(cbook)
head(acsselect)
#Distance to Farmers Market. Assumes individuals have a random address within a PUMA. As such, this variable should be used to compare geographies.
require(plyr)
require(dplyr)
library(stringr)                  
acsrecode<-function(X){
filton<-filter(cbook, Variable==X)
data.frame("level"=str_trim(as.character(unlist(unlist(strsplit(as.character(filton$Code), ".",fixed=T))))) %>% .[seq(1,length(.),2)],"label"=str_trim(unlist(unlist(strsplit(as.character(filton$Code), ".",fixed=T)))) %>% .[seq(2,length(.),2)],'variable'=X)
}

acsrecode("COW")


acsrecode.df<-function(X,INframe){
  INframe[X]<-join(data.frame("level"=as.factor(INframe[,X])),acsrecode(X), by='level',type="left",match="first")$label
  print(summary(INframe[X]))
  INframe
}

acss.code<-acsrecode.df("BROADBND",acsselect)

#make sure to change first merge to inner join

require(survey)

options( "survey.replicates.mse" = TRUE ) 

acsrepw<-svrepdesign(repweights = acsselect[204:283], weights =acsselect$PWGTP, combined.weights = TRUE, type="JK1",scale=4/80,rscales=rep(1,80),data=acsselect)

#test_res<-sample_n(acsselect,70000,weight=acsselect$PWGTP)

#mean(as.numeric(test_res$SEX))-1

wmem<-weights(acsrepw)
NFORSET=20000
bacon=list()
for(i in c(1:8,10:65,67:68,70:73,75:76,78:80)){
sampset<-which(wmem[,i]!=0)
bacon[[i]]<-sample_n(acsselect[sampset,],round(NFORSET/length(c(1:8,10:65,67:68,70:73,75:76,78:80))),weight=wmem[sampset,i],replace=FALSE)
}
bacon<-do.call("rbind",bacon)
head(bacon)
head(bacon$PUMA)
colnames(bacon)
bacon2<-bacon[,-c(6,128:202,204:283,284,286,288,290,289,290,286,288,434:513,382:433)]
colnames(bacon2)
colnames(bacon2)
summary(bacon2)
levels(as.factor(bacon2$MIGPUMA))
for(VARNAME in colnames(bacon2)){
  print(VARNAME)
  if(VARNAME%in%c("MIGPUMA","POWPUMA","PUMA","SERIALNO")==FALSE){
  if(class(bacon2[,VARNAME])=="character"){
  bacon2=acsrecode.df(VARNAME,bacon2)
  }
  }
}


library(haven)


require(rgal)
fm<-read.csv("farmersmarket.csv")
fm<-filter(fm, State%in%c("Oregon","Washington","Idaho"))
fm<-fm[which(is.na(fm$x)==FALSE),]
coordinates(fm)<-c("x","y")
gDistance(subset(fm,))
head(bacon2$PUMA)
colnames(geoshape@data)
head(geoshape@data)
colnames(bacon2)
head(bacon2)
geoshape<-spTransform(geoshape, CRS("+proj=longlat +datum=WGS84"))
proj4string(fm)<-proj4string(geoshape)
require(rgeos)
bacon2$FARMMARK<-NA
proj4string(fm)
#distancetoafarmersmmarket


for(i in 1:nrow(bacon2)){
try(bacon2$FARMMARK[i]<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==bacon2$PUMA[i]),n=1,type="random"),fm, byid=TRUE)))
}
for(i in which(is.na(bacon2$FARMMARK))
){
  try(bacon2$FARMMARK[i]<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==bacon2$PUMA[i]),n=1,type="random"),fm, byid=TRUE)))
}
colnames(geoshape@data)
test<-ddply(bacon2, .(PUMA), summarize,avgdist=mean(FARMMARK))
head(test)
colnames(test)[1]<-"PUMACE10"
require(ggplot2)
geoshape$id<-row.names(geoshape)
gfort<-fortify(geoshape)
gfort<-join(gfort,geoshape@data, by="id")

gfort<-join(gfort, test, by="PUMACE10")

ggplot(gfort)+geom_polygon(aes(x=long,y=lat, group=group, fill=avgdist))+theme_bw()

anair<-read.csv("annual_all_2014.csv")

anair<-filter(anair, State.Name%in%c("Washington","Idaho","Oregon"), Parameter.Name=="PM2.5 - Local Conditions")
head(anair)

coordinates(anair)<-c("Longitude","Latitude")
proj4string(anair)<-CRS("+proj=longlat +datum=NAD27")
anair<-spTransform(anair,proj4string(geoshape))
head(anair@data)
bacon2$PM25ARM<-NA
for(i in 1:nrow(bacon2)){
try(temp<-gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==bacon2$PUMA[i]),n=1,type="random"),anair, byid=TRUE))
    
try(bacon2$PM25ARM[i]<-anair@data$Arithmetic.Mean[which(temp==min(temp)) %>% sample(.,size=1)])
}

for(i in which(is.na(bacon2$PM25ARM))){
  try(temp<-gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==bacon2$PUMA[i]),n=1,type="random"),anair, byid=TRUE))
  
  try(bacon2$PM25ARM[i]<-anair@data$Arithmetic.Mean[which(temp==min(temp)) %>% sample(.,size=1)])
}

test2<-ddply(bacon2, .(PUMA), summarize,PM25=mean(PM25ARM))
head(test2)
colnames(test2)[1]<-"PUMACE10"
gfort<-join(gfort, test2, by="PUMACE10")
head(gfort)
ggplot(gfort)+geom_polygon(aes(x=long,y=lat,group=group, fill=PM25))+theme_bw()
anair@data$Arithmetic.Mean

##PM is a random sample from PM tests in 2014 from a site closest to the individual. This results in a mean PM measure for each individual in each county. 

head(anair)
proj4string(geoshape)
?is.projected

#
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
bacon2<-join(bacon2, temp, by="PUMA")
#Raw count of the number of Toxic Release Inventory facilities within the PUMA

for(i in 1:nrow(bacon2)){
  try(temp<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==bacon2$PUMA[i]),n=1,type="random"),tri.sp, byid=TRUE)))
  try(bacon2$TRIFACILITES[i]<-temp)
}

for(i in which(is.na(bacon2$TRIFACILITES))){
  try(temp<-min(gDistance(spsample(subset(geoshape, as.character(geoshape$PUMACE10)==bacon2$PUMA[i]),n=1,type="random"),tri.sp, byid=TRUE)))
  try(bacon2$TRIFACILITES[i]<-temp)
}

temp<-ddply(bacon2, .(PUMA), summarize,TRIFAC=mean(TRIFACILITES))

colnames(temp)[1]<-"PUMACE10"
gfort<-join(gfort, temp, by="PUMACE10")
head(gfort)
ggplot(gfort)+geom_polygon(aes(x=long,y=lat,group=group, fill=TRIFAC))+theme_bw()

##Distance to TRI Facility.

head(bacon2)

bacon2$FARMMARK.external<-bacon2$FARMMARK
bacon2$PM25ARM.external<-bacon2$PM25ARM
bacon2$TRIFACILITIES.external<-bacon2$TRIFACILITES

bacon2)
ncol(bacon2)
bacon2<-bacon2[,1:223]
colnames(bacon2)[which(colnames(bacon2)%in%c("FARMMARK","PM25ARM","TRIFACILITES"))]<-c("FARMMARK.external","PM25ARM.external","TRIFACILITIES.external")


write_sav(bacon2, "PBAF_DRAFT1.sav")



summary(lm(TRIFACILITES~HHL, data=bacon2))
#How is family income related to year of naturalization into the US?
data1<-read_sav("FILE")
data1<-bacon2
mean(data1$FINCP,na.rm=T)
sd(data1$FINCP,na.rm=T)
quantile(data1$FINCP, .5, na.rm=T)
require(ggplot2)



#Demo 1
forstudy<-dplyr::filter(data1, is.na(LAPTOP)==FALSE, is.na(FINCP)==FALSE)

summary(forstudy$)/nrow(forstudy)*100

min(forstudy$FINCP)
max(forstudy$FINCP)
median(forstudy$FINCP)
forstudy$lowINCOME<-forstudy$FINCP<=median(forstudy$FINCP)
xtabs(~forstudy$lowINCOME+as.character(forstudy$LAPTOP))

?xtabs

While the income of families in the sample range between -6500 and 910000, we read in the data dictionary that these values were rounded so as to protect personal identities. Thus, there may actually be incomes that are much greater or less than these values. The overall average income for the sample is $85150--though about half of all respondents' families earn less than $66000. 
1301/(6710+1301)
284/(7717+284)
The 90% of all respondents in the dataset own laptops. This high rate of ownership in the dataset is likely the result of the fact that respondents could report anyone in their households owned laptops. Thus, the rate of laptop ownership in the sample is far greater than found in other places owing to the nature of the question.




summary(forstudy$LAPTOP)/nrow(forstudy)*100



ggplot(data1)+geom_histogram(aes(x=FINCP), fill="blue")+theme_bw()

ggplot(filter(data1,is.na(LAPTOP)==FALSE))+geom_bar(aes(x=LAPTOP))+theme_bw()+ggtitle("Do you own a Laptop?")

ggplot(dplyr::filter(data1, is.na(LAPTOP)==FALSE, is.na(FINCP)==FALSE))+geom_bar(aes(x=LAPTOP, fill=FINCP<=median(data1$FINCP, na.rm=T)), alpha=.7)+theme_bw()+xlab("Household owns a Laptop")+ylab("Frequency")+ggtitle("Ownership of Laptops")+theme(text=element_text(size=20,family="Times"))+scale_fill_discrete("Income", labels=c("> median", "< median"))
                                                                                                                                                        
summary(acsrepsamp2$ST)

summary(data1$VEH=="No vehicles")
+scale_fill_gradient()

acsrecode

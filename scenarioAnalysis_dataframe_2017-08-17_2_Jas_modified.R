
#install.packages("foreach")
#install.packages("doParallel")

library(lubridate)
library(dplyr)
library(ggplot2)

library(zoo)
library(latticeExtra)
library(polynom)
library(car)
library(Hmisc)
library(reshape)
library(DEoptim)
library(coda)
library(dream)
library(hydromad)
library(lubridate)
library(tidyr)
library(dygraphs)
library(xts)

# ----------------------------------
# 0. Import relevant data
# ----------------------------------
#inflow
setwd("C:/Users/jthom/OneDrive - The University of Melbourne/ch2_conceptualModelling/JasModel")
Sys.getenv("TZ")
Sys.setenv(TZ="Australia/Brisbane")
Sys.getenv("TZ")
setwd("C:/Users/jthom/OneDrive - The University of Melbourne/ch2_conceptualModelling/JasModel")
WeatherData=read.table(file="Melbourne_2001-2010.txt",header=T,stringsAsFactors = FALSE) #Read rainfall data in
WeatherData <- unite(WeatherData, col=dateTime, Date:Minutes, sep=" ", remove=FALSE)
WeatherData$dateTime <- dmy_hm(WeatherData$dateTime, tz="Australia/Brisbane") #convert date to posix
WeatherData$JDay <- yday(WeatherData$Date) #create julian day
R=WeatherData$Rainfall #create rainfall vector
length(R)#check length of rainfall vector
sum(R) #calculate sum of rainfall
RainfallDate <- WeatherData$dateTime
length(RainfallDate)
#ET
ET=read.table(file="Melbourne_2000_ET.txt",header=T,stringsAsFactors = FALSE)
sum(ET$ETPM) #check total ET
ET=ET$ET.6min # Create ET vector (mm/6min)
ETtest=rep(ET,each=240) #replicate 6min ET to 6 min timestep
length(ETtest) #Check length ET == length Rainfall
sum(ETtest) # Check total ET (1214.746)
etveg=ETtest/100 # Convert ET to dm/6min
sum(etveg)
#without plants (1.3mm/day) 
etbare=rep(1.3/100,length(R))#dm/day
etbare=etbare/240#dm/6min


# ----------------------------------
# 1. Outline static parameters
# ----------------------------------
#Define Filter Profile Variables
d.pond.m=0.1 #m
d.sandyloam.m=1.5 #m
d.sand.m=0.3 #m
d.scoria.m=0.2 #m
soil.strata.m=c(0.1, 1.6, 0.9, 2.1) #m
strata.Ksu.mm.h=c(1, 1.5, 2.5, 0)#mm/hr 
medium="loamy sand"

#Define Raingarden model Variables
Af=10
Pf=14
Hf=2 #m
Ap=10 # m
Hp=14 # m 
isveg=TRUE
Ho=0 #m from base to top of saturated zone
Vstart=100*Af #Litres - Capacity is approx 240 L per m2
carea=c(1000) # meter square ratio 2.5 % 
Ksu.mm.h=c(1.5) # mm/hr Lined
outlet.rate.L.h = c(0)#L/hr medium="sand"
adj.tree.canopy.area=0
Ksf=2.5 #dm/hr 
filtertakerate=Af*Ksf

#Define tree transpiration variables
kc <- c(1,2) #Set maximum crop coefficient
decid=c(FALSE) #set deciduousness
fieldCap <- 0.367 #field capacity as a decimal 
wiltP <- 0.194 #wilting point as a decimal 
permWiltP <- 0.07 #permanent wilting point as a decimal 
t0 <- 274 #beginning of first leaf flush
t1 <- 366 #end of first flush, beginning second flush 
t2 <- 32  #end of second flush
t3 <- 117  #beginning of leaf senescence
t4 <- 153 #end of leaf senescence

# ----------------------------------
# 2. Create a grid of changing scenarios
# ----------------------------------

D0_LB_ULS <- expand.grid(code="D0_LB_ULS", 
                         outlet.rate=0, carea=c(250),
                         strata.ksu=c(1.5,2.5), ksu=0,
                         kc=c(1), decid=c(FALSE), Ho=c(0,1))
#D0_ULB_LS <- expand.grid(code="D0_ULB_LS", 
#                         outlet.rate=0, carea=c(250),
#                         strata.ksu=0, ksu=c(1.5,2.5),
#                         kc=c(1), decid=c(TRUE, FALSE))
#D0_ULB_ULS <- expand.grid(code="D0_ULB_ULS", 
#                          outlet.rate=0, carea=c(250),
#                          strata.ksu=0, ksu=c(0,1.5,2.5),
#                          kc=c(1), decid=c(TRUE, FALSE))
#D0_ULB_ULS$ksu <- D0_ULB_ULS$strata.ksu
#
#D1_LB_ULS <- expand.grid(code="D1_LB_ULS", 
#                         outlet.rate=1, carea=c(250),
#                         strata.ksu=c(0,1.5,2.5), ksu=0,
#                         kc=c(0.5,1), decid=c(TRUE, FALSE))


#scenarios <- rbind(D0_LB_ULS, D0_ULB_LS, D0_ULB_ULS, D1_LB_ULS)

#scenarios <- as.data.frame.array(scenarios)

scenarios <- D0_LB_ULS

nrow(scenarios)
scenarios$scenario <- seq(1, length=nrow(scenarios))
scenariosAll <- scenarios


# ----------------------------------
# 2. Create open results dataframes for writing model data
# ----------------------------------
results_final=NULL
results=NULL
results2=NULL
results3=NULL
results4=NULL
inflowEventsTest=NULL
cropcoeff <- NULL


# ----------------------------------
# 3. Start modelling scenarios 
# At this stage deciduousness is completed in a second for loop. 
# Would like this to be included 
# ----------------------------------


for (a in 1:nrow(scenarios)) {
  Ho.new <- scenarios$Ho[a]
  outletRate.new <- scenarios$outlet.rate[a]
  carea.new <- scenarios$carea[a]
  ksu.new <- scenarios$ksu[a]
  kc.new <- scenarios$kc[a]
  decid.new <- scenarios$decid[a]
  cropcoeff <- rep(kc.new, length(WeatherData$JDay))
  if (decid.new==TRUE) {
    for (g in 1:length(WeatherData$JDay)) {
      if(WeatherData$JDay[g] >= t0 & WeatherData$JDay[g] <= t1) { #function for kcbIni
        cropcoeff[g] <- kc.new * (WeatherData$JDay[g]-t0) * (1/(t1-t0+t2))
      } else {
        if (WeatherData$JDay[g] <= t2) { #function for kcbMid
          cropcoeff[g] <- kc.new * (WeatherData$JDay[g]*(1/(t1-t0+t2))) + kc.new * ((t1-t0)*(1/(t1-t0+t2)))
        } else {
          if (WeatherData$JDay[g] >= t3 & WeatherData$JDay[g] < t4) {
            cropcoeff[g] <- kc.new * (t4-WeatherData$JDay[g]) * (1/(t4-t3))
          } else {
            if (WeatherData$JDay[g] >= t4 & WeatherData$JDay[g] < t0) {
              cropcoeff[g] <- 0
            } else {
              cropcoeff[g] <- kc.new
            }
          }
        }
      }
    }
  }
  
  et <- etveg #=====
  
  myprofile=filter.profile(d.pond.m, d.sandyloam.m, d.sand.m, d.scoria.m, soil.strata.m, strata.Ksu.mm.h)
  filtr.prfile=myprofile$profile 
  
  runoffcoeff=0.86 #86 % of rainfall converted to runoff
  runoff=R*(carea.new) * runoffcoeff # in L - Note R is in mm and carea is in m2
  inflow=runoff #inflow into trench
  
  model_data=gardenmodel(inflow, et, Af, Pf, Hf, Ap, Hp, isveg, Ho.new, Vstart, carea.new, ksu.new, outletRate.new, filtr.prfile, 
                         adj.tree.canopy.area,medium, Ksf,cropcoeff, fieldCap, wiltP) 
  
  #Save each relevant variable from model as a vector ----
  inflow2=model_data$budget$inflow 
  adwp=model_data$budget$adwp 
  Qexf=model_data$budget$Qexf 
  et2=model_data$budget$et 
  out=model_data$budget$out 
  over=model_data$budget$over 
  pond=model_data$budget$store.pond 
  filter=model_data$budget$store.filter 
  height=model_data$budget$height.dm 
  flow2filter=model_data$budget$flow2filter 
  filtertakerate=model_data$budget$filtertakerate 
  max2filter=model_data$budget$max2filter
  LBS=model_data$budget$LBS
  LS=model_data$budget$LS
  LB=model_data$budget$LBS
  Ksu.adj=model_data$budget$Ksu.adj
  #Nconcout=model_data$budget$Nconcout 
  #Pconcout=model_data$budget$Pconcout 
  #TSSout=model_data$budget$TSSout
  #JAS ADDED
  height.dm=model_data$budget$height.dm
  pond2filterorover=model_data$budget$pond2filterorover
  
  
  
  #Create results data frame with all relevant vectors for EACH scenario   ----      
  #results=data.frame(inflow,adwp,Qexf,et,out,over,pond,filter,height,flow2filter,filtertakerate,max2filter,Nconcout,Pconcout,TSSout)
  results=data.frame(inflow2, Qexf, et2, out, over, pond, filter, height, filtertakerate, 
                     max2filter, pond2filterorover, flow2filter, cropcoeff, height.dm)
  #Add relevant scenario information
  results$Ksu=ksu.new
  results$ksu.strata=strataKsu.new
  results$Ksf=Ksf
  #results$Hp=Hpnew
  results$carea=carea.new
  results$decid=decid.new
  results$kc=kc.new
  results$outletRate=outletRate.new
  
  #Create summary dataframe for EACH scenario ----
  #Calculate summary information           
  over1=sum(results$over)/sum(results$inflow2)*100
  etsum1=sum(results$et2)/sum(results$inflow2)*100
  outsum1=sum(results$out)/sum(results$inflow2)*100
  Qexfsum1=sum(results$Qexf)/sum(results$inflow2)*100
  flow2filter1 = sum(results$flow2filter)/sum(results$inflow2)*100
  #Save as dataframe
  results2=data.frame(over1,etsum1,outsum1, flow2filter1, Qexfsum1)
  #Add relevant scenario information 
  results2$Ksu=ksu.new
  results2$ksu.strata=strataKsu.new
  results2$Ksf=Ksf
  #results2$Hp=Hpnew
  results2$carea=carea.new
  results2$decid=decid.new
  results2$kc=kc.new
  results2$outletRate=outletRate.new
  results2$Ho=Ho.new
  results2$LBS=LBS
  results2$LS=LS
  results2$LB=LB
  results2$Ksu.adj=Ksu.adj
  #results2$cover=covernew
  
  #Create data frame with ALL Results from EACH scenario ----
  results3=rbind(results3, results)
  results_final=rbind(results_final,results2)
  
  #Create event dataframe for rainfall events in EACH scenario ----
  #create dataframes for each variable with time element
  inflow2.d=data.frame(RainfallDate, inflow2)
  Qexf.d=data.frame(RainfallDate, Qexf)
  et2.d=data.frame(RainfallDate, et2)
  out.d=data.frame(RainfallDate, out)
  over.d=data.frame(RainfallDate, over)
  pond.d=data.frame(RainfallDate, pond)
  filter.d=data.frame(RainfallDate, filter)
  flow2filter.d=data.frame(RainfallDate, flow2filter)
  filtertakerate.d=data.frame(RainfallDate, filtertakerate)
  max2filter.d=data.frame(RainfallDate, max2filter)
  cropcoeff.d=data.frame(RainfallDate, cropcoeff)
  
  #create zoo objects for each variable
  inflow2.z=zoo(inflow2.d$inflow2, order.by=RainfallDate)
  Qexf.z=zoo(Qexf.d$Qexf, order.by=RainfallDate)
  et2.z=zoo(et2.d$et2, order.by=RainfallDate)
  out.z=zoo(out.d$out, order.by=RainfallDate)
  over.z=zoo(over.d$over, order.by=RainfallDate)
  pond.z=zoo(pond.d$pond, order.by=RainfallDate)
  filter.z=zoo(filter.d$filter, order.by=RainfallDate)
  flow2filter.z=zoo(flow2filter.d$flow2filter, order.by=RainfallDate)
  filtertakerate.z=zoo(filtertakerate.d$filtertakerate, order.by=RainfallDate)
  max2filter.z=zoo(max2filter.d$max2filter, order.by=RainfallDate)
  cropcoeff.z=zoo(cropcoeff.d$cropcoeff, order.by=RainfallDate)
  
  #Determine event info for each variable
  inflowEvents60=eventseq(inflow2.z, thresh = 0.0, mingap = 60, mindur = 1, inthresh = 0.0, indur = 1)
  inflow.Events=eventinfo(inflow2.z, inflowEvents60, FUN = sum)
  Qexf.Events=eventinfo(Qexf.z, inflowEvents60, FUN = sum)
  et2.Events=eventinfo(et2.z, inflowEvents60, FUN = sum)
  out.Events=eventinfo(out.z, inflowEvents60, FUN = sum)
  over.Events=eventinfo(over.z, inflowEvents60, FUN = sum)
  pond.Events=eventinfo(pond.z, inflowEvents60, FUN = sum)
  filter.Events=eventinfo(filter.z, inflowEvents60, FUN = max)
  flow2filter.Events=eventinfo(flow2filter.z, inflowEvents60, FUN = max)
  filtertakerate.Events=eventinfo(filtertakerate.z, inflowEvents60, FUN = max)
  max2filter.Events=eventinfo(max2filter.z, inflowEvents60, FUN = max)
  cropcoeff.Events=eventinfo(cropcoeff.z, inflowEvents60, FUN = max)
  
  #Calculate start, duration and number of drys days preceding event
  inflowEvents.start=inflow.Events$Time
  inflowEvents.dur = inflow.Events$Duration
  inflowEvents.depth = inflow.Events$Value
  inflowEvents.ADP = inflow.Events$PreDuration
  #Assemble into new dataframe
  inflowEvents.frame <- data.frame(mit="6 hrs", event=1:length(inflowEvents.start), start=inflowEvents.start, 
                                   duration_hrs=(inflowEvents.dur)/60, ADP_days=((inflowEvents.ADP*6)/60)/24, inflow=inflowEvents.depth, 
                                   Qexf=Qexf.Events$Value, et2=et2.Events$Value, out=out.Events$Value, over=over.Events$Value, 
                                   pond=pond.Events$Value, filter=filter.Events$Value, flow2filter=flow2filter.Events$Value, 
                                   filtertakerate=filtertakerate.Events$Value, max2filter=max2filter.Events$Value, 
                                   cropcoeff=cropcoeff.Events$Value, Ksu=ksu.new, Ksf=Ksf, Hp=Hp, carea=carea.new,
                                   decid=decid.new, kc=kc.new, outletRate=outletRate.new)
  
  #Create dataframe with ALL RAINFALL events results from ALL scenarios ----
  results4=rbind(results4, inflowEvents.frame)
  
  
}



#Check results scenarios sum to 100
results_final$check <- results_final$over1 + results_final$etsum1 + results_final$outsum1 +results_final$Qexfsum1
results_final$conditions <- "make soil.strata and strata.ksu both vectors of length 4 to match layers in filter profile"
#save to comaprison dataframe
comparisons <- rbind(comparisons, results_final)

#Save summary results, raw results, and events results.
setwd("C:/Users/jthom/OneDrive - The University of Melbourne/ch2_conceptualModelling/JasModel")
write.table(comparisons, "modelComparisons_2017-08-17.txt")
write.table(results_final,"results_daily_final_p2.txt")
write.table(results3, "raw_results_p2.txt")
write.table(results4, "events_results_p2.txt")


##########################################################
# filter.profile FUNCTION
# Creates the filter profile for use in the raingarden model
#
# d.pond.m          depth of ponding zone above filter (in m) = Hp: zero if absent
# d.sandyloam.m     depth (in m) of sandy loam layer. Enter zero if layer is missing. Sandy loam must be at surface
# d.sand.m          depth (in m) of sand layer. Enter zero if layer is missing
# d.scoria.m        depth (in m) of scoria layer. Enter zero if layer is missing. Scoria must be deepest. Total of all three depths = depth of filter Hf
# soil.strata.m     vector of depths (in m from pond surface) at which surrounding soils change in their infiltration capacity. Final value = bottom of filter
# strata.Ksu.mm.h   vector (same length as soil.strata.m) of infiltration rates for each stratum (in mm.h) (if lined system single value = 0) #function returns a
#                   data.frame with soil profile (in dm) suitable for the Ksu.profile argument of gardenmodel.R

filter.profile <- function(d.pond.m, d.sandyloam.m, d.sand.m, d.scoria.m, soil.strata.m, strata.Ksu.mm.h) {
  
  if (length(soil.strata.m) != length(strata.Ksu.mm.h)) 
    stop("soil.strata.m and strata.Ksu.mm.h vectors must be the same length", call. = FALSE)
  
  if (sum(c(d.pond.m, d.sandyloam.m, d.sand.m, d.scoria.m)) != soil.strata.m[length(soil.strata.m)])
    stop("total depth of the three media + pond must equal total depth of soil strata", call. = FALSE)
  
  if (sum(diff(soil.strata.m <= 0) > 0)) 
    stop("soil.strata.m should be a vector of increasing depths", call. = FALSE)
  
  filter.layers <- c(d.pond.m, d.sandyloam.m, d.sand.m, d.scoria.m)
  porosities <- c(1, 0.158, 0.147, 0.132)[filter.layers > 0]
  n.media <- length(porosities)
  filter.layers <- unique(cumsum(filter.layers[filter.layers > 0]))
  
  d <- unique(c(filter.layers, soil.strata.m))
  d <- d[order(d)]
  x <- data.frame(d.dm = d * 10, porosity = NA, Ksu.dm.h = NA)
  
  for (i in 1:length(d)) {
    x$porosity[i] <- porosities[which(filter.layers >= d[i])[1]]
    x$Ksu.dm.h[i] <- strata.Ksu.mm.h[which(soil.strata.m >= d[i])[1]]/100
  }
  
  # Assume pond is lined
  x$Ksu.dm.h[x$porosity == 1] <- 0
  list(profile = x, top.medium = c("loamy sand", "sand", "gravel (scoria)")[which(c(d.sandyloam.m, d.sand.m, d.scoria.m) > 0)[1]])
}

###########################################
# RAIN GARDEN MODEL
# inflow    vector of hourly inflows in L (no date or time fields)
# et        vector of hourly evapotranspiration values in mm/h (sets to zero if isveg is False)
# Af        Filter area (sq m)
# Pf        Filter perimeter (m)
# Hf        Filter depth (m)
# Ap        Ponding area (sq m)
# Hp        Ponding depth (m) 
# isveg     Is system vegetated 
# Ho        Distance from base to invert of underdrain outlet or orifice of of standpipe outlet (m)
#           If no underdrain or standpipe, set Ho to Hf + Hp
# Vstart    Vol of water in system at t1 in L
# carea     As in compile.inflow(sq m)
# Ksu.mm.h  Underlying infiltration rate: set to zero if lined bottom
# outlet.rate.L.h   Zero if no underdrain or standpipe, = Ksu if standpipe
# filtr.prfile  data.frame of depths, porosities of each filter layer, Ksus of surrounding soil layers prepare using the function filter.profile() [$profile]
# adj.tree.canopy.area  Canopy area (sq m) of trees with canopy edges < 3m from an at least partly unlined raingarden
# medium    filter medium, must be sandy loam, sand, or gravel (scoria) decides Ksf = 0.5,2.5, 36(dm/h)
#           respectively (assumes halving of K for loamy sand and sand over time) but for tim's scenario, Ksf of loamy sand = 1.5
#           must correspond to top layer in filter.profile - to be safe use filter.profile$top.medium
#Ksf        Filter media conductivity (dm/hr)

gardenmodel <- function(inflow, et, Af, Pf, Hf, Ap, Hp, isveg = TRUE, Ho, Vstart, carea, Ksu.mm.h, outlet.rate.L.h, filtr.prfile, adj.tree.canopy.area,Ksf,cropcoeff) {
  
  # Input validation
  #if (is.na(match(medium, c("loamy sand", "sand", "gravel (scoria)")))) 
  #  stop("Medium must be one of 'loamy sand', 'sand' or 'gravel (scoria)'", .call = FALSE)
  
  #if (isveg & medium == "gravel (scoria)") 
  # stop("Sorry: the calculator does not allow vegetated gravel systems", .call = FALSE)
  
  Ksu.dm.h <- Ksu.mm.h /100
  
  Ksu.dm.h <- Ksu.dm.h /10 #dm per 6 min 
  
  
  if (Ksu.dm.h == 0 & filtr.prfile$Ksu.dm.h[nrow(filtr.prfile)] > 0) 
    stop(paste("Value of Ksu.mm.h = zero, in which case so should the deepest Ksu in filter.profile (currently ", filtr.prfile$Ksu.dm.h[dim(filtr.prfile)[1]] * 
                 100, ")", sep = ""), call. = FALSE)
  
  # If system has a tub in the bottom but unsealed sides, deepest layer in filter.profile should have Ksu = 0
  lined.bottomandsides <- ifelse(sum(filtr.prfile$Ksu.dm.h) == 0 & Ksu.dm.h == 0, TRUE, FALSE)
  lined.side <- ifelse(sum(filtr.prfile$Ksu.dm.h) == 0, TRUE, FALSE)
  
  # Conversions to dm
  Ho <- Ho * 10
  Hf <- Hf * 10
  Hp <- Hp * 10
  Pf <- Pf * 10
  
  # Convert to square dm
  Afd <- Af * 100 
  
  Ho.true <- Ho
  
  # Hourly evapotranspiration
  # if (isveg) {
  #  if (!lined.bottomandsides)
  #   et <- et * (Af + 0.5 * adj.tree.canopy.area)
  #else
  # et <- et * Af
  #} else {
  # if (!lined.bottomandsides)
  #  et <- et * (0.5 * adj.tree.canopy.area)
  #else
  # et <- rep(0, nrow(inflow))
  #}
  
  # Pollutant parameters
  # if (isveg) {
  #  N1 <- c(0.59, 0.72, NA)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  # P1 <- c(0.74, 0.74, NA)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  #TSS1 <- c(0.99, 0.99, NA)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  #} else {
  # N1 <- c(-1.01, 0.32, 0.39)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  #P1 <- c(0.86, 0.86, 0.2)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  ### made up 0.2 waiting on TDF advice
  #TSS1 <- c(0.99, 0.99, 0.2)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  ### made up 0.2 waiting on TDF advice
  #  }
  
  #Define filter media hydraulic conductivity (dm/hr)
  
  #Ksf <- c(1.5, 2.5, 36)[match(medium, c("loamy sand", "sand", "gravel (scoria)"))]
  Ksf=Ksf/10 #dm/6min
  
  
  # Number of layers
  nlayers <- nrow(filtr.prfile)
  
  # Outlet rate
  outlet.rate <- outlet.rate.L.h/10 #L per 6 min
  
  # Height of filter profile in decimeters
  filtr.prfile$ht.dm <- c(filtr.prfile$d.dm[1], diff(filtr.prfile$d.dm))
  
  # Precalculate the volume (L) per decimetre (dm) depth. Used to calculate depth from volume
  filtr.prfile$voidarea <- Afd * filtr.prfile$porosity
  
  # Volume of filter layers (L)
  filtr.prfile$vol <- filtr.prfile$ht.dm * filtr.prfile$voidarea
  
  # Cumulative volume of filter layers from base to top
  filtr.prfile$cumvol <- rev(cumsum(rev(filtr.prfile$vol)))
  
  # Increase volume of pond if area of pond is bigger than filter area
  if (Ap > Af)
    filtr.prfile$vol[filtr.prfile$porosity == 1] <- filtr.prfile$vol[filtr.prfile$porosity == 1] * Ap / Af
  
  # Volume of ponding
  pondV <- sum(filtr.prfile$vol[filtr.prfile$porosity == 1])
  
  # Filter volume
  filterV <- sum(filtr.prfile$vol) - pondV
  
  # Setup filter layer parameters (height, depth to outlet etc...)
  if (Ho >= Hf + Hp) {
    # Ho adjusted as fix for field capacity
    filtr.prfile$ho.dm <- filtr.prfile$ht.dm
    # true Ho
    filtr.prfile$ho.true.dm <- filtr.prfile$ht.dm
  } else {
    # Calculate the height to the layer, is equal to layer depth is outlet not in layer        
    cumdepths <- rev(cumsum(rev(filtr.prfile$ht.dm)))
    filtr.prfile$ho.dm <- pmax(pmin(Ho - cumdepths, 0) + filtr.prfile$ht.dm, 0)
    filtr.prfile$ho.true.dm <- pmax(pmin(Ho.true - cumdepths, 0) + filtr.prfile$ht.dm, 0)
  }
  
  # Sub outlet water volume
  suboV.true <- sum(filtr.prfile$voidarea * filtr.prfile$ho.true.dm)
  
  # Initial volume of filter below outlet pipe in L volume of water in each section in Litres
  Vp0 <- max(Vstart - filterV, 0)
  Vf0 <- min(filterV, Vstart)
  
  # Inflow matrix, remove the first four columns (these are date/time vars)
  len <- length(inflow)
  
  # N removal estimates (and adwp.trigger and dry.penalty are adjustments for adwp) assumes hourly timestep
  adwp <- 0
  #adwp.trigger <- ifelse(isveg, 4 * 24, 5 * 24)
  #dry.penalty <- ifelse(isveg, 0.03, 0.02)
  
  # Setup output
  budget <- list()
  
  # Main loop
  for (i in 1:len) {
    if (i == 1) {
      Vfprev <- Vf0
      Vpprev <- Vp0
    } else {
      Vfprev <- budget$store.filter[i - 1]
      Vpprev <- budget$store.pond[i - 1]
    }
    
    # If the previous filter volume was 0 then skip all height calcs as they are = 0
    if(Vfprev == 0) {
      filterLayerVols <- 0
      filterLayerDeps <- 0
      Hwpprev <- 0
      Hwfprev <- 0
    } else {
      if (nlayers == 1) {
        # Single layer of media, volume of the filter is the volume in the single layer
        filterLayerVols <- Vfprev
        filterLayerDeps <- Vfprev / filtr.prfile$voidarea
        Hwfprev <- filterLayerDeps
        Hwpprev <- 0
      } else {
        # Determine the volume of water in each layer
        filterLayerVols <- pmax(pmin(Vfprev - filtr.prfile$cumvol, 0) + filtr.prfile$vol, 0)
        
        # Calculate the height of water (dm) in each layer from the volume and void area
        filterLayerDeps <- filterLayerVols / filtr.prfile$voidarea
        
        if (Hp > 0) {
          Hwfprev <- sum(filterLayerDeps[-1])
          Hwpprev <- filterLayerDeps[1] <- ifelse(pondV == 0, 0, Hp * Vpprev^(2/3) / pondV^(2/3))  #dm
          filterLayerVols[1] <- Vpprev
        } else {
          Hwfprev <- sum(filterLayerDeps)
          Hwpprev <- 0
        }
      }
    }
    
    # Calculate exfiltration into soil
    # No exfiltration if no water in filter or fully lined system
    if(Vfprev == 0 | lined.bottomandsides) {
      exfil <- 0
    } else {
      # Exfiltration from the base
      exfil <- Afd * Ksu.dm.h
      
      # If side is not lined then add exfiltration from sides (Have applied a constant of 0.3 for lateral infiltration)
      if(!lined.side) {
        exfil <- exfil + sum(pmin(filtr.prfile$Ksu.dm.h * filterLayerDeps * Pf*0.3, filterLayerVols))
      }
      
      #New limit on exfil-added after discussion with Chris and Rhys
      exfil <- min(Vfprev + Vpprev, exfil)
      
    }
    
    # Calculate evapotranspiration #Harry modified
    #et - dm/6min 
    #Changed 24/05/2014  
    #Vfprev-L Af*mm=L
    #if(Vfprev>59.58*Af){evap=((((cropcoeff-0.8)/(200-59.58))+0.502)*et[i])}else
    #if(Vfprev<=59.58*Af&Vfprev>24.43*Af){evap=((((0.8-0.3)/(59.58-24.43))+0.069)*et[i])}else
    #  if(Vfprev<=24.43*Af){evap=0}
    #evap=evap*Afd #evap in L
    
    #Latest model
    
    #Sand
    if(cropcoeff==1){
      evap=et[i]*Afd}else
      {
        
        if(Vfprev>=60*Af){evap=(((0.89096*Vfprev/(113*Af))+0.1752)*et[i])}else
          if(Vfprev<60*Af&Vfprev>=24*Af){evap=(((1.447*Vfprev/(113*Af))-0.12004)*et[i])}else
            if(Vfprev<24*Af){evap=0}
        evap=evap*Afd #evap in dm/6 min turned into L/6 min
      }
    
    # Volume that will either soak into filter or overflow
    pond2filterorover <- Vpprev + inflow[i]
    
    # Volume that filter can accept if enough space (given infiltration rate)
    filtertakerate <- Afd * Ksf 
    
    # Maximum volume that could flow into filter if enough space
    max2filter <- min(pond2filterorover, filtertakerate)
    # flow2filter <- min(max2filter, filterV - Vfprev + exfil + evap)
    flow2filter= max2filter
    #Harry changed outflow
    
    if (Vfprev + flow2filter - exfil -evap > (112.95*Af)) {
      outflow <-  (Vfprev + flow2filter-exfil-evap-(112.95*Af)) 
    } else {
      outflow <- 0
    }
    
    # Calculate overflow, then remove negative values
    overflow <- pond2filterorover - flow2filter - pondV
    overflow <- (overflow + abs(overflow)) / 2
    
    # Calculate pond store, then remove negative values
    store.pond <- Vpprev + inflow[i] - flow2filter - overflow
    store.pond <- (store.pond + abs(store.pond)) / 2
    
    # Calculate filter store, then remove negative values
    store.filter <- Vfprev + flow2filter - outflow - exfil - evap
    store.filter <- (store.filter + abs(store.filter)) / 2
    
    
    # If no flow then increment adwp counter
    adwp <- adwp + (inflow[i] == 0)
    
    # TSS, TN, TP outflow concentration calculations
    #if (medium == "loamy sand") {
    # Nred <- ifelse(adwp <= adwp.trigger, N1, N1 - dry.penalty * ceiling((adwp - adwp.trigger)/24))
    #} else {
    # Nred <- N1
    #}
    
    #if (outflow == 0 & exfil == 0) {
    # Nout <- 0
    #  Pout <- 0
    # Tout <- 0
    #} else {
    # Total water leaving
    # QoutT <- outflow + exfil
    
    # assume exfil and outflow have the same N concentration (given mobility of NOx)
    # Nout <- 2.2 * (1 - Nred)
    
    # weighted mean of conc in outflow and exfil: assumes v good P retention in exfiltration
    #Pout <- (outflow * 0.35 * (1 - P1) + exfil * 0.001) / QoutT
    
    # weighted mean of conc in outflow and exfil: assumes complete TSS retention in exfiltration
    #Tout <- (outflow * 150 * (1 - TSS1)) / QoutT
    #  }
    
    # Assign values calculated in this loop to output arrays
    budget$inflow[i] <- inflow[i]
    budget$adwp[i] <- adwp
    budget$Qexf[i] <- exfil
    budget$et[i] <- evap
    budget$out[i] <- outflow
    budget$over[i] <- overflow
    budget$store.pond[i] <- store.pond
    budget$store.filter[i] <- store.filter
    budget$height.dm[i] <- Hwfprev
    budget$flow2filter[i] <- flow2filter
    budget$filtertakerate[i] <- filtertakerate
    budget$max2filter[i] <- max2filter
    #budget$Nconcout[i] <- Nout
    #budget$Pconcout[i] <- Pout
    #budget$TSSout[i] <- Tout
    
  }  ## End of main loop
  
  ## Function output
  list(budget = budget, filtprof = filtr.prfile)
  
}
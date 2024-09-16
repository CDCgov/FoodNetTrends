#############
# Code for Running on Spaceballs
# chmod u+x mmwr_03132023.sh
# nohup ./mmwr_03132023.sh OR bash mmwr_02222023.sh OR ./mmwr_02222023.sh
# tail -f nohup.out

#############
# Code for Running on Aspen
# qsub filename.sh

#############
# Code for loading the conda environment and neccessary modules when running on Spaceballs or Linux
# conda activate brms
# module load R 
# module load make
# module load gcc
# module load boost
# module load jags
# module load rstan
# R


# Conda commands
## List All conda Environments
### conda env list
## Create a conda Environment
### conda create -n [name]insert name]
## Create a conda environment based on a frozen package set
### conda install --file=packages.txt will install those packages in your local environment.
## Delete a conda Environment
### conda remove -n rbartmachine --all
## Install an R package
### conda install -y -c conda-mlr conda-batchtools
Rprof(filename = "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/mmwroutput/profile_output.out") 

#############
# load libraries
# library(rlang)
# library(ggplot2)
library(brms)
library(tidybayes)
library(tidyverse)
# library(patchwork)
library(haven)
library(furrr)
library(batchtools)


#############
# Get Raw Data

# Set Travel Variables Here for Testing and Below for the Actual Dataset Desired
# run<-"Cx_13_Mar_2024"
# travel<-c("NO", "UNKNOWN", "YES")
# cidt<-c("CIDT+", "CX+", "PARASITIC")

# Create an empty dataset
bact <- data.frame(year=factor(), 
                   state=character(),
                   pathogen=character(),
                   count=character(),
                   population=character(),
                   yearn=integer())

# Load FoodnNet MMWR Data - Update this annually
# Load FoodnNet MMWR Data - Update this annually
mmwrdata<-haven::read_sas("/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623.sas7bdat")%>%
  filter(SiteID!="COEX")

# mmwrdata<-subset(mmwrdata, mmwrdata$SERO1!="NOT SPECIATED" & mmwrdata$SERO1!="")
mmwrdata$SERO2<-ifelse(mmwrdata$SERO1=="NOT SPECIATED" | mmwrdata$SERO1=="UNKNOWN" | mmwrdata$SERO1=="PARTIAL SERO" | mmwrdata$SERO1=="NOT SERO" | mmwrdata$SERO1=="", 
                       "Missing", mmwrdata$SERO1)
mmwrdata$SERO2<-ifelse(grepl("UNDET",mmwrdata$SERO2), "Missing", mmwrdata$SERO2)

# Run the fnd function to Aggregate Data by State and Year
fnd<-function(travel, cidt, datafile){
  datafile$serotypesummary<-datafile$SERO2
  fn<-datafile %>%
    setNames(tolower(names(.)))%>%
    # Exclude by detection method and travel status, other exclusions can be added but would require thought
    filter((cxcidt %in% cidt) & (travelint %in% travel))%>% 
    # Fix county names
    filter(!county %in% c("OUT OF STATE", "UNKNOWN", "99997")) %>% 
    mutate(county = if_else(county %in% c("ST. MARYS'S", "ST. MARYS"), "ST. MARY'S", county)) %>%
    mutate(county = if_else(county %in% c("PRINCE GEORGES"), "PRINCE GEORGE'S", county)) %>%
    mutate(county = if_else(county %in% c("QUEEN ANNES"), "QUEEN ANNE'S", county)) %>%
    mutate(county = if_else(county %in% c("DE BACA"), "DEBACA", county))%>%
    # Create a new variable called pathogentype since catchment is slightly different for bacterial and parasitic pathogens in some years
    mutate(pathogentype = ifelse(pathogen %in% c("CRYPTOSPORIDIUM", "CYCLOSPORA"), 
                                 "Parasitic", "Bacterial"))
  
  # Read in Census Files - Update this annually
  # Read in file for bacterial pathogens
  haven::read_sas("/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat") %>%
    setNames(tolower(names(.))) %>%
    group_by(year, state) %>%
    dplyr::summarize(population = sum(population, na.rm=TRUE)) %>%
    mutate(pathogentype = "Bacterial") %>%
    bind_rows(
      # Read in file for parasitic pathogens
      haven::read_sas("/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat") %>%
        setNames(tolower(names(.))) %>%
        group_by(year, state) %>%
        dplyr::summarize(population = sum(population, na.rm=TRUE)) %>%
        mutate(pathogentype = "Parasitic")) %>% ungroup -> census
  
  # Create a list of pathogen names excluding Listeria
  pathogens<-c("CAMPYLOBACTER", "CYCLOSPORA", "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")
  
  # Create separate datasets for each "group of pathogens that have unique exclusions or census datasets
  gtools::smartbind(as.data.frame(fn%>% filter(pathogen %in%  pathogens & pathogen!="CYCLOSPORA")), # all pathogens but Cyclospora and Listeria
                    as.data.frame(fn%>% filter(pathogen == "STEC" & stec_class=="STEC O157")%>%mutate(pathogen="STEC O157")), # make a dataset for STEC O157
                    as.data.frame(fn%>% filter(pathogen == "STEC" & (stec_class=="STEC NONO157" | stec_class== "STEC O AG UNDET"))%>%mutate(pathogen="STEC NONO157")), # make a dataset for STEC NONO157
                    as.data.frame(fn%>% filter(pathogen == "LISTERIA" & cste=="YES")))%>% # make a Listeria dataset 
    mutate(yearn=as.numeric(as.character(year)),
           year=as.factor(year))%>%
    split(.$pathogen)%>% # split the dataset by pathogen
    # Create a function and map (aka loop) over each pathogen-specific dataset to summarize disease counts by year and state
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        # Make sure that state-year combinations with 0 cases for the given pathogen are represented in the dataset
        complete(year, state = unique(fn$state), fill = list(count = 0))%>%
        as.data.frame()%>%distinct() # keep onlyunique rows
    }) %>%
    # convert from a list of datasets to a single dataset 
    bind_rows(.id = "pathogen")%>%
    mutate(year=as.numeric(as.character(year)))%>%
    # Drop year-state combinations from the dataset for years before the given state entered the FoodNet catchment
    subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | (state=="TN" & year>=2000))%>%
    # Join the aggregated FoodNet data to the census data
    merge(census%>%filter(pathogentype=="Bacterial")%>%dplyr::select(year, state, population)%>%mutate(year=as.numeric(as.character(year))), by=c("year", "state"), all.x=TRUE) %>%
    distinct()%>%
    as.data.frame()->path
  
  #######################################################################################################################################
  # Adapt the above code for Cyclospora 
  para<-if(!("CIDT+" %in% cidt)){
    para<-NA
  } else {
    para<-fn%>% 
      filter(pathogen=="CYCLOSPORA")%>%
      mutate(yearn=as.numeric(as.character(year)))%>%
      split(.$pathogen)%>%
      map(function(x) {
        x<-x%>%
          group_by(year, state) %>%
          dplyr::summarize(count = n())%>%
          ungroup()%>%
          complete(year, state = unique(fn$state), fill = list(count = 0))%>% 
          as.data.frame()%>%distinct()}) %>%
      bind_rows(.id = "pathogen")%>%
      mutate(year=as.numeric(as.character(year)))%>%
      subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | (state=="NM" & year>=2004) | 
               (state=="NY" & year>=1998) | (state=="OR") | (state=="TN" & year>=2000))%>%
      left_join(census%>%filter(pathogentype=="Parasitic")%>%mutate(year=as.numeric(as.character(year))), by=c("year", "state")) %>%
      as.data.frame()
  }
  
  #######################################################################################################################################
  # Adapt the above code to for Salmonella serotypes 
  ## Calculate Incidence for Each Combination of Serotype, State and Year
  mostcommonsero<-function(data){
    sal.rank<-data %>% 
      filter(pathogen=="SHIGELLA")%>%
      filter(serotypesummary!="NOT SEROTYPED" & serotypesummary!="")%>%
      filter(year==max(as.numeric(as.character(data$year)), na.rm=TRUE))%>%
      group_by(serotypesummary)%>%
      dplyr::summarize(count = n())%>%
      as.data.frame()
    sal.rank<-sal.rank[order(sal.rank$count, decreasing = T),]%>%filter(serotypesummary!="NOT SEROTYPED" & serotypesummary!="")
    sal.rank<-c(sal.rank$serotypesummary[1:10])
    return(sal.rank)}
  
  sal<-if(!("CIDT+" %in% cidt)){
    sal<-NA
  } else {
    path %>%
      filter(pathogen == "CAMPYLOBACTER")  %>%
      ungroup -> pop
    pop$year<-as.factor(pop$year)
    fn$year<-as.factor(fn$year)
    
    sal.rank<-mostcommonsero(fn) 
    
    fn%>% 
      filter(pathogen=="SHIGELLA" & serotypesummary %in% sal.rank)%>%
      mutate(pathogen=serotypesummary, yearn=as.numeric(as.character(year)))%>%
      mutate(year=droplevels(year))%>%
      split(.$pathogen)%>%
      map(function(x) {
        x<-x%>%
          group_by(year, state) %>%
          dplyr::summarize(count = n())%>%
          ungroup()%>%
          complete(year, state = unique(fn$state), fill = list(count = 0))%>% 
          as.data.frame()%>%distinct()
      }) %>%
      bind_rows(.id = "pathogen")%>%
      mutate(year=as.numeric(as.character(year)))%>%
      subset((state=="CA") | (state=="CO" & year>=2001) | (state=="CT") | (state=="GA") | (state=="MD" & year>=1998) | (state=="MN") | (state=="NM" & year>=2004) | 
               (state=="NY" & year>=1998) | (state=="OR") | (state=="TN" & year>=2000))%>%
      left_join(pop%>%dplyr::select(year, state, population)%>%mutate(year=as.numeric(as.character(year))), by=c("year", "state")) %>%
      as.data.frame()->sal 
  }
  ####
  # Merge back into a single dataset so it works with the furrr parallelization code
  
  bact<-if("CIDT+" %in% cidt){
    gtools::smartbind(path, para)%>%gtools::smartbind(sal) 
  } else {
    path
  }
  
  return(bact)
}
#############

# Run model
trendy<-function(datafile){ 
  # data<- data %>%filter(pathogen==target)
  travel<-if(("YES" %in% travel) & ("UNKNOWN" %in% travel)){
    travel<-"Travel Included"
  } else if (!("YES" %in% travel) & ("UNKNOWN" %in% travel)){
    travel<-"Unknown Travel Included"
  } else {
    travel<-"Excluded"
  }
  culture<-ifelse("CIDT+" %in% cidt, "CxCIDT", "Cx")
  target<-unique(datafile$pathogen)
  brm(count~ s(yearn, by=state) + state + offset(log(population)), data=datafile, family = "negbinomial",
      save_all_pars=TRUE, chains=6, iter=5001, seed=47, 
      cores=6,
      control = list(adapt_delta = 0.999, max_treedepth = 19)) ->  proposed
  # saveRDS(proposed, paste("~/FN1_24_CxCIDT/Salmonella/SplinesModel/Model", target, travel, culture, run,".Rds", sep=""))
  # proposed<-readRDS("~/2020.FoodNetSplines/mmwr_output/CAMPYLOBACTERTravel IncludedCxCIDT.Rds")
  #############
  # bayes_R2(proposed)
  # model_performance(proposed)
  #############
  
  # Draw untransformed (link-level) predictions for a new (or the original) data using add_linpred (which is an alternate spelling of add_fitted_draws) and transform them
  ## This generates a distribution of estimates for each site
  add_linpred_draws_fn<-function(data, model){
    d.prop.fdraws<-add_linpred_draws(newdata = data, model, n=NULL, transform=FALSE, value=".linpred") 
    ## Get the actual and predicted incidence by year from proposed model using the pred draws
    d.prop.fdraws$pinc<-(d.prop.fdraws$.linpred)/(d.prop.fdraws$population/100000) # predicted incidence per polity per year, dropped exp()
    d.prop.fdraws$inc<-(d.prop.fdraws$count)/(d.prop.fdraws$population/100000) # incidence per polity per year
    d.prop.fdraws$code<-paste(d.prop.fdraws$state, d.prop.fdraws$yearn)
    return(d.prop.fdraws)
  }
  
  col_list = c("CA", "CO", "CT", "GA", "MD", "MN", "NM", "NY", "OR", "TN")
  
  # Convert these draws from site-level to catchment-level estimates
  catchment<-function(linpred_draws){
    value<-prop_linpred %>% reshape2::dcast(yearn+.draw~state, value.var = ".linpred")
    value[,col_list]<-exp(value[,col_list]) #########################
    value$.value<-rowSums(value[,col_list],na.rm=TRUE)                  
    count<-prop_linpred %>% reshape2::dcast(yearn+.draw~state, value.var = "count")
    count$count<-rowSums(count[,col_list], na.rm=TRUE)
    pop<-prop_linpred %>% reshape2::dcast(yearn+.draw~state, value.var = "population")
    pop$population<-rowSums(pop[,col_list], na.rm=TRUE)
    # head(pop%>%subset(yearn>2006))
    # head(value%>%subset(yearn>2006))
    # head(count%>%subset(yearn>2006))
    
    catch<-left_join(value[c(1:2, 13)], pop[c(1:2, 13)], by=c("yearn", ".draw"))%>%
      left_join(count[c(1:2, 13)], by=c("yearn", ".draw"))
    # ggplot(value,aes(x=.value,color="red"))+ geom_density()+ geom_vline(data=count, aes(xintercept=count), color="black")+  facet_wrap(~yearn) 
    # ggplot(catch,aes(x=.value,color="red"))+ geom_density()+ geom_vline(aes(xintercept=count), color="black")+  facet_wrap(~yearn)  
    # ggplot(catch)+geom_point(aes(x=yearn, y=.value), color="red")+geom_point(aes(x=yearn, y=count))+theme_minimal()
    # ggplot(catch)+geom_boxplot(aes(group=yearn, x=yearn, y=.value), color="black")+geom_point(aes(x=yearn, y=count), color="red")+theme_minimal()
    # ggplot(catch)+geom_violin(aes(group=yearn, x=yearn, y=.value), color="black")+geom_point(aes(x=yearn, y=count), color="red")+theme_minimal()
    catch$inc<-(catch$count/(catch$population/100000))
    catch$pinc<-(catch$.value/(catch$population/100000))
    # ggplot(catch)+geom_boxplot(aes(group=yearn, x=yearn, y=pinc), color="black")+geom_point(aes(x=yearn, y=inc), color="red")+theme_minimal()
    # ggplot(catch)+geom_violin(aes(group=yearn, x=yearn, y=pinc), color="black")+geom_point(aes(x=yearn, y=inc), color="red")+theme_minimal()
    catch$inc<-(catch$count/(catch$population/100000))
    catch$pinc<-(catch$.value/(catch$population/100000))
    return(catch)
  }
  
  # Convert these draws to catchment level estimates, including equal-tailed credibility interval
  from_linpred_catchIR<-function(catchments){
    prop.fdraws<-catchments %>% group_by(yearn) %>% 
      summarise_at(.vars = c("count", "population", ".value"), 
                   .funs = list(LL=~quantile(., probs = 0.025, na.rm=TRUE),
                                median=median, UL=~quantile(., probs = 0.975, na.rm=TRUE))) %>%
      select(yearn, .value_LL, .value_UL, .value_median, population_median, count_median)
    prop.fdraws$pinc<-(prop.fdraws$.value_median)/(prop.fdraws$population_median/100000) # predicted incidence per polity per year
    prop.fdraws$pinc.LL<-(prop.fdraws$.value_LL)/(prop.fdraws$population_median/100000) # predicted incidence per polity per year
    prop.fdraws$pinc.UL<-(prop.fdraws$.value_UL)/(prop.fdraws$population_median/100000) # predicted incidence per polity per year
    prop.fdraws$inc<-(prop.fdraws$count_median)/(prop.fdraws$population_median/100000)
    
    # Calculate highest density credibility interval
    splits<-split(catchments%>%mutate(.value=log(.value))%>%select(.value), catchments$yearn)
    t<-lapply(splits, function(x) HDInterval::hdi(x)%>%exp()%>%as.data.frame())
    t<-lapply(t, function(x) cbind(x, c("lower", "upper")))
    t<-dplyr::bind_rows(t, .id = "yearn")
    t<-reshape2::dcast(t, yearn~`c("lower", "upper")`, value.var = ".value")
    colnames(t)<-c("yearn", "hdi.LL", "hdi.UL")
    t$yearn<-as.numeric(as.character(t$yearn))
    
    prop.fdraws<-left_join(prop.fdraws, t , by="yearn")
    prop.fdraws$hdi.LL<-(prop.fdraws$hdi.LL)/(prop.fdraws$population_median/100000) # predicted incidence per polity per year
    prop.fdraws$hdi.UL<-(prop.fdraws$hdi.UL)/(prop.fdraws$population_median/100000)
    
    # ggplot(prop.fdraws)+geom_line(aes(x=yearn, y=pinc), color="black")+geom_ribbon(aes(ymin = pinc.LL, ymax = pinc.UL, x=yearn), fill = "blue", alpha = .25)+geom_point(aes(x=yearn, y=inc), color="red")+theme_minimal()
    # ggplot(catchir.linpred)+geom_line(aes(x=yearn, y=pinc), color="black")+geom_ribbon(aes(ymin = pinc.LL, ymax = pinc.UL, x=yearn), fill = "blue", alpha = .5)+geom_point(aes(x=yearn, y=inc), color="red")+theme_minimal()+geom_ribbon(aes(ymin = hdi.LL, ymax = hdi.UL, x=yearn), fill = "red", alpha = .25)
    return(prop.fdraws)
  }
  
  prop_linpred<-add_linpred_draws_fn(data=(proposed$data %>% group_by(state)), model=proposed) # draws from posterior
  
  catch<-catchment(prop_linpred) # catchment-level estimates
  catchir.linpred<-from_linpred_catchIR(catch) # catchment-level estimates
  catchir.linpred$pathogen<-target
  catchir.linpred$travel<-travel
  catchir.linpred$culture<-culture
  # save estimates
  # write.table(catchir.linpred%>%select(-c(.value_LL, .value_UL)), file = paste("~/FN1_24_CxCIDT/Salmonella/SplinesModel/EstRR_Catchment", run,".csv", sep=""), append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = TRUE, row.names = TRUE)
  # write.table(siteir.linpred, file = "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/mmwroutput/DistrofEstIR_SiteCIDTTrYUN04212022.csv", append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = TRUE, row.names = TRUE)
  #############
  # Calculate RR and Percent Change for Each Year Relative to a Baseline
  ## Write a Function to Do This for Any Given Combinations of Years
  IR.comp<-function(catch, catchir.linpred, year1, year2){
    ref<-catch %>% 
      filter(yearn<=year2 & yearn >=year1) %>% group_by(.draw)%>%
      summarise_at(.vars = c("count", "population", ".value"), 
                   .funs = list(mean=mean))
    colnames(ref)<-c(".draw", "ref_count", "ref_pop", "ref_value")
    
    comb<-left_join(catch%>%select(yearn, .draw, .value, population, count), ref, by=c(".draw"))
    comb$pinc<-(comb$.value)/(comb$population/100000) # predicted incidence per polity per year
    comb$inc<-(comb$count)/(comb$population/100000)
    comb$ref.ir<-(comb$ref_count)/(comb$ref_pop /100000) # predicted incidence per polity per year
    comb$ref.est.ir<-(comb$ref_value)/(comb$ref_pop /100000) # predicted incidence per polity per year
    
    comb$rr   <-round(comb$pinc/comb$ref.est.ir, 2)
    comb$pct.change   <-round(((comb$pinc-comb$ref.est.ir)/comb$ref.est.ir)*100,2)
    # ci_hdi <- bayestestR::ci(comb$, method = "HDI")
    
    comb2<-comb%>%
      group_by(yearn) %>% 
      summarise_at(.vars = c("pinc", "inc", "ref.ir", "ref.est.ir", "ref_pop", "rr", "pct.change", ".value", "ref_value", "population", "ref_pop", "count", "ref_count"), 
                   .funs = list(LL=~quantile(., probs = 0.025, na.rm=TRUE),
                                median=median, UL=~quantile(., probs = 0.975, na.rm=TRUE))) %>%
      select(yearn, .value_median, count_median,
             ref_value_median, ref_count_median, 
             population_median, 
             ref_pop_median, 
             pinc_median, pinc_LL, pinc_UL, inc_median, 
             ref.est.ir_median, ref.est.ir_LL, ref.est.ir_UL, ref.ir_median, 
             rr_median, rr_LL,	rr_UL, 
             pct.change_median,	pct.change_LL, pct.change_UL)%>%
      dplyr::rename(inc=inc_median)
    
    return(comb2)
  }
  
  
  # Calculate for 2016-2018
  threeyear<-IR.comp(catch, catchir.linpred, 2016, 2018)
  threeyear$pathogen<-target
  threeyear$travel<-travel
  threeyear$culture<-culture  
  threeyear$startyear<-2016
  threeyear$startyear<-2018
  write.table(threeyear, "~/FN1_24_CxCIDT/Shigella/SplinesModel/EstIRR_Catchment_ShigellaSpecies.csv", append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = TRUE, row.names = TRUE)
  
  # Calculate for The Most Recent 3 Years
  recent<-IR.comp(catch, catchir.linpred, 2020, 2022)
  recent$pathogen<-target
  recent$travel<-travel
  recent$culture<-culture
  recent$startyear<-2020
  recent$startyear<-2022
  write.table(recent, "~/FN1_24_CxCIDT/Shigella/SplinesModel/EstIRR_Catchment_ShigellaSpecies.csv", append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = FALSE, row.names = TRUE)
  
  # Calculate for the Earliest Three Years of the Finalized Catchment
  finalcatchment<-IR.comp(catch, catchir.linpred, 2004, 2006)
  finalcatchment$pathogen<-target
  finalcatchment$travel<-travel
  finalcatchment$culture<-culture
  finalcatchment$startyear<-2004
  finalcatchment$startyear<-2006
  write.table(finalcatchment, "~/FN1_24_CxCIDT/Shigella/SplinesModel/EstIRR_Catchment_ShigellaSpecies.csv", append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = FALSE, row.names = TRUE)
  
  # Calculate for all Years Versus the 2006-2008 Baseline
  baseline<-IR.comp(catch, catchir.linpred, 2006, 2008)
  baseline$pathogen<-target
  baseline$travel<-travel
  baseline$culture<-culture
  baseline$startyear<-2006
  baseline$startyear<-2008
  write.table(baseline, "~/FN1_24_CxCIDT/Shigella/SplinesModel/EstIRR_Catchment_ShigellaSpecies.csv", append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = FALSE, row.names = TRUE)

  # Calculate for all Years Versus the 2006-2008 Baseline
  baseline<-IR.comp(catch, catchir.linpred, 2010, 2012)
  baseline$pathogen<-target
  baseline$travel<-travel
  baseline$culture<-culture
  baseline$startyear<-2010
  baseline$startyear<-2012
  write.table(baseline, "~/FN1_24_CxCIDT/Shigella/SplinesModel/EstIRR_Catchment_ShigellaSpecies.csv", append = TRUE, quote = TRUE, sep =",", qmethod = "double", col.names = FALSE, row.names = TRUE)
}
#############
# Pathogens with Travel and CIDT
#############
start.time <- Sys.time()
travel<-c("NO", "UNKNOWN", "YES")
cidt<-c("CIDT+", "CX+", "PARASITIC")

bact<-fnd(travel, cidt, mmwrdata)
bact<-subset(bact, pathogen=="Missing" | pathogen=="FLEXNERI"| pathogen=="SONNEI") # | pathogen=="BOYDII" | pathogen=="DYSENTERIAE")
bact$yearn<-as.numeric(as.character(bact$year))
bact$year<-as.factor(bact$year)
bact2<-split(bact, bact$pathogen)

plan(multisession, workers=10)
furrr::future_map(bact2,
                  function(x) trendy(
                    # FUN = trendy,
                    data = x),
                  .progress = FALSE,
                  .options = furrr_options(seed = TRUE))
end.time <- Sys.time()
end.time - start.time

#############
# Pathogens with Travel !="Yes" and CIDT
#############
start.time <- Sys.time()
travel<-c("NO", "UNKNOWN")
cidt<-c("CIDT+", "CX+", "PARASITIC")

bact<-fnd(travel, cidt, mmwrdata)
bact<-subset(bact, pathogen=="Missing" | pathogen=="FLEXNERI"| pathogen=="SONNEI") # | pathogen=="BOYDII" | pathogen=="DYSENTERIAE")
bact$yearn<-as.numeric(as.character(bact$year))
bact$year<-as.factor(bact$year)
bact2<-split(bact, bact$pathogen)

plan(multisession, workers=10)
furrr::future_map(bact2,
                  function(x) trendy(
                    # FUN = trendy,
                    data = x),
                  .progress = FALSE,
                  .options = furrr_options(seed = TRUE))
end.time <- Sys.time()
end.time - start.time


Rprof(NULL)
sI <- sessionInfo()
print(sI, RNG = TRUE, locale = FALSE)
toLatex(sI, file = "/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/mmwroutput/packages.bib")
##################################################################
# LIBRARY
##################################################################
library(tidyverse)
LOAD_PACKAGES <- function(pkgs){
  package_vctr <- pkgs %>% stringr::str_split('\n') %>% unlist()
  invisible(lapply(package_vctr, library, character.only = TRUE))
}

##################################################################
# REFORMATTING
##################################################################
# reformat lists
CLEAN_LIST<-function(input_string){
  ## Remove the brackets and quotes
  cleanedString <- gsub('[\\[\\]\"]', '', input_string)
  
  # Split the string by commas
  cleanedList <- strsplit(cleanedString, ",")[[1]]
  
  return(cleanedList)
}

##################################################################
# PATH Analysis
##################################################################
PATH_ANALYSIS<-function(mmwrdata,census){
  # Create a list of pathogen names excluding Listeria
  pathogens<-c("CAMPYLOBACTER", "CYCLOSPORA", 
               "SALMONELLA", "SHIGELLA", "STEC", "VIBRIO", "YERSINIA")
  
  # Create separate datasets for each "group of pathogens that have unique exclusions or census datasets
  ## all pathogens but Cyclospora and Listeria
  df1=mmwrdata %>% filter(pathogen %in% pathogens
                          & pathogen!="CYCLOSPORA")
  ## make a dataset for STEC O157
  df2=mmwrdata%>% filter(pathogen == "STEC" & 
                           stec_class=="STEC O157")%>%
    mutate(pathogen="STEC O157")
  
  ## make a dataset for STEC NONO157
  df3=mmwrdata%>% filter(stec_class=="STEC NONO157" | 
                           stec_class== "STEC O AG UNDET")%>%
    mutate(pathogen="STEC NONO157")
  
  # make a Listeria dataset 
  df4=mmwrdata%>% filter(pathogen == "LISTERIA" & cste=="YES")
  
  selectDf=gtools::smartbind(df1,df2,df3,df4)
  remove(df1,df2,df3,df4)
  
  # revise year
  selectDf = selectDf %>%
    mutate(yearn=as.numeric(as.character(year)),
           year=as.factor(year))
  
  # split the dataset by pathogen
  # Create a function and map (aka loop) over each 
  # pathogen-specific dataset to summarize disease counts 
  # by year and state
  selectDf=selectDf %>%
    # split the dataset by pathogen
    split(.$pathogen)%>% 
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        # Make sure that state-year combinations with 0 cases 
        # for the given pathogen are represented in the dataset
        complete(year, state = unique(mmwrdata$state),
                 fill = list(count = 0))%>%
        as.data.frame()%>%distinct() # keep onlyunique rows
    }) %>%
    # convert from a list of datasets to a single dataset 
    bind_rows(.id = "pathogen")
  
  # Drop year-state combinations from the dataset for years 
  # before the given state entered the FoodNet catchment
  selectDf = selectDf %>%
    mutate(year=as.numeric(as.character(year))) %>%
    subset((state=="CA") | (state=="CO" & year>=2001) | 
             (state=="CT") | (state=="GA") | 
             (state=="MD" & year>=1998) | 
             (state=="MN") | (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | 
             (state=="TN" & year>=2000))
  
  # Join the aggregated FoodNet data to the census data
  mergedDf = selectDf%>% merge(census %>%
                                 filter(pathogentype=="Bacterial") %>%
                                 dplyr::select(year, state, 
                                               population) %>%
                                 mutate(year=as.numeric(
                                   as.character(year))), 
                               by=c("year", "state"), 
                               all.x=TRUE) %>%
    distinct() %>%
    as.data.frame()
  
  saveFile=paste0(outBase,"mmwrdata.csv")
  write.csv(mergedDf,saveFile)
  
  return(mergedDf)
}

##################################################################
# CYCLOSPORA Analysis
##################################################################
CYCLOSPORA_ANALYSIS<-function(mmwrdata,census){
  para<-mmwrdata %>% 
    filter(pathogen=="CYCLOSPORA")%>%
    mutate(yearn=as.numeric(as.character(year))) %>%
    split(.$pathogen)%>%
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        complete(year, state = unique(mmwrdata$state), 
                 fill = list(count = 0))%>% 
        as.data.frame()%>%distinct()}) %>%
    bind_rows(.id = "pathogen") %>%
    mutate(year=as.numeric(as.character(year))) %>%
    subset((state=="CA") | (state=="CO" & year>=2001) | 
             (state=="CT") | (state=="GA") | 
             (state=="MD" & year>=1998) | (state=="MN") | 
             (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | 
             (state=="TN" & year>=2000)) %>%
    left_join(census%>%filter(pathogentype=="Parasitic") %>%
                mutate(year=as.numeric(as.character(year))), 
              by=c("year", "state")) %>%
    as.data.frame()
  return(para)
}

##################################################################
# SALMONELLA Analysis
##################################################################
SALMONELLA_ANALYSIS<-function(pathDf,mmwrdata,census){
  pop = pathDf %>%
    filter(pathogen == "CAMPYLOBACTER")  %>%
    ungroup
  
  # adjust year to factor
  pop$year<-as.factor(pop$year)
  mmwrdata$year<-as.factor(mmwrdata$year)
  
  # filter most common
  sal.rank = mmwrdata %>% 
    filter(pathogen=="SHIGELLA")%>%
    filter(serotypesummary!="NOT SEROTYPED" & serotypesummary!="") %>%
    filter(year==max(as.numeric(as.character(mmwrdata$year)), 
                     na.rm=TRUE))%>%
    group_by(serotypesummary)%>%
    dplyr::summarize(count = n())%>%
    as.data.frame()
  sal.rank<-sal.rank[order(sal.rank$count, decreasing = T),] %>%
    filter(serotypesummary!="NOT SEROTYPED" & serotypesummary!="")
  sal.rank<-c(sal.rank$serotypesummary[1:10])
  
  
  sal = mmwrdata%>% 
    filter(pathogen=="SHIGELLA" & serotypesummary %in% sal.rank) %>%
    mutate(pathogen=serotypesummary, 
           yearn=as.numeric(as.character(year)))%>%
    mutate(year=droplevels(year))%>%
    split(.$pathogen)%>%
    map(function(x) {
      x<-x%>%
        group_by(year, state) %>%
        dplyr::summarize(count = n())%>%
        ungroup()%>%
        complete(year, state = unique(mmwrdata$state), 
                 fill = list(count = 0))%>% 
        as.data.frame()%>%distinct()
    }) %>%
    bind_rows(.id = "pathogen") %>%
    mutate(year=as.numeric(as.character(year)))%>%
    subset((state=="CA") | (state=="CO" & year>=2001) | 
             (state=="CT") | (state=="GA") | 
             (state=="MD" & year>=1998) | (state=="MN") | 
             (state=="NM" & year>=2004) | 
             (state=="NY" & year>=1998) | (state=="OR") | 
             (state=="TN" & year>=2000))%>%
    left_join(pop%>%dplyr::select(year, state, population) %>%
                mutate(year=as.numeric(as.character(year))), 
              by=c("year", "state")) %>%
    as.data.frame()
  
  return(sal)
}

##################################################################
# BRM Modeling
##################################################################
PROPOSED_BM<-function(bact){
  proposed <- brm(
    count ~ s(yearn, by = state) + state + offset(log(population)),
    data = bact,
    family = "negbinomial",
    save_pars = save_pars(all = TRUE),
    chains = 6,
    iter = 5001,
    seed = 47,
    cores = 6,
    control = list(adapt_delta = 0.999, max_treedepth = 19)
  )
  return(proposed)
}

##################################################################
# LINPRED
##################################################################
LINPREAD_DRAW_FN<-function(data, model){
  d.prop.fdraws<-add_linpred_draws(newdata = data,
                                   model, n=NULL,
                                   transform=FALSE,
                                   value=".linpred")
  ## Get the actual and predicted incidence by year
  ## from proposed model using the pred draws
  d.prop.fdraws$pinc<-(d.prop.fdraws$.linpred)/
    (d.prop.fdraws$population/100000)
  
  # predicted incidence per polity per year, dropped exp()
  d.prop.fdraws$inc<-(d.prop.fdraws$count)/
    (d.prop.fdraws$population/100000)
  
  # incidence per polity per year
  d.prop.fdraws$code<-paste(d.prop.fdraws$state,
                            d.prop.fdraws$yearn)
  return(d.prop.fdraws)
}

##################################################################
# Catchment
##################################################################
SAVE_PLOT<-function(p,label){
  print(p)
  saveFile=paste0(outBase,label,".png")
  ggsave(saveFile,p)
}

CATCHMENT<-function(prop_linpred){
  value<-prop_linpred %>%
    reshape2::dcast(yearn+.draw~state, value.var = ".linpred")
  value[,col_list]<-exp(value[,col_list])
  value$.value<-rowSums(value[,col_list],na.rm=TRUE)
  
  count<-prop_linpred %>%
    reshape2::dcast(yearn+.draw~state, value.var = "count")
  count$count<-rowSums(count[,col_list], na.rm=TRUE)
  
  pop<-prop_linpred %>%
    reshape2::dcast(yearn+.draw~state, value.var = "population")
  pop$population<-rowSums(pop[,col_list], na.rm=TRUE)
  
  catch=left_join(value[c(1:2, 13)],
                   pop[c(1:2, 13)], 
                   by=c("yearn", ".draw")) %>%
    left_join(count[c(1:2, 13)], 
              by=c("yearn", ".draw"))
  
  p=ggplot(value,aes(x=.value,color="red"))+
    geom_density()+
    geom_vline(data=count, aes(xintercept=count),
               color="black")+  facet_wrap(~yearn)
  SAVE_PLOT(p,"catch_value")
  
  p=ggplot(catch,aes(x=.value,color="red"))+
    geom_density()+
    geom_vline(data=count, aes(xintercept=count),
               color="black")+  facet_wrap(~yearn)
  SAVE_PLOT(p,"catch")

  p=ggplot(catch)+geom_boxplot(aes(group=yearn,
                                 x=yearn, y=.value),
                             color="black")+
    geom_point(aes(x=yearn,y=count),
               color="red")+
    theme_minimal()
  SAVE_PLOT(p,"catch_boxplot")
  
  p=ggplot(catch)+
    geom_violin(aes(group=yearn, x=yearn,
                    y=.value), color="black")+
    geom_point(aes(x=yearn, y=count),
               color="red")+theme_minimal()
  SAVE_PLOT(p,"catch_vplot")
  
  catch$inc<-(catch$count/(catch$population/100000))
  catch$pinc<-(catch$.value/(catch$population/100000))
  
  p=ggplot(catch)+geom_boxplot(aes(group=yearn,
                                 x=yearn, y=pinc), color="black")+
    geom_point(aes(x=yearn, y=inc), color="red")+
    theme_minimal()
  SAVE_PLOT(p,"catch2_boxplot")
  
  p=ggplot(catch)+geom_violin(aes(group=yearn,
                                x=yearn, y=pinc),
                            color="black")+
    geom_point(aes(x=yearn, y=inc),
               color="red")+theme_minimal()
  SAVE_PLOT(p,"catch2_vplot")
  
  catch$inc<-(catch$count/(catch$population/100000))
  catch$pinc<-(catch$.value/(catch$population/100000))
  return(catch)
}

##################################################################
# Catchment IR
##################################################################
LINPRED_TO_CATCHIR<-function(catchments){
  prop.fdraws<-catchments %>%
    group_by(yearn) %>%
    summarise_at(.vars = c("count", "population", ".value"),
                 .funs = list(LL=~quantile(.,
                                           probs = 0.025,
                                           na.rm=TRUE),
                              median=median,
                              UL=~quantile(.,
                                           probs = 0.975,
                                           na.rm=TRUE))) %>%
    select(yearn, .value_LL, .value_UL,
           .value_median, population_median, count_median)
  
  # predicted incidence per polity per year
  prop.fdraws$pinc<-(prop.fdraws$.value_median)/
    (prop.fdraws$population_median/100000)
  
  # predicted incidence per polity per year
  prop.fdraws$pinc.LL<-(prop.fdraws$.value_LL)/
    (prop.fdraws$population_median/100000)
  
  # predicted incidence per polity per year
  prop.fdraws$pinc.UL<-(prop.fdraws$.value_UL)/
    (prop.fdraws$population_median/100000)
  
  prop.fdraws$inc<-(prop.fdraws$count_median)/
    (prop.fdraws$population_median/100000)
  
  # Calculate highest density credibility interval
  splits<-split(catchments %>%
                  mutate(.value=log(.value))%>%
                  select(.value), catchments$yearn)
  
  # splits
  t<-lapply(splits, function(x) HDInterval::hdi(x)%>%
              exp()%>%as.data.frame())
  t<-lapply(t, function(x) cbind(x, c("lower", "upper")))
  t<-dplyr::bind_rows(t, .id = "yearn")
  t<-reshape2::dcast(t, yearn~`c("lower", "upper")`,
                     value.var = ".value")
  colnames(t)<-c("yearn", "hdi.LL", "hdi.UL")
  t$yearn<-as.numeric(as.character(t$yearn))
  prop.fdraws<-left_join(prop.fdraws, t , by="yearn")
  
  # predicted incidence per polity per year
  prop.fdraws$hdi.LL<-(prop.fdraws$hdi.LL)/
    (prop.fdraws$population_median/100000)
  
  prop.fdraws$hdi.UL<-(prop.fdraws$hdi.UL)/
    (prop.fdraws$population_median/100000)
  
  return(prop.fdraws)
}

##################################################################
# IR COMP
##################################################################
IR_COMP<-function(catch, year1, year2,fileBase){
  ref<-catch %>%
    filter(yearn<=year2 & yearn >=year1) %>% group_by(.draw)%>%
    summarise_at(.vars = c("count", "population", ".value"),
                 .funs = list(mean=mean))
  colnames(ref)<-c(".draw", "ref_count", "ref_pop", "ref_value")
  
  comb<-left_join(catch%>%select(yearn, .draw, .value,
                                 population, count), ref,
                  by=c(".draw"))
  # predicted incidence per polity per year
  comb$pinc<-(comb$.value)/(comb$population/100000)
  comb$inc<-(comb$count)/(comb$population/100000)
  
  # predicted incidence per polity per year
  comb$ref.ir<-(comb$ref_count)/(comb$ref_pop /100000)
  
  # predicted incidence per polity per year
  comb$ref.est.ir<-(comb$ref_value)/(comb$ref_pop /100000)
  
  comb$rr   <-round(comb$pinc/comb$ref.est.ir, 2)
  comb$pct.change <-round(((comb$pinc-comb$ref.est.ir)/
                             comb$ref.est.ir)*100,2)
  
  comb2<-comb%>%
    group_by(yearn) %>%
    summarise_at(.vars = c("pinc", "inc", "ref.ir",
                           "ref.est.ir", "ref_pop",
                           "rr", "pct.change", ".value",
                           "ref_value", "population",
                           "ref_pop", "count", "ref_count"),
                 .funs = list(LL=~quantile(., probs = 0.025,
                                           na.rm=TRUE),
                              median=median,
                              UL=~quantile(., probs = 0.975,
                                           na.rm=TRUE))) %>%
    select(yearn, .value_median, count_median,
           ref_value_median, ref_count_median,
           population_median,
           ref_pop_median,
           pinc_median, pinc_LL, pinc_UL, inc_median,
           ref.est.ir_median, ref.est.ir_LL,
           ref.est.ir_UL, ref.ir_median,
           rr_median, rr_LL,	rr_UL,
           pct.change_median,	pct.change_LL, pct.change_UL)%>%
    dplyr::rename(inc=inc_median)
  
  comb2$pathogen<-target
  comb2$travel<-travel
  comb2$culture<-culture
  comb2$startyear<-year2
  
  saveFile=paste0(outBase,fileBase)
  write.table(comb2, saveFile, append = TRUE,
              quote = TRUE, sep =",",
              qmethod = "double", col.names = TRUE,
              row.names = TRUE)
}

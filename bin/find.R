#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
source("functions.R")

##############################################################
# Setup
##############################################################
# create parser object
parser <- ArgumentParser()
parser$add_argument("--debug", type="logical", 
                    default=TRUE)
parser$add_argument("--mmwrFile", type="character", 
                    help="foodnet mmwr SASS data")
parser$add_argument("--censusFile_B", type="character", 
                    help="census file for bacterial data")
parser$add_argument("--censusFile_P", type="character", 
                    help="census file for parasitic data")
parser$add_argument("--travel", type="character", 
                    help="list for travel types IE NO,UNKNOWN")
parser$add_argument("--cidt", type="character", 
                    help="cidt IE CIDT+,CX+,PARASITIC")
parser$add_argument("--projID", type="character", 
                    help="variable to name the project")
args <- parser$parse_args()


# set up args
if (args$debug == FALSE){
  # read in args
  mmwrFile=args$mmwrFile
  censusFile_B=args$censusFile_B
  censusFile_P=args$censusFile_P
  projID=args$projID
  
  # reformat list
  travel=CLEAN_LIST(args$travel)
  cidt=CLEAN_LIST(args$cidt)
  
  # set outputDir
  outDir="SplineResults"
  
} else{
  # file args
  mmwrFile="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/mmwr9623_Jan2024.sas7bdat"
  censusFile_B="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623.sas7bdat"
  censusFile_P="/scicomp/groups-pure/OID/NCEZID/DFWED/EDEB/foodnet/trends/data/cen9623_para.sas7bdat"
  projID="20240705"
  
  #reformat list
  travel=c("NO,UNKNOWN,YES")
  cidt=c("CIDT+,CX+,PARASITIC")
  travel=CLEAN_LIST(travel)
  cidt=CLEAN_LIST(cidt)
  
  # set outdir
  outDir="~/projects/ticket/Rspline_test"
}

# load packages
pkgs <- c('haven','gtools','brms')
LOAD_PACKAGES(pkgs)

# Create outputDir
dir.create(outDir)

##############################################################
# Main Code
##############################################################

# Set labels
##############################################################
if (("YES" %in% travel) || ("UNKNOWN" %in% travel)) {
  travelLabel <- "Travel Included"
} else if (!("YES" %in% travel) & ("UNKNOWN" %in% travel)) {
  travelLabel <- "Unknown Travel Included"
} else {
  travelLabel <- "Excluded"
}
culture<-ifelse("CIDT+" %in% cidt, "CxCIDT", "Cx")

# Set outfile base names
##############################################################
outBase=paste0(outDir,"/",
               projID,"_","splinesmodel_",
               paste(travelLabel,collapse=""),"_",
               paste(culture,collapse=""),"_")

# Input data
##############################################################
print("--IMPORTING")
mmwrdata<-haven::read_sas(mmwrFile)%>%
  filter(SiteID!="COEX")

# Convert to df
mmwrdata=as.data.frame(mmwrdata)

# Update SERO list
seroList=c("NOT SPECIATED","UNKNOWN","PARTIAL SERO",
           "NOT SERO","")
mmwrdata$SERO2<-ifelse(mmwrdata$SERO1 %in% seroList, 
                       "Missing", mmwrdata$SERO1)
mmwrdata$SERO2<-ifelse(grepl("UNDET",mmwrdata$SERO2), 
                       "Missing", mmwrdata$SERO2)

# relabel  SERO2 column
mmwrdata$serotypesummary<-mmwrdata$SERO2

# Update count labels
mmwrdata<-mmwrdata %>%
  setNames(tolower(names(.)))%>%
  # Exclude by detection method and travel status, other exclusions can be added but would require thought
  filter((cxcidt %in% cidt) & (travelint %in% travel))%>% 
  # Fix county names
  filter(!county %in% c("OUT OF STATE", 
                        "UNKNOWN", "99997")) %>% 
  mutate(county = if_else(county %in% c("ST. MARYS'S", "ST. MARYS"), "ST. MARY'S",
                          county)) %>%
  mutate(county = if_else(county %in% c("PRINCE GEORGES"), "PRINCE GEORGE'S", 
                          county)) %>%
  mutate(county = if_else(county %in% c("QUEEN ANNES"), "QUEEN ANNE'S", 
                          county)) %>%
  mutate(county = if_else(county %in% c("DE BACA"), "DEBACA", 
                          county))%>%
  # Create a new variable called pathogentype since catchment is slightly different for bacterial and parasitic pathogens in some years
  mutate(pathogentype = ifelse(pathogen %in% c("CRYPTOSPORIDIUM", 
                                               "CYCLOSPORA"), 
                               "Parasitic", "Bacterial"))

# Census data
##############################################################
## Read in file for bacterial pathogens
## Read in file for parasitic pathogens
print("--IMPORTING CENSUS")
census=haven::read_sas(censusFile_B) %>%
  setNames(tolower(names(.))) %>%
  group_by(year, state) %>%
  dplyr::summarize(population = sum(population, na.rm=TRUE)) %>%
  mutate(pathogentype = "Bacterial") %>%
  bind_rows(
    haven::read_sas(censusFile_P) %>%
      setNames(tolower(names(.))) %>%
      group_by(year, state) %>%
      dplyr::summarize(population = sum(population, na.rm=TRUE)) %>%
      mutate(pathogentype = "Parasitic")) %>% ungroup
census=as.data.frame(census)

# Pathogens
##############################################################
print("--RUNNING PATH")
pathDf=PATH_ANALYSIS(mmwrdata,census)

# Cyclospora, Salmonella
##############################################################
if("CIDT+" %in% cidt){
  print("--RUNNING CYCLO")
  cyloDF=CYCLOSPORA_ANALYSIS(mmwrdata,census)
  
  print("--RUNNING SALMONELLA")
  salDF=SALMONELLA_ANALYSIS(pathDf,mmwrdata,census)
  
  bact=gtools::smartbind(pathDf, cyloDF)%>%
    gtools::smartbind(salDF) 
} else{
  bact=pathDf
}

# Post Process
##############################################################
remove(mmwrdata)
print("--POST PROCESSING")
bact<-subset(bact, pathogen=="Missing" | 
               pathogen=="FLEXNERI"| pathogen=="SONNEI")
bact$yearn<-as.numeric(as.character(bact$year))
bact$year<-as.factor(bact$year)
bact<-split(bact, bact$pathogen)
target<-unique(bact$pathogen)
head(bact)

# Run model
##############################################################
print("--RUN MODEL")
brm=PROPOSED_BM(bact)
saveFile=paste0(outBase,"brm.Rds")
saveRDS(brm,saveFile)

# Draw untransformed (link-level) predictionsusing 
# add_linpred (aka add_fitted_draws) and transform them
# generates a distribution of estimates for each site
##################################################
print("--RUN POST LINPRED")
posteriorLinpred<-LINPREAD_DRAW_FN(data=(proposed$data %>%
                                           group_by(state)),
                                       model=proposed)

# Catchment
##############################################################
# catchment-level estimates
col_list = c("CA", "CO", "CT", "GA", "MD",
             "MN", "NM", "NY", "OR", "TN")

# Convert these draws from site-level to catchment-level estimates
print("--RUN CATCHMENT")
catch<-CATCHMENT(posteriorLinpred)


# Credibility intervals
##############################################################
# Convert these draws to catchment level estimates, including equal-tailed credibility interval
print("--RUN LINPRED_TO_CATCHIR")
catchir.linpred<-LINPRED_TO_CATCHIR(catch)

# redefine meta
catchir.linpred$pathogen<-target
catchir.linpred$travel<-travelLabel
catchir.linpred$culture<-culture

# save estimates
saveFile=paste0(outBase,"IRCatch.csv")
write.table(catchir.linpred%>%select(-c(.value_LL, .value_UL)),
          saveFile, append = TRUE, quote = TRUE, sep =",", 
          qmethod = "double", col.names = TRUE, row.names = TRUE)

# Calculate RR and Percent Change for Each Year Relative to a Baseline
##############################################################
print("--RUN EST IRR")

# Calculate for 2016-2018
IR_COMP(saveFile,2016,2018,"EstIRRCatch_2016_2018.csv")

# Calculate for The Most Recent 3 Years
IR_COMP(saveFile,2020,2022,"EstIRRCatch_2020_2022.csv")

# Calculate for the Earliest Three Years of the Finalized Catchment
IR_COMP(saveFile,2004,2006,"EstIRRCatch_2004_2006.csv")

# Calculate for all Years Versus the 2006-2008 Baseline
IR_COMP(saveFile,2006,2008,"EstIRRCatch_2006_2008.csv")

# Calculate for all Years Versus the 2006-2008 Baseline
IR_COMP(saveFile,2010,2012,"EstIRRCatch_2010_2012.csv")
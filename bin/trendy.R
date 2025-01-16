#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
source("functions.R")
options(warn = -1)

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
parser$add_argument("--outDir", type="character", 
                    default="./",
                    help="output directory for results (default: ./)")
args <- parser$parse_args()

# set up args
if (args$debug == FALSE){
  # read in args
  mmwrFile=args$mmwrFile
  censusFile_B=args$censusFile_B
  censusFile_P=args$censusFile_P
  projID=args$projID
  outDir
  
  # reformat list
  travel=CLEAN_LIST(args$travel)
  cidt=CLEAN_LIST(args$cidt)
  
  # set outputDir
  outDir=paste0(projID, "/","SplineResults")
  
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
  outDir=paste0(projID, "/","SplineResults")
}

dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# load packages
pkgs <- c('haven','gtools','brms','ggplot2','tidybayes','HDInterval')
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
               gsub(" ","",travelLabel),"_",
               paste(culture,collapse=""),"_")

# print variables
print("--ANALYSIS DETAILS")
print(paste0("mmwrFile: ", mmwrFile," | ",
             "censusFile_B: ",censusFile_B," | ",
             "censusFile_P: ",censusFile_P," | ",
             "projID: ",projID," | ",
             "travel: ",paste(travel,collapse = ",")," | ",
             "cidt: ",paste(cidt,collapse=",")))
# Input data
##############################################################
print("--IMPORTING")
mmwrdata<-haven::read_sas(mmwrFile)%>%
  filter(SiteID!="COEX")

# Convert to df
mmwrdata=as.data.frame(mmwrdata)

# Update SERO list: This was only there to create a “smaller” test case for building out the pipeline. We can remove it for now. Alternatively, we could add it as an option to the function if we want these values of SERO1 to be included or not. However, for the larger analysis, serotypesummary is an actual column we want to keep; we would just want to add SERO2 to that if serotypesummary was missing
# seroList=c("NOT SPECIATED","UNKNOWN","PARTIAL SERO",
#           "NOT SERO","")
# mmwrdata$SERO2<-ifelse(mmwrdata$SERO1 %in% seroList, 
#                       "Missing", mmwrdata$SERO1)
# mmwrdata$SERO2<-ifelse(grepl("UNDET",mmwrdata$SERO2), 
#                       "Missing", mmwrdata$SERO2)
# relabel  SERO2 column
# mmwrdata$serotypesummary<-mmwrdata$SERO2

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
# b.	Will the function for aggregating data work if we do not have every pathogen present in the dataset? Trying to figure out if we can apply to other datasets besides FoodNet
print("--CALCULATING INCIDENCE FOR BACTERIAL PATHOGENS")
pathDf=PATH_ANALYSIS(mmwrdata,census)

# Cyclospora, Salmonella
##############################################################
if("CIDT+" %in% cidt){
  print("--CALCULATING INCIDENCE FOR CYCLOSPORA")
  cyloDF=CYCLOSPORA_ANALYSIS(mmwrdata,census)
  
  print("--CALCULATING INCIDENCE FOR SALMONELLA SEROTYPES")
  salDF=SALMONELLA_ANALYSIS(pathDf,mmwrdata,census)
  
  bact=gtools::smartbind(pathDf, cyloDF)%>%
    gtools::smartbind(salDF) 
} else{
  bact=pathDf
}

# Post Process
##############################################################
#remove(mmwrdata)
print("--POST PROCESSING OF INCIDENCE DATA")
bact<-subset(bact, pathogen==VALUE) #, pathogen=="Missing" | pathogen=="FLEXNERI"| pathogen=="SONNEI")
    # For running as a singular model for one pathogen, need to add a way for users to choose which pathogen to run for. For the parallel version this isn't needed
target<-unique(bact$pathogen)
bact$yearn<-as.numeric(as.character(bact$year))
bact$year<-as.factor(bact$year)
bact<-split(bact, bact$pathogen)


# Run modelproposed
##############################################################
# TODO: source code line #224 filters pathogen, but this is blocked out
# after the split above (#182) the obj is split into three
# FLEXNERI, Missing, SONNEI
# unclear which is to be analyzed? for demo only analyzed FLEXNERI
print("--RUN MODEL")
# BRM modeling
results <- list() # Initialize an empty list to store results

for (pathogen_name in names(bact)) {
  print(paste0("Processing pathogen: ", pathogen_name))
  
  # Run the model
  proposed <- PROPOSED_BM(bact[[pathogen_name]])
  
  # Save the result to a file
  saveRDS(proposed, paste0(outBase, pathogen_name, "_brm.Rds"))
  
  # Store the result in the list with the pathogen name as the key
  results[[pathogen_name]] <- proposed
}

# Check the structure of the results
str(results)
# Generate posterior predictions
posteriorLinpred <- LINPREAD_DRAW_FN(
  data = proposed$data %>% 
    group_by(state, yearn, population, count) %>% 
    summarise(across(everything(), ~ unique(.)), .groups = "drop"),
  model = proposed
)

# Catchment-level estimates
catch <- CATCHMENT(posteriorLinpred)

# Convert to catchment-level IR
catchir.linpred <- LINPRED_TO_CATCHIR(catch)

# Save intermediate results
saveFile <- paste0(outBase, "IRCatch.csv")
SAFE_WRITE(catchir.linpred, saveFile)

# Run IR_COMP for different year ranges
print("--RUN EST IRR")
IR_COMP(catchir.linpred, 2016, 2018, paste0(pathogen_name, "_EstIRRCatch_2016_2018.csv"), target = pathogen_name, travel = travelLabel, culture = culture)
IR_COMP(catchir.linpred, 2020, 2022, paste0(pathogen_name, "_EstIRRCatch_2020_2022.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
IR_COMP(catchir.linpred, 2004, 2006, paste0(pathogen_name, "_EstIRRCatch_2004_2006.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
IR_COMP(catchir.linpred, 2006, 2008, paste0(pathogen_name, "_EstIRRCatch_2006_2008.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
IR_COMP(catchir.linpred, 2010, 2012, paste0(pathogen_name, "_EstIRRCatch_2010_2012.csv"),target = pathogen_name, travel = travelLabel, culture = culture)

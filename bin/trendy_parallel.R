#b.	How would one specify that they want to use trendy_parallel.R in the pipeline instead of trendy.R? Would this be something we specify in the nextflow pipeline, or do we need a pipeline for running in parallel and a pipeline for running singularly? 
#!/usr/bin/env Rscript
# Load required libraries and suppress startup messages
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("doParallel"))
source("functions.R")
options(warn = -1)

##############################################################
# Setup and Argument Parsing
##############################################################

# Create an ArgumentParser object to handle command-line arguments
parser <- ArgumentParser()
parser$add_argument("--debug", type = "logical", default = TRUE, help = "Enable debug mode")
parser$add_argument("--mmwrFile", type = "character", help = "Path to FoodNet MMWR SASS data")
parser$add_argument("--censusFile_B", type = "character", help = "Path to census data for bacterial pathogens")
parser$add_argument("--censusFile_P", type = "character", help = "Path to census data for parasitic pathogens")
parser$add_argument("--travel", type = "character", help = "Travel types (e.g., NO, UNKNOWN)")
parser$add_argument("--cidt", type = "character", help = "CIDT types (e.g., CIDT+, CX+, PARASITIC)")
parser$add_argument("--projID", type = "character", help = "Project identifier")
parser$add_argument("--outDir", type = "character", default = "./", help = "Output directory for results")

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
               projID,"_",
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
##############################################################
# Pathogen Analysis
##############################################################

# Perform PATH_ANALYSIS and include Cyclospora and Salmonella
print("--CALCULATING INCIDENCE FOR BACTERIAL PATHOGENS")
pathDf <- PATH_ANALYSIS(mmwrdata, census)

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
bact<-subset(bact) # Goal is to run models for reach value in pathogen in bact in parallel, so no need to subset
target<-unique(bact$pathogen)
bact$yearn<-as.numeric(as.character(bact$year))
bact$year<-as.factor(bact$year)
bact<-split(bact, bact$pathogen)

##############################################################
# Parallel Model Fitting
##############################################################

print("--FITTING MODELS IN PARALLEL")
cl <- makeCluster(detectCores() - 2) 
registerDoParallel(cl)

results <- foreach(
  pathogen_name = names(bact),
  .combine = 'c',
  .packages = c('brms', 'dplyr'),
  .export = c('bact', 'outBase', 'PROPOSED_BM')
) %dopar% {
  print(paste0("Processing pathogen: ", pathogen_name))
  proposed <- PROPOSED_BM(bact[[pathogen_name]])
  saveRDS(proposed, paste0(outBase, pathogen_name, "_brm.Rds"))
  setNames(list(proposed), pathogen_name)
}

stopCluster(cl)

##############################################################
# Post-Model Processing
##############################################################

# Iterate over all pathogens in `bact`
for (pathogen_name in names(bact)) {
  print(paste0("--PROCESSING PATHOGEN: ", pathogen_name))
  
  # Generate posterior predictions for the current pathogen
  print("--GENERATING POSTERIOR PREDICTIONS")
  posteriorLinpred <- LINPREAD_DRAW_FN(
    data = results[[pathogen_name]]$data %>%
      group_by(state, yearn, population, count) %>%
      summarise(across(everything(), ~ unique(.)), .groups = "drop"),
    model = results[[pathogen_name]]
  )
  
  # Calculate catchment-level estimates
  print("--CALCULATING CATCHMENT-LEVEL ESTIMATES")
  catch <- CATCHMENT(posteriorLinpred)
  
  # Calculate incidence rate comparisons
  print("--CALCULATING INCIDENCE RATE COMPARISONS")
  catchir.linpred <- LINPRED_TO_CATCHIR(catch)
   # Check dimensions to match metadata assignment
  
  # Save catchment-level estimates for the current pathogen
  save_file <- paste0(outBase, pathogen_name, "_IRCatch.csv")
  SAFE_WRITE(catchir.linpred, save_file)
  
  # Perform incidence rate comparisons for the current pathogen
  IR_COMP(catchir.linpred, 2016, 2018, paste0(pathogen_name, "_EstIRRCatch_2016_2018.csv"), target = pathogen_name, travel = travelLabel, culture = culture)
  IR_COMP(catchir.linpred, 2020, 2022, paste0(pathogen_name, "_EstIRRCatch_2020_2022.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
  IR_COMP(catchir.linpred, 2004, 2006, paste0(pathogen_name, "_EstIRRCatch_2004_2006.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
  IR_COMP(catchir.linpred, 2006, 2008, paste0(pathogen_name, "_EstIRRCatch_2006_2008.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
  IR_COMP(catchir.linpred, 2010, 2012, paste0(pathogen_name, "_EstIRRCatch_2010_2012.csv"),target = pathogen_name, travel = travelLabel, culture = culture)
  
  print(paste0("--COMPLETED PROCESSING FOR PATHOGEN: ", pathogen_name))
}

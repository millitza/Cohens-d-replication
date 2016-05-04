### GENERATED DATA ###
### STEP 0: MAKE SURE THE CURRENT FOLDER IS SET AS WORKING DIRECTORY ###
### STEP 1: LOAD NECESSARY PACKAGE AND FUNCTIONS ###
source("function2.R")
source("BIEMS data/BIEMSdataPrep.R")

### STEP 2: LOAD DATA ###
# load data
BIEMSdatafiles <- list.files("BIEMS data/",pattern=".txt")
BIEMSdataclean <- gsub(".txt","",BIEMSdatafiles)
BIEMSdatapath <- paste0("BIEMS data/",BIEMSdatafiles)
BIEMSdata <- lapply(BIEMSdatapath, read.table)
# prepare data for format function
BIEMSdata <- lapply(BIEMSdata, BIEMSdataprep)

### STEP 3: CALCULATE BAYES FACTORS FOR EACH COMBINATION OF DATA SETS UNDER 3 CONDITIONS (PHI)
# specify values for phi and original data sets
phiset <- c(.2,.35,.5,.8)
datasetorig <- "N600_d00"
datasetrepl <- c("N600_d00","N600_d02","N600_d05","N600_d08")
# calculate Bayes factors for all conditions under specified phi values and save objects in folder Workspace
for (phi in phiset[1]){
  for (origdata in datasetorig) {
    for (repldata in datasetrepl){
    
      set.seed(62878)
      
      dataset_orig <- matrix(unlist(BIEMSdata[BIEMSdataclean == origdata]),ncol=3)
      dataset_repl <- matrix(unlist(BIEMSdata[BIEMSdataclean == repldata]),ncol=3)
      
      BFobject <- BF_repl_es(data_orig=dataset_orig,data_repl=dataset_repl,ddiff=phi,iter=60000)
      nameBFobject <- paste0("Workspace/","BF_orig",origdata,"_repl",repldata,"_phi",phi,".rda")
      assign(nameBFobject,BFobject)
      save(BFobject,file=nameBFobject)
    }
  }
}

# Load all Rdata files (output files that include the results) 
Outputfiles <- list.files("Workspace/")
Outputfilespath <- paste0("Workspace/",Outputfiles)
Outputfilesname <- strsplit(Outputfiles,"\\.r")

for(f in 1:length(Outputfiles)){
  load(Outputfilespath[f],.GlobalEnv)
  nameOutput <- Outputfilesname[f][[1]][1]
  assign(nameOutput,BFobject)
}
rm(BFobject)

### STEP 4: PRESENT RESULTS (TABLE 1) ###
# extract all information from returned objects
BFobjects <- ls(pattern="BF_origN*")
for (BFobj in BFobjects){
  elements_obj <- strsplit(BFobj, "_")[1][[1]][2:6]
  phi_obj <- as.numeric(strsplit(elements_obj[5],"i")[[1]][2])
  norig_obj <- as.numeric(strsplit(elements_obj[1],"N")[[1]][2])
  nrepl_obj <- as.numeric(strsplit(elements_obj[3],"N")[[1]][2])
  d_obj <- as.numeric(strsplit(elements_obj[4],"d")[[1]][2])/10
  BF12_obj <- as.numeric(get(BFobj)$BF12)
  BF13_obj <- as.numeric(get(BFobj)$BF13)
  
  output <- c(phi_obj,norig_obj,nrepl_obj,d_obj,BF12_obj,BF13_obj)
  nameoutput <- paste0("output",BFobj)
  assign(nameoutput,output)
}

# store all elements in a table, lable and sort the table
outputobjects <- ls(pattern="outputBF*")
table <- matrix(unlist(lapply(outputobjects,get)),ncol=6,byrow=T)
colnames(table) <- c(expression(phi),"n_orig","n_repl",expression(delta),"BF_12","BF_13")
table <- table[order(table[,1],table[,2],table[,3],table[,4]),]
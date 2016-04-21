### GENERATED DATA ###
### STEP 0: MAKE SURE THE CURRENT FOLDER IS SET AS WORKING DIRECTORY ###
### STEP 1: LOAD NECESSARY PACKAGE AND FUNCTIONS ###
source("function.R")
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
phiset <- c(.2,.5,.8)
datasetorig <- c("N20_d00","N40_d00","N100_d00")
# calculate Bayes factors for all conditions under specified phi values and save objects in folder Workspace
for (phi in phiset){
  for (origdata in datasetorig) {
    for (repldata in BIEMSdataclean){
    
      set.seed(62878)
      
      dataset_orig <- matrix(unlist(BIEMSdata[BIEMSdataclean == origdata]),ncol=3)
      dataset_repl <- matrix(unlist(BIEMSdata[BIEMSdataclean == repldata]),ncol=3)
      
      BFobject <- BF_repl_es(data_orig=dataset_orig,data_repl=dataset_repl,ddiff=phi,iter=60000)
      nameBFobject <- paste0("Workspace/","BF_orig",origdata,"_repl",repldata,"_phi",phi,".rda")
      #assign(nameBFobject,BFobject)
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

#### DATA EXAMPLE ####
### LOAD ALL DATA FILES AND PERFORM ANALYSES AS DISPLAYED IN TABLE 3 ###
# onesided vs twosided
example_orig_onesided_twosided <- read.table("Data example/Original study/Data_orig_onesided_twosided.txt",header=T)
example_repl_onesided_twosided <- read.table("Data example/Replication study/Data_repl_onesided_twosided.txt",header=T)

for(phi in phiset){
  set.seed(39580)
  BFobject <- BF_repl_es(example_orig_onesided_twosided,example_repl_onesided_twosided,ddiff=phi,iter=60000)
  nameBFobject <- paste0("dataexample_onesided_twosided_phi",phi)
  assign(nameBFobject,BFobject)
}

plot_prior(dataexample_onesided_twosided_phi0.2)
plot_prior(dataexample_onesided_twosided_phi0.5)
plot_prior(dataexample_onesided_twosided_phi0.8)

# neutral vs onesided
example_orig_neutral_onesided <- read.table("Data example 2/Original study/Data_orig_neutral_onesided.txt",header=T)
example_repl_neutral_onesided <- read.table("Data example 2/Replication study/Data_repl_neutral_onesided.txt",header=T)

for(phi in phiset){
  set.seed(39580)
  BFobject <- BF_repl_es(example_orig_neutral_onesided,example_repl_neutral_onesided,ddiff=phi,iter=60000)
  nameBFobject <- paste0("dataexample_neutral_onesided_phi",phi)
  assign(nameBFobject,BFobject)
}

plot_prior(dataexample_neutral_onesided_phi0.2)
plot_prior(dataexample_neutral_onesided_phi0.5)
plot_prior(dataexample_neutral_onesided_phi0.8)

# neutral vs twosided
example_orig_neutral_twosided <- read.table("Data example 2/Original study/Data_orig_neutral_twosided.txt",header=T)
example_repl_neutral_twosided <- read.table("Data example 2/Replication study/Data_repl_neutral_twosided.txt",header=T)

for(phi in phiset){
  set.seed(39580)
  BFobject <- BF_repl_es(example_orig_neutral_twosided,example_repl_neutral_twosided,ddiff=phi,iter=60000)
  nameBFobject <- paste0("dataexample_neutral_twosided_phi",phi)
  assign(nameBFobject,BFobject)
}

plot_prior(dataexample_neutral_twosided_phi0.2)
plot_prior(dataexample_neutral_twosided_phi0.5)
plot_prior(dataexample_neutral_twosided_phi0.8)

### effect size estimates as presented in Table 2 ###
d_orig_one_two <- sampler_ESrepl(example_orig_onesided_twosided,prior=matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T),type="orig",d_diff=0)
d_orig_neut_one <- sampler_ESrepl(example_orig_neutral_onesided,prior=matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T),type="orig",d_diff=0)
d_orig_neut_two <- sampler_ESrepl(example_orig_neutral_twosided,prior=matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T),type="orig",d_diff=0)

d_repl_one_two <- sampler_ESrepl(example_repl_onesided_twosided,prior=matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T),type="orig",d_diff=0)
d_repl_neut_one <- sampler_ESrepl(example_repl_neutral_onesided,prior=matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T),type="orig",d_diff=0)
d_repl_neut_two <- sampler_ESrepl(example_repl_neutral_twosided,prior=matrix(c(0,100,0,100,.001,.001),nrow=3,ncol=2,byrow=T),type="orig",d_diff=0)
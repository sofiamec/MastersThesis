#===================================================================================================================================
# SET-UP SECTION
#===================================================================================================================================

# Setup script for assigning correct names for a given design as well as preprocessing of datasets
# This script has to be run in combination with the "Run_entire_analysis"-script

####################################################################################################################################
#===================================================================================================================================
#                      Only change these parameters for different results! 
#                            (the rest of the code should adjust)
#===================================================================================================================================
# Test-settings (Remove this in final script!)
if (onTerra==F){
  repeats = 2                                               # sets the number of runs for each case (experimental design and q)
  limitNA = 0                                                   # the lowest amount of observaitons needed to produce a results other than NA
  savePlot = T                                              # use T when plots should be saved (for many repeats)
  #loadData = F                                              # use T if it is a rerun of existing results
  effectsizes=c(1.5)#,3)                                             # q = Fold-change for downsampling
  groupSize<-c(3)#,5)#,10,30)#,30,50)                                            # m = Number of samples in each group (total nr samples = 2*m)
  sequencingDepth<-c(10000)#, 20000)#,100000)#, 5000000)#,10000,1000000,5000000,10000000)      # d = Desired sequencing depth per sample
  sequencingDepthName<-c("10 k")#, "20 k")#, "100 k")#, "500 k")# "10 k","500 k","1 M", "5 M", "10 M")         # dD = Displayed names for sequencing depths
  sequencingSaveName<-c("10k")#, "20 k")#, "100 k")#, "500 k")# "10 k","500 k","1 M", "5 M", "10 M")         # dD = Displayed names for sequencing depths
}

# Real settings
if (onTerra==T){
  repeats = 100                                              # sets the number of runs for each case (experimental design and q)
  limitNA = 10                                                   # the lowest amount of observaitons needed to produce a results other than NA
  savePlot = T                                              # use T when plots should be saved (for many repeats)
  loadData = F                                              # use T if it is a rerun of existing results
  
  # flytta eventuellt allt detta till setup-skriptet
  if (saveName =="Resistance"){
    effectsizes=c(5,10)                                     # q = Fold-change for downsampling
  } else {
    effectsizes=c(1.5,3)                                   # q = Fold-change for downsampling
  }
  
  groupSize<-c(3,5,10,30,50)                                # m = Number of samples in each group (total nr samples = 2*m)
  # sequencing depths are set later depending on dataset    # d and dD = Desired sequencing depths and how it should be displayed
}

# Extra-settings
extraL=0                                                    # unless extraDesigns are added, the length of added designs is 0                               
if (extraDesigns==T){
  # The 3 following must have equal lengths!                # Combined they give more results for trade-off curves. Here with m*d = 3M, 3M and 5M respectively
  extraSeqDepth=c(200000,500000,250000)      
  extraSeqDepthName=c("200 k","500 k","250 k")
  extraSeqSaveName=c("200k","500k","250k")
  extraGroups=c(15,6,20)
  extraL<-length(extraGroups)
}
# Strata-settings
if (runStrata==T){
  numberOfStrata = 3                                        # sets the number of groups for dividing gene abundance and variability
  strataClass<-factor(c("low","medium","high"), levels=c("low","medium","high")) #strataClass<-c("low","medium", "high")                      # should correspond to the number of stratas
}
#===================================================================================================================================
#===================================================================================================================================
####################################################################################################################################

## Loading original data
#  filtering out samples not meeting expDesign requirements
#  filtering out genes with too low counts

if(saveName == "Gut2"){
  plotName = "Human Gut II"
  Gut2Original <- read.table("../../Data/Raw_data/HumanGutII_COGcountsRaw.txt", header=T, row.names = 1)
  Gut2Intermediate = Gut2Original[,colSums(Gut2Original)>=5000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Gut2 <- remove_low_counts(Gut2Intermediate)
  Data = Gut2
  
  boldvalue2="0"
  relations<-c(3000000, 5000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000)
    sequencingDepthName<-c("10 k","100 k","500 k","1 M", "5 M")
    sequencingSaveName<-c("10k","100k","500k","1M", "5M")
  }
  
  rm(Gut2, Gut2Original, Gut2Intermediate)        # remove original and intermediate datasets
  
} else if(saveName == "Marine"){
  plotName = "Marine"
  MarineOriginal <- read.table("../../Data/Raw_data/Marine_COGcountsRaw.txt", header=T, row.names = 1)
  MarineIntermediate = MarineOriginal[,colSums(MarineOriginal)>=10000000]   # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Marine <- remove_low_counts(MarineIntermediate)
  Data = Marine
  
  boldvalue2="5e+07"                                              
  relations<-c(3000000,5000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000,10000000)
    sequencingDepthName<-c("10 k","100 k","500 k","1 M", "5 M", "10 M")
    sequencingSaveName<-c("10k","100k","500k","1M", "5M", "10M")
  }
  
  rm(Marine, MarineOriginal, MarineIntermediate)  # remove original and intermediate datasets
  
} else if(saveName == "Resistance"){
  plotName = "Resistance"
  
  ResistanceOriginal=t(read_excel("../../Data/Raw_data/GENE_QUANTIFICATIONS.raw.xlsx")[,-c(2,4)])
  colnames(ResistanceOriginal) <- ResistanceOriginal[1,]
  # Extracting Human samples:
  ResistanceOriginal<-ResistanceOriginal[,ResistanceOriginal[2,]=="Airways" | ResistanceOriginal[2,]=="Gastrointestinal" | ResistanceOriginal[2,]=="Oral" | ResistanceOriginal[2,]=="Skin" | ResistanceOriginal[2,]=="Urogenital"]
  ResistanceOriginal<-ResistanceOriginal[-2,]
  ResistanceOriginal=data.frame(row.names = row.names(ResistanceOriginal)[-1], apply(ResistanceOriginal[-1,],2,as.integer))
  ResistanceIntermediate=ResistanceOriginal[,colSums(ResistanceOriginal)>=10000000] # Filter out samples with sequencing depth below the maximum sequencing depth of the experimental design
  Resistance <- remove_low_counts(ResistanceIntermediate) # This is the dataset used in analysis where samples and genes with low counts are removed
  
  Data = Resistance
    
  runStrata=F
  boldvalue2="5e+07"                                              
  relations<-c(3000000,5000000)
  if (onTerra==T){
    sequencingDepth<-c(10000,100000,500000,1000000,5000000,10000000) #c(10000,100000,500000,1000000,5000000,10000000)
    sequencingDepthName<-c("10 k","100 k","500 k","1 M", "5 M", "10 M") #c("10k","100k","500k","1M", "5M", "10M")
    sequencingSaveName<-c("10k","100k","500k","1M", "5M", "10M")
  }
  
  rm(Resistance, ResistanceOriginal, ResistanceIntermediate)  # remove original and intermediate datasets
}

if (runStrata==T){
  DataStrata<-DESeq2_for_strata(Data,numberOfStrata)
}
#===================================================================================================================================
## Setting upp the right environment and naming designs

AllSaveDesigns=array(rep(NaN, length(effectsizes)*(length(groupSize)+extraL)*(length(sequencingDepth)+extraL)),c(length(effectsizes),(length(groupSize)+extraL),(length(sequencingDepth)+extraL))) 
AllPlotDesigns=array(rep(NaN, length(effectsizes)*(length(groupSize)+extraL)*(length(sequencingDepth)+extraL)),c(length(effectsizes),(length(groupSize)+extraL),(length(sequencingDepth)+extraL))) 

for (effect in 1:length(effectsizes)) {           # looping over q
  q=effectsizes[effect]
  
  # Creating the standard designs
  for (group in 1:length(groupSize)){             # looping over m
    m=groupSize[group]
    for (seq in 1:length(sequencingDepth)) {      # looping over d
      d=sequencingDepth[seq]
      dD=sequencingDepthName[seq]
      dS=sequencingSaveName[seq]
      
      AllSaveDesigns[effect,group,seq] <- sprintf("m%d_d%s_10q%d_f%d",m, dS, q*10, f*100)
      AllPlotDesigns[effect,group,seq] <- sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
      
      for (analysis in analyses){
        # Create folder for certain case if it doesn't exist
        if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, AllSaveDesigns[effect,group,seq]))){
          dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", AllSaveDesigns[effect,group,seq])), recursive = T)
        }
        
        if (!dir.exists(sprintf("../../Result/%s_%s/IntermediatePlots", saveName, analysis))){
          dir.create(file.path("../../Result", sprintf("%s_%s", saveName, analysis), "/IntermediatePlots"), recursive = T)
        }
      }
      
    }
  }
  
  # Creating 3 extra designs                               
  if (extraDesigns==T){
    for (i in 1:extraL) {              # looping over extra designs with fixed m and d
      m=extraGroups[i]
      d=extraSeqDepth[i]
      dD=extraSeqDepthName[i]
      dS=extraSeqSaveName[i]
      
      AllSaveDesigns[effect,group+i,seq+i] <- sprintf("m%d_d%s_10q%d_f%d", m, dS, q*10, f*100)
      AllPlotDesigns[effect,group+i,seq+i] <- sprintf("m=%d, d=%s, q=%g, f=%d%%",m,dD,q,f*100)
      
      # Create folder for certain case if it doesn't exist
      if (!dir.exists(sprintf("../../Intermediate/%s/%s", saveName, AllSaveDesigns[effect,group+i,seq+i]))){
        dir.create(file.path("../../Intermediate", sprintf("%s", saveName), sprintf("%s", AllSaveDesigns[effect,group+i,seq+i])), recursive = T)
      }
      
      if (!dir.exists(sprintf("../../Result/%s_%s/IntermediatePlots", saveName, analysis))){
        dir.create(file.path("../../Result", sprintf("%s_%s", saveName,analysis), "/IntermediatePlots"), recursive = T)
      }
    }}
}

rm(d,dD,effect, group, m, q, seq)

# PCA

######################################  Libraries  #######################################
library(plyr)
library(dplyr)
library(Hmisc)
library(stringr)
####################################### End Libraries ####################################



####################################### USER INPUT #######################################
# Location of a folder with all of your .csvs from direct_infusion_processing.R script
# Make sure no other files are in this folder. 
fileLocation <- "C:/Users/Hegeman Lab/Desktop/PCA_Testing_Input/pos"

# name for the loadings csv that is optional to generate as a last step. 
# This will be written to the current working directory which you can see by 
# running the line of code two lines below this. 
csvName <- "loadings.csv"
getwd()

# Round to value.
valueToRoundTo <- 4

# Acceptable tolerance 
tolerance <- .001

# TRUE if you want the data to be normalized (single intensity/total intensity) before plotting
# FALSE to just use raw data. Note the analysis automatically scales and centers the data. For
# more information on what that means see the documentation for prcomp()
normalization <- TRUE

# number of treatments and a list of each name
numTreatments <- 3

# number of files (samples) per treatment 
numFilesPerTreatment <- 4

# The pc to be plotted on the x and y axis respectively
xpc <- 1
ypc <- 2
####################################### END INPUT #######################################



####################################### FUNCTIONS #######################################
# function that does all the data transformation and generates a "prcomp" class object
# input files just need a column with the label 'mz' and a column with the label 'intensity'
generatePCAObject <- function(file_path, round_to_digits, match_tolerance, normalize){
  
  # get current directory for returning back to it after processing
  returnPath <- getwd()
  
  # set working directory to location with all desired csv files.
  setwd(file_path)
  
  # get file names
  fileNames <- list.files()
  
  # place holder
  mzs <- c()
  mzsLists <- list()
  iLists <- list()
  
  # get all mzs and build list of mzs and intensities in each file
  for (f in fileNames){
    
    # get masses
    temp_table <- read.csv(f)
    mzs <- c(mzs, unlist(temp_table['mz']))
    mzsLists <- append(mzsLists, temp_table['mz'])
    
    # get normalized intensities or don't normalize
    if(normalize){ totalIntensities <- sum(as.numeric(unlist(temp_table['intensity']))) }
    else { totalIntensities <- 1 }
    
    iLists <- append(iLists, temp_table['intensity']/totalIntensities)
    
  }
  
  # Round and keep unique. Rough form of getting peaks across files
  mzs4 <- unique(round(mzs, digits = round_to_digits)) 
  
  # Start a dataframe with those rounded mzs as the first rows
  pcaDF <- data.frame(mz=mzs4)
  
  # find the best matches for each value in the rounded mzs list. 
  matches <- lapply(mzsLists, find.matches, x = mzs4, tol = tolerance, maxmatch = 1)
  
  # loop over each list of intensities and generate a new column in the data frame for that file
  for (i in 1:length(iLists)){
    
    # matches uses 0 if no match was found so check for that and return a numeric 0 if that's the case, otherwise
    # return the correct intensity. 
    intensities <- lapply(matches[[i]]$matches, FUN = function(x){
      
      # indexes of 0 don't work so if that's found, return 0
      if(x==0) {0}
      
      # else return the actual intensity value.
      else {iLists[[i]][x]}
      
    })
    
    # add column to dataframe using new list of intensities. 
    pcaDF[fileNames[i]] <- unlist(intensities)
    
  }
  
  # transform for formatting
  tDF <- t(pcaDF) #### Note formatting and could use this as generic start point. 
  
  # add column names (mzs)
  colnames(tDF) <- tDF[1,]
  
  # remove mz row
  tDF <- tDF[-1, ]
  
  # run PCA analysis an center the peak areas on zero and scales them such that they all have a sd of 1
  pca <- prcomp(tDF, center = TRUE, scale=TRUE, xret=TRUE) # Things to play around with
  
  setwd(returnPath)
  
  pca
}

# Plots PCA with rough legend provided
plotPCA <- function(prcomp_object, number_of_treatments, number_of_files_per_treatment, xPC, yPC, rough_legend=TRUE){
  
  pca <- prcomp_object
  
  # Generate list for plotting symbols
  symbolList <- c()
  for(i in 1:number_of_treatments){
    symbolList <- append(symbolList, rep(i, number_of_files_per_treatment))
  }
  
  # get pca summary
  psummary <- summary(pca)
  
  # get percent variance explained by the provided pc values
  xvariance <- round(psummary$importance[2, xPC]*100, 2)
  yvariance <- round(psummary$importance[2, yPC]*100, 2)
  
  # Generate plot
  plot(pca$x[,xPC], pca$x[, yPC], 
       pch=symbolList, 
       xlab=paste0("PC", as.character(xPC), ": ", as.character(xvariance), "%"),
       ylab=paste0("PC", as.character(yPC), ": ", as.character(yvariance), "%")
  )
  
  if(rough_legend) {
    
    # generate legend
    legend("topright", 
           str_sub(row.names(pca$x)[seq(from=1, to=nrow(pca$x), by=number_of_files_per_treatment)], 1, -5),
           pch = symbolList[seq(from=1, to=nrow(pca$x), by=number_of_files_per_treatment)],
           cex = .6
    ) 
  
  }
  
}

# Plots rotation (loading) values
plotLoadings <- function(prcomp_object, xPC, yPC){
  pca <- prcomp_object
  
  # get pca summary
  psummary <- summary(pca)
  
  # get percent variance explained by the provided pc values
  xvariance <- round(psummary$importance[2, xPC]*100, 2)
  yvariance <- round(psummary$importance[2, yPC]*100, 2)
  
  # Generate plot
  plot(pca$rotation[,xPC], pca$rotation[, yPC], 
       xlab=paste0("PC", as.character(xPC), ": ", as.character(xvariance), "%"),
       ylab=paste0("PC", as.character(yPC),": ", as.character(yvariance), "%")
  )
  
}
####################################### END FUNCTIONS #######################################


# generate the PCA object (prcomp class object to be specific)
prcompObject <- generatePCAObject(file_path = fileLocation, round_to_digits = valueToRoundTo, 
                                  match_tolerance = tolerance, normalize = normalization)

#Plot the PCA. The rough legend just provides the name of the first filename in that treatment
plotPCA(prcomp_object = prcompObject, number_of_treatments = numTreatments,
        number_of_files_per_treatment = numFilesPerTreatment, xPC = xpc, 
        yPC = ypc, rough_legend = FALSE)

# Plot the loadings
plotLoadings(prcomp_object = prcompObject, xPC = xpc, yPC = ypc)

# Run this to write the loadings to a .csv file
write.csv(prcompObject$rotation, file = csvName)




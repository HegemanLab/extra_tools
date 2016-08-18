# PCA with option for PLS-DA

######################################  Libraries  #######################################
library(plyr)
library(dplyr)
library(Hmisc)
library(stringr)
library(pls)
####################################### End Libraries ####################################



####################################### USER INPUT #######################################
# Location of a folder with all of your .csvs from direct_infusion_processing.R script
# Make sure no other files are in this folder. 
fileLocation <- "C:/Users/Hegeman Lab/Desktop/Tissue_Spray_August/PCA_input_and_results/pos"

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

# number of components to use in the geration of the pls model
numComponents <- 4

# The pc to be plotted on the x and y axis of the PCA 
xpc <- 1
ypc <- 2

# The components to be plotted on the x and y axis of the PLS-DA Score Plot 
xpls <- 1
ypls <- 2

# number of components to plot in generating your optional bar plot
numBars <- 5

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

# Makes a bar plot of the top Principle Components 
makeBarPlot <- function(prcomp_object, num_bars) {
  
  # get summary of prcomp object
  psummary <- summary(prcomp_object)
  
  # Generate a barplot with the PCs in ranked order
  barplot(psummary$importance[2, 1:num_bars]*100, xlab = "Principal Components", 
          ylab = "Proportion of Variance explained", main = "Percentage of Variance by PC", 
          ylim = c(0, max(psummary$importance[2, 1:num_bars]*100)+10))
}

# Makes a plot of the top 4 or fewer PCs and shows their scatter plots side by side. 
makePairPlot <- function(prcomp_object, number_of_treatments, number_of_files_per_treatment){
  
  # Get score data
  scores <- prcomp_object$x 
  
  # generate list of symbols
  symbolList <- c()
  for(i in 1:number_of_treatments){
    symbolList <- append(symbolList, rep(i, number_of_files_per_treatment))
  }
  
  # set how many pairs to make. 
  pairs = c()
  if (ncol(scores) >= 4) {
    pairs = c(4)
  } else pairs = c(ncol(scores))
  
  pairs(scores[, 1:pairs], pch = symbolList)
  
}

# generate a PLS object. Utilizes same feature matching approach as generatePCAObject but 
# requires additional input to help supervise the modeling process. 
generatePLSObject <- function(file_path, round_to_digits, match_tolerance, normalize, 
                                 number_of_treatments, number_of_files_per_treatment, 
                                 num_components_in_model){
  
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
  
  # Make groups list
  groups <- c()
  for(i in 1:number_of_treatments){
    groups <- append(groups, rep(i, number_of_files_per_treatment))
  }
  
  plsObj <- plsr(as.matrix(groups) ~ as.matrix(tDF), ncomp = num_components_in_model, validation = "LOO")
  
  setwd(returnPath)
  
  plsObj
}

# Plots a PLS-DA Score Plot
makeScorePlot <- function(pls_object, x_comp, y_comp, number_of_treatments,
                          number_of_files_per_treatment) {
  
  # Stores maximum and minimum values for selected components (here 1 and 2) and makes them into lists
  Max.pc1 = 1.1 * (max(pls_object$score[, x_comp]))
  Min.pc1 = 1.1 * (min(pls_object$score[, x_comp]))
  Mpc1 = c(Min.pc1 * 2, Max.pc1 * 2)
  Max.pc2 = 1.1 * (max(pls_object$score[, y_comp]))
  Min.pc2 = 1.1 * (min(pls_object$score[, y_comp]))
  Mpc2 = c(Min.pc2 * 2, Max.pc2 * 2)
  
  # Empty list for symbols
  symbolList <- c()
  
  # generate symbols for each treatment
  for(i in 1:number_of_treatments){
    symbolList <- append(symbolList, rep(i, number_of_files_per_treatment))
  }
  
  # make the actual plot # labels=dn
  scoreplot(pls_object,comp=c(x_comp,y_comp), pch = symbolList, cex=0.8,xlim = c(Min.pc1, Max.pc1), ylim = c(Min.pc2, Max.pc2), main = "PLS-DA Score Plot")
  
  axis(1, at = Mpc1, pos = c(0, 0), labels = FALSE, col = "grey", lwd = 0.7)
  axis(2, at = Mpc2, pos = c(0, 0), labels = FALSE, col = "grey", lwd = 0.7)
  
  
  
}

# Generates an S plot from the provided data
makeSPlot <- function(file_path, round_to_digits, match_tolerance, normalize, 
                      number_of_treatments, number_of_files_per_treatment, 
                      num_components_in_model) {
  
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
  
  # Make groups list
  groups <- c()
  for(i in 1:number_of_treatments){
    groups <- append(groups, rep(i, number_of_files_per_treatment))
  }
  
  plsObj <- plsr(as.matrix(groups) ~ as.matrix(tDF), ncomp = num_components_in_model, validation = "LOO")
  
  setwd(returnPath)

  # Makes a model matrix from the data in pls (from the pls analysis). Model matrix seems similar to normal matrix. 
  mo <-model.matrix(plsObj)
  
  # Correlation matrix between the scores from the pls data (selected number of components), and the matrix data
  cl <- cor(mo, scores(plsObj))
  
  # Assings d the number of variables being examined
  d <-ncol(tDF)
  
  # mz are the names of all the files being examined
  mz <- colnames(tDF)
  
  # ind contains numbers from 1 to the number of variables in the data
  ind <- 1:d
  
  # dataframe containing the indexes, the variable names, and the loading data and correlation data for the first component
  d1 <-data.frame(ind,mz,loadings(plsObj)[,1], cl[,1])
  
  # Some sort of sorting function
  mat.sort <- function(mat,n){
    mat[rank(mat[,n]),] <- mat[c(1:nrow(mat)),]
    return(mat)
  }
  
  # sort dataframe d1 on column 4 (in this case the correlation) using sorting function defined above
  ds1<- mat.sort(d1, 4)
  
  # pulls out sorted list of indexes from ds1
  index <- ds1[ ,1]
  
  # Loading
  xc1<- ds1[,3]
  
  # Correlation
  yc1<- ds1[,4]
  
  # Stores maximum and minimum values for selected components (here 1 and 2) and makes them into lists
  Max.pc1 = 1.1 * (max(xc1))
  Min.pc1 = 1.1 * (min(xc1))
  Mpc1 = c(Min.pc1 * 2, Max.pc1 * 2)
  Max.pc2 = 1.1 * (max(yc1))
  Min.pc2 = 1.1 * (min(yc1))
  Mpc2 = c(Min.pc2 * 2, Max.pc2 * 2)
  
  # Plots a scatter of xc1 vs ycl then prompts to see if it's okay and then regraphs and saves the plot. 
  plot(xc1,yc1,xlim = c(Min.pc1, Max.pc1), ylim = c(Min.pc2, Max.pc2), main = "PLS-DA S-plot")
  text(xc1,yc1, labels = index, cex = 0.5, pos = 1)
  
  
}

####################################### END FUNCTIONS #######################################



######################################### WORK FLOW #########################################

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

# Make a bar plot of of the top x (numBars) PCs
makeBarPlot(prcomp_object = prcompObject, num_bars = numBars)

# Make a pair plot. Note, the scales will be different. 
makePairPlot(prcomp_object = prcompObject, number_of_treatments = numTreatments,
             number_of_files_per_treatment = numFilesPerTreatment)

# generate the PLS object (object with class 'mvr')
plsObj <- generatePLSObject(file_path = fileLocation, round_to_digits = valueToRoundTo, 
                            match_tolerance = tolerance, normalize = normalization,
                            number_of_treatments = numTreatments, 
                            number_of_files_per_treatment = numFilesPerTreatment,
                            num_components_in_model = numComponents)

# Make a plot to check the RMSEP of your current model based on the supplied number of components
plot(RMSEP(plsObj), legendpos = "topright")

# Make the score plot
makeScorePlot(pls_object = plsObj, x_comp = xpls, y_comp = ypls, 
              number_of_treatments = numTreatments, 
              number_of_files_per_treatment = numFilesPerTreatment)

# Generate the S Plot for your data
makeSPlot(file_path = fileLocation, round_to_digits = valueToRoundTo, match_tolerance = tolerance,
          normalize = normalization, number_of_treatments = numTreatments,
          number_of_files_per_treatment = numFilesPerTreatment, num_components_in_model = numComponents)





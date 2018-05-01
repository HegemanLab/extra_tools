####################################################################################
#####This is a short walkthrough for importing MS2 data for analysis using XCMS#####
####################################################################################

##IMPORTANT NOTE###

## This method removes all MS2 metadata, such as parent/precursor ion data
## DO NOT use this for any kind of untargeted analysis. 
## USE ONLY for a targeted analysis for which you are certain of the m/z and 
## retention times of your fragment(s) of interest

## The following steps happen outside of R
## 1) use msconvert to convert your file from the vendor format to .mzML format
## 2) When making the conversion select and add an MS-level filter
## 3) set the levels to 2-2. This will export only MS2 scans. The affect of keeping
##    MS1 scans in the file is unknown, but might present some issues in extracting
##    the EICs for any fragments of interest
## 4) set output to the directory of your choice

## The following is a short R protocol for converting MS2 scans to MS1 form

##First set your working directory to the directory containing your files
setwd('/your directory here')

temp <- list.files(pattern = '*.mzML') ##makes a list of the .mzML files in your directory
temp2 <- lapply(temp, xcmsRaw, includeMSn = TRUE) ## import .mzML files with MS2 scans

## you will probably get errors after this. This is normal as XCMS cannot make a profile
## matrix without MS1 data. This will not affect the conversion

temp3 <- lapply(temp2, msn2xcmsRaw) ## makes a list of new xcmsRaw objects with MS1 spectra
                                    ## filled with the MS2 scans. Note all MS2 metadata is lost

setwd('/your output directory here') ##set working directory to desired output directory

## This will write the new .mzML files to your chosen working directory
for (i in 1:length(temp)) {
  write.mzdata(temp3[[i]], filename = temp[[i]])
}

## the .mzML files may now be used in XCMS as desired
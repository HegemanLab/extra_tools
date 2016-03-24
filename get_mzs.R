# Helper script to extract a list of mz values from a list of files
library(xcms)

setwd("C:/Users/Lab/Desktop/Coding_Bits/VanKrevelen")


get_mzs <- function(in_file) {
  
  xraw <- xcmsRaw(in_file)
  
  peak_data <- findPeaks(xraw)
  
  peak_data <- peak_data@.Data
  
  peak_data[, 1]
  
}


files <- c(
  "ACM_Feb6_244-pos.mzXML",
  "ACM_Feb6_268-pos.mzXML",
  "ACM_Feb6_277-pos.mzXML",
  "ACM_Feb6_B2_255-pos.mzXML",
  "ACM_Feb6_B2_274-pos.mzXML",
  "ACM_Feb6_B3_270-pos.mzXML"
)

mzs <- lapply(files, get_mzs)

lapply(1:length(mzs), function(x){
  write.table(mzs[x], file = paste(files[x], '.csv'), col.names = "mzs", row.names = F)
})



output <- unlist(mzs)

setwd("C:/Users/Lab/Desktop/Tissue_Spray_Final/Tri/Input/")
write.table(output, file = "Tri-pos-mzs.txt", row.names = F, col.names = F)



get_both_mzs <- function(in_file) {
  
  xraw <- xcmsRaw(in_file)
  
  
  
  peak_data <- findPeaks(xraw)
  
  peak_data <- peak_data@.Data
  
  peak_data[, 1]
  
}


#Single File 
output <- get_mzs("ACM_sept16_T1R3_GL20_method1.mzXML")
setwd("C:/Users/Lab/Desktop")
write.table(output, file = "ACM_sept16_T1R3_GL20_method1.mzXML-mzs.txt", row.names = F, col.names = F)
setwd("C:/Users/Lab/Desktop/Coding_Bits/VanKrevelen")





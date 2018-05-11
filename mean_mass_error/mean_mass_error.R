# Mark Esler, Hegeman Lab. UMN
#
# Calculate mean mass error of a mass spectrometry run
# by comparing  observed and exact masses
#
# The first three columns of tabular data must be: (1) 
# standard (2) formula and (3) exact_mass
#
# All succeeding columns are samples of observed masses 
# with each row corresponding to a standard
#
# See example files

df <- read.table('input-data.tsv', header=T)

# calculate PPM mass error for standards in each sample
# append mass errors as column to df per sample
i = 1
ppm <- data.frame()
for(mass in df[,-c(0:3)]) {
  name = paste(colnames(df[3+i]),"_ppm",sep="")
  ppm <- c(abs((df$exact_mass-mass)/df$exact_mass*1e6))
  df[[name]] <- ppm
  i=i+1
}

# calculate average mass error for each standard append
# as column to df
n = (length(df)-3)/2
df$avg_ppm <- c(rowMeans(df[,(3+n):(3+2*n)]))

print(paste("Average PPM mass error:",mean(df$avg_ppm)))

write.csv(df, file='output-data.csv')

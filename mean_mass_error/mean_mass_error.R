# Hegeman Lab UMN
# Mark Esler, 2018
#
# Calculate mean mass error of a mass spectrometry run
# by comparing observed and known masses of standards.
#
# The first three columns of the input file are named:
# (1) standard (2) formula & (3) known_mass.
# All succeeding columns are samples.
# Sample columns can be given any unique name.
# Each row corresponds to a standard.
# See example input and output files.

df <- read.table('input-data.tsv',header=T)
#df <- read.table('output-data.tsv', header=T)

# Calculate PPM mass error for standards in each sample.
# Append mass errors as a column to df per sample.
i = 1
ppm <- data.frame()
for(mass in df[,-c(1:3)]) {
  ppm <- c(abs((df$known_mass-mass)/df$known_mass*1e6))
  name = paste(colnames(df[3+i]),"_ppm",sep="")
  df[[name]] <- ppm
  i=i+1
}

# Calculate mean mass error for each standard.
# Append as column to df.
n = (length(df)-3)/2
df$avg_ppm <- c(rowMeans(df[,(3+n):(3+2*n)]))

print(paste("Average PPM mass error:",mean(df$avg_ppm)))

write.table(df,file='output-data.tsv',sep="\t",row.names=F)

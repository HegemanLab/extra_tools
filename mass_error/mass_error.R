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

# Calculate PPM mass error for standards in each sample.
# Append mass errors as a column to df per sample.
i = 1
ppm <- data.frame()
# copy without columns 1 thru 3
mass_df <- df[,-c(1:3)]
n <- ncol(mass_df)
for(mass in mass_df) {
  ppm <- abs((df$known_mass-mass)/df$known_mass*1e6)
  name = paste(colnames(df[3+i]),"_ppm",sep="")
  df[[name]] <- ppm
  i=i+1
}

# Calculate mean mass error for each standard.
# Append as column to df.
avg_ppm <- rowMeans(df[,(3+n+1):(3+n+n)])
var_pm <- rowSums( df[,(3+n+1):(3+n+n)]^2 ) / (n - 1 )
df$avg_ppm <- avg_ppm
df$var_pm <- var_pm
df$n <- n

print(paste("Average PPM mass error:",mean(df$avg_ppm)))
print(paste("Overall SD per million:",sqrt(sum(df$var_pm)/nrow(mass_df))))

#install.packages("ggplot2",repos = "http://cran.us.r-project.org")
library(ggplot2)
write.table(df,file='output-data.tsv',sep="\t",row.names=F)
p <- ggplot(df, aes(known_mass, avg_ppm)) + geom_point() + expand_limits(y=0)
ggsave("output-data.png")

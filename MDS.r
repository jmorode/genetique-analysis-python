# MDS.r
#input: matrix pairwise distance between genotypes
#output: MDS plot
# --------------------------------------------------------------------------

# Input file
inputfile <- "projects/analysis_jan_2025_18loci/outputs/pairwise_differences/pairwise_differences_everything_all.csv"

# Read csv and prepare dataframe
df <- data.frame(read.csv(inputfile, sep=";", stringsAsFactors=T))
row.names(df) <- df$X
df[is.na(df)] <- 0

# Remove Year from Population 
df$Population <- as.character(df$Population)
df$POP <- sapply(strsplit(df$Population, " - "), function(x) x[1])
df$Population <- NULL
df$POP <- factor(df$POP)

subdf <- df
subdf$X <- NULL
subdf$POP <- NULL
subdf$Sample <- NULL


# Apply cmd
mds <- cmdscale(subdf, k=2, eig=TRUE)

# Plot result
print(levels(df$POP))
pchs = c(15,15,15,15,15,15,15, 4)
cols= c('green', 'blue', 'red', 'orange', 'purple', 'brown', 'cyan', 'gray')
plot(mds$points[,1], mds$points[,2], col = cols[df$POP], pch = pchs[df$POP], xlab="Dim. 1", ylab="Dim. 2")
legend('topleft', col=cols, legend=levels(df$POP), pch = pchs, cex = 0.75, bty = 'n')







args <- commandArgs(trailingOnly = TRUE)

exp_id = args[1]
ref_file= args[2]
output_dir = args[3]

print(exp_id)
print(ref_file)
print(output_dir)

library("FlowSorted.Blood.450k")
library('IlluminaHumanMethylation450kmanifest')
library('minfi')
library('BiocParallel')

options(warn = 0)
idat_files = as.matrix(read.table(ref_file, header=FALSE, sep='', dec='.'))

rgset = read.metharray(idat_files, verbose=TRUE, force=TRUE)
rgset.ssNoob = preprocessNoob(rgset, dyeMethod="single")
grset = mapToGenome(rgset.ssNoob)

# calculate beta matrix 
b_matrix = getBeta(rgset.ssNoob)

# returns dataframe, median Meth / median unmeth / predicted sex
sex_QC = minfiQC(rgset.ssNoob, fixOutliers=TRUE)

# returns dataframe
cell_counts = estimateCellCounts(rgset, processMethod="preprocessNoob")


gz1 <- gzfile(paste(output_dir, exp_id, '_methmatrix.gz', sep=''), "w")
write.csv(b_matrix, gz1)
close(gz1)

gz2 <- gzfile(paste(output_dir, exp_id, '_qc.gz', sep=''), "w")
write.csv(sex_QC$qc, gz2)
close(gz2)

gz3 <- gzfile(paste(output_dir, exp_id, '_cell_counts.gz', sep=''), "w")
write.csv(cell_counts, gz3)
close(gz3)

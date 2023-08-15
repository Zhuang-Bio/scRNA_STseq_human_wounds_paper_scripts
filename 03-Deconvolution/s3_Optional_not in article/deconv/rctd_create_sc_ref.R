library(RCTD)
library(Matrix)
library(optparse)

# parse args
parser <- OptionParser()
option_list = list(
  make_option(c("-i","--input"), action="store", help="single cell prefix", type="character"),
  make_option(c("-o","--outdir"), action="store", help="output directory", type="character")
)
opt = parse_args(OptionParser(option_list=option_list))
#print(opt)
      
#Rscript rctd_create_sc_ref.R  -i ../../deconv/inputs/sc_data/s1_subsampled -o ../../deconv/results/rctd/sc_ref 

# s1_subsampled_counts.mtx
# s1_subsampled_features.csv
# s1_subsampled_meta.csv

# Setup direcotries etc.
dir.create(opt$outdir, showWarnings = F)

#Check for sc files
extensions = c("_counts.mtx","_features.csv", "_meta.csv")
for (e in extensions){
  if(!file.exists(paste0(opt$input, e))){
    stop(sprintf("ERROR: No file %s", paste0(opt$input, e)))
  }
}

scC = readMM(paste0(opt$input, "_counts.mtx"))
# convert to dgCMatrix
scC = scC*1

genes = read.table(paste0(opt$input, "_features.csv"), header = F)
rownames(scC) = genes[,1]

scMeta = read.csv(paste0(opt$input,  "_meta.csv"))
colnames(scC) = scMeta$barcode


# create the reference object.

umi = scMeta$nUMI
names(umi) = scMeta$barcode
cl = factor(scMeta$cluster)
names(cl) = scMeta$barcode

reference = Reference(scC, cl, umi)

saveRDS(reference, file = file.path(opt$outdir, "reference.rds"))


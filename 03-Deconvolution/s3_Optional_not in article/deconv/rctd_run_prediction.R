library(RCTD)
library(Matrix)
library(optparse)

# parse args
parser <- OptionParser()
option_list = list(
  make_option(c("-i","--indir"), action="store", help="input directory", type="character"),
  make_option(c("-o","--outdir"), action="store", help="output directory", type="character"),
  make_option(c("-r","--reference"), action="store", help="single cell reference", type="character"),
  make_option(c("-s","--sample"), action="store", help="ST sample to run on", type="character")
)
opt = parse_args(OptionParser(option_list=option_list))
print(opt)
      
#Rscript rctd_run_prediction.R -i ../../deconv/inputs/st_data/ -o ../../deconv/results/rctd/ -r ../../deconv/results/rctd/sc_ref -s Donor1_Wound7

# Setup direcotries etc.
dir.create(opt$outdir, showWarnings = F)

sample = opt$sample
outdir = file.path(opt$outdir, sample)
dir.create(outdir, showWarnings = F)

sc_path = file.path(opt$reference, "reference.rds")
if(!file.exists(sc_path)){
  stop(sprintf("ERROR: No sc ref file %s",sc_path))
}

if (!file.exists(file.path(opt$indir,sample,"counts.mtx"))) {
  stop(sprintf("ERROR: No sc ref file %s",file.path(opt$indir,sample,"counts.mtx")))
}


# read ST data
stC = readMM(file.path(opt$indir,sample,"counts.mtx"))
stC = stC*1
# read gene names
genes = read.csv(file.path(opt$indir, "st_features.csv"), header = F)
rownames(stC) = genes[,1]
# coordinates
coord = read.csv(file.path(opt$indir,sample,"coordinates.csv"))
colnames(stC) = coord$X


# load SC ref
print("loading sc ref")
reference = readRDS( file.path(opt$reference, "reference.rds"))


# create ST object
print("reading st data")
nUMI = colSums(stC)
cc = data.frame(xcoord=coord$imagerow, ycoord=coord$imagecol)
rownames(cc) = coord$X
puck <- SpatialRNA(cc, stC, nUMI)

# run RCTD
print("running rctd")
myRCTD <- create.RCTD(puck, reference)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

results <- myRCTD@results

write.csv(as.matrix(results$weights), file = file.path(outdir,"RCTD_res.csv"))
saveRDS(myRCTD, file = file.path(outdir, "RCTDobj.rds"))


#Error in h(simpleError(msg, call)) :
#  error in evaluating the argument 'x' in selecting a method for function 'as.matrix': object 'results' not found
#Calls: write.csv ... is.data.frame -> as.matrix -> .handleSimpleError -> h
#Execution halted


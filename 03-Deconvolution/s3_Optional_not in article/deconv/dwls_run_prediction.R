library(Giotto)
library(data.table)
library(optparse)

# parse args
parser <- OptionParser()
option_list = list(
  make_option(c("-i","--indir"), action="store", help="input directory", type="character"),
  make_option(c("-o","--outdir"), action="store", help="output directory", type="character"),
  make_option(c("-r","--reference"), action="store", help="single cell reference Sig matrix", type="character"),
  make_option(c("-s","--sample"), action="store", help="ST sample to run on", type="character"),
  make_option(c("-p","--python"), action="store", help="path to conda python", type="character")
)
opt = parse_args(OptionParser(option_list=option_list))
#print(opt)
      
#Rscript dwls_run_prediction.R -i ../../deconv/inputs/st_data/ -o ../../deconv/results/dwls/ -r ../../deconv/inputs/sc_data/degs_fc0.5_pval0.01_sig.rds -s Donor1_Wound7 -p /Users/asbj/miniconda3/envs/R_giotto/bin/python3

#python_path = "/Users/asbj/miniconda3/envs/R_giotto/bin/python3"
instrs = createGiottoInstructions(python_path = opt$python,
                                  show_plot = F, return_plot = T, save_plot = T,
                                  dpi = 300, height = 9, width = 9)


# Setup direcotries etc.
dir.create(opt$outdir, showWarnings = F)

sample = opt$sample
outdir = file.path(opt$outdir, sample)
dir.create(outdir, showWarnings = F)

indir = file.path(opt$indir, sample)
raw_path = readLines(file.path(indir,"rawdata_path.csv"))

if(!file.exists(opt$reference)){
  stop(sprintf("ERROR: No sc ref file %s",opt$reference))
}

if (!file.exists(raw_path)) {
  stop(sprintf("ERROR: No sc ref file %s",raw_path))
}


# read ST data
print("Reading ST data...")
barcodes = read.csv(file.path(indir, "barcodes.csv"), header = F)
raw_matrix<-get10Xmatrix_h5(file.path(raw_path, "filtered_feature_bc_matrix.h5"), gene_ids = "symbols")
# gives a list with dgCMatrix in [[1]]
raw_matrix = raw_matrix[[1]]

spatial_results=data.table::fread(file.path(raw_path, "spatial", "tissue_positions_list.csv"))
colnames(spatial_results) = c('barcode', 'in_tissue', 'array_row', 'array_col', 'col_pxl', 'row_pxl')

m1 = match(barcodes[,1], spatial_results$barcode)
spatial_results = spatial_results[m1,]

m2 = match(barcodes[,1], colnames(raw_matrix))
raw_matrix = raw_matrix[,m2]

st_data <- createGiottoObject(raw_exprs = raw_matrix,spatial_locs = spatial_results[,.(row_pxl,-col_pxl)],
                                   instructions = instrs,
                                   cell_metadata = spatial_results[,.(in_tissue, array_row, array_col)])

# OBS! Need cluster column in st data to run DWLS, so run with the giotto pipeline:
print("Processing ST data...")
st_data <- normalizeGiotto(gobject = st_data)
st_data <- calculateHVG(gobject = st_data)
gene_metadata = fDataDT(st_data)
featgenes = gene_metadata[hvg == 'yes']$gene_ID
st_data <- runPCA(gobject = st_data, genes_to_use = featgenes, scale_unit = F)
#signPCA(st_data, genes_to_use = featgenes, scale_unit = F)
st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 15)
st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)

# load SC ref
print("loading sc ref...")
Sig = readRDS(opt$reference)



# run DWLS
print("Running DWLS....")
st_data <- runDWLSDeconv(st_data,sign_matrix = Sig, n_cell = 20)

write.csv(st_data@spatial_enrichment$DWLS, file = file.path(outdir, "dwls_proportions.csv"))



###### HELPERS #######

# Find script directory
# -------------------------------
get_script_dir <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  fileArgName <- "--file="
  scriptPath <- sub(fileArgName, "", cmdArgs[grep(fileArgName, cmdArgs)])
  
  if (length(scriptPath) > 0) {
    return(dirname(normalizePath(scriptPath)))
  }
  
  return(getwd())
}

# set up script to create .bat file
# download to local directory
# R.home to get where is R and then create .bat file
# source setup from Dora's github page set up script
# then have .bat file on desktop
# then download from mageck and 
 
# user.id=Sys.info()[["user"]]
# wd = paste0('C:/Users/',user.id, '/Box/Dora/MCD_CAGE_UI')
# setwd(wd)
# print('Hello World!')
script_dir <- get_script_dir()

###### START HERE #######

library(utils)
cat("Select counts data file \n")
cnt.path = choose.files(caption = "Select counts data file")
cnts = read.csv(cnt.path)

cat("First few lines of Counts Matrix \n")
head(cnts)

cat("Dimensions of Counts Matrix: \n")
dim(cnts)

###### AUTODETECT METADATA #######
# create metadata and allow user editing

# look for key words in counts matrix
# TODO: check if Z first, then beta
# TODO: If it is not a mageck result, then a counts file w/ column names

# zb.names = grep('z|beta', colnames(cnts), value =T)
zb.names = grep('z', colnames(cnts), value =T)
cnts = cnts[, zb.names]
samps = strsplit(zb.names, split = '\\.|_')
meta = as.data.frame(do.call(rbind, samps),
                    stringsAsFactors = FALSE)

colnames(meta) = c("cell line", "treatment", "treatment type")
meta$sample = zb.names

# TODO: try and find out treatment type and assign to treatment type
cat('Metadata from samples looks like: \n')
head(meta)
cat('Edit metadata file \n')
cat('Columns for cell line, treatment, treatment type are required')
meta = edit(meta)

# there has to be a cell line and treatment column and treatment type column
# stop if not present and provide

# now create treatment - control
cell_col = "cell line"
type_col = "treatment type"
sample_col = "sample"

lsplit = split(meta, meta$`cell line`)

res = lapply(lsplit, function(l){
  
  
  ctrl  = l[l$`treatment type` == "control", sample_col]
  treat = l[l$`treatment type` == "treatment", sample_col]
  
  # check that these two are subtractable
  if (length(ctrl) != 1)
    stop("Each cell line must have exactly one control")
  
  
  out = cnts[, treat, drop = FALSE] - cnts[, ctrl,  drop = FALSE]
  colnames(out) = paste0(treat, "_vs_", ctrl)
  
  out
  
})

cdiff = do.call(cbind, res)


###### Choose Return Directory ############
try(file.choose(), silent = TRUE)
res.dir = choose.dir(default = wd)

timestamp = format(Sys.time(), "%Y%m%d_%H%M")

fname = tail(unlist(strsplit(cnt.path, split = '\\\\')), 1)
fname = sub("\\.[^.]+$", "", fname)
meta.fname = file.path(
  res.dir, 
  paste0("sample_metadata_", fname, '_',timestamp, ".csv"))
write.csv(df,meta.fname)

# produce excel file 
# md5sum


############# Run Pipeline ###################

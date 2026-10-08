#  Packages and Python configuration

library(reticulate)
library(RaceID)
library(Matrix)


python_path <- Sys.getenv("PYTHON_PATH")

if (python_path == "") {
  stop(
    "PYTHON_PATH is not set. ",
    "Set it to your local Python executable."
  )
}

reticulate::use_python(
  python_path,
  required = TRUE
)

modules_available <-
  reticulate::py_module_available("leidenalg") &&
  reticulate::py_module_available("igraph")

if (!modules_available) {
  warning(
    "Python modules 'leidenalg' and/or 'igraph' ",
    "are not available."
  )
}

if (!reticulate::py_available()) {
  stop("Python is not available through reticulate.")
}


#  Local directories

data_dir <- Sys.getenv("DATA_DIR")
output_dir <- Sys.getenv("OUTPUT_DIR")

if (data_dir == "") {
  stop(
    "DATA_DIR is not set. ",
    "Set it to the directory containing your input data."
  )
}

if (output_dir == "") {
  output_dir <- file.path(
    getwd(),
    "results"
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



# Function for loading count matrices


load_counts <- function(
    sample_dir,
    min_counts = 1000
) {
  
  mtx_file <- file.path(
    sample_dir,
    "result_dev.mtx"
  )
  
  genes_file <- file.path(
    sample_dir,
    "result_dev.genes.txt"
  )
  
  barcodes_file <- file.path(
    sample_dir,
    "result_dev.barcodes.txt"
  )
  
  
  required_files <- c(
    mtx_file,
    genes_file,
    barcodes_file
  )
  
  
  missing_files <- required_files[
    !file.exists(required_files)
  ]
  
  
  if (length(missing_files) > 0) {
    
    stop(
      "Missing input file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  }
  
  
  x <- readMM(mtx_file)
  
  f <- read.csv(
    genes_file,
    sep = "\t",
    header = FALSE
  )
  
  b <- read.csv(
    barcodes_file,
    sep = "\t",
    header = FALSE
  )
  
  
  xM <- as(
    x,
    "dgCMatrix"
  )
  
  
  dimnames(xM) <- list(
    as.character(b$V1),
    as.character(f$V1)
  )
  
  
  # Transpose so that genes are rows
  # and cells are columns.
  xM <- t(xM)
  
  
  # Filter cells by total counts.
  cs <- colSums(xM)
  
  xM <- xM[
    ,
    cs > min_counts,
    drop = FALSE
  ]
  
  
  return(xM)
}


# Sample configuration


sample_dirs <- list(
  
  sample_01 = file.path(
    data_dir,
    "sample_01"
  ),
  
  sample_02 = file.path(
    data_dir,
    "sample_02"
  ),
  
  sample_03 = file.path(
    data_dir,
    "sample_03"
  ),
  
  sample_04 = file.path(
    data_dir,
    "sample_04"
  ),
  
  sample_05 = file.path(
    data_dir,
    "sample_05"
  ),
  
  sample_06 = file.path(
    data_dir,
    "sample_06"
  ),
  
  sample_07 = file.path(
    data_dir,
    "sample_07"
  ),
  
  sample_08 = file.path(
    data_dir,
    "sample_08"
  ),
  
  sample_09 = file.path(
    data_dir,
    "sample_09"
  ),
  
  sample_10 = file.path(
    data_dir,
    "sample_10"
  ),
  
  sample_11 = file.path(
    data_dir,
    "sample_11"
  ),
  
  sample_12 = file.path(
    data_dir,
    "sample_12"
  ),
  
  sample_13 = file.path(
    data_dir,
    "sample_13"
  ),
  
  sample_14 = file.path(
    data_dir,
    "sample_14"
  ),
  
  sample_15 = file.path(
    data_dir,
    "sample_15"
  ),
  
  sample_16 = file.path(
    data_dir,
    "sample_16"
  )
)



# Load count matrices


sample_01 <- load_counts(
  sample_dirs$sample_01
)

sample_02 <- load_counts(
  sample_dirs$sample_02
)

sample_03 <- load_counts(
  sample_dirs$sample_03
)

sample_04 <- load_counts(
  sample_dirs$sample_04
)

sample_05 <- load_counts(
  sample_dirs$sample_05
)

sample_06 <- load_counts(
  sample_dirs$sample_06
)

sample_07 <- load_counts(
  sample_dirs$sample_07
)

sample_08 <- load_counts(
  sample_dirs$sample_08
)

sample_09 <- load_counts(
  sample_dirs$sample_09
)

sample_11 <- load_counts(
  sample_dirs$sample_11
)

sample_12 <- load_counts(
  sample_dirs$sample_12
)

sample_13 <- load_counts(
  sample_dirs$sample_13
)

sample_14 <- load_counts(
  sample_dirs$sample_14
)

sample_15 <- load_counts(
  sample_dirs$sample_15
)

sample_16 <- load_counts(
  sample_dirs$sample_16
)



# Assign anonymized cell/sample labels


colnames(sample_11) <- paste(
  colnames(sample_11),
  "NC11",
  sep = "/"
)

colnames(sample_15) <- paste(
  colnames(sample_15),
  "NC10",
  sep = "/"
)

colnames(sample_14) <- paste(
  colnames(sample_14),
  "NC12",
  sep = "/"
)

colnames(sample_13) <- paste(
  colnames(sample_13),
  "NC9",
  sep = "/"
)

colnames(sample_05) <- paste(
  colnames(sample_05),
  "C1",
  sep = "/"
)

colnames(sample_02) <- paste(
  colnames(sample_02),
  "C4",
  sep = "/"
)

colnames(sample_10) <- paste(
  colnames(sample_10),
  "P1",
  sep = "/"
)

colnames(sample_04) <- paste(
  colnames(sample_04),
  "C3",
  sep = "/"
)

colnames(sample_06) <- paste(
  colnames(sample_06),
  "NC6",
  sep = "/"
)

colnames(sample_03) <- paste(
  colnames(sample_03),
  "NC8",
  sep = "/"
)

colnames(sample_07) <- paste(
  colnames(sample_07),
  "NC7",
  sep = "/"
)

colnames(sample_01) <- paste(
  colnames(sample_01),
  "NC5",
  sep = "/"
)

colnames(sample_09) <- paste(
  colnames(sample_09),
  "NC2",
  sep = "/"
)

colnames(sample_08) <- paste(
  colnames(sample_08),
  "NC3",
  sep = "/"
)

colnames(sample_12) <- paste(
  colnames(sample_12),
  "NC4",
  sep = "/"
)

colnames(sample_16) <- paste(
  colnames(sample_16),
  "NC1",
  sep = "/"
)

colnames(sample_10) <- paste(
  colnames(sample_10),
  "C2",
  sep = "/"
)



# Combine datasets


prdata <- cbind(
  sample_14,
  sample_13,
  sample_16,
  sample_15,
  sample_05,
  sample_02,
  sample_03,
  sample_07,
  sample_04,
  sample_06,
  sample_01,
  sample_10,
  sample_12,
  sample_08,
  sample_09,
  sample_16,
  sample_10
)




# Initial clustering


sc <- SCseq(
  prdata
)

sc <- filterdata(
  sc,
  mintotal = 1000,
  CGenes = rownames(sc@expdata)[
    grep(
      "^(MT|RP(L|S)|GM\\D|GYPA)",
      rownames(sc@expdata)
    )
  ]
)

expData <- getExpData(sc)

res <- pruneKnn(
  expData,
  large = TRUE,
  regNB = TRUE,
  knn = 25,
  seed = 12345,
  no_cores = 32,
  do.prune = FALSE
)

cl <- graphCluster(
  res,
  pvalue = 0.01,
  use.leiden = TRUE,
  leiden.resolution = 1.5
)

sc <- updateSC(
  sc,
  res = res,
  cl = cl,
  flo = 0.1
)

sc <- compumap(
  sc
)

plotmap(
  sc,
  um = TRUE,
  cex = 0.5
)



# Remove low-quality cells / doublets


C <- names(sc@cpart)[
  sc@cpart %in% c(
    7,
    13,
    14
  )
]

TF <- !colnames(sc@ndata) %in% C

prdata_clean <- prdata[
  ,
  TF,
  drop = FALSE
]



# Re-run clustering on cleaned data


sc <- SCseq(
  prdata_clean
)

sc <- filterdata(
  sc,
  mintotal = 1000,
  CGenes = rownames(sc@expdata)[
    grep(
      "^(MT|RP(L|S)|GM\\D|GYPA)",
      rownames(sc@expdata)
    )
  ]
)

expData <- getExpData(
  sc
)

res <- pruneKnn(
  expData,
  large = TRUE,
  regNB = TRUE,
  knn = 25,
  seed = 12345,
  no_cores = 32,
  do.prune = FALSE,
  pcaComp = 40
)

plotPC(
  res
)

cl <- graphCluster(
  res,
  pvalue = 0.01,
  use.leiden = TRUE,
  leiden.resolution = 1.5
)

sc <- updateSC(
  sc,
  res = res,
  cl = cl,
  flo = 0.1
)

sc <- compumap(
  sc
)

plotmap(
  sc,
  um = TRUE,
  cex = 0.5
)



# Second low-quality cell / doublet filtering


C <- names(sc@cpart)[
  sc@cpart %in% c(
    20,
    42
  )
]

TF <- !colnames(sc@ndata) %in% C

prdata_clean2 <- prdata_clean[
  ,
  TF,
  drop = FALSE
]



#  Final clustering on second cleaned dataset

sc <- SCseq(
  prdata_clean2
)

sc <- filterdata(
  sc,
  mintotal = 1000,
  CGenes = rownames(sc@expdata)[
    grep(
      "^(MT|RP(L|S)|GM\\D|GYPA)",
      rownames(sc@expdata)
    )
  ]
)

expData <- getExpData(
  sc
)

res <- pruneKnn(
  expData,
  large = TRUE,
  regNB = TRUE,
  knn = 25,
  seed = 12345,
  no_cores = 32,
  do.prune = FALSE,
  pcaComp = 40
)

plotPC(
  res
)

cl <- graphCluster(
  res,
  pvalue = 0.01,
  use.leiden = TRUE,
  leiden.resolution = 1.5
)

sc <- updateSC(
  sc,
  res = res,
  cl = cl,
  flo = 0.1
)

sc <- compumap(
  sc
)

plotmap(
  sc,
  um = TRUE,
  cex = 0.5
)



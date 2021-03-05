calculateVAF <- function(folder.path) {
  require(data.table)
  require(tidyverse)
  require(tictoc)
  tic()
  # find all xls files
  xls_paths <- list.files(folder.path, pattern = ".xls$", full.names = TRUE)
  cat(
    paste(length(xls_paths), "files found:"), 
    basename(xls_paths), 
    sep = "\n"
  )
  # import xls files
  xls_files <- lapply(xls_paths, fread)
  # change file names in the list
  names(xls_files) <- gsub("(.*)\\.GATK.*$", "\\1", basename(xls_paths))
  # change column 76 in each dataframe
  xls_files <- lapply(xls_files, function(x) {
    names(x)[76] <- "info_detail"
    return(x)
  })
  # bind datasets
  d_whole <- rbindlist(xls_files, idcol = "file")
  # filter by column Func
  d_whole <- d_whole[Func %in% c("exonic", "splicing", "exonic;splicing"), ]
  cat("\nKeep Func: exonic, splicing, exonic;splicing", "\n")
  cat("Remain ", d_whole[, .N], " rows", "\n")
  # filter by column ExonicFunc
  d_whole <- d_whole[ExonicFunc %in% c(
    ".", "frameshift insertion", "frameshift deletion", "missense SNV",
    "unknown", "stopgain", "stoploss"), ]
  cat("\nKeep ExonicFunc: ., frameshift insertion/deletion, missense SNV, unknown, stopgain, stoploss", "\n")
  cat("Remain ", d_whole[, .N], " rows", "\n")
  # elements in column FORMAT
  v_format_names <- strsplit(d_whole$FORMAT[1], ":")[[1]]
  # create new columns for each element in column FORMAT
  d_whole[, (v_format_names) := tstrsplit(info_detail, ":", fixed = TRUE)]
  # create new columns for each element in column AD
  d_whole[, c("ADREF", "ADALT", "ADALT2") := tstrsplit(AD, ",", fixed = TRUE)]
  # convert chr to num
  v_num_vars <- c("ADREF", "ADALT", "ADALT2", "DP")
  d_whole[, (v_num_vars) := lapply(.SD, as.numeric), .SDcols = v_num_vars]
  # calculate VAF
  d_whole[, VAF := 1 - ADREF / DP]
  # print unique GTs
  cat("\nUnique GTs: \n")
  print(d_whole[, .N, by = GT])
  cat("\n")
  # tictoc
  toc()
  cat("\n")
  return(d_whole)
}

folder_path <- "Y:/numerabilis/gongwei/SMZ1.WES.Novogene.Jan.2021/03.Variant_and_Annotation/03.Variant_and_Annotation/Mutation"
indel_path <- file.path(folder_path, "INDEL/Annotation")
snp_path <- file.path(folder_path, "SNP/Annotation")

d_indel <- calculateVAF(indel_path)
d_snp <- calculateVAF(snp_path)
# Combine d_indel and d_snp
d_indel_snp <- rbindlist(
  list(indel = d_indel, snp = d_snp), 
  idcol = "mutation"
)
# Subset columns
v_select_vars <- c(names(d_indel_snp)[1:12], "ExonicFunc", "FORMAT",
"info_detail", names(d_indel_snp)[92:100])
d_indel_snp_sub <- d_indel_snp[, .SD, .SDcols = v_select_vars]
fwrite(d_indel_snp_sub, "indel_snp_VAR_2021_3_5.csv")

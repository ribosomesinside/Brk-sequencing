library(data.table)

#sample list from asf of sequencing file
asf_samplefile <- "RN23010new.csv"

#relative fastq file directory as it appears on your machine
rel_fastq_dir <- "./RN23010"

#absolute fastq file directory as it appears on nemo (i.e. using pwd -P)
abs_fastq_dir <- "/nemo/stp/sequencing/inputs/instruments/fastq/230811_A01366_0434_BHHFLWDSX7/fastq/RN23010"

#samplesheet output csv file
outfile <- "./samplesheetRN23010.csv"

samples_dt <- fread(file = asf_samplefile,
                    skip = 1,
                    select = c(1, 2, 3),
                    col.names = c("sampleID", "sampleName", "date"),
                    header = TRUE)

#subset here if there are more samples on your submission document then in your sequencing results,
#I will usually subset by date.

fastq_files <- list.files(path = rel_fastq_dir,
                          pattern = "*.fastq.gz",
                          recursive = TRUE,
                          full.names = TRUE)

fastq_dt <- data.table(filepath = fastq_files)
fastq_dt[, filename := basename(filepath)]
fastq_dt[, c("sampleID", "sampleNumber", "lane", "read") := tstrsplit(filename, "_")[1:4]]

samplesheet_dt <- merge.data.table(samples_dt, fastq_dt,
                                   by = "sampleID")

samplesheet_dt[, NEWfilepath := gsub(pattern = rel_fastq_dir,
                                     replacement = abs_fastq_dir,
                                     x = filepath)]

dt <- samplesheet_dt[, .(fastq_1 = NEWfilepath[grepl("*_R1_001.fastq.gz", filename)],
                         fastq_2 = NEWfilepath[grepl("*_R2_001.fastq.gz", filename)]),
                     by = c("sampleNumber", "sampleID", "sampleName")]

dt[, sample := sampleName]
dt[, strandedness := "reverse"]

samplesheet_out <- dt[, c("sample", "fastq_1", "fastq_2", "strandedness")]

fwrite(samplesheet_out, file = "./samplesheetRN23010.csv")
#################################
species <- c("Sexi")

samplesheetR <- function(x) {
  
  date <- format(Sys.Date(), "%Y%m%d")
  
  dt <- samplesheet_dt[species == x]
  
  dt <- dt[, .(fastq_1 = NEWfilepath[grepl("*_R1_001.fastq.gz", filename)],
                             fastq_2 = NEWfilepath[grepl("*_R2_001.fastq.gz", filename)]),
                         by = c("sampleNumber", "sampleID", "sampleName")]
  
  dt[, sample := sampleName]
  dt[, strandedness := "reverse"]
  
  samplesheet_out <- dt[, c("sample", "fastq_1", "fastq_2", "strandedness")]
  
  return(samplesheet_out)
  
}

test <- samplesheetR(species[1])

date <- format(Sys.Date(), "%Y%m%d")

samplesheet_out <- samplesheetR(species[1])

fwrite(samplesheet_out, file = outfile)

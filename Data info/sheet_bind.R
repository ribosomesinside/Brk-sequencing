library(data.table)

sheet1 <-fread("samplesheetRN23010.csv")


sheetn<-sheet1[, sample:= paste0(sample, "_1")]

sheet2 <-fread("samplesheetRN23010_23.csv")

samplesheetRN23010_full<-rbind(sheetn, sheet2)

fwrite(samplesheetRN23010_full, file = "./samplesheetRN23010_full.csv")

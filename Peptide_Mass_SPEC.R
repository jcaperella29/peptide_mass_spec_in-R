library("doParallel")

library("MSnbase")


## downloads the data
library("rpx")
px1 <- PXDataset("PXD000001")
mzf <- pxget(px1, 7)

## reads the data
ms <- openMSfile(mzf)
hd <- header(ms)

## a set of spectra of interest: MS1 spectra eluted
## between 30 and 35 minutes retention time
ms1 <- which(hd$msLevel == 1)
rtsel <- hd$retentionTime[ms1] / 60 > 30 &
  hd$retentionTime[ms1] / 60 < 35

## the map
M <- MSmap(ms, ms1[rtsel], 521, 523, .005, hd, zeroIsNA = TRUE)


plot3D(M)

i <- ms1[which(rtsel)][1]
j <- ms1[which(rtsel)][2]
M2 <- MSmap(ms, i:j, 100, 1000, 1, hd)
plot3D(M2)


sp <- itraqdata[["X1"]]

plot(sp, reporters = iTRAQ4, full = TRUE)



# find path to a mzXML file
quantFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "mzXML$")
## find path to a mzIdentML file
identFile <- dir(system.file(package = "MSnbase", dir = "extdata"),
                 full.name = TRUE, pattern = "dummyiTRAQ.mzid")
## create basic MSnExp
msexp <- readMSData(quantFile, verbose = FALSE)
head(fData(msexp), n = 2)


msexp <- addIdentificationData(msexp, id = identFile)
head(fData(msexp), n = 2)


itraqdata2 <- pickPeaks(itraqdata, verbose=FALSE)
i <- 14
s <- as.character(fData(itraqdata2)[i, "PeptideSequence"])

plot(itraqdata2[[i]], s, main = s)

msexp <- removeNoId(msexp)


ms_seq <- fData(msexp)$sequence

itraqdata2_fet<-fData(itraqdata2)
itraqdata2_seq<-itraqdata2_fet$PeptideSequence

#comparing two sequences
unique_to_ms_seq<-setdiff(ms_seq,itraqdata2_seq)

unique_to_intraqdata2_seq<-setdiff(itraqdata2_seq,ms_seq)


#outputting text files of unique values


write.table(unique_to_ms_seq,"C:/Users/ccape/Downloads/unique_seq_1.txt")



write.table(unique_to_intraqdata2_seq,"C:/Users/ccape/Downloads/unique_seq_2.txt")

#original code came from https://www.bioconductor.org/packages/release/bioc/vignettes/MSnbase/inst/doc/v01-MSnbase-demo.html#11_Spectra_comparison. 
#I added lines to find unique/differentially expressed protiens  and output them as text-files.


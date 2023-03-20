library(DepecheR)
#Here, the neighbor smoothing for the histological inflammation scores for figure
#5, Kokkinou et al can be reproduced. 
metaData <- readRDS("Data/metaData.rds")
pcaData <- readRDS("Data/pcaData.rds")
umapData <- readRDS("Data/umapData.rds")

table(metaData$CHASTea, metaData$donor1, useNA = "always")
#       D1  D10  D11  D16  D19   D6 <NA>
#1    3776    0 1810    0 2662    0    0
#3    4766    0    0 7267    0    0    0
#7       0    0    0 3851    0 4731    0
#9       0    0    0    0 2131    0    0
#10      0    0    0    0    0 6799    0
#11      0 3153 1471    0    0    0    0
#12      0 2300    0    0    0    0    0
#<NA>    0    0    0    0    0    0    0

#It is hard to make a perfect analysis here, but we could possibly categorize
#1 and 3 as low, 7 as mid and 9-12 as high. In that case, we need to have as
#many from each of them to start with. We therefore need to sample in two steps. 
#First, each subcategoy needs to have equal numbers from each donor, so that
#the sampling procedure will pick as many from each. Then, from this selection, 
#we will pick the min from each group. 
roughChas <- sapply(as.character(metaData$CHASTea), switch, 
                    "1" = 1, "3" = 1, "7" = 2, "9" = 3, "10" = 3, "11" = 3, "12" = 3)
splitMetaData <- split(metaData, f = roughChas)
splitNRows <- split(seq_len(nrow(metaData)), f = roughChas)
roughChasSampling1 <- lapply(seq_along(splitMetaData),
                            function(x){
                              locMeta <- splitMetaData[[x]]
                              locDonInflam <- 
                                paste0(locMeta$donor1, locMeta$Tissue)
                              minNRows <- min(unlist(lapply(split(locMeta, locDonInflam), nrow)))
                              locRows <- splitNRows[[x]]
                              neighRows <- unlist(lapply(split(locRows, locDonInflam), 
                                                         sample, minNRows))
                            })

#Did we get any duplicates?
length(which(duplicated(unlist(roughChasSampling1))))
#0, nice. 


minSamplingGroup <- min(unlist(lapply(roughChasSampling1, length)))


neighRows <- unlist(lapply(roughChasSampling1, sample, minSamplingGroup))

#So with that, we are over to the analysis!
smoothCHAS <- neighSmooth(metaData$CHASTea, 
                         euclidSpaceData = pcaData,
                         neighRows = neighRows)
dir.create("Results")
contours <- dContours(umapData)

dColorPlot(smoothCHAS, xYData = umapData,
           colorScale = c("blue", "white", "red"), densCont = contours,
           plotName = "Results/smooth_CHASTea_bluered")

#We also create a filtered version
nDonorVector <- nUniqueNeighDons(metaData$donor, 
                                 euclidSpaceData = pcaData,
                                 neighRows = neighRows)

dColorPlot(nDonorVector, xYData = umapData, 
           plotName = "Results/nDonors_among_neighbors")

#Here, we check where it is suitable to put a gate on the inflammation vector. 
pdf("Results/Inflammation_histogram.pdf")
hist(smoothCHAS, breaks = 200)
abline(v = 4, col = "blue")
abline(v = 9, col = "red")
dev.off()

#OK, so we now want to cut out the data that has values above 9 and more than 6 donors
inMany <- which(nDonorVector > 5)

inflamInMany <- which(nDonorVector > 5 & smoothCHAS > 9)

dColorPlot(smoothCHAS[inflamInMany], controlData = c(0, 12), colorScale = c("blue", "white", "red"),
           xYData = umapData[inflamInMany,], densCont = contours,
           plotName = "Results/smooth_CHASTea_above_9_dens_above_5_dons_bluered", dotSize = 4)

#For this, we need to find a subset of cells that are present in many samples and
#Contains primarily non-inflammatory cells
nonInflamInMany <- which(nDonorVector > 5 & smoothCHAS < 4)
dColorPlot(smoothCHAS[nonInflamInMany], controlData = c(0, 12),colorScale = c("blue", "white", "red"),
           xYData = umapData[nonInflamInMany,], densCont = contours,
           plotName = "Results/smooth_CHASTea_below_4_dens_above_5_dons_bluered", dotSize = 3)

#Now, we add this information to the metadata
metaData$neighbor_smooth_CHASTea <- smoothCHAS
metaData$Inflam <- FALSE
metaData$Inflam[inflamInMany] <- TRUE
metaData$NonInflam <- FALSE
metaData$NonInflam[nonInflamInMany] <- TRUE
metaData$above5samplesInNeighbors <- FALSE
metaData$above5samplesInNeighbors[which(nDonorVector > 5)] <- TRUE

#Now, we export the new complete metadata daetaset
write.csv(metaData, "Complete_metadata.csv")

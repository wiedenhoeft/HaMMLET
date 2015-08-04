# extract simulated data from R's serialization object and write them to a file, the first column represents the true copy number

dir.create("breakpoint_detection_and_merging/")
setwd("breakpoint_detection_and_merging/")
download.file("http://www.cbs.dtu.dk/~hanni/aCGH/20chromosome.simulated.data.RData", "20chromosome.simulated.data.RData")
load("20chromosome.simulated.data.RData") 
for (sample in ls(simulated.data)){
	filename <- paste(sample,".csv", sep="")
	result <- simulated.data[[sample]]$sample[,c("copynumber", "log2ratios")]
	write.table(result, filename, sep="\t", row.names=FALSE, col.names=FALSE)
}

dir.create("../spatial_resolution_study/")
setwd("../spatial_resolution_study/")
download.file("http://www.cbs.dtu.dk/~hanni/aCGH/simulated.below5.data.RData", "simulated.below5.data.RData")
download.file("http://www.cbs.dtu.dk/~hanni/aCGH/simulated.5to10.data.RData", "simulated.5to10.data.RData")
download.file("http://www.cbs.dtu.dk/~hanni/aCGH/simulated.10to20.data.RData", "simulated.10to20.data.RData")
download.file("http://www.cbs.dtu.dk/~hanni/aCGH/simulated.above20.data.RData", "simulated.above20.data.RData")
for (simType in c("below5", "5to10", "10to20", "above20")){
	for (s in 1:100){
		load(paste("simulated.", simType, ".data.RData", sep=""))
		sample = paste("sample", s, sep="")
		result = simulated.data[[sample]]$sample[,c("copynumber", "log2ratios")]
		filename = paste(simType, "_sample_", s, ".csv", sep="")
		write.table(result, filename, sep="\t", row.names=FALSE, col.names=FALSE)
		}
}

dir.create("../testing/")
setwd("../testing/")
download.file("http://www.cbs.dtu.dk/~hanni/aCGH/Heterogeneous.simulated.data.RData", "Heterogeneous.simulated.data.RData")
load("Heterogeneous.simulated.data.RData")
for (dataset in 1:500){
	for (sample in 1:20){
		copynumber = simulated.data[[paste("dataset", dataset, sep="")]]$copynumbers[,sample]
		log2ratios = simulated.data[[paste("dataset", dataset, sep="")]]$datamatrix[,sample]
		result = cbind(copynumber, log2ratios)
		filename = paste("dataset_", dataset, "_sample_", sample, ".csv", sep="")
		write.table(result, filename, sep="\t", row.names=FALSE, col.names=FALSE)
	}
}
inputFile<-system.file("data", "simulinfile", package = "stepwise")
breaks<-c(548, 735, 832)
WinHalfWidth<-30
permReps<-10
b<-phylpro(inputFile, breaks, WinHalfWidth, permReps)
summary.phylpro(b)

inputFile<- system.file("data", "simulinfile", package = "stepwise")
breaks<-c(531,616,810)
WinHalfWidth<-30
permReps<-10
a<-maxchi(inputFile, breaks, WinHalfWidth, permReps)
summary.maxchi(a)

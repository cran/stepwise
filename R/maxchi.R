"maxchi" <-
    function(input_file, breaks, winHalfWidth, permReps) {
    input_file<-c("Rinterface", input_file)
    out<- .Call("Rmaxchi", input_file = as.character(input_file), 
       breaks = as.integer(breaks), winHalfWidth = as.integer(WinHalfWidth), 
       permReps = as.integer(permReps), PACKAGE="MaxChi")
    if (length(out$winlocs) > 0) class(out) = "maxchi"
    out
}

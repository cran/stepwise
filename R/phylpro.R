"phylpro" <-
    function(input_file, breaks, winHalfWidth, permReps) {
    input_file<-c("Rinterface", input_file)
    out<-.Call("Rphylpro", input_file = as.character(input_file), 
      breaks = as.integer(breaks), winHalfWidth = as.integer(WinHalfWidth), 
      permReps = as.integer(permReps), PACKAGE="Phylpro")    
    if (length(out$winlocs) > 0) class(out)<-"phylpro"
    out
}

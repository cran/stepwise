"summary.phylpro" <-
    function(object, ...) {
    attach(object[c("corrs","winlocs","target.seqs","quants")]);
    sigSites <- 1
    siteWinlocs <- -1
    siteChisqs <- -1.0
    if (length(winlocs)<0) {
        cat("\nNo significant Minimum Correlation statistics found.\n\n")
        return()
    }
    siteChisqs[1] <- corrs[1];
    siteWinlocs[1] <- winlocs[1];
    pairs <- target.seqs[1]
    for (i in 2:length(winlocs)) {
        if(winlocs[i] != siteWinlocs[sigSites]) { # found a new site 
            sigSites <- sigSites+1;       
            siteChisqs[sigSites]=corrs[i];
            siteWinlocs[sigSites]=winlocs[i]; 
            pairs[sigSites] <- target.seqs[i]
        } 
        else {
            if(corrs[i] < siteChisqs[sigSites]) {
                siteChisqs[sigSites] <- corrs[i]
                pairs[sigSites] <- target.seqs[i]
            } 
            else if(corrs[i] == siteChisqs[sigSites]) {
            # concatenate this pair on end of current list of pairs 
                pairs[sigSites] <- paste(pairs[sigSites],"\n\t\t\t   ", target.seqs[i], sep="")
            }
        } # end of else
    } # end of for loop
    cat("--------------------------------------------------\n");
    cat("There were", sigSites,"site-specific Minimum Correlation statistics significant at the\n");
    cat("10 percent level (10th percentile = ", sprintf("%5.3lf", quants[1]), ", 5th percentile = ",sprintf("%5.3lf", quants[2]),"):\n\n", sep="");
    cat ("Number Location  MinCor   targets\n")
    for (i in 1:length(siteWinlocs)) {
        if(siteChisqs[i]<quants[1]) star<-"*" else star<-" "
        cat(sprintf("%6d  %7d  %5.3lf%s   %s\n",i, as.integer(siteWinlocs[i]), siteChisqs[i], star,  pairs[i]));
    }
    cat("--------------------------------------------------\n");
    cat("Notes - \"Location\" is the polymorphic site just before the proposed breakpoint.\n");
    cat("      - MinCor statistics significant at the 5 percent level indicated by a * \n\n");
    
    a <- list(siteWinlocs = siteWinlocs, siteChisqs = siteChisqs, pairs = pairs)
    detach(object[c("corrs","winlocs","target.seqs","quants")]);
    invisible(a)
}

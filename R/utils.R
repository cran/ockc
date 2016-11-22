#
#  Copyright (C) 2005-2008 Friedrich Leisch

MClapply <- function(X, FUN, multicore=TRUE, ...)
{
    if(inherits(multicore, "cluster"))
        parLapply(multicore, X, FUN)
    else if(multicore)
        mclapply(X, FUN, ...)
    else
        lapply(X, FUN, ...)
}


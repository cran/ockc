ockc <- function(x, k, family=kccaFamily("kmeans"), order=NULL,
                 control=NULL, save.data=FALSE, multicore=FALSE, ...)
{
    MYCALL <- match.call()
    control <- as(control, "flexclustControl")

    x <- as(x, "matrix")
    x <- family@preproc(x)
    n <- nrow(x)
    p <- ncol(x)
    h <- as.integer(control@min.size)

    k <- as.integer(k)
    if(any(k < 2)) stop("The number of clusters has to be at least 2.")

    x <- as.matrix(x)
    
    if(n < max(k*h)) stop("It is not possible to partition ", n, " objects into ",
                          k, " clusters of minimum size ", h, ".")

    ss <- as.integer(control@subsampling)
    if(ss < 1) stop("subsampling must be a positive integer")

    
    solutions	= list()    
    proximity 	= family@dist(x, x)
    # proximity 	= as.matrix(dist(x))
    
    ## Step 1: seriation
    if(is.null(order)) {
        rNSsuc <- requireNamespace("seriation", quietly=TRUE)
        if(rNSsuc){
            order = seriation::get_order(seriation::seriate(as.dist(proximity), ...))
        } else {
            stop("Argument 'order' or package 'seriation' necessary.")
        }
    }
    proximity = proximity[order, order]
    x <- x[order,]
    xm <- ModelEnvMatrix(designMatrix=x)

    ## Step 2: partitioning
    # if(control@subsampling > 1)
    # {
    #     ## <FIXME>
    #     subsample 	= sort(rep(1:ceiling(n/ss), ss))[1:n]
    #     sub_x 		= family@allcent(x, subsample)
    #     z <- ockc(sub_x, k=k, family=family, control=control, order=1:nrow(sub_x), ...)
    #     for(s in 1:length(k)) {
    #         ockc@solutions[[s]]@clustering = ockc@solutions[[s]]@clustering*control@subsampling
    #         ssd 	= 0
    #         start 	= 1
    #         for(cluster in 2:ockc@solutions[[s]]@k+1){
    #             end 	= ockc@solutions[[s]]@clustering[cluster-1]
    #             ssd 	= ssd + sum(proximity[start:end, start:end])/(2*length(start:end))
    #             start 	= end + 1
    #         }
    #         ockc@solutions[[s]]@ssd = ssd
    #     }
    #     ## </FIXME>
    # }
    # else {
        # function for calculation of cumulated distances
        triang.j <- function(j, di, proximity) {
          result <- diag((t(di) %*% proximity %*% di))/ seq(1, n-j+1)
          return(result)
        }

        # function to calculate optimal solutions by dynamic programming
        getMinSSDPart <- function(n, maxk, h, ssd) {
          
          # initialize dp
          # curRes is a list of lists, which saves in sublist k all possible
          # current values of the summed distances (cv), all possible total
          # length of all clusters in iteration k (l) and the corresponding lengths
          # of cluster k (p).
          curRes <- list()
          # vector of endpoints of first segment (also the length p of the first
          # cluster)
          l1 <- h:n
          # curRes[[1]] <- list(cv=ssd[1,l1], l=l1, p=l1)
          curRes[[1]] <- list(cv=ssd[[1]][l1], l=l1, p=l1)
          
          for(ck in 2:maxk){
            # number of solutions after last iteration.
            nlold <- length(curRes[[ck-1]]$l)
            # number of solutions in current iteration.
            nlnew <- n-ck*h+1
            # matrix for all possible new summed distances, initialized with
            # Inf for easier retrieval of best solutions.
            # a[i,j] is the summed distance for a total length of (ck*h)+i-1, with
            # cluster size of j+h-1 for the last cluster
            # nlold-h columns because for the last h total lengths of the last
            # iteration an additional cluster is not possible
            a <- matrix(Inf, nlnew, nlold-h)
            for(i in 1:(nlold-h)){
              nstart <- curRes[[ck-1]]$l[i] + 1
              a[i:nlnew,i] <- curRes[[ck-1]]$cv[i] +
                ssd[[nstart]][h:(n-nstart+1)]
                # ssd[nstart,(nstart+h-1):n]
            }
            # which column of a has lowest cv:
            wmina <- max.col(-a, ties.method="first")
            # calculat length of last cluster and total length after current
            # iteration
            p <- h:(n-(ck-1)*h) - wmina +1
            l <- (ck*h):n

            # remeber only the best solutions
            w <- sapply(1:nlnew, function(x) a[x,wmina[x]])

            curRes[[ck]] <- list(cv=w, l=l, p=p)
          }
          return(curRes)
        }

        di <- upper.tri(proximity, diag=TRUE)
        ssd <- MClapply(1:n, function(x) triang.j(x, di=di[x:n,x:n, drop=FALSE], proximity=proximity[x:n, x:n, drop=FALSE]), multicore=multicore)

        z <- new("ockc", order = as.integer(order))
	
        dpRes <- getMinSSDPart(as.integer(n), max(k), as.integer(h), ssd)
        for(s in 1:length(k))
        {
            # retrieve results
            y <- numeric(k[s])
            ind <- length(dpRes[[k[s]]]$l)
            y[k[s]] <- dpRes[[k[s]]]$l[ind]
            for(j in (k[s]-1):1){
              y[j] <- dpRes[[j+1]]$l[ind] - dpRes[[j+1]]$p[ind]
              ind <- which(dpRes[[j]]$l == y[j])
            }

            cluster <- rep(1:k[s], diff(c(0,y[1:k[s]])))
            centers <- family@allcent(x, cluster)

            maxdist <- avdist <- double(k[s])
            for(r in 1:k[s]){
                maxdist[r] <- max(proximity[cluster==r, cluster==r])
                avdist[r] <- mean(proximity[cluster==r, cluster==r])
            }
            clusinfo <- data.frame(size=as.integer(table(cluster)),
                                   max_dist=maxdist,
                                   av_dist=avdist)
                
            z@models[[s]] <- new("kccasimple",
                                 k=k[s],
                                 cluster=cluster,
                                 iter=as(1, "integer"),
                                 converged=TRUE,
                                 call=MYCALL,
                                 control=control,
                                 centers=centers,
                                 family=family,
                                 clusinfo=clusinfo
                                 )
            if(save.data) z@data <- xm
        }
    # }
    
    z@call 	= MYCALL
    z@order 	= as.integer(order)
    return(z)
}

bootockc <- function(x, k, nboot=100, order=NULL, correct=TRUE, seed=NULL,
                          multicore=TRUE, verbose=FALSE, ...)
{
    MYCALL <- match.call()
    
    # if(multicore) require(parallel)
    if(!is.null(seed)) set.seed(seed)

    if(is.null(order)){
        order <- ockc(x=x, k=k, verbose=verbose)@order
    }

    nk <- length(k)
    nx <- nrow(x)

    index1 <- matrix(integer(1), nrow=nx, ncol=nboot)
    index2 <- index1
    
    ## empirical experiments show parallization does not pay for this
    ## (sample is too fast)
    for(b in 1:nboot){
        index1[,b] <- sort(sample(1:nx, nx, replace=TRUE))
        index2[,b] <- sort(sample(1:nx, nx, replace=TRUE))
    }

    BFUN <- function(b){
        if(verbose &! multicore){
            if((b %% 100) == 0)
                cat("\n")
            if((b %% 10) == 0)
                cat(b, "")
        }

        s1 <- ockc(x[index1[,b],,drop=FALSE], k=k, verbose=FALSE, order=order[index1[,b]])
        s2 <- ockc(x[index2[,b],,drop=FALSE], k=k, verbose=FALSE, order=order[index2[,b]])

        clust1 <- clust2 <- matrix(integer(1), nrow=nx, ncol=nk)
        cent1 <- cent2 <- list()
        rand <- double(nk)
        
        for(l in 1:nk)
        {
            cl1 <- getModel(s1, l)
            cl2 <- getModel(s2, l)

            # clust1[,l] <- unlist(mapply(rep, 1:k[l], cl1@clusinfo$size))
            # clust2[,l] <- unlist(mapply(rep, 1:k[l], cl2@clusinfo$size))
            clust1[,l] <- clusters(cl1, newdata=x)
            clust2[,l] <- clusters(cl2, newdata=x)

            cent1[[l]] <- cl1@centers
            cent2[[l]] <- cl2@centers

            rand[l] <- randIndex(table(clust1[,l], clust2[,l]),
                              correct=correct)
        }
        list(cent1=cent1, cent2=cent2, clust1=clust1, clust2=clust2,
             rand=rand)
        
    }
    
    ## empirical experiments show parallization does not pay for the 
    ## following (element extraction from list is too fast)
    z <- MClapply(as.list(1:nboot), BFUN, multicore=multicore)

    clust1 <- unlist(lapply(z, function(x) x$clust1))
    clust2 <- unlist(lapply(z, function(x) x$clust2))
    dim(clust1) <- dim(clust2) <- c(nx, nk, nboot)

    cent1 <- cent2 <- list()
    for(l in 1:nk){
        cent1[[l]] <- unlist(lapply(z, function(x) x$cent1[[l]]))
        cent2[[l]] <- unlist(lapply(z, function(x) x$cent2[[l]]))
        dim(cent1[[l]]) <- dim(cent2[[l]]) <- c(k[l], ncol(x), nboot)
    }

    # rand <- t(matrix(unlist(lapply(z, function(x) x$rand)), ncol=nboot))
    # rand <- t(do.call(cbind, lapply(z, function(x) x$rand)))
    rand <- t(sapply(z, function(x) x$rand))
    colnames(rand) <- k
    
    if(verbose) cat("\n")

    new("bootFlexclust", k=as.integer(k), centers1=cent1, centers2=cent2,
        cluster1=clust1, cluster2=clust2, index1=index1, index2=index2,
        rand=rand, call=MYCALL)
}

setClass("ockc",
         contains = "stepFlexclust",
         representation(order = "integer"))

setMethod("show", "ockc",
function(object)
{
    cat("ockc object of family",
        sQuote(object@models[[1]]@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\n")
    
    z <- data.frame(distsum = sapply(object@models,
                                     function(x) info(x, "distsum")))
    
    print(z, na.string="")
})

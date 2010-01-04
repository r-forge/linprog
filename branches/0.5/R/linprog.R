solveLP <- function( cvec, bvec, Amat, maximum=FALSE, maxiter=1000,
                     zero=1e-10, lpSolve=FALSE, verbose = 0 )
{

   result <- list()  # list for results that will be returned
   result$status <- 0

   nVar <- length(cvec)  # number of variables
   nCon <- length(bvec)  # number of constraints

   if( !all.equal( dim( Amat ), c( nCon, nVar ) ) == TRUE ) {
      stop( paste( "Matrix A must have as many rows as constraints (=elements",
         "of vector b) and as many columns as variables (=elements of vector c).\n" ) )
   }

   ## Labels
   if( is.null(names(cvec))) {
      clab <- as.character(1:nVar)
   } else {
      clab <- names(cvec)
   }
   if( is.null(names(bvec))) {
      blab <- as.character(1:nCon)
   } else {
      blab <- names(bvec)
   }

   ## lpSolve
   if( lpSolve ) {
      library( lpSolve )
      if( maximum ) {
         direction <- "max"
      } else {
         direction <- "min"
      }
      lpres <- lp (direction = direction, cvec, Amat, rep("<=",length(bvec)), bvec )

      solution  <- lpres$solution
      objval    <- lpres$objval
      Tab       <- NULL
      iter1     <- NULL
      iter2     <- NULL
      allvar    <- NULL
      basvar    <- NULL
      con       <- NULL

   } else {
      ## Simplex algorithm
      iter1 <- 0
      iter2 <- 0

      if(maximum) cvec <- -cvec

      ## Slack Variables
      for(i in 1:nCon) clab <- c( clab, paste("S", blab[i] ) )
      cvec <- c( cvec, rep( 0, nCon ) )

      ## Tableau ( Basic Variables, Slack,Variables, P0, Z-C )
      Tab <- rbind( cbind( Amat, diag( 1, nCon, nCon ), bvec ),
                 c( cvec, 0 ) )
      rownames(Tab) <- c( blab, "Z-C" )
      colnames(Tab) <- c( clab, "P0" )
      if( verbose >= 3 ) {
         print("initial Tableau")
         print(Tab)
      }

      ## searching for feasible solution for starting
      # basis: Zero Solution ( Basic Variables = Slack Variables )
      basis <- c( (nVar+1) : (nVar+nCon) )
      if(min(Tab[ 1:nCon, nVar+nCon+1]) < 0 ) {
         Tab2 <- Tab
         Tab2 <- rbind( Tab2, rep(0, ncol(Tab2) ) )
         rownames(Tab2)[nCon+2] <- "M Z-C"     # additional artificial 'Z-C' row
         basis2 <- basis
         nArt   <- 0       # number of artificial variables
         for(i in 1:nCon) {
            if(Tab[ i, nVar+nCon+1] < 0 ) {
               Tab2[ i, ] <- -Tab2[ i, ]
               Tab2 <- cbind( Tab2[ , 1:(nVar+nCon+nArt) ], rep(0,nCon+2),
                              Tab2[ , (nVar+nCon+nArt+1)] )
               nArt <- nArt + 1
               colnames(Tab2)[ nVar+nCon+nArt ] <- paste("M", rownames(Tab2)[i] )
               Tab2[ i, nVar+nCon+nArt ] <- 1
               Tab2[ nCon+2, nVar+nCon+nArt ] <- 1
               # put artificial variables in basis
               rownames(Tab2)[ i ] <- paste("M", rownames(Tab2)[i] )
               basis2[i] <- nVar+nCon+nArt
            }
         }
         for(i in 1:nCon) {    # artificial objective function (Z-C)
            if(Tab[ i, nVar+nCon+1] < 0 ) {
               Tab2[nCon+2, 1:(nVar+nCon+nArt)] <- Tab2[nCon+2, 1:(nVar+nCon+nArt)] -
                                                   Tab2[ i , 1:(nVar+nCon+nArt)]
            }
         }
         for(i in 1:nCon) {    # value of artificial objective function
            Tab2[nCon+2, nVar+nCon+nArt+1 ] <- Tab2[nCon+2, nVar+nCon+nArt+1 ] -
                         Tab2[i, nVar+nCon+nArt+1] * Tab2[nCon+2, basis[i] ]
         }
         colnames(Tab2)[ nVar+nCon+nArt+1 ] <- "P0"
         if( verbose >= 3 ) {
            print("initial Tableau for Phase 1")
            print(Tab2)
         }

         ## Simplex algorithm (Phase 1)
         while( min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ] ) < 0 & iter1 < maxiter) {
            iter1 <- iter1 + 1
            ## Pivot
            Tab[ abs(Tab) < zero ] <- 0
            pcolumn <- which.min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ]) # Pivot column
            rwerte  <- Tab2[ 1:nCon , nVar+nCon+nArt+1 ] / Tab2[ 1:nCon , pcolumn ] # R-values
            rwerte[ Tab2[1:nCon, pcolumn ] <= 0 ] <- max(rwerte)+1
            prow  <- which.min( rwerte )    # Pivot row
            if( verbose >=2 ) {
               cat( paste( "\nPivot Column:", as.character(pcolumn),
                           "(",colnames(Tab2)[pcolumn],")\n" ) )
               cat( paste( "Pivot Row:", as.character( prow ),
                  "(", rownames(Tab2)[prow], ")\n\n" ) )
            }

            ## New Basis
            basis[prow] <- pcolumn
            rownames(Tab2)[prow] <- colnames(Tab2)[pcolumn]

            ## New Tableau
            Tab2[ prow, ] <- Tab2[ prow, ] / Tab2[ prow, pcolumn ]
            for( i in 1:(nCon+2) ) {
               if( i != prow ) {
                  Tab2[ i, ] <- Tab2[ i, ] - Tab2[ prow, ] * Tab2[ i, pcolumn ]
               }
            }
            if( verbose >= 4 ) print(Tab2)
         }
         if(iter1 >= maxiter ) warning("Simplex algorithm (phase 1) did reach optimum.")
         Tab <- cbind( Tab2[ 1:(nCon+1), 1:(nCon+nVar) ], Tab2[ 1:(nCon+1), nVar+nCon+nArt+1 ] )
         if( verbose >= 3 ) {
            print("New starting Tableau for Phase II")
            print(Tab)
         }
      }
      ## Simplex algorithm (Phase 2)
      while( min( Tab[ nCon+1, 1:(nVar+nCon) ] ) < 0  & iter2 < maxiter) {
         iter2 <- iter2 + 1
         ## Pivot
         Tab[ abs(Tab) < zero ] <- 0
         pcolumn <- which.min( Tab[ nCon+1, 1:(nVar+nCon) ]) # Pivot column
         rwerte  <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon , pcolumn ]     # R-values
         rwerte[ Tab[1:nCon, pcolumn ] <= 0 ] <- max(rwerte)+1
         prow  <- which.min( rwerte )    # Pivot row
         if( verbose >= 2 ) {
            cat( paste( "\nPivot Column:", as.character(pcolumn),
                        "(",colnames(Tab)[pcolumn],")\n" ) )
            cat( paste( "Pivot Row:", as.character( prow ) ,
               "(",rownames(Tab)[prow],")\n\n") )
         }

         ## New Basis
         basis[prow] <- pcolumn
         rownames(Tab)[prow] <- colnames(Tab)[pcolumn]

         ## New Tableau
         Tab[ prow, ] <- Tab[ prow, ] / Tab[ prow, pcolumn ]
         for( i in 1:(nCon+1) ) {
            if( i != prow ) {
               Tab[ i, ] <- Tab[ i, ] - Tab[ prow, ] * Tab[ i, pcolumn ]
            }
         }
         if( verbose >= 4 ) print(Tab)
      }
      if(iter2 >= maxiter ) warning("Simplex algorithm (phase 2) did reach optimum.")
      ## Results: Basic Variables
      basvar <- matrix( NA, nCon, 1 )
      colnames(basvar) <- c("opt")
      rownames(basvar) <- rep("a",nCon)
      for( i in 1:nCon ) {
         rownames(basvar)[i] <- clab[sort(basis)[i]]
         basvar[i,1] <- Tab[ which(basis==sort(basis)[i]), nVar+nCon+1 ]
      }

      ## Results: All Variables (Including Slack Variables)
      allvar <- matrix( NA, nVar+nCon, 6 )
      colnames(allvar) <- c("opt", "c", "min c","max c","marg.", "marg.reg." )
      rownames(allvar) <- clab
      for( i in 1:(nVar+nCon) ) {
         if(i %in% basis ) {
            allvar[i,1] <- Tab[ which(basis==i), nVar+nCon+1 ]
            ## Stability of Basic Variables
            quot <- Tab[ nCon+1, 1:(nVar+nCon) ] / Tab[ which(basis==i), 1:(nVar+nCon) ]
            if(maximum) {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  op <- options()
                  options(warn=-1)
                  allvar[i,3] <- -cvec[ i ] - min(quot[quot>0 & !is.na(quot)])
                  options(op)
               }
               if(min(quot[!is.na(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ) {
                     allvar[i,4] <- -cvec[ i ] - max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar[i,4] <- Inf
                  }
               } else {
                  allvar[i,4] <- Inf
               }
            } else {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  allvar[i,4] <- cvec[ i ] + min(quot[quot>0 & !is.na(quot)])
               }
               if(min(quot[!is.na(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ){
                     allvar[i,3] <- cvec[ i ] + max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar[i,3] <- -Inf
                  }
               } else {
                  allvar[i,3] <- -Inf
               }
            }
         } else {
            allvar[i,1] <- 0
            if( i <= nVar ) {
               if(maximum) {
                  allvar[i,3] <- -Inf
                  allvar[i,4] <- Tab[ nCon+1, i ] - cvec[i]
               } else {
                  allvar[i,3] <- -Tab[ nCon+1, i ] + cvec[i]
                  allvar[i,4] <- Inf
               }
            }
         }
         if(maximum) {
            allvar[i,2] <- -cvec[i]
         } else {
            allvar[i,2] <- cvec[i]
         }
         # marginal contribution to objective function (Shadow prices)
         if( !( ( i %in% basis ) & ( i <= nVar ) ) ) {
            allvar[ i, 5] <- Tab[ nCon+1, i ] * (-1)^maximum
            quot <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon, i ]
            op <- options()
            options(warn=-1)
            if(min(quot[quot>0 & !is.na(quot)]) != Inf ) {
               allvar[i,6] <- min(quot[quot>0 & !is.na(quot)])
            } else {
               allvar[i,6] <- NA
            }
            options(op)
         }
      }
      allvar[ (allvar[,3] >  1e16), 3 ] <-  Inf
      allvar[ (allvar[,3] < -1e16), 3 ] <- -Inf

      ## Results: Constraints
      con <- matrix( NA, nCon, 5 )
      colnames(con) <- c("max", "actual", "diff", "dual price","dual.reg")
      rownames(con) <- blab
      con[ , 1 ] <- bvec
      for(i in 1: nCon) {
         if( (i+nVar) %in% basis ) {
            con[ i, 2 ] <- bvec[i] - Tab[ which((i+nVar)==basis), nVar+nCon+1 ]
         } else {
            con[ i, 2 ] <- bvec[i]
         }
         con[ i, 4 ] <- -allvar[ i+nVar , 5 ]
         con[ i, 5 ] <-  allvar[ i+nVar , 6 ]
      }
      con[ , 3 ] <- con[ , 1 ] - con[ , 2 ]

      objval   <- -Tab[ nCon+1, nCon+nVar+1 ] * (-1)^maximum
      basvar   <- round( basvar, digits=10 )
      con      <- round( con, digits=10 )
      allvar   <- round( allvar, digits=10 )
      solution <- allvar[ 1 : nVar, 1 ]
   }

   ## List of Results
   result$opt      <- round( objval, digits=10 )
   result$iter1    <- ifelse( iter1 < maxiter, iter1, -iter1 )
   result$iter2    <- ifelse( iter2 < maxiter, iter2, -iter2 )
   result$solution <- solution
   result$basvar   <- basvar
   result$con      <- con
   result$allvar   <- allvar
   result$maximum  <- maximum
   result$Tab      <- Tab
   result$lpSolve  <- lpSolve
   class(result)   <- "solveLP"
   result
}

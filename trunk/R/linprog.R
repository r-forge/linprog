solveLP <- function( cvec, bvec, Amat, maximum=FALSE,
               const.dir=ifelse( rep( maximum, length(bvec) ), rep("<=",length(bvec)),
                 rep(">=",length(bvec))),
               maxiter=1000, maxiterLpSolve=20, zero=1e-9,
               tol=1e-6, dualtol = tol, lpSolve=FALSE, solve.dual=FALSE, verbose=0 )
{

   result <- list()  # list for results that will be returned
   result$status <- 0

   rdigits <- -round( log10( zero ) )

   nVar <- length(cvec)  # number of variables
   nCon <- length(bvec)  # number of constraints

   if( !all.equal( dim( Amat ), c( nCon, nVar ) ) == TRUE ) {
      stop( paste( "Matrix A must have as many rows as constraints (=elements",
         "of vector b) and as many columns as variables (=elements of vector c).\n" ) )
   }
   if( length( const.dir ) != nCon ) {
      stop( paste( "'const.dir' must have as the elements as constraints",
         "(=elements of vector b).\n" ) )
   }
   if( sum( const.dir == ">=" | const.dir == ">" | const.dir == "=" | const.dir == "<=" |
            const.dir == "<" ) < nCon ) {
       stop( "'const.dir' may only contain '>=', '>', '=', '<=' or '<'.\n" )
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
   const.dir2 <- rep( 0, nCon )
   const.dir2[ const.dir == ">=" | const.dir == ">" ] <-  1
   const.dir2[ const.dir == "<=" | const.dir == "<" ] <- -1

   ## lpSolve
   if( lpSolve ) {
      library( lpSolve )
      if( maximum ) {
         direction <- "max"
      } else {
         direction <- "min"
      }
      lpres <- lp( direction = direction, cvec, Amat, const.dir, bvec )
      if( lpres$status == 0 ) {
         if( min( lpres$solution ) < -tol ) {
            lpres$status <- 7
         } else if( max( round( bvec - c( Amat %*% lpres$solution ),
            digits=rdigits ) * ( -1 )^maximum   ) > tol ) {
            lpres$status <- 3
         }
      }

      result$lpStatus <- c( result$lpStatus, lpres$status )

      if( lpres$status != 0 ) {
         iter <- 0
         while( lpres$status != 0 && iter <= maxiterLpSolve ) {
            iter <- iter + 1
            if( round( iter / 2 ) != iter / 2 ) {
               colPerm <- order( cvec == 0 ) # FALSE sorts before TRUE
            } else {
               colPerm <- sample( 1:length( cvec ) )
            }
            if( iter >= 2 ) {
               rowPerm <- sample( 1:length( bvec ) )
            } else {
               rowPerm <- c( 1:length( bvec ) )
            }
            lpres <- lp( direction = direction, cvec[ colPerm ], Amat[ rowPerm, colPerm ],
               const.dir, bvec[ rowPerm ] )
            reColPerm <- apply( as.matrix(c( 1:length( colPerm ))), 1, match, colPerm )
            lpres$solution <- lpres$solution[ reColPerm ]
            if( lpres$status == 0 ) {
               if( min( lpres$solution ) < -tol ) {
                  lpres$status <- 7
               } else if( max( round( bvec - c( Amat %*% lpres$solution ),
                  digits=rdigits ) * ( -1 )^maximum ) > tol ) {
                  lpres$status <- 3
               }
            }
            result$lpStatus <- c( result$lpStatus, lpres$status )
         }
         if( lpres$status == 0 ) {
            warning( paste( "lpSolve returned solution after", as.character( iter ),
               "permutation(s) of columns/rows (Error status:",
               paste( as.character( result$lpStatus ), collapse=" " ), ")" ) )
         } else {
            result$status <- 1
         }
      }
      if( result$status == 0 ) {
         result$solution            <- lpres$solution
         names( result$solution )   <- clab
         result$opt                 <- lpres$objval

         result$con <- data.frame( actual=NA, dir=const.dir, bvec=bvec, free=NA )
         result$con$actual <- round( c( Amat %*% result$solution ), digits=rdigits )
         names( result$con$actual ) <- blab
         result$con$free   <- round( result$con$bvec - result$con$actual, digits=rdigits )
         result$con$free[ const.dir2 == 1 ] <- -result$con$free[ const.dir2 == 1 ]
         result$con$free[ const.dir2 == 0 ] <- -abs( result$con$free[ const.dir2 == 0 ] )
      }
      if( result$status == 0 && solve.dual ) {
         if( sum( const.dir2 == 0 ) > 0 ) {
            stop( paste("At the moment the dual problem can not be solved with",
                        "equality constraints" ) )
         }
         if( maximum ) {
            direction <- "min"
            const.dir.dual <- rep(">=",nVar)
         } else {
            direction <- "max"
            const.dir.dual <- rep("<=",nVar)
         }
         dualres <- lp( direction, bvec * const.dir2 * (-1)^maximum,
            t( Amat * const.dir2 ) * (-1)^maximum, const.dir.dual, cvec )
         if( dualres$status == 0 ) {
            if( min( dualres$solution ) < -dualtol ) {
               dualres$status <- 7
            } else if( max( round( cvec - c( ( t( Amat * const.dir2 ) *
               (-1)^maximum ) %*% dualres$solution ), digits=rdigits ) *
               ( -1 )^(!maximum) ) > dualtol ) {
               dualres$status <- 3
            }
         }
         result$dualStatus <- c( result$dualStatus, dualres$status )

         if( dualres$status == 0 ) {
            result$con$dual <- dualres$solution
         } else {
            iter <- 0
            while( dualres$status != 0 && iter <= maxiterLpSolve ) {
               iter <- iter + 1
               if( round( iter / 2 ) != iter / 2 ) {
                  colPerm <- order( bvec == 0 ) # FALSE sorts before TRUE
               } else {
                  colPerm <- sample( 1:length( bvec ) )
               }
               if( iter >= 2 ) {
                  rowPerm <- sample( 1:length( cvec ) )
               } else {
                  rowPerm <- c( 1:length( cvec ) )
               }
               dualres <- lp( direction, ( bvec * const.dir2 * (-1)^maximum )[ colPerm ],
                  t( Amat * const.dir2 )[ rowPerm, colPerm ] * (-1)^maximum, const.dir.dual,
                  cvec[ rowPerm ] )
               if( dualres$status == 0 ) {
                  if( min( dualres$solution ) < -dualtol ) {
                     dualres$status <- 7
                  } else if( max( round( cvec[ rowPerm ] - c( ( t( Amat *
                     const.dir2 )[ rowPerm, colPerm ] * (-1)^maximum ) %*%
                     dualres$solution ), digits=rdigits ) *
                     ( -1 )^(!maximum) ) > dualtol ) {
                     dualres$status <- 3
                  }
               }
               reColPerm <- apply( as.matrix(c( 1:length( colPerm ))), 1, match, colPerm )
               dualres$solution <- dualres$solution[ reColPerm ]

               result$dualStatus <- c( result$dualStatus, dualres$status )
            }
            if( dualres$status == 0 ) {
               warning( paste( "lpSolve for the dual problem returned solution after",
                  as.character( iter ), "permutation(s) of columns/rows (Error status:",
                  paste( as.character( result$dualStatus ), collapse=" " ), ")" ) )
               # reperm <- apply( c( 1:length( perm ) ), 1, match, perm )
               result$con$dual <- dualres$solution
            } else {
               result$status <- 2
            }
         }
      }
   } else {
      ## Simples algorithm
      iter1 <- 0
      iter2 <- 0

      ## Slack Variables
      for(i in 1:nCon) clab <- c( clab, paste("S", blab[i] ) )
      cvec2 <- c( cvec, rep( 0, nCon ) )

      ## Tableau ( Basic Variables, Slack,Variables, P0, Z-C )
      Tab <- rbind( cbind( -Amat * const.dir2, diag( 1, nCon, nCon ), -bvec * const.dir2 ),
                 c( cvec2 * (-1)^maximum, 0 ) )
      rownames(Tab) <- c( blab, "Z-C" )
      colnames(Tab) <- c( clab, "P0" )
      if( verbose >= 3 ) {
         print("initial Tableau")
         print(Tab)
      }

      ## searching for feasible solution for starting
      # basis: Zero Solution ( Basic Variables = Slack Variables )
      basis <- c( (nVar+1) : (nVar+nCon) )
      if( min(Tab[ 1:nCon, nVar+nCon+1]) < 0 ) {
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
         while( min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ] ) < -zero & iter1 < maxiter) {
            iter1 <- iter1 + 1
           ## Pivot
            Tab[ abs(Tab) < zero ] <- 0
#            pcolumn <- which.min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ]) # Pivot column
            decval <- array( NA, nVar+nCon )
            for( pcolumnt in 1:(nVar+nCon+nArt) ) {
               if( Tab2[ nCon+2, pcolumnt ] < 0 ) {
                  rwerte  <- Tab2[ 1:nCon, nVar+nCon+nArt+1 ] / Tab2[ 1:nCon , pcolumnt ]
                       # R-values
                  rwerte[ Tab2[1:nCon, pcolumnt ] <= 0 ] <- max(rwerte,na.rm=TRUE)+1
                  prow  <- which.min( rwerte )    # Pivot row
                  if( length( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] ) >= 1 ) {
                     decval[ pcolumnt ] <- Tab2[ nCon+2, pcolumnt ] *
                           min( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] )
                  }
               }
            }
            if( min( decval, na.rm=TRUE ) < -zero ) {
               pcolumn <- which.min( decval ) # Pivot column
            } else {
               pcolumn <- which.min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ]) # Pivot column
            }
            rwerte  <- Tab2[ 1:nCon , nVar+nCon+nArt+1 ] / Tab2[ 1:nCon , pcolumn ] # R-values
            rwerte[ Tab2[1:nCon, pcolumn ] <= 0 ] <- max(rwerte, na.rm=TRUE)+1
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

         if(iter1 >= maxiter ) warning("Simplex algorithm (phase 1) did not reach optimum.")
         Tab <- cbind( Tab2[ 1:(nCon+1), 1:(nCon+nVar) ],
            Tab2[ 1:(nCon+1), nVar+nCon+nArt+1 ] )
         if( verbose >= 3 ) {
            print("New starting Tableau for Phase II")
            print(Tab)
         }
      }
      ## Simplex algorithm (Phase 2)
      while( min( Tab[ nCon+1, 1:(nVar+nCon) ] ) < -zero  & iter2 < maxiter ) {
         iter2 <- iter2 + 1
         ## Pivot
         Tab[ abs(Tab) < zero ] <- 0
#         pcolumn <- which.min( Tab[ nCon+1, 1:(nVar+nCon) ]) # Pivot column
         decval <- array( NA, nVar+nCon )
         for( pcolumnt in 1:(nVar+nCon) ) {
            if( Tab[ nCon+1, pcolumnt ] < 0 ) {
               rwerte  <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon , pcolumnt ]
                    # R-values
               rwerte[ Tab[1:nCon, pcolumnt ] <= 0 ] <- max(rwerte,na.rm=TRUE)+1
               prow  <- which.min( rwerte )    # Pivot row
               if( length( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] ) >= 1 ) {
                  decval[ pcolumnt ] <- Tab[ nCon+1, pcolumnt ] *
                        min( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] )
               }
            }
         }
         if( min( decval, na.rm=TRUE ) < -zero ) {
            pcolumn <- which.min( decval ) # Pivot column
         } else {
            pcolumn <- which.min( Tab[ nCon+1, 1:(nVar+nCon) ]) # Pivot column
         }
         rwerte  <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon , pcolumn ]     # R-values
         rwerte[ Tab[1:nCon, pcolumn ] <= 0 ] <- max(rwerte,na.rm=TRUE)+1
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
      if(iter2 >= maxiter ) warning("Simplex algorithm (phase 2) did not reach optimum.")
      ## Results: Basic Variables
      basvar <- matrix( NA, nCon, 1 )
      colnames(basvar) <- c("opt")
      rownames(basvar) <- rep("a",nCon)
      for( i in 1:nCon ) {
         rownames(basvar)[i] <- clab[sort(basis)[i]]
         basvar[i,1] <- Tab[ which(basis==sort(basis)[i]), nVar+nCon+1 ]
      }

      ## Results: All Variables (Including Slack Variables)
      allvar <- data.frame( opt=rep( NA, nVar+nCon ), cvec=cvec2, min.c=NA,
                              max.c=NA, marg=NA, marg.reg=NA )
      rownames(allvar) <- clab
      for( i in 1:(nVar+nCon) ) {
         if(i %in% basis ) {
            allvar$opt[ i ] <- Tab[ which(basis==i), nVar+nCon+1 ]
            ## Stability of Basic Variables
            quot <- Tab[ nCon+1, 1:(nVar+nCon) ] / Tab[ which(basis==i), 1:(nVar+nCon) ]
            if( maximum ) {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  op <- options()
                  options(warn=-1)
                  allvar$min.c[ i ] <- cvec2[ i ] - min(quot[quot>0 & !is.na(quot)])
                  options(op)
               }
               if(min(quot[!is.na(quot) & is.finite(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ) {
                     allvar$max.c[ i ] <- cvec2[ i ] - max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar$max.c[ i ] <- Inf
                  }
               } else {
                  allvar$max.c[ i ] <- Inf
               }
            } else {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  op <- options()
                  options(warn=-1)
                  allvar$max.c[ i ] <- cvec2[ i ] + min(quot[quot>0 & !is.na(quot)])
                  options(op)
               }
               if(min(quot[!is.na(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ){
                     allvar$min.c[ i ] <- -cvec2[ i ] + max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar$min.c[ i ] <- NA
                  }
               } else {
                  allvar$min.c[ i ] <- NA
               }
            }
         } else {
             allvar$opt[ i ] <- 0
             if( i <= nVar ) {
                if( maximum ) {
                   allvar$min.c[ i ] <- -Inf
                   allvar$max.c[ i ] <- Tab[ nCon+1, i ] + cvec2[i]
                } else {
                   allvar$min.c[ i ] <- 99#-Tab[ nCon+1, i ] - cvec2[i]
                   allvar$max.c[ i ] <- 77#Inf
                }
             }
         }
         allvar$cvec[ i ] <- cvec2[ i ]

         # marginal contribution to objective function (Shadow prices)
         if( !( ( i %in% basis ) & ( i <= nVar ) ) ) {
            allvar$marg[ i ] <- Tab[ nCon+1, i ] * (-1)^maximum
            if( !( i %in% basis ) & ( i > nVar ) ) {
               if( maximum ) {
                  allvar$max.c[ i ] <- Tab[ nCon+1, i ] #* (-1)^maximum
                  allvar$min.c[ i ] <- -Inf
               } else {
                  allvar$min.c[ i ] <- -Tab[ nCon+1, i ]
                  allvar$max.c[ i ] <-  Inf
               }
            }
            quot <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon, i ]
            op <- options()
            options(warn=-1)
            if( !( i %in% basis) ) {
               allvar$marg.reg[ i ] <- min(quot[quot>0 & !is.na(quot)])
            } else {
               allvar$marg.reg[ i ] <- NA
            }
            options(op)
         }
      }
      allvar$min.c[ allvar$min.c >  1e16 ] <-  Inf
      allvar$min.c[ allvar$min.c < -1e16 ] <- -Inf

      ## Results: Constraints
      con <- data.frame( actual=NA, dir=const.dir, bvec=bvec, free=NA, dual=NA, dual.reg=NA )
      names( con$actual ) <- blab
      for(i in 1: nCon) {
         if( (i+nVar) %in% basis ) {
            con$actual[ i ] <- round( bvec[i] + Tab[ which((i+nVar)==basis),
                                 nVar+nCon+1 ] * const.dir2[ i ], digits=rdigits )
         } else {
            con$actual[ i ] <- round( bvec[i], digits=rdigits )
         }
         if( -allvar$opt[ i+nVar ] == 0 ) {
            con$dual[ i ]     <- allvar$marg[ i+nVar ] * (-1)^maximum
            con$dual.reg[ i ] <- allvar$marg.reg[ i+nVar ]
         } else {
            con$dual[ i ]     <- 0
            con$dual.reg[ i ] <- allvar$opt[ i+nVar ]
         }
      }
      con$free <- round( con$bvec - con$actual, digits=rdigits )
      con$free[ const.dir2 == 1 ] <- -con$free[ const.dir2 == 1 ]
      con$free[ const.dir2 == 0 ] <- -abs( con$free[ const.dir2 == 0 ] )

      if( solve.dual ) {
         if( sum( const.dir2 == 0 ) > 0 ) {
            stop( paste( "At the moment the dual problem can not be solved",
               "with equality constraints" ) )
         }
         con$dual.p <- con$dual
         dualres <- solveLP( bevec, cvec, t(Amat), rep(">=",length(cvec)), cvec )
         con$dual <- dualres$solution
      }

      result$opt      <- round( -Tab[ nCon+1, nCon+nVar+1 ], digits=rdigits ) * (-1)^maximum
      result$iter1    <- iter1
      result$iter2    <- iter2
      result$allvar   <- round( allvar, digits=rdigits )
      result$basvar   <- round( basvar, digits=rdigits )
      result$solution <- allvar$opt[ 1 : nVar ]
      names( result$solution ) <- clab[ 1: nVar ]
      result$con      <- con
      if( verbose >= 1 ) result$Tab <- Tab
      if( iter1 >= maxiter ) result$status <- 4
      if( iter2 >= maxiter ) result$status <- 5
   }

   if( result$status == 0 ) {
      if( min ( result$con$free ) < - tol ) {
         result$status <- 3
      }
   }

   ## List of Results
   result$maximum    <- maximum
   result$lpSolve    <- lpSolve
   result$solve.dual <- solve.dual
   result$maxiter    <- maxiter
   class(result)     <- "solveLP"
   result
}

## print the (summary) results
summary.solveLP <- function(object,...) {
   cat("\n\nResults of Linear Programming / Linear Optimization\n")

   cat("\nObjective function")
   if( object$maximum ) {
      cat(" (Maximum): ")
   } else {
      cat(" (Minimum): ")
   }
   cat( object$opt, "\n" )

   cat("\nSolution\n")
   object$solution <- as.matrix(object$solution)
   colnames( object$solution ) <- c("opt")
   print( object$solution )
   cat("\n")
}

## print the results
print.solveLP <- function( x, digits=6,... ) {
   object <- x

   save.digits <- unlist(options(digits=digits))
   on.exit(options(digits=save.digits))

   cat("\n\nResults of Linear Programming / Linear Optimization\n")
   if( object$lpSolve ) cat("(using lpSolve)\n")

   if( object$status %in% c( 0, 4, 5 ) ) {
      cat("\nObjective function")
      if( object$maximum ) {
         cat(" (Maximum): ")
      } else {
         cat(" (Minimum): ")
      }
      cat( object$opt, "\n" )

      if( !is.null( object$iter1 ) ) {
         cat("\nIterations in phase 1: ")
         cat( object$iter1 )
         if( object$iter1 >= object$maxiter ) {
            cat(" (equals 'maxiter' !!!)")
         }
         cat("\nIterations in phase 2: ")
         cat( object$iter2 )
         if( object$iter2 >= object$maxiter ) {
            cat(" (equals 'maxiter' !!!)")
         }
      }
      cat("\nSolution\n")
      object$solution <- as.matrix(object$solution)
      colnames( object$solution ) <- c("opt")
      print( object$solution )

      if( !is.null( object$basvar ) ) {
         cat("\nBasic Variables\n")
         print( object$basvar )
      }

      cat("\nConstraints\n")
      print( object$con )

      if( !is.null( object$allvar ) ) {
         cat("\nAll Variables (including slack variables)\n")
         print( object$allvar )
      }
   } else if( object$status == 1 ) {
      if( sum( object$lpStatus == 2 ) == length( object$lpStatus ) ) {
         cat( "lpSolve did not find a feasible solution (Error: status 2)",
            " after ", length( object$lpStatus ) - 1, " permutation(s) of rows" )
      } else {
         cat( "lpSolve returned only error codes: ", object$lpStatus )
      }
   }
   if( object$status == 2 ) {
      if( sum( object$dualStatus == 2 ) == length( object$lpStatus ) ) {
         cat( "lpSolve did not find a feasible solution for the dual problem",
            "(Error: status 2) after ", length( object$lpStatus ) - 1,
            " permutation(s) of rows" )
      } else {
         cat( "lpSolve for the dual problem did not succeed. It returned",
            "following error codes: ", object$lpStatus )
      }
   } else if( object$status == 3 ) {
      print( result$con[ 1: 4 ] )
      cat( "The Constraints are violated. This is most likely due to rounding errors" )
   } else if( object$status == 4 ) {
      cat( "Simplex algorithm Phase 1 did not succed" )
   } else if( object$status == 5 ) {
      cat( "Simplex algorithm Phase 2 did not succed" )
   }

   cat("\n")
}

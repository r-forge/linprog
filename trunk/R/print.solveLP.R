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
      print( object$con[ 1: 4 ] )
      cat( "The Constraints are violated. This is most likely due to rounding errors" )
   } else if( object$status == 4 ) {
      cat( "Simplex algorithm Phase 1 did not succed" )
   } else if( object$status == 5 ) {
      cat( "Simplex algorithm Phase 2 did not succed" )
   }

   cat("\n")
   invisible( x )
}

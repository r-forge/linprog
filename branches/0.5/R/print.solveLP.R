## print the results
print.solveLP <- function( x, digits=6,... ) {
   object <- x

   save.digits <- unlist(options(digits=digits))
   on.exit(options(digits=save.digits))

   cat("\n\nResults of Linear Programming / Linear Optimization\n")
   cat("\nObjective function")
   if( object$maximum ) {
      cat(" (Maximum):\n")
   } else {
      cat(" (Minimum):\n")
   }
   print( object$opt )
   cat("\nIterations in phase 1: ")
   cat(abs(object$iter1))
   if( !object$lpSolve ) {
      if( object$iter1 < 0 ) cat(" (equals 'maxiter' !!!)")
   }
   cat("\nIterations in phase 2: ")
   cat(abs(object$iter2))
   if( !object$lpSolve ) {
      if( object$iter2 < 0 ) cat(" (equals 'maxiter' !!!)")
   }
   cat("\n\nBasic Variables\n")
   print( object$basvar )
   cat("\nConstraints\n")
   print( object$con )
   cat("\nAll Variables (including slack variables)\n")
   print( object$allvar )
   cat("\n")
}

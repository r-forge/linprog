library( linprog )

## Example 1
## Steinhauser, Langbehn and Peters (1992)
cvec <- c(1800, 600, 600)  # gross margins
names(cvec) <- c("Cows","Bulls","Pigs")
bvec <- c(40, 90, 2500)  # endowment
names(bvec) <- c("Land","Stable","Labor")
Amat <- rbind( c(  0.7,   0.35,   0 ),
               c(  1.5,   1,      3 ),
               c( 50,    12.5,   20 ) )
result1a <- solveLP( cvec, bvec, Amat, TRUE, verbose = 1 )
print( result1a )
# print summary results
summary( result1a )
# print all elements of the returned object
print.default( result1a )

# estimation with verbose = TRUE
result1b <- solveLP( cvec, bvec, Amat, TRUE, verbose = 4 )
all.equal( result1a, result1b )

# estimation with lpSolve
result1c <- solveLP( cvec, bvec, Amat, TRUE, lpSolve = TRUE, verbose = 4 )
print( result1c )
# print summary results
summary( result1c )
# print all elements of the returned object
print.default( result1c )

# using argument const.dir
const.dir <- c( ">=", ">=", ">=" )
result1d <- solveLP( cvec, -bvec, -Amat, maximum = TRUE, verbose = 1,
   const.dir = const.dir )
print( result1d )
all.equal( result1a[-8], result1d[-8] )

# using argument const.dir and lpSolve
result1e <-solveLP( cvec, -bvec, -Amat, maximum = TRUE, verbose = 1,
   const.dir = const.dir, lpSolve = TRUE )
print( result1e )
all.equal( result1c[-4], result1e[-4] )


## Example 2
## example 1.1.3 of Witte, Deppe and Born (1975)
cvec <- c(2.5, 2 )  # prices of feed
names(cvec) <- c("Feed1","Feed2")
bvec <- c( -10, -1.5, 12)
names(bvec) <- c("Protein","Fat","Fibre")
Amat <- rbind( c(-1.6,-2.4 ),
               c(-0.5,-0.2 ),
               c( 2.0, 2.0 ) )
result2a <- solveLP( cvec, bvec, Amat, verbose = 1 )
print( result2a )
# print summary results
summary( result2a )
# print all elements of the returned object
print.default( result2a )

# estimation with verbose = TRUE
result2b <- solveLP( cvec, bvec, Amat, verbose = 4 )
all.equal( result1a, result1b )

# estimation with lpSolve
result2c <- solveLP( cvec, bvec, Amat, lpSolve = TRUE, verbose = 4 )
print( result2c )
# print summary results
summary( result2c )
# print all elements of the returned object
print.default( result2c )

# using argument const.dir
const.dir <- c( ">=", ">=", "<=" )
result2d <- solveLP( cvec, abs( bvec ), abs( Amat ), verbose = 1,
   const.dir = const.dir )
print( result2d )
all.equal( result2a[-8], result2d[-8] )

# using argument const.dir and lpSolve
result2e <- solveLP( cvec, abs( bvec ), abs( Amat ), verbose = 1,
   const.dir = const.dir, lpSolve = TRUE )
print( result2e )
all.equal( result2c[-4], result2e[-4] )



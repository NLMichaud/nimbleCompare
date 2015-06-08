# nimbleCompare
A comparison of NIMBLE run times using Daniel's linear Gaussian model with the correlated parametrization.  

*  nimCompare.Rmd has all of the code needed to compare NIMBLE to the 'SMC' R package, POMP,  LibBi, and Biips.
*  The pomp, Biips, LibBi, and smc folders have code which sets up the model in those specific packages
*  The filters folder has the necessary particle filters

The "pomp" and "SMC" packages should be installed.

The comparison to LibBi requires [LibBi](http://libbi.org/getting-started.html) to be installed, and uses  the "Rbi" package which can be found [here](https://github.com/libbi/RBi), although it may not install correctly on the most recent version of R.

The comparison to Biips requires the RBiips package which can be found [here](https://alea.bordeaux.inria.fr/biips/doku.php?id=Download).

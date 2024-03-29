this is the Modern Fortran (Fortran 2003/2008) implementation of the Laplace-Transform Analytic Element Method (LT-AEM) for computing the head and velocity due to superposition of elements consisting of combinations of circles, points, ellipses, and line segments. The method also has several particle tracking routines implemented.

This software has been described in the following publications
----------------------------------------------------------
1) Bakker, M. and K.L. Kuhlman, 2011. "Computational issues and applications of line-elements to model subsurface flow governed by the modified Helmholtz equation", Advances in Water Resources. http://dx.doi.org/10.1016/j.advwatres.2011.02.008

2) Kuhlman, K.L. and S.P. Neuman, 2009. "Laplace-transform analytic-element method for transient porous-media flow", Journal of Engineering Mathematics, 64(2), 113-130. http://dx.doi.org/10.1007/s10665-008-9251-1

3) Kuhlman, K.L., 2008. "Laplace transform analytic element method for transient groundwater flow simulation", PhD Dissertation, University of Arizona. http://sabio.arizona.edu/record=b5975534~S9

4) Kuhlman, K.L., and S.P. Neuman, 2006. "Recent Advances in LT-AEM Theory and Application to Transient Groundwater Flow", Copenhagen, Denmark, Conference on Computational Methods in Water Resources XVI http://dx.doi.org/10.4122/1.1000000586

Dependencies and Compilation Instructions
============================================================
The software was primarliy developed in Linux using freely available Fortran compilers.  LAPACK is used for matrix operations (e.g., QR algorithm), which can be downloaded freely from NETLIB at: http://www.netlib.org/lapack/

LT-AEM and LAPACK itself depend on (and the LAPACK source provides) BLAS for lower-level matrix operations.  System blas and laplack libraries should be available with -lblas and -llapack or static libraries for blas and lapack should be linked or copied to the source directory as blas.a and lapack.a, or the Makefile needs to be changed to point to them.

A not-so recent version of gfortran is required (>= 5.0), because of the use of the ieee_arithmetic intrinsic Fortran2003 module.

arb (http://arblib.org) is an optional dependency. Specifying ARB=1 as an argument to make will comple wrapper code to use the arbitrary precision arb library for K and I bessel functions, instead of the Amos fortran code. The arb library should work for any argument (very large or small) and will always give double-precision accurate results, but the library is much slower than amos (double precision only). Arb must be a recent version, since the scaled bessel functions aren't in all versions of the library. Arb then requires flint, mpfr, and mpc.

With the two (or optionally three) libraries in place, compilation should work with just the simple "make" command from the source directory.  This copies the compiled executable into the "run" sub directory, where there are many input files (unfortunately of various vintages -- some might not work with the latest version of the code).


Author Contact Info
============================================================ 
If you need assistance or have questions, contact me (klkuhlm _at_ sandia _dot_ gov) regarding the code. My website includes more reference and publications: http://kris.kuhlmans.net  I may be able to provide you with an executable for your platform, if you have issues compiling it.  If you run into a bug or missing feature, I might be able to accommodate you, but this depends on the amount of work required to fix/implement the feature of concern.

# manipulatecube

A Fortran tool used by me and my group to work with volume data files in the 
'cube' format.

(c) Jochen Autschbach, 2022

Files dlasrt.f, lsame.f, xerbla.f come from LAPACK, with their own license
(BSD) and copyright. See below for further details.

See https://github.com/jautschbach/mathematica-notebooks for some Mathematica
notebooks that can be used to visualize this kind of volume data.

If you have GNU Fortran (gfortran) and the `make` utility, simply type
`make` on the command line to build the executable, after downloading this code.

There is a sorting routine (DLASRT) used in the code that comes from Lapack, so 
you need to have Lapack libraries installed. Needed at link time. Alternatively,
the relevant Fortran source files (dlasrt.f, lsame.f, xerbla.f) are included in 
this directory. Simply omment/uncomment the relevant lines in the Makefile as 
indicated, to compile the Lapack routines along with the manipulatecube sources. Lapack routines were downloaded from https://netlib.org on 2022-09-14. The License is 
https://netlib.org/lapack/LICENSE.txt and included here in file LICENSE.lapack

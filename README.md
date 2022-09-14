# manipulatecube

A Fortran tool used by me and my group to work with volume data files in the 
'cube' format.

See https://github.com/jautschbach/mathematica-notebooks for some Mathematica
notebooks that can be used to visualize this kind of volume data.

If you have GNU Fortran (gfortran) and the `make` utility, simply type
`make` on the command line to build the executable, after downloading this code.

There is a sorting routine (DLASRT) used in the code that comes from lapack, so you need to have lapack libraries installed. Needed at link time. Alternatively, download the source code for DLASRT + dependencies here: https://netlib.org/lapack/explore-html/df/ddf/dlasrt_8f.html Then add the Fortran source files (dlasrt.f, lsame.f, xerbla.f) to this directory, and comment/uncomment the relevant lines in the Makefile as indicated. 

# manipulatecube

A Fortran tool used by me and my group to work with volume data files in the 
'cube' format

See https://github.com/jautschbach/mathematica-notebooks for some Mathematica
notebooks to visualize this kind of volume data

If you have GNU Fortran (gfortran) and the `make` utility, simply type
`make` on the command line to build the executable.

There is a sorting routine used in the code that comes from lapack, so you need to have lapack libraries installed. Needed at compile time. 

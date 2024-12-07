# manipulatecube

A Fortran tool used by me and my group to work with volume data files in the 
'cube' format.

(c) Jochen Autschbach, 2022

Files dlasrt.f, lsame.f, and xerbla.f come from LAPACK, with their own license
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

The output of `manipulatecube` is also in cube format, either `results.cube`
or a pair of cube files (`results1.cube` and `results2.cube`), depending
on the command. Suppose you have input cube files `inp1.cube` and `inp2.cube`.
Then, `manipulatecube` can perform the following operations:

1. Multiply cube by a constant factor:
```
 manipulatecube inp1.cube fac <factor>
```		
Note that the output of this command also prints a numerical integration
of the resulting cube data, which is useful when one wants to find an isosurface
that corresponds to a prescribed fraction of the total orbital probability.
This is described in the README file in the `example` directory.

2. Fix a cube file if the grid vectors are not in the order x, y, z or
if it goes in a negative direction along one of the coordinates:
```
 manipulatecube inp1.cube fix
```
Some volume data visualizers have trouble if the grid is not properly
ordered or has negative steps. If the input grid is not rectangular, 
`manipulatecube` will unfortunately not be able to fix it. 

3. Add or subtract the data in two cube files:
```
 manipulatecube inp1.cube add|sub inp2.cube
```
Evidently, the grid must be the same for the two cubes, and they should be
for the same molecule. `manipulatecube` does some rudimentary checks and rejects
any command involving two cubes if the cube file headers appear to be
different.

4. Multiply the data from two cubes with each other point-by-point:
```
 manipulatecube inp1.cube mul inp2.cube
```
Above, the same file name may be used for `inp1.cube` and `inp2.cube`, in which 
case the result is the square of the original cube data. This is useful
to create volume data for orbital densities. 

5. Take a linear combination of the orbitals respresented by two cubes:
```
 manipulatecube inp1.cube mix inp2.cube
```
Here we have two resulting cubes. Assuming that the two input cubes are
for two different orbitals *f* and *g*, say,the results are cube files for the orbitals (*f*+*g*)/sqrt(2) and 
(*f*-*g*)/sqrt(2), respectively. If the input orbitals are orthonormal, then
the output cubes will also correspond to orthonormal orbitals. 

6. A generalization of the `mix` command is
```
 manipulatecube inp1.cube rot inp2.cube <angle>
```
This will perform a 'rotation' (in orbital space) of the two input orbitals. The
output cubes correspond to
 *c* `inp1.cube` + *s* `inp2.cube` and -*s* `inp1.cube` + *c* `inp2.cube`, respectively, where *c* = cos(angle), *s* = sin(angle), and the angle is in degrees. Apart from a negative phase for the second orbital, the results for the 
angle being 45 degrees is equivalent to the command under item 5. Note that
the `rot` command does not correspond to a rotation of an orbital. Rather,
the command is equivalent to applying a rotation to a pair of vectors in the
space spanned by said vectors.  


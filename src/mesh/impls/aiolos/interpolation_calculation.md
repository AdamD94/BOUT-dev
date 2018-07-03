Notes on the interp_to calculation
==================================


In the case of a equi-spaced grid in interpolation direction, the
stencils can e.g. by copied from the main code. To also support
forwards and backwards stencils, and to be less error prone, the
stencil values are calculated in get_interp_vals(), defined in
stencils.py

This includes an inversion of a 4x4 matrix. In the more general case
of an abitrary, non-uniform mesh, the stencils are potentially
different for each point of the mesh.
Inverting the matrix for each point in x and y, seems however rather
complicated.
Thus the matrix is inverted, using maxima.
The calculation done in maxima are:

```
A : matrix([1,1,1,1],[a,b,c,d],[2*a**2,2*b**2,2*c**2,2*d**2],[6*a**3,6*b**3,6*c**3,6*d**3]);
Ai : invert(A);
C: matrix(
 [1], 
 [0], 
 [0], 
 [0]
);
r: Ai . C;
simplify(r);

B : matrix([1,1,1,1],[a,b,-b,c],[2*a**2,2*b**2,2*b**2,2*c**2],[6*a**3,6*b**3,-6*b**3,6*c**3]);
Bi : invert(B);
r2 : Bi . C;
simplify(r2);

factor(ratsimp(r));
factor(ratsimp(r2));
```

The first case (A, Ai, r) is the more general case.
`a` ... `d` are the distances of the input points with respect to the
output point.
The second case is assuming that the output point is in the center
between the second and third point.
The computational cost for the first calculation is 6*4
multiplication, 3*4 substractions and 4 fractions.
The second case has only 4 multiplications and 2 substractions less,
thus implementing the more general case should be sufficient.

The output of the above commands is stored in interpolation_calculation.html

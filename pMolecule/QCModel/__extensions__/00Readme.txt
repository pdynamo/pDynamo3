----
CXC:
----

Version of cxc taken from libxc-1.2.0.

Some automatically generated files are required (funcs_XX.h and funcs_YY.c).
Also some other minor modifications for compilation to work (including changing
the names of all funcs_*.c files to funcs_*.i).

---------------------------------
Gaussian Basis Integral Notation:
---------------------------------

BxEyNz = x bases, y electrons, z nuclei (or point charges).

The presence of x bases implies neither x different bases nor centers.

---------
Problems:
---------

Maximum angular momentum for a Gaussian basis is "G". However, it looks like there are problems for integrals between G functions,
in which case they should be avoided for the moment (e.g. in QZVP basis).

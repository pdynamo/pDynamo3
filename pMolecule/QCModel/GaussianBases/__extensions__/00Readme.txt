---------------------------------
Gaussian Basis Integral Notation:
---------------------------------

The following notation is used for the different methods:

        Basis, Operator, Basis, Properties

Basis                   f<n>, g<n>, h<n>        a different symbol implies a different basis
                                                the number of basis functions, n, is required
                        m<n>, n<n>              points with charges
                        p<n>, q<n>              points without charges

Operator                A                       anti-Coulomb
                        C                       Coulomb
                        D                       dipole
                        O                       overlap (= identity)
                        Q                       quadrupole

                        X                       a generic operator - only used in file names as a catch-all

Properties              i                       raw integrals
                        r<n>                    raw derivative integrals with respect to Cartesian coordinates of basis centers
                                                the order of the derivative, n, is required
                                                multiple orders are possible

                        E                       an energy
                        I and R<n>              the raw integrals contracted with the appropriate densities, for properties (e.g.
                                                dipoles) and center coordinate derivatives.
                        V                       potentials

Notes:

  - Lower and upper case property symbols correspond to raw and contracted or full quantities, respectively.
  - Multiple operators are possible if the functions calculate them automatically or as alternatives.
  - Points have delta function or very peaked Gaussian s-function distributions.

Examples:

        f1Cg2R1         electron-fit derivatives
        f1Op1i          the values of the basis functions at a set of grid points
        f2Cf2i          two-electron integrals
        m1Cn1ER1        the Coulomb energy and coordinate first derivatives between two sets of point charges

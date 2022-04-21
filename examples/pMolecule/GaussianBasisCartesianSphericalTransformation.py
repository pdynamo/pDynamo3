"""Test Gaussian basis C->S (X) and S->C (Y) transformations."""

import math, os, os.path

from Definitions                     import yamlPath
from pCore                           import logFile             , \
                                            TestScriptExit_Fail , \
                                            YAMLUnpickle
from pMolecule.QCModel.GaussianBases import GaussianBasis
from pScientific.Arrays              import Array               , \
                                            ArrayPrint2D

#===================================================================================================================================
# . Options.
#===================================================================================================================================
# . Maximum and minimum angular momenta (S P D F G H I).
_LMax = 6
_LMin = 0

# . Tolerances.
_Precision = 1.0e-10

#===================================================================================================================================
# . Transformation matrices.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Read in reference data.
xReference = YAMLUnpickle ( os.path.join ( yamlPath, "C2STransformations_0to6.yaml" ) )

# . Generate the transformations and print them.
failures = 0
for l in range ( _LMin, _LMax + 1 ):

    # . Get the transformations and print them.
    X           = GaussianBasis.CartesianToSphericalTransformation ( l, l    )
    Y           = GaussianBasis.SphericalToCartesianTransformation ( l, l, X )
    itemsPerRow = max ( 2*l+1, 5 )
    ArrayPrint2D ( X, itemsPerRow = itemsPerRow, title = "L={:d} X (C->S)".format ( l ) )
    ArrayPrint2D ( Y, itemsPerRow = itemsPerRow, title = "L={:d} Y (S->C)".format ( l ) )

    # . Generate the identity (Y^T * X) and check it.
    I = Array.WithExtents ( X.columns, X.columns )
    I.MatrixMultiply ( Y, X, xTranspose = True )
    for i in range ( I.rows ): I[i,i] -= 1.0
    deviation = I.AbsoluteMaximum ( )
    if deviation > _Precision:
        failures += 1
        logFile.Paragraph ( "Excessive deviation ({:f}) from the identity for the product Y^T*X for L={:d}.".format ( deviation, l ) )
        ArrayPrint2D ( I, itemsPerRow = itemsPerRow, title = "Invalid Identity Y^T*X" )

    # . Compare X versus the reference values.
    X0        = xReference[l]
    deviation = 0.0
    for r in range ( X.rows ):
        for c in range ( X.columns ):
            deviation = max ( deviation, math.fabs ( X[r,c] - X0[r][c] ) )
    if deviation > _Precision:
        failures += 1
        logFile.Paragraph ( "Excessive deviation ({:f}) between calculated and reference X values for L={:d}.".format ( deviation, l ) )

# . Footer.
logFile.Footer ( )
if failures > 0: TestScriptExit_Fail ( )

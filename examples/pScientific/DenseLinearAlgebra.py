"""Tests for dense linear algebra."""

from pCore                     import CPUTime                 , \
                                      logFile                 , \
                                      TestScriptExit_Fail
from pScientific.Arrays        import Array                   , \
                                      StorageType
from pScientific.LinearAlgebra import EigenPairs              , \
                                      LinearEquations         , \
                                      LinearLeastSquaresBySVD , \
                                      MachineConstants
from pScientific.RandomNumbers import NormalDeviateGenerator  , \
                                      RandomNumberGenerator

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Matrix sizes.
_Extents = ( 1, 5, 10, 25, 100, 500, 1000 )

# . Precision.
_Tolerance = 1.0e-10

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def _CheckEigenPairs ( s, e, v ):
    n  = s.rows
    sF = _SymmetricToSquare ( s )
    t  = Array.WithExtents  ( n, n )
    t.MatrixMultiply  ( sF, v )
    sF.MatrixMultiply ( v, t, xTranspose = True )
    for i in range ( n ):
        sF[i,i] -= e[i]
    return sF.iterator.AbsoluteMaximum ( )

def _CheckSquareLinearEquations ( a, b, c ):
    a.VectorMultiply ( c, b, alpha = 1.0, beta = -1.0 )
    return b.iterator.AbsoluteMaximum ( )

def _CheckSymmetricLinearEquations ( a, b, c ):
    aF = _SymmetricToSquare ( a )
    return _CheckSquareLinearEquations ( aF, b, c )

def _EigenPairs ( extent, cpuTimer, ndg ):
    s = Array.WithExtent  ( extent, storageType = StorageType.Symmetric ) ; s.Set ( 0.0 )
    e = Array.WithExtent  ( extent         ) ; e.Set ( 0.0 )
    v = Array.WithExtents ( extent, extent ) ; v.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( i+1 ): s[i,j] = ndg.NextDeviate ( )
        s[i,i] *= 2.0
    tStart = cpuTimer.Current ( )
    try:
        results = EigenPairs ( s, e, v, preserveInput = True )
        error   = _CheckEigenPairs ( s, e, v )
        success = True
    except Exception as e:
        print ( e )
        error   = 0.0
        success = False
    return ( cpuTimer.Current ( ) - tStart, success, error )

def _SquareMatrixLinearEquations ( extent, cpuTimer, ndg ):
    a = Array.WithExtents ( extent, extent ) ; a.Set ( 0.0 )
    b = Array.WithExtent  ( extent         ) ; b.Set ( 0.0 )
    c = Array.WithExtent  ( extent         ) ; c.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( extent ):
            a[i,j] = ndg.NextDeviate ( )
        b[i] = ndg.NextDeviate ( )
    tStart = cpuTimer.Current ( )
    try:
        results = LinearEquations ( a, b, preserveInput = True, solution = c )
        error   = _CheckSquareLinearEquations ( a, b, c )
        success = True
    except Exception as e:
        print ( e )
        error   = 0.0
        success = False
    return ( cpuTimer.Current ( ) - tStart, success, error )

def _SquareMatrixLinearEquationsBySVD ( extent, cpuTimer, ndg ):
    a = Array.WithExtents ( extent, extent ) ; a.Set ( 0.0 )
    b = Array.WithExtent  ( extent         ) ; b.Set ( 0.0 )
    c = Array.WithExtent  ( extent         ) ; c.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( extent ):
            a[i,j] = ndg.NextDeviate ( )
        b[i] = ndg.NextDeviate ( )
    tStart = cpuTimer.Current ( )
    try:
        results = LinearLeastSquaresBySVD ( a, b, preserveInput = True, solution = c )
        error   = _CheckSquareLinearEquations ( a, b, c )
        success = True
    except Exception as e:
        print ( e )
        error   = 0.0
        success = False
    return ( cpuTimer.Current ( ) - tStart, success, error )

def _SymmetricMatrixLinearEquations ( extent, cpuTimer, ndg ):
    a = Array.WithExtent ( extent , storageType = StorageType.Symmetric ); a.Set ( 0.0 )
    b = Array.WithExtent ( extent ) ; b.Set ( 0.0 )
    c = Array.WithExtent ( extent ) ; c.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( i+1 ):
            a[i,j] = ndg.NextDeviate ( )
        b[i] = ndg.NextDeviate ( )
    tStart = cpuTimer.Current ( )
    try:
        results = LinearEquations ( a, b, preserveInput = True, solution = c )
        error   = _CheckSymmetricLinearEquations ( a, b, c )
        success = True
    except Exception as e:
        print ( e )
        error   = 0.0
        success = False
    return ( cpuTimer.Current ( ) - tStart, success, error )

def _SymmetricToSquare ( s ):
    n = s.rows
    f = Array.WithExtents ( n, n )
    f.Set ( -999999999999999.0 )
    for i in range ( n ):
        for j in range ( i ):
            f[i,j] = s[i,j]
            f[j,i] = s[i,j]
        f[i,i] = s[i,i]
    return f

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK     = True
cpuTimer = CPUTime ( )
rng      = RandomNumberGenerator.WithSeed ( 314159 )
ndg      = NormalDeviateGenerator.WithRandomNumberGenerator ( rng, mu = 0.0, sigma = 5.0 )

# . Print machine constants.
constants = MachineConstants ( )

# . Run the tests.
numberFailed     = 0
numberInaccurate = 0
for ( function, label ) in ( ( _EigenPairs                       , "EigenPairs"    ) ,
                             ( _SquareMatrixLinearEquations      , "Square LEs"    ) ,
                             ( _SquareMatrixLinearEquationsBySVD , "SVD LEs"       ) ,
                             ( _SymmetricMatrixLinearEquations   , "Symmetric LEs" ) ):
    results = [ function ( extent, cpuTimer, ndg ) for extent in _Extents ]
    table   = logFile.GetTable ( columns = [ 10, 20, 20, 20 ] )
    table.Start   ( )
    table.Title   ( label + " Results" )
    table.Heading ( "Extent" )
    table.Heading ( "Times", columnSpan = 2 )
    table.Heading ( "Error" )
    for ( extent, ( time, success, error ) ) in zip ( _Extents, results ):
        table.Entry ( "{:d}"  .format ( extent ) )
        table.Entry ( "{:.3f}".format ( time   ) )
        table.Entry ( CPUTime.TimeToString ( time ) )
        if success:
            if error >= _Tolerance:
                tag = "{:.2g} *".format ( error )
                numberInaccurate += 1
            else:
                tag = "{:.2g}  ".format ( error )
        else:
            tag = "Failed  "
            numberFailed += 1        
        table.Entry ( tag )
    table.Stop ( )

# . Footer.
logFile.Footer ( )
if ( numberFailed > 0 ) or ( numberInaccurate > 0 ): TestScriptExit_Fail ( )

"""Check the pseudo-inverse."""
"""Tests for dense linear algebra."""

from pCore                     import CPUTime                 , \
                                      logFile                 , \
                                      TestScriptExit_Fail
from pScientific.Arrays        import Array                   , \
                                      ArrayPrint2D            , \
                                      StorageType
from pScientific.LinearAlgebra import EigenPairs              , \
                                      LinearEquations         , \
                                      LinearLeastSquaresBySVD , \
                                      MachineConstants        , \
                                      MatrixPseudoInverse
from pScientific.RandomNumbers import NormalDeviateGenerator  , \
                                      RandomNumberGenerator

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Matrix shapes.
_Shapes = ( (    1,    1 ), (    5,   5 ) ,
            (   10,   10 ), (   25,  25 ) ,
            (  100,  100 ), (  500, 500 ) ,
            ( 1000, 1000 ),
            (    1,    5 ), (    5,   1 ) ,
            (    5,   10 ), (   10,   5 ) ,
            (   10,   25 ), (   25,  10 ) ,
            (   25,  100 ), (  100,  25 ) ,
            (  100,  500 ), (  500, 100 ) ,
            (  500, 1000 ), ( 1000, 500 ) )

# . Options.
_MaximumExtent = 1001
_Tolerance     = 1.0e-10

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

def _CheckPseudoInverse ( a, b ):
    m = a.rows
    n = a.columns
    if ( m == n ) or ( m < n ):
        c = Array.WithExtents ( m, m ) # . A * B.
        c.Set ( 0.0 )
        c.MatrixMultiply ( a, b )
        s = 0
        for i in range ( m ):
            if ( c[i,i] - 1.0 ) <= _Tolerance:
                c[i,i] = 0.0
                s += 1
        error = c.iterator.AbsoluteMaximum ( )
    else: error = 0.0
    if ( m == n ) or ( m > n ):
        c = Array.WithExtents ( n, n ) # . B * A.
        c.Set ( 0.0 )
        c.MatrixMultiply ( b, a )
        t = 0
        for i in range ( n ):
            if ( c[i,i] - 1.0 ) <= _Tolerance:
                c[i,i] = 0.0
                t += 1
        error = max ( error, c.iterator.AbsoluteMaximum ( ) )
    return error

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

def _PseudoInverse ( extent0, extent1, cpuTimer, ndg ):
    a = Array.WithExtents ( extent0, extent1 ) ; a.Set ( 0.0 )
    b = Array.WithExtents ( extent1, extent0 ) ; b.Set ( 0.0 )
    for i in range ( extent0 ):
        for j in range ( extent1 ): a[i,j] = ndg.NextDeviate ( )
    tStart = cpuTimer.Current ( )
    try:
        results = MatrixPseudoInverse ( a, b, preserveInput = True )
        error   = _CheckPseudoInverse ( a, b )
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
for ( function, label, isSquare ) in ( ( _EigenPairs                       , "EigenPairs"    , True  ) ,
                                       ( _PseudoInverse                    , "Pseudo-Inverse", False ) ,
                                       ( _SquareMatrixLinearEquations      , "Square LEs"    , True  ) ,
                                       ( _SquareMatrixLinearEquationsBySVD , "SVD LEs"       , True  ) ,
                                       ( _SymmetricMatrixLinearEquations   , "Symmetric LEs" , True  ) ):
    results = []
    for ( e0, e1 ) in _Shapes:
        if ( e0 > _MaximumExtent ) or ( e1 > _MaximumExtent ): continue
        if isSquare:
            if e0 == e1:
                results.append ( ( e0, e1, function ( e0, cpuTimer, ndg ) ) )
        else:
            results.append ( ( e0, e1, function ( e0, e1, cpuTimer, ndg ) ) )
    table   = logFile.GetTable ( columns = [ 4, 2, 4, 20, 20, 20 ] )
    table.Start   ( )
    table.Title   ( label + " Results" )
    table.Heading ( "Shape", columnSpan = 3 )
    table.Heading ( "Times", columnSpan = 2 )
    table.Heading ( "Error" )
    for ( extent0, extent1, ( time, success, error ) ) in results:
        table.Entry ( "{:d}"  .format ( extent0 ) )
        table.Entry ( "x" )
        table.Entry ( "{:d}"  .format ( extent1 ) )
        table.Entry ( "{:.3f}".format ( time    ) )
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

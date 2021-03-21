"""Linear algebra (matrix) benchmarks."""

from pCore                     import CPUTime                , \
                                      logFile
from pScientific.Arrays        import Array                  , \
                                      StorageType
from pScientific.LinearAlgebra import EigenPairs
from pScientific.RandomNumbers import NormalDeviateGenerator , \
                                      RandomNumberGenerator

# . Matrix sizes.
eExtents = ( 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000 )
mExtents = ( 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000 )

# . Functions.
def EigenValues ( extent, cpuTimer, ndg ):
    s = Array.WithExtent  ( extent, storageType = StorageType.Symmetric ) ; s.Set ( 0.0 )
    e = Array.WithExtent  ( extent         ) ; e.Set ( 0.0 )
    v = Array.WithExtents ( extent, extent ) ; v.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( i+1 ): s[i,j] = ndg.NextDeviate ( )
        s[i,i] *= 2.0
    tStart = cpuTimer.Current ( )
    EigenPairs ( s, e, v )
    return ( cpuTimer.Current ( ) - tStart )

def MatrixMultiply ( extent, cpuTimer, ndg ):
    a = Array.WithShape ( [ extent, extent ] ); a.Set ( 0.0 )
    b = Array.WithShape ( [ extent, extent ] ); b.Set ( 0.0 )
    c = Array.WithShape ( [ extent, extent ] ); c.Set ( 0.0 )
    for i in range ( extent ):
        for j in range ( extent ):
            a[i,j] = ndg.NextDeviate ( )
            b[i,j] = ndg.NextDeviate ( )
    tStart = cpuTimer.Current ( )
    c.MatrixMultiply ( a, b )
    return ( cpuTimer.Current ( ) - tStart )

# . Header.
logFile.Header ( )

# . Initialization.
cpuTimer = CPUTime ( )
rng      = RandomNumberGenerator.WithSeed ( 314159 )
ndg      = NormalDeviateGenerator.WithRandomNumberGenerator ( rng, mu = 0.0, sigma = 5.0 )

# . Calculation.
for ( extents, function, tag ) in ( ( eExtents, EigenValues   , "EigenValue"      ) ,
                                    ( mExtents, MatrixMultiply, "Matrix Multiply" ) ):
    times = [ function ( extent, cpuTimer, ndg ) for extent in extents ]
    table = logFile.GetTable ( columns = [ 10, 20, 20 ] )
    table.Start   ( )
    table.Title   ( tag + " Timings" )
    table.Heading ( "Extent" )
    table.Heading ( "Times", columnSpan = 2 )
    for ( extent, time ) in zip ( extents, times ):
        table.Entry ( "{:d}"  .format ( extent ) )
        table.Entry ( "{:.3f}".format ( time   ) )
        table.Entry ( CPUTime.TimeToString ( time ) )
    table.Stop ( )

# . Footer.
logFile.Footer ( )

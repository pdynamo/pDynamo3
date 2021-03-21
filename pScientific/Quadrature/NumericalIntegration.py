"""Functions for numerical integration."""

import heapq

from math import fabs

# . Simple low-precision integration algorithms for well-behaved functions.

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
# . Based on implementations from wikipedia page but with addition of a priority queue.
# . The method is OK but should not be used for functions which have zero-valued regions
# . as 
def AdaptiveSimpsonsRule ( F, a, b, tolerance, maximumIterations ):
    """Calculate integral of F from a to b."""
    # . Number of function evaluations is 3 + 2 * iterations.
    # . Initialization.
    h0         = fabs ( a - b )
    epsilon    = tolerance / h0
    integral   = 0.0
    intervals  = []
    iterations = 0
    state      = {}
    # . The first interval covers the whole range.
    p1    = 0.5 * ( a + b )
    f0    = F ( a  )
    f1    = F ( p1 )
    f2    = F ( b  )
    whole = ( f0 + 4.0 * f1 + f2 ) * h0 / 6.0
    heapq.heappush ( intervals, ( 0.0, h0, a, p1, b, f0, f1, f2, whole ) )
    # . Loop over intervals on the heap.
    while ( iterations < maximumIterations ) and ( len ( intervals ) > 0 ):
        ( error, h, p0, p2, p4, f0, f2, f4, whole ) = heapq.heappop ( intervals )
        p1       = 0.5 * ( p0 + p2 )
        p3       = 0.5 * ( p2 + p4 )
        f1       = F ( p1 )
        f3       = F ( p3 )
        left     = ( f0 + 4.0 * f1 + f2 ) * h / 12.0
        right    = ( f2 + 4.0 * f3 + f4 ) * h / 12.0
        lrw      = left + right - whole
        error    = fabs ( lrw )
        estimate = ( left + right + lrw / 15.0 )
        if error < 15.0 * h * epsilon:
            integral += estimate
        else:
            h *= 0.5
            heapq.heappush ( intervals, ( - error, h, p0, p1, p2, f0, f1, f2, left  ) )
            heapq.heappush ( intervals, ( - error, h, p2, p3, p4, f2, f3, f4, right ) )
        iterations += 1
    # . Accumulate unfinished parts of the integral.
    state["Iterations"] = iterations
    if len ( intervals ) > 0:
        state["Is Converged"] = False
        maximumError = 0.0
        for ( error, h, p0, p2, p4, f0, f2, f4, whole ) in intervals:
            integral    += whole
            maximumError = max ( maximumError, fabs ( error ) * h0 / ( 30.0 * h ) )
        state["Error"] = maximumError
    else:
        state["Is Converged"] = True
    # . Finish up.
    if a > b: integral *= -1.0
    return ( integral, state )

def TrapezoidalRule ( F, a, b, n ):
    """Calculate integral of F from a to b."""
    integral = 0.0
    if n > 0:
        h        = ( b - a ) / float ( n )
        integral = 0.5 * ( F ( a ) + F ( b ) )
        x        = a
        for i in range ( 1, n ):
            x        += h
            integral += F ( x )
        integral *= h
    return integral

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    from math import cos, sin, tan
    analytic           = 1.0 - cos ( 1.0 )
    ( numeric, state ) = AdaptiveSimpsonsRule ( sin, 0.0, 1.0, 1.0e-15, 1000 )
    print ( "Test Integration Results:" )
    print ( "Analytic   = {:20.10f}".format ( analytic ) )
    print ( "Numerical  = {:20.10f}".format ( numeric  ) )
    print ( "Difference = {:20.10f}".format ( fabs ( analytic - numeric ) ) )
    print ( "\nState: {!r}".format ( state ) )


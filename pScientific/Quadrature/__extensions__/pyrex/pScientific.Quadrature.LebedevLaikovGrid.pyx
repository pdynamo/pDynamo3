"""Functions for Lebedev-Laikov grids."""

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def LebedevLaikovGrid_GetGridPoints ( CInteger  gridAngularMomentum ):
    """Get grid points given an angular momentum."""
    cdef CInteger     numberOfPoints
    cdef Coordinates3 gridPoints
    cdef RealArray1D  weights
    gridPoints     = None
    weights        = None
    numberOfPoints = LebedevLaikov_Number_Of_Points ( gridAngularMomentum )
    if numberOfPoints > 0:
        gridPoints = Coordinates3.WithExtent ( numberOfPoints )
        weights    = RealArray1D.WithExtent  ( numberOfPoints )
        LebedevLaikov_GridPointsWeights ( numberOfPoints, gridPoints.cObject, weights.cObject )
    return ( gridPoints, weights )

"""ESP charge fitting."""

import math

from pCore                     import logFile       , \
                                      LogFileActive
from pScientific               import PeriodicTable , \
                                      Units
from pScientific.Arrays        import Array
from pScientific.Geometry3     import Coordinates3  , \
                                      Vector3 
from pScientific.LinearAlgebra import LinearLeastSquaresBySVD
from pScientific.Quadrature    import LebedevLaikovGrid_GetGridPoints

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Default values - RESP parameters.
_DefaultAResp = 0.0005 # . Weak constraints.
_DefaultBResp = 0.1

# . Default values - RESP iterator.
_DefaultFTolerance        = 1.0e-10
_DefaultLogFrequency      = 1
_DefaultMaximumIterations = 250
_DefaultQTolerance        = 1.0e-8

#===================================================================================================================================
# . Perform the charge fitting.
#===================================================================================================================================
def ESPChargeFitting ( system, aRESP = _DefaultAResp, bRESP = _DefaultBResp, doRESP = False, log = logFile ):
    """Perform an ESP charge fit."""

    # . Check for a system with a qcModel.
    if system.qcModel is None: raise ValueError ( "System does not have a QC model." )

    # . Initialization.
    natoms = len ( system.qcState.qcAtoms )
    ndim   = natoms + 1
    qtotal = system.qcModel.charge

    # . Get the grid points for the QC atoms only and convert to bohrs.
    gridPoints = GenerateVanDerWaalsSurface ( system, log = log, qcAtomsOnly = True, scalingFactors = [ 1.4, 1.6, 1.8, 2.0 ] )
    gridPoints.Scale ( Units.Length_Angstroms_To_Bohrs )

    # . Allocate space - one larger than necessary for fInteraction.
    fInteraction = Array.WithExtents ( ndim, gridPoints.rows )
    fInteraction.Set ( 0.0 )

    # . Get the observed potentials at the grid points.
    phi = system.qcModel.GridPointPotentials ( system.scratch, gridPoints )

    # . Get the interaction terms for each atom at the grid points.
    coordinates3 = system.scratch.qcCoordinates3AU
    GetInteractionTerms ( coordinates3, gridPoints, fInteraction )

    # . Get the A matrix and the B vector.
    A = Array.WithExtents ( ndim, ndim ) ; A.Set ( 0.0 )
    B = Array.WithExtent  ( ndim )       ; B.Set ( 0.0 )
    A.MatrixMultiply ( fInteraction, fInteraction, yTranspose = True )
    fInteraction.VectorMultiply ( phi, B )

    # . Add the total charge constraint terms.
    for i in range ( natoms ):
        A[i,ndim-1] = 1.0
        A[ndim-1,i] = 1.0
    B[ndim-1] = qtotal

    # . Get |phi^2|.
    phi2 = phi.DotSelf ( )

    # . Iterative solution.
    if doRESP:
        ( isConverged, sos, solution, condition, rank ) = RESPIterator ( natoms, A, B, fInteraction, phi, aRESP, bRESP, log = log )

    # . Non-iterative solution.
    else:

        # . Solve the equations.
        results  = LinearLeastSquaresBySVD ( A, B )
        solution = results["Solution"]

        # . Determine the sum of squares.
        fInteraction.VectorMultiply ( solution, phi, alpha = -1.0, beta = 1.0, transpose = True )
        sos = phi.DotSelf ( )

    # . Determine the error measure.
    error = math.sqrt ( sos / phi2 )

    # . Do some printing.
    if LogFileActive ( log ):
        items = [ ( "Charges"   , "{:d}"  .format ( natoms                      ) ) ,
                  ( "Rank"      , "{:d}"  .format ( results["Rank"]             ) ) ,
                  ( "Error"     , "{:.5g}".format ( error                       ) ) ,
                  ( "Condition" , "{:.5g}".format ( results["Condition Number"] ) ) ]
        if doRESP: items.append ( ( "Converged", repr ( isConverged ) ) )
        log.SummaryOfItems ( items, title = "ESP Charge Fitting Summary" )

    # . Return the charges.
    charges = Array.WithExtent ( natoms )
    for i in range ( natoms ):
        charges[i] = solution[i]
    return charges 

#===================================================================================================================================
# . Generate grid points on the van der Waals surface for a system.
# . Maybe could simplify this by scaling a single complete surface?
#===================================================================================================================================
def GenerateVanDerWaalsSurface ( system, gridAngularMomentum = 21, log = logFile, qcAtomsOnly = False, scalingFactors = [ 1.0 ] ):
    """Generate a superposition of van der Waals surfaces represented by grid points."""

    # . QC atoms only (including link atoms).
    if qcAtomsOnly:
        atomicNumbers = system.qcState.atomicNumbers
        coordinates3  = system.scratch.qcCoordinates3
    # . All atoms.
    else:
        atomicNumbers = [ atom.atomicNumber for atom in system.atoms ]
        coordinates3  = system.coordinates3

    # . Find the number of atoms.
    natoms = len ( atomicNumbers )

    # . Set radii.
    radii = Array.WithExtent ( natoms )
    for ( i, n ) in enumerate ( atomicNumbers ):
        radii[i] = PeriodicTable[n].vdwRadius

    # . Get the grid points for a single center.
    ( basicgrid, weights ) = LebedevLaikovGrid_GetGridPoints ( gridAngularMomentum )

    # . Allocate space for the possible number of grid points.
    npossible  = natoms * basicgrid.rows * len ( scalingFactors )
    gridPoints = Coordinates3.WithExtent ( npossible )
    gridPoints.Set ( 0.0 )

    # . Initialization.
    nfound      = 0
    atomgrid    = Coordinates3.WithExtent ( basicgrid.rows )
    translation = Vector3.Null ( )
    atomgrid.Set ( 0.0 )

    # . Loop over scaling factors.
    for factor in scalingFactors:

        # . Loop over points.
        for i in range ( natoms ):

            # . Get the radius.
            iradius = factor * radii[i]

            # . Get the translation.
            translation[0] = coordinates3[i,0]
            translation[1] = coordinates3[i,1]
            translation[2] = coordinates3[i,2]

            # . Get the scaled grid centered at the point.
            basicgrid.CopyTo   ( atomgrid    )
            atomgrid.Scale     ( iradius     )
            atomgrid.Translate ( translation )

            # . Remove points that are within the current scaled radii of other points.
            for p in range ( atomgrid.rows ):
                QOK = True
                x   = atomgrid[p,0]
                y   = atomgrid[p,1]
                z   = atomgrid[p,2]
                for j in range ( natoms ):
                    if j != i:
                        dx       = coordinates3[j,0] - x
                        dy       = coordinates3[j,1] - y
                        dz       = coordinates3[j,2] - z
                        jradius2 = ( factor * radii[j] )**2
                        if ( dx**2 + dy**2 + dz**2 ) <= jradius2:
                            QOK = False
                            break
                if QOK:
                    gridPoints[nfound,0] = x
                    gridPoints[nfound,1] = y
                    gridPoints[nfound,2] = z
                    nfound += 1

    # . Reduce the array size if necessary.
    if nfound < npossible:
        newpoints = Coordinates3.WithExtent ( nfound )
        for p in range ( nfound ):
            newpoints[p,0] = gridPoints[p,0]
            newpoints[p,1] = gridPoints[p,1]
            newpoints[p,2] = gridPoints[p,2]
        gridPoints = newpoints

    # . Do some printing.
    if LogFileActive ( log ):
        items = [ ( "Atoms"           , "{:d}".format ( natoms                 ) ) ,
                  ( "Surfaces"        , "{:d}".format ( len ( scalingFactors ) ) ) ,
                  ( "Found Points"    , "{:d}".format ( nfound                 ) ) ,
                  ( "Possible Points" , "{:d}".format ( npossible              ) ) ]
        log.SummaryOfItems ( items, title = "van der Waals Surface Generation Summary" )

    # . Finish up.
    return gridPoints

#===================================================================================================================================
# . Get the interaction terms between the model charges and grid points.
#===================================================================================================================================
def GetInteractionTerms ( coordinates3, gridPoints, fInteraction ):
    """Get interaction terms."""
    # . To start with use simple 1/r terms. Should be generalized to allow non-delta function (e.g. Gaussian) densities for atoms.
    # . In which case would require radii or Gaussian widths as well.
    fInteraction.Set ( 0.0 )
    for i in range ( coordinates3.rows ):
        x = coordinates3[i,0]
        y = coordinates3[i,1]
        z = coordinates3[i,2]
        for g in range ( gridPoints.rows ):
            dx = gridPoints[g,0] - x
            dy = gridPoints[g,1] - y
            dz = gridPoints[g,2] - z
            r  = math.sqrt ( dx**2 + dy**2 + dz**2 )
            if r != 0.0:
                fInteraction[i,g] = 1.0 / r

#===================================================================================================================================
# . Get the RESP constraints.
#===================================================================================================================================
def GetRESPConstraints ( natoms, aRESP, bRESP, charges, A ):
    """Add in RESP constraints to the diagonal elements of the A matrix."""
    for i in range ( natoms ):
        q = charges[i]
        A[i,i] += aRESP * q / math.sqrt ( q**2 + bRESP**2 )

#===================================================================================================================================
# . RESP iterator function.
#===================================================================================================================================
def RESPIterator ( natoms                                        ,
                   A                                             ,
                   rhs                                           ,
                   fInteraction                                  ,
                   phi                                           ,
                   aRESP                                         ,
                   bRESP                                         ,
                   ftolerance        = _DefaultFTolerance        ,
                   log               = logFile                   ,
                   logFrequency      = _DefaultLogFrequency      ,
                   maximumIterations = _DefaultMaximumIterations ,
                   qtolerance        = _DefaultQTolerance        ):

    """Solve the RESP equations by simple iteration."""

    # . Allocate space.
    n       = len ( rhs )
    Atemp   = Array.WithExtents ( n, n )
    Btemp   = Array.WithExtent  ( n )
    phitemp = Array.WithExtent  ( len ( phi ) )
    q       = Array.WithExtent  ( n, initializer = 0.0 ) # . With initial values.
    qold    = Array.WithExtent  ( n )

    # . Check for printing.
    doPrinting = ( logFrequency > 0 ) and ( logFrequency < maximumIterations ) and LogFileActive ( log )
    if doPrinting:
        table = log.GetTable ( columns = [ 10, 20, 20, 20, 20, 10 ] )
        table.Start ( )
        table.Title ( "RESP Solver" )
        table.Heading ( "Iteration"   )
        table.Heading ( "Function"    )
        table.Heading ( "Change in F" )
        table.Heading ( "Change in Q" )
        table.Heading ( "Condition"   )
        table.Heading ( "Rank"        )

    # . Determine the sum of squares with initial zero charges.
    f = phi.DotSelf ( )

    # . Initialization.
    niterations = 0
    isConverged  = False

    # . Perform the iterations.
    while ( niterations < maximumIterations ) and ( not isConverged ):

        # . Save old data.
        fold = f
        q.CopyTo ( qold )

        # . Fill new A and RHS.
        A.CopyTo   ( Atemp )
        rhs.CopyTo ( Btemp )

        # . Add in the constraints.
        GetRESPConstraints ( natoms, aRESP, bRESP, q, Atemp )

        # . Solve the equations.
        results = LinearLeastSquaresBySVD ( Atemp, Btemp )
        results["Solution"].CopyTo ( q )

        # . Determine the sum of squares.
        phi.CopyTo ( phitemp )
        fInteraction.VectorMultiply ( q, phitemp, alpha = -1.0, beta = 1.0, transpose = True )
        f = phitemp.DotSelf ( )

        # . Check for convergence.
        qold.Add ( q, scale = -1.0 )
        fdifference = math.fabs ( fold - f )
        qdifference = qold.AbsoluteMaximum ( )
        isConverged  = ( fdifference < ftolerance ) and ( qdifference < qtolerance )

        # . Printing.
        if doPrinting and ( niterations % logFrequency == 0 ):
            table.Entry ( "{:d}"  .format ( niterations                 ) )
            table.Entry ( "{:.6g}".format ( f                           ) )
            table.Entry ( "{:.6g}".format ( fdifference                 ) )
            table.Entry ( "{:.6g}".format ( qdifference                 ) )
            table.Entry ( "{:.6g}".format ( results["Condition Number"] ) )
            table.Entry ( "{:d}"  .format ( results["Rank"]             ) )

        # . End of loop.
        niterations += 1

    # . Stop printing.
    if doPrinting:
        table.Stop ( )
        if isConverged: log.Paragraph ( "RESP procedure converged."              )
        else:           log.Paragraph ( "Warning: RESP procedure not converged." )

    # . Finish up.
    return ( isConverged, f, q, condition, rank )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

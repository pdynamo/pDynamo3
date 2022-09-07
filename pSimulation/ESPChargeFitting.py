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

#
# . The function being minimized is:
#
#   F^2 = sum_g ( phi^ref_g - phi^model_g )^2    where the sum is over grid points, g.
#
#   phi^ref   is the reference potential calculated from the QC model.
#   phi^model is the potential from the RESP model and is equal to:
#
#   phi^model_g = sum_i q_i / r_ig      where the sum is over atoms, i, each with a charge q_i.
#
#   Minimizing gives:
#
#   dF/dq_i = - 2 sum_g ( phi^ref_g - phi^model_g ) * d phi^model_g/dq_i = 0
#
#   or:
#
#   sum_g phi^ref_g / r_ig = sum_gj q_j / r_ig * r_jg
#
#   This can be rewritten as A*Q = B where:
#
#   A is the matrix whose elements are A_ij = sum_g 1 / r_ig * r_jg
#   B is the vector whose elements are B_i  = sum_g phi^ref_g / r_ig
#
#   If R is the matrix of 1/r_ig then A = R*R^T, B = R * phi^ref and phi^model = R^T * Q.
#   Once phi^model is known then the original function (the sum of squares) is easily
#   calculated as the square of the residual ( phi^ref_g - phi^model_g ).
#
# . The above summary omits constraints although linear constraints of the form C^T * Q + l = 0 are easily added.
#
#   In the function below, the optional constraints argument is a tuple containing the vectors C and the values l.
#   If absent, a total charge constraint is added by default.
#

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Default values - RESP parameters suggested by Bayly et al (J. Phys. Chem. 1993, 97, 10269-10280).
_DefaultAResp = 0.0005 # . Weak constraints. Use 0.001 for strong constraints.
_DefaultBResp = 0.1

# . Default values - RESP iterator.
_DefaultFTolerance        = 1.0e-10
_DefaultLogFrequency      = 1
_DefaultMaximumIterations = 250
_DefaultQTolerance        = 1.0e-8

#===================================================================================================================================
# . Perform the charge fitting.
#===================================================================================================================================
def ESPChargeFitting ( system                          ,
                       aRESP           = _DefaultAResp ,
                       bRESP           = _DefaultBResp ,
                       constraints     = None          , # . A tuple of ( vectors, values ).
                       doRESP          = False         ,
                       fitBases        = None          ,
                       fitCoefficients = None          ,
                       includeNuclear  = True          ,
                       log             = logFile       ):
    """Perform an ESP charge fit."""
    # . Check for a system with a qcModel.
    if system.qcModel is None:
        raise Exception ( "System does not have a QC model." )
    # . Check for constraints.
    nAtoms = len ( system.qcState.qcAtoms )
    if constraints is None: # . A constraint on total charge is always added.
        # . Get the total charge.
        qTotal = float ( system.electronicState.charge )
        if not includeNuclear: qTotal -= sum ( system.qcState.nuclearCharges )
        nConstraints = 1
        values       = [ qTotal ]
        vectors      = Array.WithExtents ( nAtoms, 1 )
        vectors.Set ( 1.0 )
    elif len ( constraints ) > 0:
        ( vectors, values ) = constraints
        nConstraints = vectors.columns
        if ( len ( values ) != nConstraints ) or ( vectors.rows != nAtoms ):
            raise Exception ( "Invalid constraint dimensions." )
    else:
        nConstraints = 0
    # . Get the grid points for the QC atoms only and convert to bohrs.
    gridPoints = GenerateVanDerWaalsSurface ( system, log = log, qcAtomsOnly = True, scalingFactors = [ 1.4, 1.6, 1.8, 2.0 ] )
    gridPoints.Scale ( Units.Length_Angstroms_To_Bohrs )
    # . Allocate space - one larger than necessary for fInteraction.
    nDim         = nAtoms + nConstraints
    fInteraction = Array.WithExtents ( nDim, gridPoints.rows )
    fInteraction.Set ( 0.0 )
    # . Get the observed potentials at the grid points due to the full density.
    if ( fitBases is None ) and ( fitCoefficients is None ):
        phi = system.qcModel.GridPointPotentials ( system, gridPoints, includeNuclear = includeNuclear )
    # . Get the observed potentials at the grid points due to the fitted density.
    else:
        phi = Array.WithExtent ( gridPoints.shape[0] )
        phi.Set ( 0.0 )
        system.qcModel.gridPointEvaluator.f1Cp1V ( system, fitCoefficients, gridPoints, phi, fitBases = fitBases )
        if includeNuclear: system.qcModel.gridPointEvaluator.m1Cp1V ( system, gridPoints, phi )
    # . Get the interaction terms for each atom at the grid points.
    coordinates3 = system.scratch.qcCoordinates3AU
    GetInteractionTerms ( coordinates3, gridPoints, fInteraction )
    # . Get the A matrix and the B vector.
    A = Array.WithExtents ( nDim, nDim ) ; A.Set ( 0.0 )
    B = Array.WithExtent  ( nDim       ) ; B.Set ( 0.0 )
    A.MatrixMultiply ( fInteraction, fInteraction, yTranspose = True )
    fInteraction.VectorMultiply ( phi, B )
    # . Add the constraint terms.
    for c in range ( nConstraints ):
        for i in range ( nAtoms ):
            A[i,nAtoms+c] = vectors[i,c]
            A[nAtoms+c,i] = vectors[i,c]
        B[nAtoms+c] = values[c]
    # . Get |phi^2|.
    phi2 = phi.DotSelf ( )
    # . Iterative solution.
    if doRESP:
        ( isConverged, sos, solution, results ) = RESPIterator ( nAtoms, A, B, fInteraction, phi, aRESP, bRESP, log = log )
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
        items = [ ( "Charges"     , "{:d}"  .format ( nAtoms                      ) ) ,
                  ( "Constraints" , "{:d}"  .format ( nConstraints                ) ) ,
                  ( "Rank"        , "{:d}"  .format ( results["Rank"]             ) ) ,
                  ( "Error"       , "{:.5g}".format ( error                       ) ) ,
                  ( "Condition"   , "{:.5g}".format ( results["Condition Number"] ) ) ]
        if doRESP: items.append ( ( "Converged", repr ( isConverged ) ) )
        log.SummaryOfItems ( items, title = "ESP Charge Fitting Summary" )
    # . Return the charges.
    charges = Array.WithExtent ( nAtoms )
    for i in range ( nAtoms ):
        charges[i] = solution[i]
    # . Add in nuclear charges (electron charges are already negative as potentials are negative).
    if not includeNuclear: charges.Add ( system.qcState.nuclearCharges )
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
    nAtoms = len ( atomicNumbers )
    # . Set radii.
    radii = Array.WithExtent ( nAtoms )
    for ( i, n ) in enumerate ( atomicNumbers ):
        radii[i] = PeriodicTable[n].vdwRadius
    # . Get the grid points for a single center.
    ( basicGrid, weights ) = LebedevLaikovGrid_GetGridPoints ( gridAngularMomentum )
    # . Allocate space for the possible number of grid points.
    nPossible  = nAtoms * basicGrid.rows * len ( scalingFactors )
    gridPoints = Coordinates3.WithExtent ( nPossible )
    gridPoints.Set ( 0.0 )
    # . Initialization.
    nFound      = 0
    atomGrid    = Coordinates3.WithExtent ( basicGrid.rows )
    translation = Vector3.Null ( )
    atomGrid.Set ( 0.0 )
    # . Loop over scaling factors.
    for factor in scalingFactors:
        # . Loop over points.
        for i in range ( nAtoms ):
            # . Get the radius.
            iRadius = factor * radii[i]
            # . Get the translation.
            translation[0] = coordinates3[i,0]
            translation[1] = coordinates3[i,1]
            translation[2] = coordinates3[i,2]
            # . Get the scaled grid centered at the point.
            basicGrid.CopyTo   ( atomGrid    )
            atomGrid.Scale     ( iRadius     )
            atomGrid.Translate ( translation )
            # . Remove points that are within the current scaled radii of other points.
            for p in range ( atomGrid.rows ):
                isOK = True
                x    = atomGrid[p,0]
                y    = atomGrid[p,1]
                z    = atomGrid[p,2]
                for j in range ( nAtoms ):
                    if j != i:
                        dX       = coordinates3[j,0] - x
                        dY       = coordinates3[j,1] - y
                        dZ       = coordinates3[j,2] - z
                        jRadius2 = ( factor * radii[j] )**2
                        if ( dX**2 + dY**2 + dZ**2 ) <= jRadius2:
                            isOK = False
                            break
                if isOK:
                    gridPoints[nFound,0] = x
                    gridPoints[nFound,1] = y
                    gridPoints[nFound,2] = z
                    nFound += 1
    # . Reduce the array size if necessary.
    if nFound < nPossible:
        newPoints = Coordinates3.WithExtent ( nFound )
        for p in range ( nFound ):
            newPoints[p,0] = gridPoints[p,0]
            newPoints[p,1] = gridPoints[p,1]
            newPoints[p,2] = gridPoints[p,2]
        gridPoints = newPoints
    # . Do some printing.
    if LogFileActive ( log ):
        items = [ ( "Atoms"           , "{:d}".format ( nAtoms                 ) ) ,
                  ( "Surfaces"        , "{:d}".format ( len ( scalingFactors ) ) ) ,
                  ( "Found Points"    , "{:d}".format ( nFound                 ) ) ,
                  ( "Possible Points" , "{:d}".format ( nPossible              ) ) ]
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
def GetRESPConstraints ( nAtoms, aRESP, bRESP, charges, A ):
    """Add in RESP constraints to the diagonal elements of the A matrix."""
    for i in range ( nAtoms ):
        q = charges[i]
        A[i,i] += aRESP * q / math.sqrt ( q**2 + bRESP**2 )

#===================================================================================================================================
# . RESP iterator function.
#===================================================================================================================================
def RESPIterator ( nAtoms                                        ,
                   A                                             ,
                   rhs                                           ,
                   fInteraction                                  ,
                   phi                                           ,
                   aRESP                                         ,
                   bRESP                                         ,
                   fTolerance        = _DefaultFTolerance        ,
                   log               = logFile                   ,
                   logFrequency      = _DefaultLogFrequency      ,
                   maximumIterations = _DefaultMaximumIterations ,
                   qTolerance        = _DefaultQTolerance        ):

    """Solve the RESP equations by simple iteration."""
    # . Allocate space.
    n       = len ( rhs )
    ATemp   = Array.WithExtents ( n, n )
    BTemp   = Array.WithExtent  ( n )
    phiTemp = Array.WithExtent  ( len ( phi ) )
    q       = Array.WithExtent  ( n ) ; q.Set ( 0.0 ) # . With initial values.
    qOld    = Array.WithExtent  ( n )
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
    nIterations = 0
    isConverged = False
    # . Perform the iterations.
    while ( nIterations < maximumIterations ) and ( not isConverged ):
        # . Save old data.
        fOld = f
        q.CopyTo ( qOld )
        # . Fill new A and RHS.
        A.CopyTo   ( ATemp )
        rhs.CopyTo ( BTemp )
        # . Add in the constraints.
        GetRESPConstraints ( nAtoms, aRESP, bRESP, q, ATemp )
        # . Solve the equations.
        results = LinearLeastSquaresBySVD ( ATemp, BTemp )
        results["Solution"].CopyTo ( q )
        # . Determine the sum of squares.
        phi.CopyTo ( phiTemp )
        fInteraction.VectorMultiply ( q, phiTemp, alpha = -1.0, beta = 1.0, transpose = True )
        f = phiTemp.DotSelf ( )
        # . Check for convergence.
        qOld.Add ( q, scale = -1.0 )
        fDifference = math.fabs ( fOld - f )
        qDifference = qOld.AbsoluteMaximum ( )
        isConverged = ( fDifference < fTolerance ) and ( qDifference < qTolerance )
        # . Printing.
        if doPrinting and ( nIterations % logFrequency == 0 ):
            table.Entry ( "{:d}"  .format ( nIterations                 ) )
            table.Entry ( "{:.6g}".format ( f                           ) )
            table.Entry ( "{:.6g}".format ( fDifference                 ) )
            table.Entry ( "{:.6g}".format ( qDifference                 ) )
            table.Entry ( "{:.6g}".format ( results["Condition Number"] ) )
            table.Entry ( "{:d}"  .format ( results["Rank"]             ) )
        # . End of loop.
        nIterations += 1
    # . Stop printing.
    if doPrinting:
        table.Stop ( )
        if isConverged: log.Paragraph ( "RESP procedure converged."              )
        else:           log.Paragraph ( "Warning: RESP procedure not converged." )
    # . Finish up.
    return ( isConverged, f, q, results )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

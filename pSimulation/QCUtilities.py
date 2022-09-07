"""QC model utilities."""

# . Miscellaneous utilities including density and wavefunction analysis.

import math

from pCore                           import Align                  , \
                                            Clone                  , \
                                            logFile                , \
                                            LogFileActive
from pMolecule.QCModel               import FockConstruction_MakeCoefficientsFromFitIntegrals , \
                                            FockConstruction_MakeFromFitIntegralsCoulomb      , \
                                            FockConstruction_MakeFromTEIsCoulomb              , \
                                            OccupancyType                                     , \
                                            QCModelError                                      , \
                                            QCOrbitals
from pMolecule.QCModel.GaussianBases import GaussianBasisContainer , \
                                            GaussianBasisOperator  , \
                                            GaussianBasisType
from pScientific                     import Units
from pScientific.Arrays              import Array                  , \
                                            ArrayPrint             , \
                                            ArrayPrint2D           , \
                                            StorageType
from pScientific.Geometry3           import Matrix33               , \
                                            Vector3
from pScientific.LinearAlgebra       import Determinant            , \
                                            MatrixPseudoInverse

#
# . Center dependence of multipoles (Rc wrt R0):
#
#   dipole_C     = dipole_O     - Q * Rc
#   quadrupole_C = quadrupole_O - Rc x dipole_O^T - dipole_O x Rc^T + Q * Rc x Rc^T
#
#   x = outer product
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Tolerances for transition dipole check.
_DeviationTolerance = 1.0e-06
_OverlapTolerance   = 1.0e-03

#===================================================================================================================================
# . Determine multipole properties using density fit approaches.
#===================================================================================================================================
def DensityFitMultipoles ( target                     ,
                           center           = None    ,
                           fitBasis         = None    ,
                           fitOperator      = None    ,
                           log              = logFile ,
                           printBasisValues = False   ,
                           testFitFunction  = False   ):
    """Calculate QC multipoles based on density fitting."""
    # . A QC calculation should have been done on the target before calling this method.
    # . Initialization.
    qcModel      = target.qcModel
    state        = getattr ( target, qcModel.__class__._stateName )
    scratch      = target.scratch
    # . Fit basis and operator.
    evaluator = qcModel.integralEvaluator
    if fitBasis is None:
        fitBases = state.fitBases
        fitBasis = qcModel.fitBasis
    else:
        fitBases = GaussianBasisContainer.FromParameterDirectory ( fitBasis                              ,
                                                                   state.atomicNumbers                   ,
                                                                   basisType = GaussianBasisType.Density )
    if fitOperator is None:
        operator = qcModel.fitOperator
    else:
        operator = fitOperator
    # . Evaluate integrals.
    evaluator.f1Xg2i      ( target                                   ,
                            attribute       = "propertyFitIntegrals" ,
                            fitBases        = fitBases               ,
                            operator        = operator               ,
                            reportTag       = "Property Fit"         )
    evaluator.f1Xf1i_f1Oi ( target                                   ,
                            attribute       = "propertyFitMatrix"    , 
                            fitBases        = fitBases               ,
                            operator        = operator               ,
                            withConstraints = True                   )
    # . Do the fit.
    n               = scratch.propertyFitMatrix.rows # . With fit constraints.
    fitCoefficients = Array.WithExtent ( n )
    fitVector       = Array.WithExtent ( n )
    FockConstruction_MakeCoefficientsFromFitIntegrals ( scratch.onePDMP.density      ,
                                                        scratch.propertyFitIntegrals ,
                                                        scratch.propertyFitMatrix    ,
                                                        scratch.onePDMP.totalCharge  ,
                                                        fitCoefficients              ,
                                                        fitVector                    )
    # . Do the fit function.
    if testFitFunction:
        # . Get the TEIs with the appropriate operator.
        teis = None
        if operator is GaussianBasisOperator.Coulomb:
            teis = scratch.Get ( "twoElectronIntegrals" )
        # . Recalculate TEIs as operator is not Coulomb or Coulomb TEIs were not found.
        if teis is None:
            evaluator.f2Xf2i ( target                                ,
                               attribute = "fitTwoElectronIntegrals" ,
                               operator  = operator                  ,
                               reportTag = "Fit"                     )
            teis = scratch.Get ( "fitTwoElectronIntegrals" )
        # . Construct fit function.
        dTotal = scratch.onePDMP.density
        fTotal = Array.WithExtent ( dTotal.rows, storageType = StorageType.Symmetric )
        fTotal.Set ( 0.0 )
        eC = FockConstruction_MakeFromTEIsCoulomb ( dTotal ,
                                                    teis   ,
                                                    fTotal )
        fTotal.Set ( 0.0 )
        eF = FockConstruction_MakeFromFitIntegralsCoulomb ( dTotal                       ,
                                                            scratch.propertyFitIntegrals ,
                                                            scratch.propertyFitMatrix    ,
                                                            scratch.onePDMP.totalCharge  ,
                                                            fitCoefficients              ,
                                                            fTotal                       )
        fitFunction  = 2.0 * ( eC - eF )
        fitReference = 2.0 * eC
        if fitFunction >= 0.0:
            fitErrorAbs = math.sqrt ( fitFunction )
            fitErrorRel = math.sqrt ( fitFunction / math.fabs ( fitReference ) )
        else:
            fitErrorAbs = -1.0
            fitErrorRel = -1.0
    # . Multipoles.
    # . Integrals.
    overlap                          = evaluator.f1Oi ( target, fitBases = fitBases                  )
    ( dX, dY, dZ )                   = evaluator.f1Di ( target, fitBases = fitBases, center = center )
    ( qXX, qYY, qZZ, qXY, qXZ, qYZ ) = evaluator.f1Qi ( target, fitBases = fitBases, center = center )
    # . Atomic multipoles.
    # . Nuclear contributions.
    coordinates3 = scratch.qcCoordinates3AU
    cX = cY = cZ = 0.0
    if center is not None:
        cX = center[0] ; cY = center[1] ; cZ = center[2]
    # . Charges.
    charges = Clone ( state.nuclearCharges )
    # . Dipoles.
    dipoles = Array.WithExtents ( len ( charges ), 3 )
    for ( a, q ) in enumerate ( charges ):
        dipoles[a,0] = q * ( coordinates3[a,0] - cX )
        dipoles[a,1] = q * ( coordinates3[a,1] - cY )
        dipoles[a,2] = q * ( coordinates3[a,2] - cZ )
    # . Quadrupoles.
    quadrupoles = Array.WithExtents ( len ( charges ), 6 )
    for ( a, q ) in enumerate ( charges ):
        x = coordinates3[a,0] - cX
        y = coordinates3[a,1] - cY
        z = coordinates3[a,2] - cZ
        quadrupoles[a,0] = q * x * x
        quadrupoles[a,1] = q * y * y
        quadrupoles[a,2] = q * z * z
        quadrupoles[a,3] = q * x * y
        quadrupoles[a,4] = q * x * z
        quadrupoles[a,5] = q * y * z
    # . Electronic contributions.
    functionCenters = fitBases.functionCenters
    for ( f, a ) in enumerate ( functionCenters ):
        c                 = fitCoefficients[f]
        charges[a]       -= c * overlap[f]
        dipoles[a,0]     -= c * dX[f]
        dipoles[a,1]     -= c * dY[f]
        dipoles[a,2]     -= c * dZ[f]
        quadrupoles[a,0] -= c * qXX[f]
        quadrupoles[a,1] -= c * qYY[f]
        quadrupoles[a,2] -= c * qZZ[f]
        quadrupoles[a,3] -= c * qXY[f]
        quadrupoles[a,4] -= c * qXZ[f]
        quadrupoles[a,5] -= c * qYZ[f]
    # . Total multipoles.
    charge = sum ( charges )
    # . Dipole.
    dipole    = Vector3.Null ( )
    dipole[0] = sum ( dipoles[:,0] ) * Units.Dipole_Atomic_Units_To_Debyes
    dipole[1] = sum ( dipoles[:,1] ) * Units.Dipole_Atomic_Units_To_Debyes
    dipole[2] = sum ( dipoles[:,2] ) * Units.Dipole_Atomic_Units_To_Debyes
    # . Quadrupole.
    quadrupole      = Matrix33.Null ( )
    quadrupole[0,0] = sum ( quadrupoles[:,0] ) * Units.Quadrupole_Atomic_Units_To_Buckinghams
    quadrupole[1,1] = sum ( quadrupoles[:,1] ) * Units.Quadrupole_Atomic_Units_To_Buckinghams
    quadrupole[2,2] = sum ( quadrupoles[:,2] ) * Units.Quadrupole_Atomic_Units_To_Buckinghams
    quadrupole[0,1] = sum ( quadrupoles[:,3] ) * Units.Quadrupole_Atomic_Units_To_Buckinghams
    quadrupole[0,2] = sum ( quadrupoles[:,4] ) * Units.Quadrupole_Atomic_Units_To_Buckinghams
    quadrupole[1,2] = sum ( quadrupoles[:,5] ) * Units.Quadrupole_Atomic_Units_To_Buckinghams
    quadrupole[1,0] = quadrupole[0,1]
    quadrupole[2,0] = quadrupole[0,2]
    quadrupole[2,1] = quadrupole[1,2]
    # . Printing.
    if LogFileActive ( log ):

        # . Header.
        log.Heading ( "Density Fit Multipoles using Basis \"{:s}\" and Operator \"{:s}\"".format ( fitBasis, operator.name.lower ( ) ), includeBlankLine = True )

        # . Basis values.
        if printBasisValues:
            fitLabels = fitBases.functionLabels
            table     = log.GetTable ( columns = [ 10 ] + 12 * [ 12 ] )
            table.Start   ( )
            table.Title   ( "Basis Set Values" )
            table.Heading ( "Atom"                        )
            table.Heading ( "Function"                    )
            table.Heading ( "Fit Value"                   )
            table.Heading ( "Overlap"                     )
            table.Heading ( "Dipole"     , columnSpan = 3 )
            table.Heading ( "Quadrupole" , columnSpan = 6 )
            for tag in ( "", "", "", "", "X", "Y", "Z", "XX", "YY", "ZZ", "XY", "XZ", "YZ" ):
                table.Heading ( tag )
            for f in range ( len ( fitCoefficients ) - 1 ):
                a = functionCenters[f]
                table.Entry ( target.atoms[a].path, align = Align.Left )
                table.Entry (         fitLabels[f], align = Align.Left )
                table.Entry ( "{:.3f}".format ( fitCoefficients[f] ) )
                table.Entry ( "{:.3f}".format ( overlap        [f] ) )
                table.Entry ( "{:.3f}".format ( dX             [f] ) )
                table.Entry ( "{:.3f}".format ( dY             [f] ) )
                table.Entry ( "{:.3f}".format ( dZ             [f] ) )
                table.Entry ( "{:.3f}".format ( qXX            [f] ) )
                table.Entry ( "{:.3f}".format ( qYY            [f] ) )
                table.Entry ( "{:.3f}".format ( qZZ            [f] ) )
                table.Entry ( "{:.3f}".format ( qXY            [f] ) )
                table.Entry ( "{:.3f}".format ( qXZ            [f] ) )
                table.Entry ( "{:.3f}".format ( qYZ            [f] ) )
            table.Stop ( )

        # . Atomic values.
        table = log.GetTable ( columns = [ 10 ] + 10 * [ 15 ] )
        table.Start   ( )
        table.Title   ( "Atomic Multipoles (Atomic Units)" )
        table.Heading ( "Atom"                        )
        table.Heading ( "Charge"                      )
        table.Heading ( "Dipole"     , columnSpan = 3 )
        table.Heading ( "Quadrupole" , columnSpan = 6 )
        for tag in ( "", "", "X", "Y", "Z", "XX", "YY", "ZZ", "XY", "XZ", "YZ" ):
            table.Heading ( tag )
        for a in range ( len ( charges ) ):
            table.Entry ( target.atoms[a].path, align = Align.Left )
            table.Entry ( "{:.3f}".format ( charges    [a  ] ) )
            table.Entry ( "{:.3f}".format ( dipoles    [a,0] ) )
            table.Entry ( "{:.3f}".format ( dipoles    [a,1] ) )
            table.Entry ( "{:.3f}".format ( dipoles    [a,2] ) )
            table.Entry ( "{:.3f}".format ( quadrupoles[a,0] ) )
            table.Entry ( "{:.3f}".format ( quadrupoles[a,1] ) )
            table.Entry ( "{:.3f}".format ( quadrupoles[a,2] ) )
            table.Entry ( "{:.3f}".format ( quadrupoles[a,3] ) )
            table.Entry ( "{:.3f}".format ( quadrupoles[a,4] ) )
            table.Entry ( "{:.3f}".format ( quadrupoles[a,5] ) )
        table.Stop ( )

        # . System values.
        log.Paragraph ( "Total System Charge (e) = {:.3f}".format ( charge ) )
        ArrayPrint    ( dipole    , itemFormat = "{:.3f}", log = log, title = "Dipole (Debyes)" )
        ArrayPrint2D  ( quadrupole, itemFormat = "{:.3f}", log = log, title = "Quadrupole (Buckinghams)" )

        # . Fit function.
        if testFitFunction:
            log.SummaryOfItems ( [ ( "Fit Basis"              , fitBasis                         ) , 
                                   ( "Fit Function"           , "{:.5f}".format ( fitFunction  ) ) ,
                                   ( "Fit Error Absolute"     , "{:.5f}".format ( fitErrorAbs  ) ) ,
                                   ( "Fit Error Relative (%)" , "{:.1%}".format ( fitErrorRel  ) ) ,
                                   ( "Fit Operator"           , operator.name.lower ( )          ) ,
                                   ( "Fit Reference"          , "{:.5f}".format ( fitReference ) ) ] , title = "Fit Function" )
    # . Finish up.
    results = { "Charges"          : charges         ,
                "Dipole"           : dipole          ,
                "Fit Bases"        : fitBases        ,
                "Fit Coefficients" : fitCoefficients ,
                "Fit Operator"     : fitOperator     ,
                "Quadrupole"       : quadrupole      }
    if testFitFunction:
        results["Fit Function"          ] = fitFunction
        results["Fit Reference"         ] = fitReference
        results["Fit Error Absolute"    ] = fitErrorAbs
        results["Fit Error Relative (%)"] = fitErrorRel
    return results

#===================================================================================================================================
# . Determine the overlap between two single determinant wavefunctions with mutually non-orthogonal orbitals.
#===================================================================================================================================
#
# . The overlap is:
#
#   <Psi_1|Psi_2> = Det [<phi^occ_1|phi^occ_2>_alpha] x Det [<phi^occ_1|phi^occ_2>_beta]
#
#   It is, of course, zero if the number of occupied alpha and beta orbitals are different
#   for the two wavefunctions.
#
#   Note that to be completely general the orbitals should be scaled by their respective occupancies.
#
#   This method can be used for any one-electron operator (not just the overlap).
#

def DetermineWavefunctionOverlap ( wavefunction1, wavefunction2, overlap = None ):
    """Determine the overlap between two wavefunctions."""
    # . The wavefunction arguments are either single sets or tuples of two sets (alpha/beta) of QCOrbitals instances.
    # . Gather orbitals.
    ( noSpin, C1P, C1Q, C2P, C2Q ) = _GatherOrbitals ( wavefunction1, wavefunction2 )
    # . Calculate the overlap.
    detP = _Overlap ( C1P, C2P, S = overlap )
    if noSpin: detQ = detP
    else:      detQ = _Overlap ( C1Q, C2Q, S = overlap )
    # . Finish up.
    return ( detP * detQ )

def _GatherOrbitals ( wavefunction1, wavefunction2 ):
    """Gather orbitals from two wavefunctions."""
    data   = []
    noSpin = True
    for Psi in ( wavefunction1, wavefunction2 ):
        if   isinstance ( Psi, QCOrbitals ):
            data.append ( ( Psi, Psi ) )
        elif isinstance ( Psi, ( list, tuple ) ) and ( len ( Psi ) == 2 ) and all ( [ isinstance ( C, QCOrbitals ) for C in Psi ] ):
            data.append ( Psi )
            noSpin = False
        else: raise QCModelError ( "Invalid wavefunction input." )
    return ( noSpin, data[0][0], data[0][1], data[1][0], data[1][1] )

def _Overlap ( o1, o2, S = None ):
    """Overlap between two sets of orbitals."""
    d = 0.0
    m = o1.occupancyHandler.numberOccupied
    if m == o2.occupancyHandler.numberOccupied:
        c1      = o1.orbitals[:,0:m]
        c2      = o2.orbitals[:,0:m]
        overlap = Array.WithExtents ( m, m )
        if S is None:
            Sc2 = c2
        else:
            Sc2 = Array.WithExtents ( S.extent, m )
            S.PostMatrixMultiply ( c2, Sc2 )
        overlap.MatrixMultiply ( c1, Sc2, xTranspose = True )
        d = Determinant ( overlap )
    return d

#===================================================================================================================================
# . Enforce the charge restraints on the densities in target.
# . This is only useful if a previous QC calculation (normally without restraints) has been performed
# . so that the densities and the overlap matrix already exist.
# . All information is taken about the restraints is taken from target.
# . This method assumes of course that the restraints are linear in the densities!
# . Note that application of the constraints means that the densities may no longer be "valid"
# . and satisfy the correct one-PDM conditions.
#===================================================================================================================================
def EnforceChargeRestraints ( target ):
    """Enforce charge restraints on the densities in target."""
    # . Basic checks.
    if ( target.chargeRestraintModel                is None ) or \
       ( target.chargeRestraintModelState           is None ) or \
       ( target.chargeRestraintModelState.evaluator is None ) or \
       ( target.qcModel                        is None ) or \
       ( len ( target.chargeRestraintModel ) == 0 ):
        raise QCModelError ( "The target has missing charge restraints or QC model." )
    if ( target.electronicState.isSpinRestricted ) and \
       ( target.chargeRestraintModel.spinRestraints > 0 ):
        raise QCModelError ( "Spin restraints require an unrestricted QC model." )
    if ( target.chargeRestraintModel.chargeRestraints > 0 ) and \
       ( target.scratch.Get ( "onePDMP", None ) is None   ):
        raise QCModelError ( "There are charge restraints but no total charge density matrix." )
    if ( target.chargeRestraintModel.spinRestraints > 0 ) and \
       ( target.scratch.Get ( "onePDMQ", None ) is None ):
        raise QCModelError ( "There are spin restraints but no spin density matrix." )
    # . Gather data.
    scratch    = target.scratch
    evaluator  = target.chargeRestraintModelState.evaluator
    restraints = target.chargeRestraintModel.restraints
    n  = target.qcState.orbitalBases.numberOfFunctions
    nT = ( n * ( n + 1 ) ) // 2
    Qt = sum ( target.qcState.nuclearCharges ) - target.electronicState.charge
    Qs = target.electronicState.multiplicity - 1
    S  = target.scratch.Get ( "overlapMatrix", None )
    if S is None:
        S = Array.WithExtent ( n, storageType = StorageType.Symmetric )
        S.Set ( 0.0 )
        for i in range ( n ): S[i,i] = 1.0
    # . Gather charge and spin restraints separately.
    toApply = []
    for ( ( spin, q0 ) ) in ( ( False, Qt ) ,
                              ( True , Qs ) ):
        WQ = []
        for key in sorted ( restraints.keys ( ) ):
            restraint = restraints[key]
            if restraint.isSpin is spin:
                ( w, c ) = evaluator.ChargeRestraintMatrix ( target, restraint )
                WQ.append ( ( w, restraint.target - c ) )
        if len ( WQ ) > 0:
            Q = Array.WithExtents (     len ( WQ ) + 1 )
            W = Array.WithExtents ( nT, len ( WQ ) + 1 )
            Q[0] = q0
            S.ScaleOffDiagonal ( 2.0 )
            S.iterator.CopyTo ( W[:,0] )
            S.ScaleOffDiagonal ( 0.5 )
            for ( c, ( w, q ) ) in enumerate ( WQ ):
                Q[c+1] = q
                w.ScaleOffDiagonal ( 2.0 )
                w.iterator.CopyTo ( W[:,c+1] )
            nVectors         = W.OrthogonalizeColumns ( constants = Q )
            linearDependence = ( W.columns - nVectors ) # . For output.
            if linearDependence > 0:
                Q = Q[  0:nVectors]
                W = W[:,0:nVectors]
            if nVectors > 1: toApply.append ( ( spin, W, Q ) )
    # . Apply the restraints to the target densities (and for completeness also the Fock matrices!).
    for ( spin, W, Q ) in toApply:
        if spin: PDM = scratch.Get ( "onePDMQ", None )
        else:    PDM = scratch.Get ( "onePDMP", None )
        P = PDM.density
        F = PDM.fock
        for c in range ( len ( Q ) ):
            fFactor =      - F.iterator.Dot ( W[:,c] )
            pFactor = Q[c] - P.iterator.Dot ( W[:,c] )
            F.iterator.Add ( W[:,c], scale = fFactor )
            P.iterator.Add ( W[:,c], scale = pFactor )

#===================================================================================================================================
# . Make excited states by swapping occupied and virtual orbitals and their energies.
#===================================================================================================================================#
# . The swaps to perform are given in swapsP and swapsQ.
#
#   swapsP - the orbitals in a spin-restricted state or the alpha orbitals in a spin-unrestricted state.
#   swapsQ - the beta orbitals in a spin-unrestricted state.
#
#   The swaps are lists of index pairs of the form ( occupied index, virtual index ). The indices are
#   not absolute indices but relative indices so that the occupied index is relative to the homo, and
#   the virtual index to the lumo. For example:
#
#   (  0, 0 ) = homo   <-> lumo
#   ( -1, 0 ) = homo-1 <-> lumo
#   ( -2, 2 ) = homo-2 <-> lumo+2
#
# . The occupancies of the swapped orbital set remain unchanged.
#
def MakeExcitedStateByOrbitalSwapping ( target, swapsP = None, swapsQ = None ):
    """Create an excited state by swapping occupied and virtual orbitals in target from a previous QC calculation."""
     # . Swap input.
    isSpinRestricted = target.electronicState.isSpinRestricted
    if (       isSpinRestricted   and ( ( swapsP is None ) or  ( swapsQ is not None ) ) ) or \
       ( ( not isSpinRestricted ) and ( ( swapsP is None ) and ( swapsQ is     None ) ) ):
        raise QCModelError ( "Invalid swaps input." )
    # . Orbitals and densities.
    scratch   = target.scratch
    onePDMP   = scratch.Get ( "onePDMP"  , None )
    orbitalsP = scratch.Get ( "orbitalsP", None )
    isOK      = ( onePDMP is not None ) and ( orbitalsP is not None )
    toProcess = [ ( orbitalsP, onePDMP, swapsP ) ]
    if not isSpinRestricted:
        onePDMQ   = scratch.Get ( "onePDMQ"  , None )
        orbitalsQ = scratch.Get ( "orbitalsQ", None )
        isOK      = isOK and ( onePDMQ is not None ) and ( orbitalsQ is not None )
        toProcess.append ( ( orbitalsQ, onePDMQ, swapsQ ) )
    if not isOK: raise QCModelError ( "Invalid wavefunction to make excited state." )
    # . Swap the orbitals and update the corresponding density.
    for ( orbitals, onePDM, swaps ) in toProcess:
        if swaps is not None:
            numberOccupied = orbitals.occupancyHandler.numberOccupied
            numberOrbitals = orbitals.numberOrbitals
            homo           = numberOccupied - 1
            lumo           = numberOccupied
            energies       = orbitals.energies
            mos            = orbitals.orbitals
            occupancies    = orbitals.occupancies
            for ( o, v ) in swaps:
                occupied = homo - abs ( o )
                virtual  = lumo + abs ( v )
                if ( occupied >= 0 ) and ( virtual < numberOrbitals ):
                    ( energies[occupied], energies[virtual] ) = ( energies[virtual], energies[occupied] )
                    mos[:,occupied].Swap ( mos[:,virtual] )
                else: raise QCModelError ( "Swap orbital indices ({:d} <-> {:d}) out of range [0,{:d}].".format ( occupied, virtual, numberOrbitals ) )
            onePDM.density.MakeFromEigenSystem ( numberOccupied ,
                                                 occupancies    ,
                                                 mos            )
    
#===================================================================================================================================
# . Set up a MOM excited state calculation by modifying in-place and existing QC system.
#===================================================================================================================================#
def SetUpMOMExcitedStateCalculation ( target, momElectronicState ):
    """Set up a MOM excited state calculation from an existing QC system."""
    # . The system must have a QC model and its densities and orbitals must exist.
    isOK = ( target.electronicState is not None ) and \
           ( target.qcModel         is not None ) and \
           ( target.scratch         is not None )
    if isOK:
        scratch   = target.scratch
        onePDMP   = scratch.Pop ( "onePDMP"   )
        onePDMQ   = scratch.Pop ( "onePDMQ"   )
        orbitalsP = scratch.Pop ( "orbitalsP" )
        orbitalsQ = scratch.Pop ( "orbitalsQ" )
        isOK      = ( onePDMP   is not None ) and \
                    ( orbitalsP is not None )
        if isOK and ( not target.electronicState.isSpinRestricted ):
            isOK  = ( onePDMQ   is not None ) and \
                    ( orbitalsQ is not None )
    if not isOK: raise QCModelError ( "Invalid system for setting up MOM calculation." )
    # . The new electronic state must use MOM occupancy handlers.
    if momElectronicState.occupancyType is not OccupancyType.MOM:
        raise QCModelError ( "The input electronic state does not use MOM occupancy handlers." )
    # . The old and new electronic states must also be compatible.
    if momElectronicState.isSpinRestricted and ( not target.electronicState.isSpinRestricted ):
        raise QCModelError ( "Unable to setup a restricted excited state calculation from an unrestricted target." )
    # . The set up is already done if the target already uses MOM occupancy handlers.
    if target.electronicState.occupancyType is not OccupancyType.MOM:
        # . Change the electronic state in target which also clears its scratch.
        target.electronicState = momElectronicState
        # . Rebuild the target's (guess) densities and orbitals.
        target.qcModel.SetUpDensities ( target )
        target.qcModel.SetUpOrbitals  ( target )
        # . Copy back the original densities, orbitals and orbital energies.
        # . Copying the densities is not strictly necessary but is done for completeness.
        if target.electronicState.isSpinRestricted and ( not momElectronicState.isSpinRestricted ):
            onePDMQ   = onePDMP
            orbitalsQ = orbitalsP
        onePDMP.density  .CopyTo  ( scratch.onePDMP.density    )
        orbitalsP.energies.CopyTo ( scratch.orbitalsP.energies )
        orbitalsP.orbitals.CopyTo ( scratch.orbitalsP.orbitals )
        if not momElectronicState.isSpinRestricted:
            onePDMQ.density  .CopyTo  ( scratch.onePDMQ.density    )
            orbitalsQ.energies.CopyTo ( scratch.orbitalsQ.energies )
            orbitalsQ.orbitals.CopyTo ( scratch.orbitalsQ.orbitals )

#===================================================================================================================================
# . Transition dipole moment and oscillator strength between two single determinant wavefunctions with mutually non-orthogonal
# . orbitals.
# . deltaE is the energy difference in hartrees between the states, and the wavefunction input
# . is the same as for DetermineWavefunctionOverlap.
# . The position rather than the velocity formulation of the transition dipole moment is used.
#===================================================================================================================================#
def TransitionDipoleMoment ( target, deltaE, wavefunction1, wavefunction2, center = None, doCheck = False ):
    """Determine the transition dipole moment and oscillator strength between two wavefunctions.
    
    The transition dipole moment is returned in Debyes and the oscillator strength is dimensionless.
    """
    # . Dipole integrals.
    S = target.scratch.Get ( "overlapMatrix", None )
    ( dX, dY, dZ ) = target.qcModel.integralEvaluator.f1Df1i ( target, center = center )
    # . Gather orbitals.
    ( noSpin, C1P, C1Q, C2P, C2Q ) = _GatherOrbitals ( wavefunction1, wavefunction2 )
    # . Transition dipole moment.
    tDipole = Vector3.Null ( )
    tDipole.Set ( 0.0 )
    detP = _TransitionDipole ( C1P, C2P, dX, dY, dZ, tDipole, S = S )
    if noSpin: detQ = 2.0 * detP # . 2 to account for double occupancy. 
    else:      detQ = _TransitionDipole ( C1Q, C2Q, dX, dY, dZ, tDipole, S = S )
    O12 = detP * detQ
    tDipole.Scale ( O12 )
    # . If the overlap is one then check the TDs versus trace ( P * Dx ), etc.
    if doCheck and ( math.fabs ( O12 - 1.0 ) <= _OverlapTolerance ):
        # . Make Pt.
        Pa = Array.WithExtent ( dX.extent, storageType = StorageType.Symmetric )
        Pb = Array.WithExtent ( dX.extent, storageType = StorageType.Symmetric )
        Pa.MakeFromEigenSystem ( C1P.occupancyHandler.numberOccupied, C1P.occupancies, C1P.orbitals )
        Pb.MakeFromEigenSystem ( C1Q.occupancyHandler.numberOccupied, C1Q.occupancies, C1Q.orbitals )
        Pa.Add ( Pb )
        # . Calculate the electronic dipole the standard way.
        mu    = Vector3.Null ( )
        mu[0] = Pa.TraceOfProduct ( dX )
        mu[1] = Pa.TraceOfProduct ( dY )
        mu[2] = Pa.TraceOfProduct ( dZ )
        mu.Add ( tDipole, scale = -1.0 )
        deviation = mu.AbsoluteMaximum ( )
        if deviation > _DeviationTolerance:
            raise QCModelError ( "Failed transition dipole consistency check with deviation of {:.5g}.".format ( deviation ) )
    # . Return the dimensionless oscillator strength and the transition dipole moment in Debyes.
    f = 2.0 * tDipole.DotSelf ( ) * deltaE / 3.0
    tDipole.Scale ( - Units.Dipole_Atomic_Units_To_Debyes ) # . Negative as electrons.
    return ( O12, tDipole, f )

def _TransitionDipole ( o1, o2, dX, dY, dZ, tDipole, S = None ):
    """Transition dipole between two sets of orbitals."""
    d = 0.0
    m = o1.occupancyHandler.numberOccupied
    if m == o2.occupancyHandler.numberOccupied:
        c1     = o1.orbitals[:,0:m]
        c2     = o2.orbitals[:,0:m]
        Mc2    = Array.WithExtents ( dX.extent, m )
        M12    = Array.WithExtents ( m        , m )
        S12    = Array.WithExtents ( m        , m )
        S12inv = Array.WithExtents ( m        , m )
        if S is None:
            Sc2 = c2
        else:
            Sc2 = Mc2
            S.PostMatrixMultiply ( c2, Sc2 )
        S12.MatrixMultiply ( c1, Sc2, xTranspose = True )
        MatrixPseudoInverse ( S12, S12inv )
        d = Determinant ( S12 )
        for ( c, M ) in enumerate ( ( dX, dY, dZ ) ):
            M.PostMatrixMultiply ( c2, Mc2 )
            M12.MatrixMultiply ( c1, Mc2, xTranspose = True )
            tDipole[c] += S12inv.TraceOfProduct ( M12 )
    return d

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

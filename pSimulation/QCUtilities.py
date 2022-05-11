"""QC model utilities."""

# . Miscellaneous utilities including density and wavefunction analysis.

import math

from pCore                           import Align                  , \
                                            Clone                  , \
                                            logFile                , \
                                            LogFileActive
from pMolecule.QCModel               import FockConstruction_MakeCoefficientsFromFitIntegrals , \
                                            FockConstruction_MakeFromFitIntegralsCoulomb      , \
                                            FockConstruction_MakeFromTEIsCoulomb
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

#
# . Center dependence of multipoles (Rc wrt R0):
#
#   dipole_C     = dipole_O     - Q * Rc
#   quadrupole_C = quadrupole_O - Rc x dipole_O^T - dipole_O x Rc^T + Q * Rc x Rc^T
#
#   x = outer product
#

#===================================================================================================================================
# . Methods.
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
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

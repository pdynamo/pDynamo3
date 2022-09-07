"""Test for calculating the energies of various systems with charge and spin restraints."""

import os, os.path

from Definitions        import dataPath                        , \
                               referenceDataPath
from pBabel             import ImportSystem          
from pCore              import Align                           , \
                               DataType                        , \
                               logFile                         , \
                               TestScriptExit_Fail             , \
                               YAMLUnpickle               
from pMolecule          import EnergyModelPriority             , \
                               RestraintEnergyModel            , \
                               SystemGeometryObjectiveFunction
from pMolecule.QCModel  import ChargeModel                     , \
                               ChargeRestraint                 , \
                               ChargeRestraintModel            , \
                               DFTGridAccuracy                 , \
                               DFTGridIntegrator               , \
                               DIISSCFConverger                , \
                               ElectronicState                 , \
                               QCModelDFT                      , \
                               QCModelMNDO
from pScientific.Arrays import Array
from pSimulation        import EnforceChargeRestraints

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_DoDFT                = False
_DoMulliken           = False
_ForceUnrestricted    = False
_MaximumGradientAtoms = 10
_TestGradients        = False

# . Restraint options - a larger force constant is OK for DFT and HF but fails for many of the MNDO cases.
_ForceConstant        = 1.0 # . Hartrees.
_RestraintEnergyModel = RestraintEnergyModel.Harmonic ( 0.0, _ForceConstant )

# . Paths.
xyzPath = os.path.join ( dataPath, "xyz" )

# . QC options.
# . Use a density tolerance of 10^-10 to get properly converged gradients at the expense of some energy failures!
if _TestGradients: _DensityTolerance = 1.0e-10
else:              _DensityTolerance = 1.0e-08 
_Iterations       = 500
_MaximumMemory    = 6.0

# . QC models.
converger = DIISSCFConverger.WithOptions ( densityTolerance  = _DensityTolerance ,
                                           maximumIterations = _Iterations       )
_QCModels = [ ( "MNDO"         , QCModelMNDO.WithOptions ( converger     = converger      ,
                                                           hamiltonian   = "am1"          ), 100, ChargeModel.MNDO    ) ,
              ( "HF/def2-sv(p)", QCModelDFT.WithOptions  ( converger     = converger      ,
                                                           functional    = "hf"           ,
                                                           maximumMemory = _MaximumMemory ,
                                                           orbitalBasis  = "def2-sv(p)"   ),  25, ChargeModel.Loewdin ) ]
if _DoDFT:
    _DFTGridAccuracy  = DFTGridAccuracy.Medium
    gridIntegrator    = DFTGridIntegrator.WithOptions ( accuracy = _DFTGridAccuracy ,
                                                        inCore   = True             )
    _QCModels.append ( ( "BLYP/dgauss-dzvp", QCModelDFT.WithOptions ( converger      = converger           ,
                                                                      fitBasis       = "dgauss-a1-dftjfit" ,
                                                                      functional     = "blyp"              ,
                                                                      gridIntegrator = gridIntegrator      ,
                                                                      maximumMemory  = _MaximumMemory      ,
                                                                      orbitalBasis   = "dgauss-dzvp"       ),  25, ChargeModel.Loewdin ) )
if _DoMulliken:
    _QCModels.append ( ( "HF/3-21g", QCModelDFT.WithOptions  ( converger     = converger           ,
                                                               functional    = "hf"                ,
                                                               maximumMemory = _MaximumMemory      ,
                                                               orbitalBasis  = "3-21g"             ),  50, ChargeModel.Mulliken ) )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Gather systems with charge constraints.
_SystemData = YAMLUnpickle ( os.path.join ( referenceDataPath, "ChargeRestraintTestSystems.yaml" ) )
systems     = []
for ( name, value ) in _SystemData.items ( ):
    charge         = value["Charge"      ]
    multiplicity   = value["Multiplicity"]
    path           = os.path.join ( xyzPath, value["Structure File"] )
    # . Charge restraints.
    restraints = [ ChargeRestraint.WithOptions ( energyModel = _RestraintEnergyModel ,
                                                 indices     = Array.FromIterable ( constraint["Indices"], dataType = DataType.Integer ) ,
                                                 isSpin      = constraint["IsSpin"]  ,
                                                 target      = constraint["Target"]  ) for constraint in value["Restraints"] ]
    chargeRestraintModel = ChargeRestraintModel ( )
    for ( isSpin, tag ) in ( ( False, "Charge" ), ( True, "Spin" ) ):
        n = 0
        for restraint in restraints:
            if restraint.isSpin is isSpin:
                chargeRestraintModel["{:s}{:d}".format ( tag, n )] = restraint
                n += 1
    # . System.
    hasSpin                = ( multiplicity != 1 ) or any ( [ item.isSpin for item in restraints ] )
    system                 = ImportSystem ( path, log = None )
    system.label           = name
    system.electronicState = ElectronicState.WithOptions ( charge           = charge       , 
                                                           isSpinRestricted = not ( hasSpin or _ForceUnrestricted ) ,
                                                           multiplicity     = multiplicity )
    systems.append ( ( name, system, chargeRestraintModel ) )

# . Loop over systems.
eFailures                = []
eTests                   = 0
maximumGradientDeviation = 0.0
for ( name, system, chargeRestraintModel ) in systems:
    hasSpin = ( not system.electronicState.isSpinRestricted )
    # . Loop over qcModels.
    for ( qcTag, qcModel, qcMaximumAtoms, rChargeModel ) in _QCModels:
        if len ( system.atoms ) > qcMaximumAtoms: continue
        system.DefineQCModel ( qcModel )
        results = []
        # . Do unrestrained and restrained calculations.
        energy0 = None
        for ( doRestrained, rTag ) in ( ( ( False, "Unrestrained" ) ,
                                          ( True , "Restrained"   ) ) ):
            if doRestrained:
                chargeRestraintModel.chargeModel = rChargeModel
                system.AddEnergyModel ( "chargeRestraintModel", chargeRestraintModel, priority = EnergyModelPriority.QCAddOns, valueClass = ChargeRestraintModel )
                EnforceChargeRestraints ( system )
                eTests += 1
            system.Summary ( )
            try:
                energy  = system.Energy ( )
                if energy0 is None: energy0 = energy
                deltaE  = energy - energy0
                charges = system.qcModel.AtomicCharges ( system, chargeModel = rChargeModel )
                dipole  = system.DipoleMoment ( )
                if hasSpin: spins = system.qcModel.AtomicSpins ( system, chargeModel = rChargeModel )
                else:       spins = None
                cValues = chargeRestraintModel.Values ( charges, spins )
                results.append ( ( rTag, deltaE, charges, spins, dipole.Norm2 ( ), cValues ) )
                if _TestGradients and doRestrained and ( len ( system.atoms ) < _MaximumGradientAtoms ):
                    of = SystemGeometryObjectiveFunction.FromSystem ( system )
                    gradientDeviation        = of.TestGradients ( )
                    maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )
            except Exception as error:
                if doRestrained: eFailures.append ( ( name, qcTag ) )
                print ( error.args[0] )
        # . Output the results.
        if len ( results ) > 0:
            atomNames = [ atom.path for atom in system.atoms ]
            cTags     = sorted ( chargeRestraintModel.restraints.keys ( ) )
            n         = len ( atomNames )
            columns   = [ 10, 10 ] + len ( results ) * [ 20 ]
            table     = logFile.GetTable ( columns = columns )
            table.Start   ( )
            table.Title   ( "Results for {:s} with QC Model {:s}/{:s}".format ( name, qcTag, rChargeModel.name ) )
            table.Heading ( "Property"       , columnSpan = 2 )
            table.Heading ( "Restraint Model", columnSpan = len ( results ) )
            table.Heading ( "" )
            table.Heading ( "" )
            for record in results: table.Heading ( record[0] )
            table.Entry   ( "Relative Energy", align = Align.Left, columnSpan = 2 )
            for record in results: table.Entry ( "{:.3f}".format ( record[1] ) )
            for i in range ( n ):
                if i == 0: table.Entry ( "Charge", align = Align.Left )
                else:      table.Entry ( "" )
                table.Entry ( atomNames[i] )
                for record in results: table.Entry ( "{:.3f}".format ( record[2][i] ) )
            if hasSpin:
                for i in range ( n ):
                    if i == 0: table.Entry ( "Spin", align = Align.Left )
                    else:      table.Entry ( "" )
                    table.Entry ( atomNames[i] )
                    for record in results: table.Entry ( "{:.3f}".format ( record[3][i] ) )
            table.Entry ( "Dipole", align = Align.Left )
            table.Entry ( "" )
            for record in results: table.Entry ( "{:.3f}".format ( record[4] ) )
            for ( i, cTag ) in enumerate ( cTags ):
                if i == 0: table.Entry ( "Restraints", align = Align.Left )
                else:      table.Entry ( "" )
                table.Entry ( "{:s}".format ( cTag ) )
                for record in results: table.Entry ( "{:.3f}".format ( record[5][cTag] ) )
            table.Stop ( )

# . Print out failures.
if len ( eFailures ) > 0:
    logFile.Paragraph ( "Only {:d} restrained energy calculations out of {:d} completed.".format ( eTests - len ( eFailures ), eTests ) )
    table     = logFile.GetTable ( columns = [ 50, 20 ] )
    table.Start   ( )
    table.Title   ( "Restrained Energy Failures" )
    table.Heading ( "System"   )
    table.Heading ( "QC Model" )
    name0 = ""
    for ( name, qcTag ) in eFailures:
        if name == name0:
            table.Entry ( "" )
        else:
            table.Entry ( name, align = Align.Left )
            name0 = name
        table.Entry ( qcTag )
    table.Stop ( )
else:
    logFile.Paragraph ( "All {:d} restrained energy calculations completed.".format ( eTests ) )
if _TestGradients:
    logFile.Paragraph ( "Maximum gradient deviation of completed gradient calculations = {:.3f}.".format ( maximumGradientDeviation ) )

# . Finish up.
logFile.Footer ( )
if len ( eFailures ) > 0: TestScriptExit_Fail ( )

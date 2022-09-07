"""Test for calculating the excited states of various simple molecules using the MOM method."""

# . The MOM method is sensitive (!) and often does not converge easily.
# . This can also happen during geometry optimizations. Those below all
# . converge except for the formaldehyde cases (minimizations are not done by default).
# . The accuracy of the gradients can be improved by using a density tolerance
# . of 1.0e-10. However, this leads to a couple of cases not converging.

import glob, math, os, os.path, sys

from Definitions          import dataPath                                 , \
                                 outPath                                  , \
                                 referenceDataPath
from pBabel               import ExportSystem                             , \
                                 ImportSystem          
from pCore                import Align                                    , \
                                 Clone                                    , \
                                 DataType                                 , \
                                 logFile                                  , \
                                 TestScriptExit_Fail                      , \
                                 YAMLUnpickle                                
from pMolecule            import SystemGeometryObjectiveFunction
from pMolecule.QCModel    import DIISSCFConverger                         , \
                                 ElectronicState                          , \
                                 OccupancyType                            , \
                                 QCModelDFT                               , \
                                 QCModelMNDO
from pScientific          import Units
from pScientific.Arrays   import Array                                    , \
                                 StorageType
from pScientific.Symmetry import Find3DGraphPointGroup                    , \
                                 IdentifyIrreducibleRepresentations       , \
                                 PrintIrreducibleRepresentations 
from pSimulation          import ConjugateGradientMinimize_SystemGeometry , \
                                 DetermineWavefunctionOverlap             , \
                                 MakeExcitedStateByOrbitalSwapping        , \
                                 SetUpMOMExcitedStateCalculation          , \
                                 TransitionDipoleMoment

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Tolerances.
_OrbitalDegeneracyTolerance =  1.0e-3 # . To distinguish different orbitals in symmetry determination (in Hartrees).

#===================================================================================================================================
# . Methods.
#===================================================================================================================================
def FindPointGroup ( system, log = logFile ):
    """Find the point group of the structure."""
    masses  = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
    numbers = [ atom.atomicNumber for atom in system.atoms ]
    report  = Find3DGraphPointGroup ( numbers                                ,
                                      system.coordinates3                    ,
                                      doCharacterSymmetryOperations = False  ,
                                      log                           = log    ,
                                      weights                       = masses )
    return report

#-----------------------------------------------------------------------------------------------------------------------------------
def OrbitalCharacterFunction ( system, useOrbitalsP = True ):
    """The orbital character function."""
    def GetItemCharacters ( operation ):
        """Get item characters."""
        return system.qcModel.OrbitalCharacters ( system, operation.transformationMatrix, operation.mapping, useOrbitalsP = useOrbitalsP )
    return GetItemCharacters

def OrbitalSymmetries ( system, pgReport, log = logFile ):
    """Identify orbital symmetries."""
    orbitalIRs = {}
    for ( attribute, flag ) in ( ( "orbitalsP", True ), ( "orbitalsQ", False ) ):
        orbitals = system.scratch.Get ( attribute, None )
        if orbitals is not None:
            energies = orbitals.energies
            GetItemCharacters   = OrbitalCharacterFunction ( system, useOrbitalsP = flag )
            ( iRs, characters ) = IdentifyIrreducibleRepresentations ( pgReport                    ,
                                                                       energies                    ,
                                                                       GetItemCharacters           ,
                                                                       _OrbitalDegeneracyTolerance ,
                                                                       maximumIRs = 100           ) # . d orbitals.
            iRs = [ iR.lower ( ) for iR in iRs ]
            PrintIrreducibleRepresentations ( pgReport              ,
                                              energies              ,
                                              iRs                   ,
                                              characters            ,
                                              itemName   = "{:s} Orbital".format ( attribute[-1] ),
                                              itemFormat = "{:.5f}" ,
                                              log        = log      ,
                                              valueName  = "Energy" )
            orbitalIRs[attribute[-1]] = [ ( iR, e ) for ( iR, e ) in zip ( iRs, energies ) ]
    return orbitalIRs

#===================================================================================================================================
# . System data.
#===================================================================================================================================
_SystemData = YAMLUnpickle ( os.path.join ( referenceDataPath, "MOMTestSystems.yaml" ) )

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_DensityTolerance  = 1.0e-08
_ForceUnrestricted = True
_Iterations        = 250
_MaximumAtoms      = 8
_TestGradients     = True
_TestMinimization  = False

# . Paths.
molPath = os.path.join ( os.path.join ( dataPath, "mol" ) )
if _TestMinimization:
    outPath = os.path.join ( outPath, "mom" )
    if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . QC models.
_OrbitalBasis = "def2-sv(p)"
_Converger    = DIISSCFConverger.WithOptions ( densityTolerance  = _DensityTolerance ,
                                               maximumIterations = _Iterations       )
_QCModels     = ( ( "MNDO", QCModelMNDO.WithOptions ( converger         = _Converger        ,
                                                      hamiltonian       = "am1"             ) ) ,
                  ( "HF"  , QCModelDFT.WithOptions  ( converger         = _Converger        ,
                                                      functional        = "hf"              ,
                                                      maximumMemory     = 6.0               ,
                                                      orbitalBasis      = _OrbitalBasis     ) ) )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initial pass over systems to order them.
systems = {}
toOrder = []
for ( name, value ) in _SystemData.items ( ):
    # . Set up system.
    charge            = value["Charge"       ]
    multiplicity      = value["Multiplicity" ]
    swaps             = value["Orbital Swaps"]
    isSpinRestricted  = ( multiplicity == 1 ) and ( not _ForceUnrestricted )
    system            = ImportSystem ( os.path.join ( molPath, value["Structure File"] ), log = None )
    system.label      = name
    electronicStateGS = ElectronicState.WithOptions ( charge           = charge                 ,
                                                      isSpinRestricted = isSpinRestricted       ,
                                                      multiplicity     = multiplicity           ,
                                                      occupancyType    = OccupancyType.Cardinal )
    electronicStateES = ElectronicState.WithOptions ( charge           = charge                 ,
                                                      isSpinRestricted = isSpinRestricted       ,
                                                      multiplicity     = multiplicity           ,
                                                      occupancyType    = OccupancyType.MOM      )
    # . Find point group.
    pgReport   = FindPointGroup ( system, log = None )
    pointGroup = pgReport.get ( "Point Group", None )
    if pointGroup is None:
        report = None
        logFile.Paragraph ( "System \"{:s}\" read with undefined point group.".format ( name ) )
    else:
        report = pgReport
        logFile.Paragraph ( "System \"{:s}\" read with {:s} point group.".format ( name, pointGroup.label ) )
    # . Save data.
    systems[name] = ( system, electronicStateGS, electronicStateES, report, swaps )
    toOrder.append ( ( len ( system.atoms ), sum ( [ atom.mass for atom in system.atoms ] ), name ) )

# . Loop over ordered systems.
eFailures                = []
eTests                   = 0
gFailures                = 0
gTests                   = 0
maximumGradientDeviation = 0.0
mFailures                = 0
mTests                   = 0
for ( _, _, name ) in sorted ( toOrder ):
    ( system, electronicStateGS, electronicStateES, pgReport, swaps ) = systems[name]
    coordinates3 = Clone ( system.coordinates3 ) # . Save the original coordinates.

    # . Loop over QC models.
    for ( qcTag, qcModel ) in _QCModels:

        # . Ground state calculation.
        hasSpin = ( not electronicStateGS.isSpinRestricted )
        system.electronicState = electronicStateGS
        system.DefineQCModel ( qcModel )
        system.Summary ( )
        energy  = system.Energy        ( )
        charges = system.AtomicCharges ( )
        dipole  = system.DipoleMoment  ( )
        if hasSpin: spins = system.qcModel.AtomicSpins ( system )
        else:       spins = None
        if pgReport is not None: orbitalIRs = OrbitalSymmetries ( system, pgReport, log = None )
        else:                    orbitalIRs = None
        results1 = [ ( "Ground State", energy, charges, spins, dipole.Norm2 ( ), orbitalIRs ) ]
        if hasSpin:
            wavefunctions = [ ( Clone ( system.scratch.orbitalsP ) ,
                                Clone ( system.scratch.orbitalsQ ) ) ]
        else:
            wavefunctions = [   Clone ( system.scratch.orbitalsP )   ]

        # . Excited state calculations.
        SetUpMOMExcitedStateCalculation ( system, electronicStateES )

        # . Loop over swaps making sure swapping is done on the GS orbitals!
        for ( i, swap ) in enumerate ( swaps ):
            esTag  = "Excited State {:d}".format ( i )
            swapsP = swap.get ( "SwapsP", None )
            swapsQ = swap.get ( "SwapsQ", None )
            if hasSpin:
                wavefunctions[0][0].orbitals.CopyTo ( system.scratch.orbitalsP.orbitals )
                wavefunctions[0][1].orbitals.CopyTo ( system.scratch.orbitalsQ.orbitals )
            else:
                wavefunctions[0].orbitals.CopyTo ( system.scratch.orbitalsP.orbitals )
            MakeExcitedStateByOrbitalSwapping ( system, swapsP = swapsP, swapsQ = swapsQ )
            eTests += 1
            try:
                energy  = system.Energy        ( )
                charges = system.AtomicCharges ( )
                dipole  = system.DipoleMoment  ( )
                if system.electronicState.isSpinRestricted: spins = None
                else:                                       spins = system.qcModel.AtomicSpins ( system )
                if pgReport is not None: orbitalIRs = OrbitalSymmetries ( system, pgReport, log = None )
                else:                    orbitalIRs = None
                results1.append ( ( esTag, energy, charges, spins, dipole.Norm2 ( ), orbitalIRs ) )
                if hasSpin: wavefunctions.append ( ( Clone ( system.scratch.orbitalsP ) ,
                                                     Clone ( system.scratch.orbitalsQ ) ) )
                else:       wavefunctions.append (   Clone ( system.scratch.orbitalsP )   )
                if len ( system.atoms ) <= _MaximumAtoms:
                    if _TestGradients:
                        gTests += 1
                        try:
                            of = SystemGeometryObjectiveFunction.FromSystem ( system )
                            gradientDeviation    = of.TestGradients ( delta = 1.0e-4 )
                            maximumGradientDeviation = max ( maximumGradientDeviation, gradientDeviation )
                        except:
                            gFailures += 1
                    if _TestMinimization:
                        mTests += 1
                        try:
                            cgReport = ConjugateGradientMinimize_SystemGeometry ( system                      ,
                                                                                  logFrequency         =  100 ,
                                                                                  maximumIterations    = 2000 ,
                                                                                  rmsGradientTolerance =  0.5 )
                            if cgReport["Converged"]:
                                energyM = system.Energy ( doGradients = True )
                                masses  = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
                                system.coordinates3.Superimpose ( coordinates3, weights = masses )
                                rms = system.coordinates3.RootMeanSquareDeviation ( coordinates3, weights = masses )
                                logFile.Paragraph ( "{:s} Energy difference and RMS coordinate deviation after minimization for {:s}/{:s} = {:.1f} and {:.3f}.".format ( qcTag, name, esTag, energyM - energy, rms ) )
                                ExportSystem ( os.path.join ( outPath, "{:s}_{:s}_es{:d}.xyz".format ( name.lower ( ).replace ( " ", "_" ), qcTag.lower ( ), i ) ), system )
                            else:
                                mFailures += 1
                        except:
                            mFailures += 1
                        coordinates3.CopyTo ( system.coordinates3 ) # . Copy back the original coordinates.
            except:
                eFailures.append ( ( name, qcTag, esTag ) )

        # . Determine overlaps, transition dipole moments and oscillator strengths.
        # . All systems neutral so center argument to TransitionDipoleMoment method is ignored.
        overlap  = system.scratch.Get ( "overlapMatrix", None )
        results2 = []
        for i in range ( len ( wavefunctions ) ):
            name1 = results1[i][0]
            e1    = results1[i][1]
            for j in range ( i, len ( wavefunctions ) ):
                name2 = results1[j][0]
                dE12  = ( results1[j][1] - e1 ) / Units.Energy_Electron_Volts_To_Kilojoules_Per_Mole
                if dE12 < 0.0:
                    ( name1, name2 ) = ( name2, name1 )
                    dE12 *= -1.0
                #O12 = DetermineWavefunctionOverlap ( wavefunctions[i], wavefunctions[j], overlap = overlap )
                ( O12, tDipole, fOsc ) = TransitionDipoleMoment ( system, dE12, wavefunctions[i], wavefunctions[j], doCheck = False )
                results2.append ( ( name1, name2, dE12, O12, tDipole[0], tDipole[1], tDipole[2], fOsc ) )

        # . Output the results for one state.
        # . Energies and properties.
        atomNames = [ atom.path for atom in system.atoms ]
        n         = len ( atomNames )
        columns   = [ 10, 10 ] + len ( results1 ) * [ 20 ]
        table     = logFile.GetTable ( columns = columns )
        table.Start   ( )
        table.Title   ( "{:s} Single State Results for {:s}".format ( qcTag, name ) )
        table.Heading ( "Property", columnSpan = 2 )
        table.Heading ( "Model"   , columnSpan = len ( results1 ) )
        table.Heading ( "" )
        table.Heading ( "" )
        for record in results1: table.Heading ( record[0] )
        table.Entry   ( "Energy" )
        table.Entry   ( "" )
        for record in results1: table.Entry ( "{:.3f}".format ( record[1] ) )
        for i in range ( n ):
            if i == 0: table.Entry ( "Charge" )
            else:      table.Entry ( "" )
            table.Entry ( atomNames[i] )
            for record in results1: table.Entry ( "{:.3f}".format ( record[2][i] ) )
        if not system.electronicState.isSpinRestricted:
            for i in range ( n ):
                if i == 0: table.Entry ( "Spin" )
                else:      table.Entry ( "" )
                table.Entry ( atomNames[i] )
                for record in results1: table.Entry ( "{:.3f}".format ( record[3][i] ) )
        table.Entry ( "Dipole" )
        table.Entry ( "" )
        for record in results1: table.Entry ( "{:.3f}".format ( record[4] ) )
        table.Stop ( )

        # . Orbital characters.
        toProcess = [ ( tag, IRs ) for ( tag, _, _, _, _, IRs ) in results1 if IRs is not None ]
        if len ( toProcess ) > 0:
            columns   = [ 5, 5 ] + len ( toProcess ) * [ 10, 10 ]
            table     = logFile.GetTable ( columns = columns )
            table.Start   ( )
            table.Title   ( "{:s} Orbital Irreducible Representations and Energies for {:s}".format ( qcTag, name ) )
            table.Heading ( "Orbitals", columnSpan = 2 )
            for ( tag, _ ) in toProcess: table.Heading ( tag, columnSpan = 2 )
            table.Heading ( "" )
            table.Heading ( "" )
            for i in range ( len ( toProcess ) ):
                table.Heading ( "IR"     )
                table.Heading ( "Energy" )
            tags = [ "P" ]
            if hasSpin: tags.append ( "Q" )
            for tag in tags:
                data = [ IRs.get ( tag, None ) for ( _, IRs ) in toProcess ]
                if all ( [ datum is not None for datum in data ] ):
                    n = len ( data[0] )
                    for i in range ( n ):
                        if i == 0: table.Entry ( tag )
                        else:      table.Entry ( ""  )
                        table.Entry ( "{:d}".format ( i ) )
                        for IRs in data:
                            ( IR, e ) = IRs[i]
                            table.Entry ( "{:s}".format   ( IR ) )
                            table.Entry ( "{:.5f}".format ( e  ) )
            table.Stop ( )

        # . Output the results for two states.
        columns   = [ 20, 20 ] + 6 * [ 10 ]
        table     = logFile.GetTable ( columns = columns )
        table.Start   ( )
        table.Title   ( "{:s} Two State Results for {:s}".format ( qcTag, name ) )
        table.Heading ( "States", columnSpan = 2 )
        table.Heading ( "deltaE (eV)"  )
        table.Heading ( "Overlap" )
        table.Heading ( "Transition Dipole (D)", columnSpan = 3 )
        table.Heading ( "fOsc"    )
        for ( name1, name2, dE12, O12, tX, tY, tZ, fOsc ) in sorted ( results2 ):
            table.Entry ( name1, align = Align.Left )
            table.Entry ( name2, align = Align.Left )
            table.Entry ( "{:.5f}".format ( dE12 ) )
            table.Entry ( "{:.3f}".format ( O12  ) )
            table.Entry ( "{:.3f}".format ( tX   ) )
            table.Entry ( "{:.3f}".format ( tY   ) )
            table.Entry ( "{:.3f}".format ( tZ   ) )
            table.Entry ( "{:.3f}".format ( fOsc ) )
        table.Stop ( )

# . Print out failures.
if len ( eFailures ) > 0:
    logFile.Paragraph ( "Only {:d} excited state energy calculations out of {:d} completed.".format ( eTests - len ( eFailures ), eTests ) )
    table     = logFile.GetTable ( columns = [ 20, 10, 20 ] )
    table.Start   ( )
    table.Title   ( "Energy Failures" )
    table.Heading ( "System" )
    table.Heading ( "Method" )
    table.Heading ( "State"  )
    name0 = ""
    for ( name, qcTag, state ) in eFailures:
        if name == name0:
            table.Entry ( "" )
        else:
            table.Entry ( name, align = Align.Left )
            name0 = name
        table.Entry ( qcTag )
        table.Entry ( state )
    table.Stop ( )
else:
    logFile.Paragraph ( "All {:d} excited state energy calculations completed.".format ( eTests ) )
if _TestGradients:
    if gFailures > 0: logFile.Paragraph ( "Only {:d} excited state gradient calculations out of {:d} completed.".format ( gTests - gFailures, gTests ) )
    else:             logFile.Paragraph ( "All {:d} excited state gradient calculations completed.".format ( gTests ) )
    logFile.Paragraph ( "Maximum gradient deviation of completed gradient calculations = {:.3f}.".format ( maximumGradientDeviation ) )
if _TestMinimization:
    if mFailures > 0: logFile.Paragraph ( "Only {:d} excited state minimization calculations out of {:d} completed.".format ( mTests - mFailures, mTests ) )
    else:             logFile.Paragraph ( "All {:d} excited state minimization calculations completed.".format ( mTests ) )

# . Finish up.
logFile.Footer ( )
if ( len ( eFailures ) + gFailures + mFailures ) > 0: TestScriptExit_Fail ( )

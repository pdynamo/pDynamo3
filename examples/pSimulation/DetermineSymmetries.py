"""Find CIP labels for several bALA configurations."""

import os, os.path

from Definitions          import structuresPath
from pBabel               import ImportSystem
from pCore                import logFile                                , \
                                 TestScriptExit_Fail
from pMolecule.QCModel    import CIDiagonalization                      , \
                                 CIMethod                               , \
                                 DIISSCFConverger                       , \
                                 QCModelDFT                             , \
                                 QCModelMNDO                            , \
                                 QCModelMNDOCI
from pScientific.Arrays   import Array                                  , \
                                 ArrayPrint
from pScientific.Symmetry import Find3DGraphPointGroup                  , \
                                 IdentifyIrreducibleRepresentations     , \
                                 PrintIrreducibleRepresentations 
from pSimulation          import NormalModes_IrreducibleRepresentations , \
                                 NormalModes_SystemGeometry

# . The test is relatively sensitive to tolerances.
# . The higher virtual orbitals of DFT/HF wavefunctions are also often not identified (and so excluded from pass/fail test).

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Tolerances.
_OrbitalDegeneracyTolerance    =  1.0e-3 # . Hartrees.
_NormalModeDegeneracyTolerance =  5.0    # . cm^-1.
_StateDegeneracyTolerance      = 20.0    # . kJ/mol.

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
def OrbitalCharacterFunction ( system ):
    """The orbital character function."""
    def GetItemCharacters ( operation ):
        """Get item characters."""
        return system.qcModel.OrbitalCharacters ( system, operation.transformationMatrix, operation.mapping )
    return GetItemCharacters

def OrbitalSymmetries ( system, pgReport, log = logFile ):
    """Identify orbital symmetries."""
    GetItemCharacters   = OrbitalCharacterFunction ( system )
    ( iRs, characters ) = IdentifyIrreducibleRepresentations ( pgReport                          ,
                                                               system.scratch.orbitalsP.energies ,
                                                               GetItemCharacters                 ,
                                                               _OrbitalDegeneracyTolerance       ,
                                                               maximumIRs = 6                    ) # . Cartesian d orbitals.
    iRs = [ iR.lower ( ) for iR in iRs ]
    PrintIrreducibleRepresentations ( pgReport                          ,
                                      system.scratch.orbitalsP.energies ,
                                      iRs                               ,
                                      characters                        ,
                                      itemName   = "Orbital"            ,
                                      itemFormat = "{:.5f}"             ,
                                      log        = logFile              ,
                                      valueName  = "Energy"             )
    return iRs

#-----------------------------------------------------------------------------------------------------------------------------------
def StateCharacterFunction ( system ):
    """The state character function."""
    def GetItemCharacters ( operation ):
        """Get item characters."""
        return system.qcModel.StateCharacters ( system, operation.transformationMatrix, operation.mapping )
    return GetItemCharacters

def StateSymmetries ( system, pgReport, log = logFile ):
    """Identify state symmetries."""
    GetItemCharacters   = StateCharacterFunction ( system )
    ( iRs, characters ) = IdentifyIrreducibleRepresentations ( pgReport                     ,
                                                               system.scratch.ci.ciEnergies ,
                                                               GetItemCharacters            ,
                                                               _StateDegeneracyTolerance    )
    PrintIrreducibleRepresentations ( pgReport                     ,
                                      system.scratch.ci.ciEnergies ,
                                      iRs                          ,
                                      characters                   ,
                                      itemName   = "State"         ,
                                      itemFormat = "{:.3f}"        ,
                                      log        = logFile         ,
                                      valueName  = "Energy"        )
    return iRs

#===================================================================================================================================
# . Systems and QC models.
#===================================================================================================================================
# . Systems.
_XYZFiles = ( "C2", "C2v", "C3h", "D6h", "Oh_a" )

# . SCF converger.
convergerCI = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-14, maximumIterations = 250 )
convergerHF = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-12, maximumIterations = 250 )

# . CI space - electrons, then orbitals.
# . Care must be taken here as symmetry is easily broken.
# . E.g. if degenerate orbitals of the same IR are split between different spaces,
#        if not all degenerate states of a given IR are not included in the number of states resolved,
#        with some other choices of active space and CI method. 
_ActiveSpace   = { "C2"   : ( 6, 6 ) ,
                   "C2v"  : ( 6, 6 ) ,
                   "C3h"  : ( 8, 8 ) ,
                   "D6h"  : ( 8, 7 ) ,
                   "Oh_a" : ( 4, 6 ) }

# . Options for QC models.
_HFOptions     = { "converger"           : convergerHF             ,
                   "orbitalBasis"        : "631gs"                 }
_MNDOCIOptions = { "ciDiagonalization"   : CIDiagonalization.Dense ,
                   "ciMethod"            : CIMethod.Full           ,
                   "converger"           : convergerCI             ,
                   "multiplicity"        : 3                       ,
                   "minimalMultiplicity" : 3                       ,
                   "numberOfStates"      : 9                       ,
                   "requiredRoot"        : 1                       }
_MNDOOptions   = { "converger"           : convergerHF             ,
                   "hamiltonian"         : "pm6"                   }

# . QC models.
_QCModels = ( ( QCModelMNDO   , _MNDOOptions   , True , True , 1000, False ) ,
              ( QCModelMNDOCI , _MNDOCIOptions , False, False, 1000, False ) ,
              ( QCModelDFT    , _HFOptions     , True , True ,    5, True  ) )

#===================================================================================================================================
# . Calculation.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
xyzRoot  = os.path.join ( structuresPath, "pointGroupExamples", "xyz" )

# . Loop over paths and QC models.
nBadHF    = 0
nFailedIR = 0
nFailedPG = 0
for xyzFile in _XYZFiles:

    # . Read system.
    xyzPath       = os.path.join ( xyzRoot, xyzFile + ".xyz" )
    system        = ImportSystem ( xyzPath )
    numberOfAtoms = len ( system.atoms )

    # . Find point group.
    pgReport   = FindPointGroup ( system, log = None )
    pointGroup = pgReport.get ( "Point Group", None )
    if pointGroup is None: outGroup = None
    else:                  outGroup = pointGroup.label
    if outGroup is None:
        nFailedPG += 1
    else:
        logFile.Paragraph ( "Point group for system \"{:s}\" is {:s}.".format ( xyzFile, outGroup ) )

        #. Do QC calculations.
        for ( qcModelClass, qcModelOptions, doModes, doOrbitals, maximumAtoms, checkVirtuals ) in _QCModels:
            isCI = ( qcModelClass == QCModelMNDOCI )
            if isCI:
                qcModelOptions = qcModelOptions.copy ( )
                qcModelOptions["activeElectrons"] = _ActiveSpace[xyzFile][0]
                qcModelOptions["activeOrbitals" ] = _ActiveSpace[xyzFile][1]
            qcModel = qcModelClass.WithOptions ( **qcModelOptions )
            system.DefineQCModel ( qcModel )
            system.Summary ( )
            system.Energy ( doGradients = True )
            if isCI:
                stateIRs = StateSymmetries ( system, pgReport )
                nFailedIR += stateIRs.count ( "?" )
            if doModes and ( numberOfAtoms <= maximumAtoms ):
                NormalModes_SystemGeometry ( system )
                modeIRs    = NormalModes_IrreducibleRepresentations ( system, degeneracyTolerance = _NormalModeDegeneracyTolerance, results = pgReport )
                nFailedIR += modeIRs.count ( "?" )
            if doOrbitals:
                orbitalIRs = OrbitalSymmetries ( system, pgReport )
                if checkVirtuals:
                    orbitals   = system.scratch.orbitalsP
                    n          = orbitals.numberOrbitals
                    m          = orbitals.occupancyHandler.numberOccupied
                    nFailedIR += orbitalIRs[0:m].count ( "?" )
                    nBadHF    += orbitalIRs[m:n].count ( "?" )
                else:
                    nFailedIR += orbitalIRs.count ( "?" )

# . Print unassigned virtuals.
if nBadHF > 0: logFile.Paragraph ( "There were {:d} unassigned virtual orbitals.".format ( nBadHF ) )

# . Success/failure.
isOK = ( nFailedIR == 0 ) and ( nFailedPG == 0 )
if not isOK:
    string = "There were"
    if nFailedIR > 0:
        string += " {:d} failed irreducible representation".format ( nFailedIR )
        if nFailedPG > 0: string += " and"
    if nFailedPG > 0:
        string += " {:d} failed point group".format ( nFailedPG )
    string += " assignments."
    logFile.Paragraph ( string )
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )

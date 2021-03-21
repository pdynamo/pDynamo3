"""Test for calculating dihydrogen dissociation curves."""

import glob, math, os, os.path

from Definitions       import dataPath         , \
                              _FullVerificationSummary
from pBabel            import ImportSystem
from pCore             import logFile          , \
                              TestDataSet      , \
                              TestReal         , \
                              TestScriptExit_Fail
from pMolecule.QCModel import CIMethod         , \
                              ElectronicState  , \
                              QCModelMNDO      , \
                              QCModelMNDOCI
from pScientific       import Units

# . UHF should give a lower energy than RHF. Problem with initial guess?

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Distance parameters.
_NumberIncrements = 59
_XIncrement       = 0.1
_XStart           = 0.2

# . QC models.
# . MNDO is used because it gives good results.
_UnrestrictedElectronicState = ElectronicState.WithOptions ( isSpinRestricted = False )
_RestrictedQCModel           = QCModelMNDO.WithOptions   ( hamiltonian      = "mndo" )
_UnrestrictedQCModel         = QCModelMNDO.WithOptions   ( hamiltonian      = "mndo" )
_S0CIQCModel                 = QCModelMNDOCI.WithOptions ( hamiltonian      = "mndo" ,
                                                           ciMethod         = CIMethod.Full ,
                                                           activeElectrons  = 2      ,
                                                           activeOrbitals   = 2      ,
                                                           multiplicity     = 1      ,
                                                           requiredRoot     = 0      )
_S1CIQCModel                 = QCModelMNDOCI.WithOptions ( hamiltonian      = "mndo" ,
                                                           ciMethod         = CIMethod.Full ,
                                                           activeElectrons  = 2      ,
                                                           activeOrbitals   = 2      ,
                                                           multiplicity     = 1      ,
                                                           requiredRoot     = 1      )
_T1CIQCModel                 = QCModelMNDOCI.WithOptions ( hamiltonian      = "mndo" ,
                                                           ciMethod         = CIMethod.Full ,
                                                           activeElectrons  = 2      ,
                                                           activeOrbitals   = 2      ,
                                                           multiplicity     = 3      ,
                                                           requiredRoot     = 0      )
_QCModels = ( ( "RHF"  , _RestrictedQCModel  , None ) ,
              ( "UHF"  , _UnrestrictedQCModel, _UnrestrictedElectronicState ) ,
              ( "CI S0", _S0CIQCModel        , None ) ,
              ( "CI S1", _S1CIQCModel        , None ) ,
              ( "CI T1", _T1CIQCModel        , None ) )

# . Reference values.
# . The experimental values for the ground state are ~ 4.52 eV with a minimum of ~ 0.74 A.
_DissociationEnergy = 4.52
_MinimumDistance    = 0.74

# . Tolerances.
_DistanceTolerance = 0.1
_EnergyTolerance   = 0.1

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPath, "xyz" )

# . Initialize the distances.
distances = [ _XStart + float ( i ) * _XIncrement for i in range ( _NumberIncrements ) ]

# . Loop over the QC models.
results = {}
for ( label, qcModel, electronicState ) in _QCModels:

    # . Define the system.
    molecule = ImportSystem ( os.path.join ( dataPath, "dihydrogen.xyz" ) )
    if electronicState is not None: molecule.electronicState = electronicState
    molecule.DefineQCModel ( qcModel    )
    molecule.Summary       ( log = None )

    # . Initialize the coordinates.
    molecule.coordinates3.Set ( 0.0 )

    # . Calculate energies at different distances.
    energies = []
    for distance in distances:
        molecule.coordinates3[0,1] = distance
        energies.append ( molecule.Energy ( log = None ) / Units.Energy_Electron_Volts_To_Kilojoules_Per_Mole )
    minimumEnergy = min ( energies )
    for i in range ( len ( energies ) ):
        energies[i] -= minimumEnergy

    # . Save the results.
    results[label] = energies

# . Print the distances.
table = logFile.GetTable ( columns = [ 10 ] + len ( results ) * [ 20 ] )
table.Start ( )
table.Title ( "Dissociation Curves (Angstroms/eV)" )
table.Heading ( "Distance" )
labels = list ( results.keys ( ) )
labels.sort ( )
for label in labels:
    table.Heading ( label )
for i in range ( len ( distances ) ):
    table.Entry ( "{:4.1f}".format ( distances[i] ) )
    for label in labels:
        table.Entry ( "{:20.2f}".format ( results[label][i] ) )
table.Stop ( )

# . Get the observed data.
energies           = results["CI S0"]
dissociationEnergy = energies[-1]
minimumDifference  = dissociationEnergy
minimumDistance    = distances[0]
minimumEnergy      = min ( energies )
for ( i, ( distance, energy ) ) in enumerate ( zip ( distances, energies ) ):
    energyDifference = math.fabs ( energy - minimumEnergy )
    if ( energyDifference < minimumDifference ):
        minimumDifference = energyDifference
        minimumDistance   = distance
observed = { "Dissociation Energy" : dissociationEnergy, "Minimum Distance" : minimumDistance }

# . Generate reference data.
referenceData = TestDataSet.WithOptions ( label = "Dihydrogen Dissociation" )
referenceData.AddDatum ( TestReal.WithOptions ( label = "Dissociation Energy", value = _DissociationEnergy, parent = referenceData, absoluteErrorTolerance = _EnergyTolerance  , toleranceFormat = "{:.2f}", valueFormat = "{:.2f}" ) )
referenceData.AddDatum ( TestReal.WithOptions ( label = "Minimum Distance"   , value = _MinimumDistance   , parent = referenceData, absoluteErrorTolerance = _DistanceTolerance, toleranceFormat = "{:.2f}", valueFormat = "{:.2f}" ) )

# . Footer.
results = referenceData.VerifyAgainst ( observed )
results.Summary ( fullSummary = _FullVerificationSummary )
isOK    = results.WasSuccessful ( )
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )

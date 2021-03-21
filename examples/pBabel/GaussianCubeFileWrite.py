"""Test for writing Gaussian cube files."""

import glob, os, os.path

from Definitions        import dataPathM              , \
                               outPath
from pBabel             import GaussianCubeFileWriter , \
                               ImportSystem
from pCore              import logFile                , \
                               TestScriptExit_Fail
from pMolecule.QCModel  import ElectronicState        , \
                               QCModelDFT             , \
                               QCModelMNDO
from pScientific.Arrays import ArrayPrint

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The destination for results.
_Destination = "gaussianCubeFiles"

# . The file extension.
_Extension = ".cube"

# . Grid spacing.
_GridSpacing = 0.3

# . Energy models - use defaults.
_QCModels = [ ( "pm6", QCModelMNDO.WithOptions ( hamiltonian = "pm6" ) ) ,
              ( "hf" , QCModelDFT.WithDefaults ( )                     ) ]

# . The system.
_SystemLabel = "F-:CH3Cl complex"
_SystemPath  = "fch3cl"

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Initialization.
isOK = True
logFile.Header ( )

# . Output setup.
dataPath = os.path.join ( dataPathM, "xyz"        )
outPath  = os.path.join ( outPath  , _Destination )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Loop over energy models.
for ( qcLabel, qcModel ) in _QCModels:

    # . Read the system.
    system                 = ImportSystem ( os.path.join ( dataPath, _SystemPath + ".xyz" ) )
    system.label           = "{:s}".format ( _SystemLabel )
    system.electronicState = ElectronicState.WithOptions ( charge = -1 )
    system.DefineQCModel ( qcModel )
    system.Summary ( )
    system.Energy  ( )

    # . Orbital data.
    orbitalsP = system.scratch.orbitalsP
    energies  = orbitalsP.energies
    LUMO      = orbitalsP.occupancyHandler.numberOccupied
    HOMO      = LUMO - 1
    if energies is not None: ArrayPrint ( energies, itemFormat = "{:.5f}", title = "Orbital Energies" )
    logFile.Paragraph ( "HOMO and LUMO indices = {:d} and {:d}.".format ( HOMO, LUMO ) )

    # . Write out the cube files.
    path = os.path.join ( outPath, _SystemPath + "_" + qcLabel + "_" + "{:s}" + _Extension )
    GaussianCubeFileWriter.PathFromSystemDensity   ( path.format ( "density"   ), system,                 gridSpacing = _GridSpacing )
    GaussianCubeFileWriter.PathFromSystemOrbitals  ( path.format ( "orbitals"  ), system, [ HOMO, LUMO ], gridSpacing = _GridSpacing )
    GaussianCubeFileWriter.PathFromSystemPotential ( path.format ( "potential" ), system,                 gridSpacing = _GridSpacing )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )

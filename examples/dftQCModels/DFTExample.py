"""Modification of Example 5 that includes some DFT QC Models."""

import os, os.path

from pBabel            import ImportSystem
from pCore             import Align                           , \
                              logFile                         , \
                              TestScript_InputDataPath
from pMolecule         import SystemGeometryObjectiveFunction
from pMolecule.QCModel import DIISSCFConverger                , \
                              ElectronicState                 , \
                              QCModelDFT                      , \
                              QCModelMNDO                     , \
                              DFTGridAccuracy                 , \
                              DFTGridIntegrator

# . The DFT QC Model is QCModelDFT. Its constructor takes the following arguments:
#
#   converger        - the SCF converger (currently only DIIS).
#   fitBasis     - the basis set for Coulomb fitting. The default is "demon".
#   functional       - the DFT functional. Many possibilities with "hf" the default.
#   gridIntegrator   - the grid integrator (currently only DFTGridIntegrator).
#   orbitalBasis     - the orbital basis set. The default is "321g".
#
#   DFTGridIntegrator options:
#
#   accuracy         - instance of DFTGridAccuracy enum, either VeryLow, Low, Medium (default), High or VeryHigh.
#   inCore           - recalculate or not (default) the basis function values at the grid points.
#                      True takes a lot of storage!
#
# . Good combinations of orbital and density basis sets are (in increasing order of cost):
#
#   321g/demon; 631gs/ahlrichs; svp/weigend.
#
# . To perform QC/MM calculations, proceed in the usual way. Any of the NB models that can be used with
#   QCModelMNDO can also be used with QCModelDFT.
#

# . Start.
logFile.Header ( )

# . Options.
_TestGradients = False

# . Define the SCF converger.
converger      = DIISSCFConverger.WithOptions  ( densityTolerance = 1.0e-10, maximumIterations = 250 )
gridIntegrator = DFTGridIntegrator.WithOptions ( accuracy = DFTGridAccuracy.Medium, inCore = True )

# . Define the energy models.
_EnergyModels = [ ( "am1"                 , QCModelMNDO.WithOptions ( hamiltonian = "am1" ) ),
                  ( "lda/321g/demon"      , QCModelDFT.WithOptions  ( converger = converger, fitBasis = "demon"   , functional = "lda"  , gridIntegrator = gridIntegrator, orbitalBasis = "321g"  ) ) ,
                  ( "blyp/321g/demon"     , QCModelDFT.WithOptions  ( converger = converger, fitBasis = "demon"   , functional = "blyp" , gridIntegrator = gridIntegrator, orbitalBasis = "321g"  ) ) ,
                  ( "b3lyp/321g"          , QCModelDFT.WithOptions  ( converger = converger,                        functional = "b3lyp", gridIntegrator = gridIntegrator, orbitalBasis = "321g"  ) ) ,
                  ( "hf/321g"             , QCModelDFT.WithOptions  ( converger = converger,                        functional = "hf"   ,                                  orbitalBasis = "321g"  ) ) ,
                  ( "lda/631gs/ahlrichs"  , QCModelDFT.WithOptions  ( converger = converger, fitBasis = "ahlrichs", functional = "lda"  , gridIntegrator = gridIntegrator, orbitalBasis = "631gs" ) ) ,
                  ( "blyp/631gs/ahlrichs" , QCModelDFT.WithOptions  ( converger = converger, fitBasis = "ahlrichs", functional = "blyp" , gridIntegrator = gridIntegrator, orbitalBasis = "631gs" ) ) ,
                  ( "b3lyp/631gs"         , QCModelDFT.WithOptions  ( converger = converger,                        functional = "b3lyp", gridIntegrator = gridIntegrator, orbitalBasis = "631gs" ) ) ,
                  ( "hf/631gs"            , QCModelDFT.WithOptions  ( converger = converger,                        functional = "hf"   ,                                  orbitalBasis = "631gs" ) ) ,
                  ( "lda/svp/weigend"     , QCModelDFT.WithOptions  ( converger = converger, fitBasis = "weigend" , functional = "lda"  , gridIntegrator = gridIntegrator, orbitalBasis = "svp"   ) ) ,
                  ( "blyp/svp/weigend"    , QCModelDFT.WithOptions  ( converger = converger, fitBasis = "weigend" , functional = "blyp" , gridIntegrator = gridIntegrator, orbitalBasis = "svp"   ) ) ,
                  ( "b3lyp/svp"           , QCModelDFT.WithOptions  ( converger = converger,                        functional = "b3lyp", gridIntegrator = gridIntegrator, orbitalBasis = "svp"   ) ) ,
                  ( "hf/svp"              , QCModelDFT.WithOptions  ( converger = converger,                        functional = "hf"   ,                                  orbitalBasis = "svp"   ) ) ]

# . Get the fileName.
_Names  = ( ( "chloride", -1 ), ( "dinitrogen", 0 ), ( "water", 0 ), ( "methane", 0 ), ( "formaldehyde", 0 ), ( "ethene", 0 ), ( "ozone", 0 ), ( "fch3cl", -1 ), ( "nitromethane", 0 ), ( "glycine", 0 ) )
_Source = os.path.join ( TestScript_InputDataPath ( "pMolecule" ), "xyz" )

# . Loop over the energy models.
dG      = None
nQ      = 0
results = []
for ( name, charge ) in _Names:
    fileName = os.path.join ( _Source, "{:s}.xyz".format ( name ) )
    for ( label, model ) in _EnergyModels:
        system = ImportSystem ( fileName )
        system.electronicState = ElectronicState.WithOptions ( charge = charge )
        system.DefineQCModel ( model )
        system.Summary ( )
        energy  = system.Energy ( )
        charges = system.AtomicCharges ( )
        dipole  = system.DipoleMoment  ( )
        nQ      = max ( nQ, len ( system.atoms ) )
        if _TestGradients:
            of = SystemGeometryObjectiveFunction.FromSystem ( system )
            dG = of.TestGradients ( )
        results.append ( ( name, label, energy, list ( charges ), list ( dipole ), dipole.Norm2 ( ), dG ) )

# . Output the results.
nColumns = nQ + 7
if _TestGradients: nColumns += 1
table    = logFile.GetTable ( columns = nColumns * [ 20 ] )
table.Start  ( )
table.Title  ( "Energy Model Results" )
table.Heading ( "System" )
table.Heading ( "Model"  )
table.Heading ( "Energy" )
table.Heading ( "Charges"       , columnSpan = nQ )
table.Heading ( "Dipole Vector" , columnSpan =  3 )
table.Heading ( "Dipole Norm"   )
if _TestGradients: table.Heading ( "Gradient Error" )
for ( name, label, energy, charges, dipole, dipoleNorm, dG ) in results:
    table.Entry ( name , align = Align.Left )
    table.Entry ( label, align = Align.Left )
    table.Entry ( "{:.1f}".format ( energy ) )
    for charge in charges: table.Entry ( "{:.3f}".format ( charge ) )
    for i in range ( len ( charges ), nQ ): table.Entry ( "-" )
    for value in dipole: table.Entry ( "{:.3f}".format ( value ) )
    table.Entry ( "{:.3f}".format ( dipoleNorm ) )
    if _TestGradients: table.Entry ( "{:.5f}".format ( dG ) )
table.Stop ( )

# . Stop.
logFile.Footer ( )

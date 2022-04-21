"""Modification of Example 5 that uses a selection of different DFT QC Models."""

# . These currently work fine except for cases with G-functions and higher (problem with Rys polynomials)!

import os, os.path

from pBabel                          import ImportSystem
from pCore                           import Align                           , \
                                            logFile                         , \
                                            TestScript_InputDataPath
from pMolecule                       import SystemGeometryObjectiveFunction
from pMolecule.QCModel               import ChargeModel                     , \
                                            DIISSCFConverger                , \
                                            ElectronicState                 , \
                                            QCModelDFT                      , \
                                            QCModelMNDO                     , \
                                            DFTGridAccuracy                 , \
                                            DFTGridIntegrator
from pMolecule.QCModel.GaussianBases import GaussianBasisOperator
from pScientific.Arrays              import ArrayPrint2D
from pSimulation                     import DensityFitMultipoles

#
# . The DFT QC Model is QCModelDFT. Its constructor takes the following arguments:
#
#   converger      - the SCF converger (currently only DIIS).
#   fitBasis       - the basis set for density fitting. The default is "dgauss-a1-dftjfit".
#   fitOperator    - the operator for density fitting. The default is "Coulomb".
#   functional     - the DFT functional. Many possibilities with "hf" the default.
#   gridIntegrator - the grid integrator (currently only DFTGridIntegrator).
#   orbitalBasis   - the orbital basis set. The default is "dgauss-dzvp".
#
#   DFTGridIntegrator options:
#
#   accuracy         - instance of DFTGridAccuracy enum, either VeryLow, Low, Medium (default), High or VeryHigh.
#   inCore           - recalculate or not (default) the basis function values at the grid points.
#                      True takes a lot of storage!
#
# . Good combinations of orbital and density basis sets are (in increasing order of cost):
#
#   dgauss-dzvp/dgauss-a1-dftjfit; def2-sv(p)/def2-sv(p)-rifit; def2-svp/def2-tzvp-rifit; def2-svp/def2-tzvp-rifit.
#
# . To perform QC/MM calculations, proceed in the usual way. Any of the NB models that can be used with
#   QCModelMNDO can also be used with QCModelDFT.
#

# . Start.
logFile.Header ( )

# . Options.
_DefaultFitBasis = { "sto3g"       : "dgauss-a1-dftjfit" ,
                     "dgauss-dzvp" : "dgauss-a1-dftjfit" ,
                     "def2-sv(p)"  : "def2-sv(p)-rifit"  ,
                     "def2-svp"    : "def2-svp-rifit"    }
_MaximumGAtoms   = 6
_TestGradients   = False

# . Define the QC models.
converger      = DIISSCFConverger.WithOptions  ( densityTolerance = 1.0e-8, maximumIterations = 250 )
gridIntegrator = DFTGridIntegrator.WithOptions ( accuracy = DFTGridAccuracy.Medium, inCore = True )
_EnergyModels  = []
for ( fitBasis, fitOperator, functional, orbitalBasis ) in ( ( "dgauss-a1-dftjfit" , GaussianBasisOperator.AntiCoulomb, "lda"  , "dgauss-dzvp" ) ,
                                                             ( "dgauss-a1-dftjfit" , GaussianBasisOperator.AntiCoulomb, "blyp" , "dgauss-dzvp" ) ,
                                                             ( "dgauss-a1-dftjfit" , GaussianBasisOperator.Coulomb    , "lda"  , "dgauss-dzvp" ) ,  
                                                             ( "dgauss-a1-dftjfit" , GaussianBasisOperator.Coulomb    , "blyp" , "dgauss-dzvp" ) ,  
                                                             ( "dgauss-a1-dftjfit" , GaussianBasisOperator.Overlap    , "lda"  , "dgauss-dzvp" ) ,  
                                                             ( "dgauss-a1-dftjfit" , GaussianBasisOperator.Overlap    , "blyp" , "dgauss-dzvp" ) ,  
                                                             ( None                , None                             , "b3lyp", "dgauss-dzvp" ) ,  
                                                             ( None                , None                             , "hf"   , "dgauss-dzvp" ) ,  
                                                             ( "def2-sv(p)-rifit"  , GaussianBasisOperator.AntiCoulomb, "lda"  , "def2-sv(p)"  ) ,  
                                                             ( "def2-sv(p)-rifit"  , GaussianBasisOperator.AntiCoulomb, "blyp" , "def2-sv(p)"  ) ,  
                                                             ( "def2-sv(p)-rifit"  , GaussianBasisOperator.Coulomb    , "lda"  , "def2-sv(p)"  ) ,  
                                                             ( "def2-sv(p)-rifit"  , GaussianBasisOperator.Coulomb    , "blyp" , "def2-sv(p)"  ) ,  
                                                             ( "def2-sv(p)-rifit"  , GaussianBasisOperator.Overlap    , "lda"  , "def2-sv(p)"  ) ,  
                                                             ( "def2-sv(p)-rifit"  , GaussianBasisOperator.Overlap    , "blyp" , "def2-sv(p)"  ) ,  
                                                             ( None                , None                             , "b3lyp", "def2-sv(p)"  ) ,  
                                                             ( None                , None                             , "hf"   , "def2-sv(p)"  ) ,  
                                                             ( "def2-svp-rifit"    , GaussianBasisOperator.AntiCoulomb, "lda"  , "def2-svp"    ) , 
                                                             ( "def2-svp-rifit"    , GaussianBasisOperator.AntiCoulomb, "blyp" , "def2-svp"    ) , 
                                                             ( "def2-svp-rifit"    , GaussianBasisOperator.Coulomb    , "lda"  , "def2-svp"    ) , 
                                                             ( "def2-svp-rifit"    , GaussianBasisOperator.Coulomb    , "blyp" , "def2-svp"    ) , 
                                                             ( "def2-svp-rifit"    , GaussianBasisOperator.Overlap    , "lda"  , "def2-svp"    ) , 
                                                             ( "def2-svp-rifit"    , GaussianBasisOperator.Overlap    , "blyp" , "def2-svp"    ) , 
                                                             ( None                , None                             , "b3lyp", "def2-svp"    ) , 
                                                             ( None                , None                             , "hf"   , "def2-svp"    ) ):
    label = "{:s}/{:s}".format ( functional, orbitalBasis )
    if fitBasis is not None:
        label += "/{:s}".format ( fitBasis )
        if   fitOperator is GaussianBasisOperator.AntiCoulomb: label += "/A"
        elif fitOperator is GaussianBasisOperator.Coulomb    : label += "/C"
        else                                                 : label += "/O"
    _EnergyModels.append ( ( label, QCModelDFT.WithOptions ( converger      = converger      ,
                                                             fitBasis       = fitBasis       ,
                                                             fitOperator    = fitOperator    ,
                                                             functional     = functional     ,
                                                             gridIntegrator = gridIntegrator ,
                                                             orbitalBasis   = orbitalBasis   ) ) )

# . Get the fileName.
_Names  = ( ( "chloride"     , -1 ) ,
            ( "dihydrogen"   ,  0 ) ,
            ( "dinitrogen"   ,  0 ) ,
            ( "water"        ,  0 ) ,
            ( "methane"      ,  0 ) ,
            ( "formaldehyde" ,  0 ) ,
            ( "ethene"       ,  0 ) ,
            ( "ozone"        ,  0 ) ,
            ( "fch3cl"       , -1 ) ,
            ( "nitromethane" ,  0 ) ,
            ( "glycine"      ,  0 ) )
_Source = os.path.join ( TestScript_InputDataPath ( "pMolecule" ), "xyz" )

# . Loop over the energy models.
lModel         = 0
lName          = 0
nEnergyFails   = 0
nGradientFails = 0
nQ             = 0
results        = []
for ( name, charge ) in _Names:
    fileName = os.path.join ( _Source, "{:s}.xyz".format ( name ) )
    for ( label, model ) in _EnergyModels:
        # . Set up system.
        system = ImportSystem ( fileName, log = None )
        system.electronicState = ElectronicState.WithOptions ( charge = charge )
        system.DefineQCModel ( model )
        system.coordinates3.ToPrincipalAxes ( weights = system.qcState.nuclearCharges ) # . To avoid setting a multipole center.
        system.Summary ( )
        # . Energy and properties.
        try:
            energy= system.Energy ( )
        except:
            logFile.Paragraph ( "Energy failure for system {:s} with model {:s}.".format ( name, label ) )
            nEnergyFails += 1
            continue
        charges    = system.qcModel.AtomicCharges    ( system, chargeModel = ChargeModel.Loewdin )
        dipole     = system.qcModel.DipoleMoment     ( system )
        quadrupole = system.qcModel.QuadrupoleMoment ( system )
        multipoles = list ( dipole  ) + [ dipole.Norm2  ( ) ] + [ quadrupole [0,0], quadrupole [1,1], quadrupole [2,2], quadrupole [0,1], quadrupole [0,2], quadrupole [1,2] ]
        properties = [ [ "Density", list ( charges ), multipoles ] ]
        if system.qcModel.fitBasis is None: fitBasis = _DefaultFitBasis[system.qcModel.orbitalBasis]
        else:                               fitBasis = None

        for ( operator, tag ) in ( ( GaussianBasisOperator.AntiCoulomb , "Anti-Coulomb Fit" ) ,
                                   ( GaussianBasisOperator.Coulomb     , "Coulomb Fit"      ) ,
                                   ( GaussianBasisOperator.Overlap     , "Overlap Fit"      ) ):
            items      = DensityFitMultipoles ( system, fitBasis = fitBasis, fitOperator = operator, log = None )
            charges    = items["Charges"   ]
            dipole     = items["Dipole"    ]
            quadrupole = items["Quadrupole"]
            multipoles = list ( dipole ) + [ dipole.Norm2 ( ) ] + [ quadrupole [0,0], quadrupole [1,1], quadrupole [2,2], quadrupole [0,1], quadrupole [0,2], quadrupole [1,2] ]
            properties.append ( [ tag, list ( charges ), multipoles ] )
        # . Test gradients.
        dG = None
        if _TestGradients and ( len ( system.atoms ) <= _MaximumGAtoms ):
            of = SystemGeometryObjectiveFunction.FromSystem ( system )
            try:
                dG = of.TestGradients ( )
            except:
                logFile.Paragraph ( "Gradient failure for system {:s} with model {:s}.".format ( name, label ) )
                nGradientFails += 1
                continue
        # . Save data for output.
        lModel = max ( lModel, len ( label ) )
        lName  = max ( lName , len ( name  ) )
        nQ     = max ( nQ, len ( system.atoms ) )
        results.append ( ( name, label, energy, dG, properties ) )

# . Output the results.
nBefore = 4
if _TestGradients: nBefore += 1
nColumns = nBefore + nQ + 10
table = logFile.GetTable ( columns = [ lName+2, lModel+2 ] + ( nBefore - 3 ) * [ 12 ] + [ 18 ] + ( nColumns - nBefore ) * [ 10 ] )
table.Start  ( )
table.Title  ( "Energy Model Results" )
table.Heading ( "System"       )
table.Heading ( "QC Model"     )
table.Heading ( "Energy"       )
if _TestGradients: table.Heading ( "Gradient Error" )
table.Heading ( "Charge Model" )
table.Heading ( "Charges"      , columnSpan = nQ )
table.Heading ( "Dipole"       , columnSpan =  4 )
table.Heading ( "Quadrupole"   , columnSpan =  6 )
for h in ( nBefore + nQ )* [ "" ] + [ "X", "Y", "Z", "Norm", "XX", "YY", "ZZ", "XY", "XZ", "ZZ" ]:
    table.Heading ( h )
for ( name, label, energy, dG, properties ) in results:
    table.Entry ( name , align = Align.Left )
    table.Entry ( label, align = Align.Left )
    table.Entry ( "{:.1f}".format ( energy ) )
    if _TestGradients:
        if dG is None: table.Entry ( "-" )
        else:          table.Entry ( "{:.5f}".format ( dG ) )
    for ( p, ( model, charges, multipoles ) ) in enumerate ( properties ):
        if p > 0:
            for i in range ( nBefore-1 ): table.Entry ( "" )
        table.Entry ( "  " + model, align = Align.Left )
        for charge in charges: table.Entry ( "{:.3f}".format ( charge ) )
        for i in range ( len ( charges ), nQ ): table.Entry ( "-" )
        for multipole in multipoles: table.Entry ( "{:.3f}".format ( multipole ) )
table.Stop ( )

# . Other output.
if nEnergyFails   > 0: logFile.Paragraph ( "There were {:d} energy failures.".format   ( nEnergyFails   ) )
if nGradientFails > 0: logFile.Paragraph ( "There were {:d} gradient failures.".format ( nGradientFails ) )

# . Stop.
logFile.Footer ( )

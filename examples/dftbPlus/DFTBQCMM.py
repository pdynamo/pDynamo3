"""DFTB+ QC calculations."""

import os, os.path

from Definitions       import dataPath                                 , \
                              dataPathM                                , \
                              outPath
from pBabel            import ExportSystem                             , \
                              ImportSystem
from pCore             import logFile                                  , \
                              Selection
from pMolecule         import SystemGeometryObjectiveFunction
from pMolecule.MMModel import MMModelOPLS
from pMolecule.NBModel import NBModelDFTB
from pMolecule.QCModel import ElectronicState                          , \
                              QCModelDFTB
from pSimulation       import ConjugateGradientMinimize_SystemGeometry , \
                              ModifyOption                             , \
                              NormalModes_SystemGeometry

# . Header.
logFile.Header ( )

# . DFTB+ settings for the 3ob parameter set used in this example (in skfPath).
eiHO   = "  ThirdOrderFull = Yes\n  HCorrection = Damping { Exponent = 4.0 }\n  HubbardDerivs = { \n    H = -0.1857 \n    O = -0.1575 \n }"
eiCHON = "  ThirdOrderFull = Yes\n  HCorrection = Damping { Exponent = 4.0 }\n  HubbardDerivs = { \n    H = -0.1857 \n    O = -0.1575 \n    C = -0.1492 \n    N = -0.1535 \n }"

# . Loop over molecules and QC model options (QC/MM only works with SCC).
for ( mLabel, qIndices ) in ( ( "bAla_c7eq", [ 10, 11, 12, 13 ] ), ( "waterDimer_cs", [ 0, 1, 2 ] ) ):
    for ( useSCC, qLabel ) in ( ( True , "scc"   ), ):

        # . Select settings.
        if "bAla" in mLabel:
           ei = eiCHON
        else:
           ei = eiHO

        # . Define various models.
        electronicState = ElectronicState.WithOptions ( charge = 0, multiplicity = 1 )
        mmModel         = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
        nbModel         = NBModelDFTB.WithDefaults ( )
        qcModel         = QCModelDFTB.WithOptions ( deleteJobFiles = False                            ,
                                                    randomScratch  = True                             ,
                                                    scratch        = outPath                          ,
                                                    skfPath        = os.path.join ( dataPath, "skf" ) ,
                                                    useSCC         = useSCC                           )

        # . Define the system.
        system                 = ImportSystem ( os.path.join ( dataPathM, "mol", mLabel + ".mol" ) )
        system.electronicState = electronicState
        system.DefineMMModel ( mmModel )
        system.DefineQCModel ( qcModel, qcSelection = Selection.FromIterable ( qIndices ) )
        system.DefineNBModel ( nbModel )
        system.Summary ( )

        # . Calculate an energy.
        system.Energy ( doGradients = True )
        system.nbModel.StatisticsSummary ( system )

        # . Test gradients.
        of = SystemGeometryObjectiveFunction.FromSystem ( system )
        of.TestGradients ( )

        # . Optimization.
        eStart = system.Energy ( )
        ConjugateGradientMinimize_SystemGeometry ( system                      ,
                                                   logFrequency         =  100 ,
                                                   maximumIterations    = 2000 ,
                                                   rmsGradientTolerance =  0.1 )
        eStop = system.Energy ( doGradients = True )
        logFile.Paragraph ( "Energy difference after minimization = {:.1f}.".format ( eStop - eStart ) )

        # . Normal modes.
        NormalModes_SystemGeometry ( system, modify = ModifyOption.Project )

        # . Finish up.
        ExportSystem ( os.path.join ( outPath, "{:s}_{:s}_qcmm.xyz".format ( mLabel, qLabel ) ), system )

# . Footer.
logFile.Footer ( )

"""DFTB+ QC calculations."""

import os, os.path

from Definitions       import dataPath                                 , \
                              dataPathM                                , \
                              outPath
from pBabel            import ExportSystem                             , \
                              ImportSystem
from pCore             import logFile
from pMolecule         import SystemGeometryObjectiveFunction
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

# . Loop over molecules and QC model options (only +/- SCC in this case).
for mLabel in ( "bAla_c7eq", "water", "waterDimer_cs" ):
    for ( useSCC, qLabel ) in ( ( False, "noSCC" ) ,
                                ( True , "scc"   ) , ):

        # . Select settings.
        ei = None
        if useSCC:
           if "bAla" in mLabel:
              ei = eiCHON
           else:
              ei = eiHO
        # . Set up the system.
        system                 = ImportSystem ( os.path.join ( dataPathM, "mol", mLabel + ".mol" ) )
        system.electronicState = ElectronicState.WithOptions ( charge = 0, multiplicity = 1 )
        qcModel                = QCModelDFTB.WithOptions ( deleteJobFiles = False                            ,
                                                           randomScratch  = True                             ,
                                                           scratch        = outPath                          ,
                                                           skfPath        = os.path.join ( dataPath, "skf" ) ,
                                                           useSCC         = useSCC                           ,
                                                           extendedInput  = ei   )
        system.DefineQCModel ( qcModel )
        system.Summary ( )

        # . Energy and gradients.
        e  = system.Energy ( doGradients = True )
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
        ExportSystem ( os.path.join ( outPath, "{:s}_{:s}.xyz".format ( mLabel, qLabel ) ), system )

# . Footer.
logFile.Footer ( )

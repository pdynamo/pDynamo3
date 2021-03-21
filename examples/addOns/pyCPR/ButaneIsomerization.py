"""pyCPR test case."""

import os, os.path

from addOns.pyCPR      import ConjugatePeakRefinementOptimizePath
from Definitions       import dataPath                 , \
                              outPath
from pBabel            import ExportTrajectory         , \
                              ImportSystem
from pCore             import logFile                  , \
                              TestScriptExit_Fail
from pMolecule.QCModel import QCModelMNDO
from pSimulation       import FromCoordinateFiles      , \
                              ToCoordinateFiles

# . Header.
logFile.Header ( )

# . Paths.
molIn  = os.path.join ( dataPath, "mol"         )
xyzIn  = os.path.join ( dataPath, "xyz"         )
resOut = os.path.join ( outPath , "cpr.restart" )
xyzOut = os.path.join ( outPath , "xyz"         )
if not os.path.exists ( xyzOut ): os.mkdir ( xyzOut )

# . Define the molecule and its energy model.
system  = ImportSystem ( os.path.join ( molIn, "butane.mol" ) )
qcModel = QCModelMNDO.WithOptions ( hamiltonian = "rm1" )
system.DefineQCModel ( qcModel )
system.Summary ( )

# . Generate initial trajectory (at least 2 structures).
xyzPaths = [ os.path.join ( xyzIn, "butane-anti-opt.xyz"   ) ,
             os.path.join ( xyzIn, "butane-gauche-opt.xyz" ) ]

# . Directory to save the path structures to.
trjPath    = os.path.join ( outPath, "cpr.ptGeo" )
FromCoordinateFiles ( xyzPaths, trjPath, system )
trajectory = ExportTrajectory ( trjPath, system, append = True )

# . Run CPR on the system starting from the initial trajectory.
ConjugatePeakRefinementOptimizePath ( system                              ,
                                      trajectory                          ,
                                      breakIfTauReached          = False  ,
				      finalUnrefinableRefinement = False  ,
                                      interpolationIncreasement  = 0      ,
                                      outputRestartFile          = resOut ,
				      rmsGradientTolerance       = 0.001  )

# Convert the final path to XYZ format.
ToCoordinateFiles ( trjPath, xyzOut, system, outFormat = "xyz", padding = 3 )

# . Footer.
isOK = True
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )

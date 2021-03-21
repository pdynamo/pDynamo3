"""Example 12."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the molecule and its energy model.
molecule              = ImportSystem       ( os.path.join ( molPath    , "cyclohexane_chair.mol"  ) )
molecule.coordinates3 = ImportCoordinates3 ( os.path.join ( scratchPath, "cyclohexane_saddle.xyz" ) )
molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
molecule.Summary ( )

# . Calculate an energy.
molecule.Energy ( )

# . Optimization.
trajectory = ExportTrajectory ( os.path.join ( scratchPath, "cyclohexane_sdpath.ptGeo" ), molecule )
SteepestDescentPath_SystemGeometry ( molecule                        ,
                                     functionStep      =   2.0       ,
                                     logFrequency      =  10         ,
                                     maximumIterations = 400         ,
                                     pathStep          =   0.025     ,
                                     saveFrequency     =  10         ,
                                     trajectory        = trajectory  ,
                                     useMassWeighting  = True        )

# . Footer.
logFile.Footer ( )

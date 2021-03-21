"""Example 14."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the molecule and its energy model.
molecule              = ImportSystem       ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.coordinates3 = ImportCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_chair.xyz" ) )
molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
molecule.Summary ( )

# . Calculate the normal modes.
NormalModes_SystemGeometry ( molecule, modify = ModifyOption.Project )

# . Generate a trajectory for one of the modes.
trajectory = ExportTrajectory ( os.path.join ( scratchPath, "cyclohexane_chair_mode7.ptGeo" ), molecule )
NormalModesTrajectory_SystemGeometry ( molecule            ,
                                       trajectory          ,
                                       mode        =  7    ,
                                       cycles      = 10    ,
                                       frames      = 21    ,
                                       temperature = 600.0 )

# . Footer.
logFile.Footer ( )

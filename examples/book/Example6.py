"""Example 6."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Generate the molecule.
molecule = ImportSystem ( os.path.join ( xyzPath, "bala_c7eq.xyz" ) )
molecule.DefineQCModel ( QCModelMNDO.WithDefaults ( ) )
molecule.Summary ( )

# . Create an objective function for the molecule.
of = SystemGeometryObjectiveFunction.FromSystem ( molecule )

# . Test the gradients.
of.TestGradients ( )

# . Footer.
logFile.Footer ( )

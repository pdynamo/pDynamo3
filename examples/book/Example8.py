"""Example 8."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the MM, NB and QC models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )
qcModel = QCModelMNDO.WithDefaults ( )

# . Define the molecule.
molecule = ImportSystem ( os.path.join ( molPath, "waterDimer_cs.mol" ) )

# . Define the selection for the first molecule.
firstWater = Selection.FromIterable ( [ 0, 1, 2 ] )

# . Define the energy model.
molecule.DefineMMModel ( mmModel )
molecule.DefineQCModel ( qcModel, qcSelection = firstWater )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Calculate an energy.
molecule.Energy ( doGradients = True )
molecule.nbModel.StatisticsSummary ( molecule )

# . Footer.
logFile.Footer ( )

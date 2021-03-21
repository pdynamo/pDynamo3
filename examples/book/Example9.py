"""Example 9."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the MM, NB and QC models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )
qcModel = QCModelMNDO.WithDefaults ( )

# . Define the molecule.
molecule = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )

# . Define the selection for the first molecule.
methylGroup = Selection.FromIterable ( [ 10, 11, 12, 13 ] )

# . Define the energy model.
molecule.DefineMMModel ( mmModel )
molecule.DefineQCModel ( qcModel, qcSelection = methylGroup )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Calculate an energy.
molecule.Energy ( doGradients = True )
molecule.nbModel.StatisticsSummary ( molecule )

# . Footer.
logFile.Footer ( )

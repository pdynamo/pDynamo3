"""Simple PM6-ML test."""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RESULTS (kJ/mol)
# PM6: pDynamo -414.12  MOPAC/cuby -416.31
# D3:  pDynamo  -80.73  MOPAC/cuby  -78.94
# ML:  pDynamo -387.98  MOPAC/cuby -387.98
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os.path

from pCore              import Selection, TestScript_InputDataPath, NotInstalledError, logFile
from pBabel             import ImportSystem 

from pMolecule          import EnergyModelPriority
from pMolecule.QCModel  import QCModelMNDO
from pScientific        import Units

try:
    from addOns.PM6ML       import QCDeltaMLModel
    import numpy, torch, dftd3
except:
    raise NotInstalledError( "The QC Delta-ML model requires the numpy, torchmd-net, simple-dftd3 and dftd3-python python modules installed." )
    # See instalation notes at $PDYNAMO3_HOME/addOns/PM6ML/__init__.py

# . Start.
logFile.Header ( )

# . The input data paths.
dataPath = TestScript_InputDataPath ( "book" )
molPath  = os.path.join ( dataPath, "mol" )

# . Define the molecule.
molecule  = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )

# . Define the energy model with DeltaML correction.
molecule.DefineQCModel  ( QCModelMNDO.WithOptions ( hamiltonian = "pm6" ) )
eml = QCDeltaMLModel ( ml_model_file = "PM6-ML_correction_seed8_best.ckpt" )
molecule.AddEnergyModel ('QC Delta-ML correction', eml, priority = EnergyModelPriority.QCAddOns )

# . Calculate an energy.
molecule.Energy( doGradients=True )

# . Stop.
logFile.Footer ( )


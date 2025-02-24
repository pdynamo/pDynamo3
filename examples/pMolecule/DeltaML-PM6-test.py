"""Simple PM6-ML test."""

import os.path

from pCore              import Selection, TestScript_InputDataPath
from pBabel             import ImportSystem 

from pMolecule          import EnergyModelPriority
from pMolecule.QCModel  import QCDeltaMLModel , QCModelMNDO

# . The input data paths.
dataPath = TestScript_InputDataPath ( "book" )
molPath  = os.path.join ( dataPath, "mol" )

# . Define the molecule.
molecule  = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )

# . Define the energy model with DeltaML correction.
molecule.DefineQCModel  ( QCModelMNDO.WithOptions ( hamiltonian = "pm6" ) )
eml = QCDeltaMLModel("/home/rezac/Calculations/PM6-ML/MOPAC-ML/models/PM6-ML_correction_seed8_best.ckpt")
molecule.AddEnergyModel ('QC Delta-ML correction', eml, priority = EnergyModelPriority.QCAddOns )

# . Calculate an energy.
molecule.Energy( doGradients=True )


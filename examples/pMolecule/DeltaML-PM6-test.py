"""Simple PM6-ML test."""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RESULTS (kJ/mol)
# PM6: pDynamo -414.12  MOPAC/cuby -416.31
# D3:  pDynamo  -80.73  MOPAC/cuby  -78.94
# ML:  pDynamo  -92.73  MOPAC/cuby  -92.73
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import os.path

from pCore              import Selection, TestScript_InputDataPath
from pBabel             import ImportSystem 

from pMolecule          import EnergyModelPriority
from pMolecule.QCModel  import QCDeltaMLModel , QCModelMNDO
from pScientific        import Units

import numpy as np

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

# Print gradient
#print(molecule.__dict__)
qcGradients3 = molecule._scratch.qcGradients3AU # It is actually in kJ/mol/A
print("#" * 80)
print("Gradient")
rmsgrad2 = qcGradients3.RootMeanSquare()
print(f"RMS kJ/mol/A: {rmsgrad2}")
print(f"RMS kcal/mol/A: {rmsgrad2 * 0.2390057361376673}")
print("Gradient in kcal/mol/A")
print(np.array(qcGradients3).reshape((-1,3)) * 0.2390057361376673)
print("#" * 80)




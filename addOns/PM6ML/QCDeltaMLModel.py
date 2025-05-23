"""Delta-ML (Machine Learning) correction to QC PM6 energy and gradients."""


from   pCore              import logFile               
from   pScientific        import Units , PeriodicTable
from   pMolecule.QCModel  import QCModelError
from   pMolecule          import EnergyClosurePriority , EnergyModelState , EnergyModel


# NumPy (dftd3 expects numpy arrays on input)
try:
    import numpy as np
    numpy_available = True
except:
    numpy_available = False

# Torch and TorchMD-NET
try:
    import torch
    from torchmdnet.models.model import load_model
    torchmd_available = True
except:
    torchmd_available = False

# Dispersion correction
try:
    from dftd3.interface import RationalDampingParam, DispersionModel
    dftd3_available = True
except:
    dftd3_available = False


#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . State class.
#===================================================================================================================================
class QCDeltaMLModelState ( EnergyModelState ):
    """QC DeltaML Model State."""

    _attributable = dict ( EnergyModelState._attributable )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCDeltaMLModel ( EnergyModel ):
    """QCDeltaML model."""

    # . Defaults.
    _attributable = dict ( EnergyModel._attributable )
    _classLabel   = "QC DeltaML Model"
    _stateName    = "qcDeltaMLState"
    _stateObject  = QCDeltaMLModelState
    _summarizable = dict ( EnergyModel._summarizable )

    # Constants
    Z_TO_ATYPE = {
         1: 10, # H
         3: 14, # Li
         6:  3, # C
         7: 17, # N
         8: 21, # O
         9:  9, # F
        11: 19, # Na
        12: 15, # Mg
        15: 23, # P
        16: 26, # S
        17:  7, # Cl
        19: 13, # K
        20:  5, # Ca
        35:  1, # Br
        53: 12, # I
    }

    def __init__(self, ml_model_file=None):
        """Initialize the model and load ML correction data."""
        # Load ML model
        if not ml_model_file:
            raise QCModelError("The QC Delta-ML model requires the PM6-ML correction weights.")
        else:    
            self.model = load_model(ml_model_file, derivative=True)
            self.device = torch.get_default_device()
            self.model.to(self.device)

    def BuildModel ( self, target ):
        """Build the model."""
        # Python modules check
        if not numpy_available:
            raise QCModelError("The QC Delta-ML model requires the numpy python module installed.")
        if not torchmd_available:
            raise QCModelError("The QC Delta-ML model requires the torchmd-net python module installed.")
        if not dftd3_available:
            raise QCModelError("The QC Delta-ML model requires the simple-dftd3 and dftd3-python python modules installed.")
        # Method checks
        if target.qcState is None:
            raise QCModelError ( "This QC Delta-ML model requires a pre-defined QC model." )
        elif target.qcModel.hamiltonian != 'pm6':     
            raise QCModelError ( "This QC Delta-ML correction is only valid for the PM6 hamiltonian." )
        else:
            state = super ( QCDeltaMLModel, self ).BuildModel ( target )

        # Check for parametrized elements
        for an in target.qcState.atomicNumbers:
            if not an in self.__class__.Z_TO_ATYPE:
                raise QCModelError(f"The QC Delta-ML model does not have parameters for element {PeriodicTable.Symbol(an)} .")

    def Energy ( self, target ):
        """Calculate ML and D3 corrections to the quantum chemical energy."""
        # Coordinates and atomic numbers
        coordinates3  = target.scratch.qcCoordinates3AU
        atomicNumbers = target.qcState.atomicNumbers 

        # Compute D3 correction
        disp = DispersionModel(
            numbers=np.array(atomicNumbers), positions=np.array(coordinates3) # Coordinates read in AU
        )
        disp_res = disp.get_dispersion(
            RationalDampingParam(s6=1.0, s8=0.3908, a1=0.566, a2=3.128), grad=target.scratch.doGradients
        )
        eD3_au = disp_res.get("energy") # Energy in hartrees

        # Compute ML correction
        # Convert atomic numbers to atom types
        atomtypes = []
        for i, an in enumerate(atomicNumbers):
            if not an in self.__class__.Z_TO_ATYPE:
                raise QCModelError(f"The Delta-ML correction is not available for element {PeriodicTable.Symbol(an)} .")
            atomtypes.append(self.__class__.Z_TO_ATYPE[an])
        # Convert coordinates to Angstrom, reshape the array to 3*natom
        geom = np.array(coordinates3) * Units.Length_Bohrs_To_Angstroms
        geom.shape = (geom.size//3, 3)
        # Prepare data for torch
        types = torch.tensor(atomtypes, dtype=torch.long)
        types = types.to(self.device)
        pos = torch.tensor(geom, dtype=torch.float32)
        pos = pos.to(self.device)
        # Run the calculation
        ml_energy, ml_forces = self.model.forward(types, pos)
        eML = ml_energy.item() # Returns energy in kJ/mol

        # Add energy corrections.
        target.scratch.energyTerms["QC Delta-ML correction"]      = eML
        target.scratch.energyTerms["QC Dispersion D3 correction"] = eD3_au * Units.Energy_Hartrees_To_Kilojoules_Per_Mole

        doGradients  = target.scratch.doGradients
        if doGradients:
            n = len ( target.qcState.atomicNumbers )
            # D3 gradient, to be retrieved from disp_res computed above
            gradients3_d3 = coordinates3.WithExtent(n)
            # Get gradient, it is in a.u.
            dftd3_grad = disp_res.get("gradient")
            # Copy to gradients3_d3
            for i in range(n):
                for j in range(3):
                    gradients3_d3[i][j] = dftd3_grad[i][j]

            # ML corection gradient, available in ml_forces
            gradients3_ml = coordinates3.WithExtent(n)
            for i in range(n):
                for j in range(3):
                    # Not gradient but forces
                    # in kJ/mol
                    gradients3_ml[i][j] = (ml_forces[i][j].item() * -1.0
                        / Units.Energy_Hartrees_To_Kilojoules_Per_Mole 
                        / Units.Length_Angstroms_To_Bohrs
                        )

            # Add gradient corrections to gradients3
            qcGradients3 = target.scratch.qcGradients3AU # QM gradient in a.u.
            # qcGradients3.Set(0.0) # For testing corrections only
            gradients3_d3.ScatterAdd(1.0, qcGradients3)
            gradients3_ml.ScatterAdd(1.0, qcGradients3)

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ):
            results = self.Energy ( target )
        return [ ( EnergyClosurePriority.QCGradients, a, "QC Delta-ML and D3 correction" ) ]

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCDeltaMLModel, self ).SummaryItems ( )
        return items

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

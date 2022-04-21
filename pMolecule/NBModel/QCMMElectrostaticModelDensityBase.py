"""Base class for QC/MM electrostatic density models."""

from   pMolecule.QCModel      import FockClosurePriority
from   pScientific            import Units
from   pScientific.Arrays     import Array                  , \
                                     StorageType
from  .QCMMElectrostaticModel import QCMMElectrostaticModel
from ..EnergyModel            import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelDensityBase ( QCMMElectrostaticModel ):
    """Base class for density QC/MM electrostatic models."""

# . NB model after QC model ... so can call QCModel .pyx from here.
# . NB model can also ask for appropriate QC/MM model from QCModel?
# . E.g. Full/ABFS, etc.

    _attributable = dict ( QCMMElectrostaticModel._attributable )
    _classLabel   = "Density QC/MM Electrostatic Model"
    _attributable.update ( { "evaluator" : None } )

    def BuildModel ( self, target ):
        """Build the model."""
        state = super ( QCMMElectrostaticModelDensityBase, self ).BuildModel ( target )
        target.qcState.AddFockModel ( "QC/MM Electrostatic", self )
        return state

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.QCBPPotentials ( target )
        def b ( ): self.QCBPGradients  ( target )
        closures = super ( QCMMElectrostaticModelDensityBase, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, a, "QC/BP Electrostatic Potentials" ) ,
                            ( EnergyClosurePriority.QCGradients, b, "QC/BP Electrostatic Gradients"  ) ] )
        return closures

    def EnergyInitialize ( self, target ):
        """Energy initialization"""
        super ( QCMMElectrostaticModelDensityBase, self ).EnergyInitialize ( target )
        # . Potentials array - QC/MM and QC/BP not separated for the moment.
        potentials = target.scratch.Get ( "qcmmPotentials", None )
        if potentials is None:
            n = len ( target.qcState.orbitalBases )
            potentials = Array.WithExtent ( n, storageType = StorageType.Symmetric )
            target.scratch.qcmmPotentials = potentials
        potentials.Set ( 0.0 )

    def Fock ( self, target ):
        """Energy and Fock matrix contributions."""
        scratch    = target.scratch
        dTotal     = scratch.onePDMP.density
        fTotal     = scratch.onePDMP.fock
        potentials = scratch.qcmmPotentials
        eQCMM      = dTotal.TraceOfProduct ( potentials ) # . Potentials have -1 factor for electrons.
        fTotal.Add ( potentials )
        scratch.energyTerms["QC/MM Electrostatic"] = ( eQCMM * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        return eQCMM

    def FockClosures ( self, target ):
        """Fock closures."""
        def a ( ):
            return self.Fock ( target )
        return [ ( FockClosurePriority.Low, a ) ]

    def QCBPGradients ( self, target ):
        """Calculate the QC/BP electrostatic gradients."""
        pass

    def QCBPPotentials ( self, target ):
        """Calculate the QC/BP electrostatic potentials."""
        pass

    def UnbuildModel ( self, target ):
        """Unbuild the model."""
        target.qcState.AddFockModel ( "QC/MM Electrostatic", None )
        super ( QCMMElectrostaticModelDensityBase, self ).UnbuildModel ( target )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

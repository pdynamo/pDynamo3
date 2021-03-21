"""Base class for QC/MM electrostatic multipole models."""

from   pMolecule.QCModel      import FockClosurePriority
from   pScientific.Arrays     import Array
from  .QCMMElectrostaticModel import QCMMElectrostaticModel
from ..EnergyModel            import EnergyClosurePriority

# . Multipole order - O to L.
# . Spherical multipoles - 1, 3, 5,  7, ... 2L+1         : cumulative sum is (L+1)^2.
# . Cartesian multipoles - 1, 3, 6, 10, ... (L+1)(L+2)/2 : cumulative sum is (L+1)(L+2)(L+3)/6.

# . For the moment, the multipole and potential arrays are flattened nMult (Cartesian) x nQC row-major matrices.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelMultipoleBase ( QCMMElectrostaticModel ):
    """Base class for multipole QC/MM electrostatic models."""

    _attributable = dict ( QCMMElectrostaticModel._attributable )
    _summarizable = dict ( QCMMElectrostaticModel._summarizable )
    _classLabel   = "Multipole QC/MM Electrostatic Model"
    _attributable.update ( { "multipoleOrder" : 0                 } )
    _summarizable.update ( { "multipoleOrder" : "Multipole Order" } )

    def BuildModel ( self, target ):
        """Build the model."""
        state = super ( QCMMElectrostaticModelMultipoleBase, self ).BuildModel ( target )
        target.qcState.AddFockModel ( "QC/MM Electrostatic", self )
        return state

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        # . "Weighted density".
        # . Note that this is not really a weighted density.
        # . It is the coordinate derivative term of the multipoles which involve the
        # . derivatives of the overlap matrix. These terms are conveniently added in
        # . here to avoid calculating the overlap derivatives twice.
        closures = super ( QCMMElectrostaticModelMultipoleBase, self ).EnergyClosures ( target )
        if hasattr ( target.qcModel.multipoleEvaluator, "WeightedDensity" ):
           def a ( ): self.GetWeightedDensity ( target )
           closures.append ( ( EnergyClosurePriority.QCPreGradients, a, "QC/MM Weighted Density" ) )
        return closures

    def EnergyInitialize ( self, target ):
        """Energy initialization"""
        super ( QCMMElectrostaticModelMultipoleBase, self ).EnergyInitialize ( target )
        # . Multipole and potential arrays.
        l = self.multipoleOrder
        n = len ( target.qcState.qcAtoms ) * ( ( l + 1 ) * ( l + 2 ) * ( l + 3 ) ) // 6
        multipoles = target.scratch.Get ( "qcmmMultipoles", None )
        potentials = target.scratch.Get ( "qcmmPotentials", None )
        target.scratch.qcmmMultipoleOrder = self.multipoleOrder
        if multipoles is None:
            multipoles = Array.WithExtent ( n )
            target.scratch.qcmmMultipoles = multipoles
        if potentials is None:
            potentials = Array.WithExtent ( n )
            target.scratch.qcmmPotentials = potentials
        multipoles.Set ( 0.0 )
        potentials.Set ( 0.0 )

    def Fock ( self, target ):
        """Energy and Fock matrix contributions."""
        return 0.0

    def FockClosures ( self, target ):
        """Fock closures."""
        def a ( ):
            return self.Fock ( target )
        return [ ( FockClosurePriority.Low, a ) ]

    def UnbuildModel ( self, target ):
        """Unbuild the model."""
        target.qcState.AddFockModel ( "QC/MM Electrostatic", None )
        super ( QCMMElectrostaticModelMultipoleBase, self ).UnbuildModel ( target )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

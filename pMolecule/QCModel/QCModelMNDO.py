"""pDynamo's in-built MNDO QC model."""

from   math                    import pow
from   pCore                   import Clone
from   pScientific.Arrays      import Array                          , \
                                      StorageType
from  .GaussianBases           import GaussianBasisIntegralEvaluator
from  .MNDOIntegralEvaluator   import MNDOIntegralEvaluator   
from  .MNDOMultipoleEvaluator  import MNDOMultipoleEvaluator
from  .MNDOParametersContainer import MNDOParametersContainer 
from  .QCDefinitions           import ChargeModel                    , \
                                      FockClosurePriority
from  .QCModelBase             import QCModelBase                    , \
                                      QCModelBaseState
from  .QCModelError            import QCModelError            
from ..EnergyModel             import EnergyClosurePriority   

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelMNDOState ( QCModelBaseState ):
    """A QC MNDO model state."""

    _attributable = dict ( QCModelBaseState._attributable )
    _attributable.update ( { "mndoParameters" : None } )

    def __setstate__ ( self, state ):
        """Set the state of the object."""
        super ( QCModelMNDOState, self ).__setstate__ ( state )
        if self.mndoParameters is not None: self.mndoParameters.ScaleBasesExponents ( self.orbitalBases )

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelMNDOState, self ).SummaryItems ( )
        items.append ( ( "MNDO Parameter Entries", "{:d}".format ( len ( self.mndoParameters.uniqueEntries ) ) ) )
        return items

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelMNDO ( QCModelBase ):
    """A MNDO QC model."""

    _attributable = dict ( QCModelBase._attributable )
    _chargeModels = { ChargeModel.MNDO : MNDOMultipoleEvaluator }
    _classLabel   = "MNDO QC Model"
    _stateObject  = QCModelMNDOState
    _summarizable = dict ( QCModelBase._summarizable )
    _attributable.update ( { "hamiltonian"        : "am1"                  ,
                             "integralEvaluator"  : MNDOIntegralEvaluator  ,
                             "multipoleEvaluator" : MNDOMultipoleEvaluator } )
    _summarizable.update ( { "hamiltonian"        :"Hamiltonian"           } )

    def EnergyClosureGradients ( self, target ):
        """Gradient energy closure."""
        def a ( ): self.integralEvaluator.ResonanceGradients          ( target )
        def b ( ): self.integralEvaluator.ElectronNuclearTEIGradients ( target )
        return [ ( EnergyClosurePriority.QCGradients, a, "QC Resonance Gradients"            ) ,
                 ( EnergyClosurePriority.QCGradients, b, "QC Electron-Nuclear/TEI Gradients" ) ]

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.integralEvaluator.CoreCoreEnergy              ( target )
        def b ( ): self.integralEvaluator.ElectronNuclearTEIIntegrals ( target )
        def c ( ): self.integralEvaluator.ResonanceIntegrals          ( target )
        closures = super ( QCModelMNDO, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, a, "QC Nuclear"                        ) ,
                            ( EnergyClosurePriority.QCIntegrals, b, "QC Electron-Nuclear/TEI Integrals" ) ,
                            ( EnergyClosurePriority.QCIntegrals, c, "QC Resonance Integrals"            ) ] )
        closures.extend ( self.EnergyClosureGradients ( target ) )
        return closures

    def GetNonOrthogonalDensity ( self, target, spinType = None ):
        """Get the density in the non-orthogonal basis."""
        densityA = super ( QCModelMNDO, self ).GetNonOrthogonalDensity ( target, spinType = spinType )
        self.GetOrthogonalizer ( target, doInverse = False, doLoewdin = False )
        orthogonalizer = target.scratch.orthogonalizer
        density = Array.WithExtent ( orthogonalizer.shape[0], storageType = StorageType.Symmetric )
        densityA.Transform ( orthogonalizer, density, useTranspose = True )
        return density

    def GetNonOrthogonalOrbitals ( self, target, indices = None, spinType = None ):
        """Get the orbitals in the non-orthogonal basis."""
        orbitalsA = super ( QCModelMNDO, self ).GetNonOrthogonalOrbitals ( target, indices = indices, spinType = spinType )
        self.GetOrthogonalizer ( target, doInverse = False, doLoewdin = False )
        orthogonalizer = target.scratch.orthogonalizer
        orbitals       = Array.WithShape ( [ orthogonalizer.shape[0], orbitalsA.shape[1] ] )
        orbitals.MatrixMultiply ( orthogonalizer, orbitalsA )
        return orbitals

    def GetOrthogonalizer ( self, target, doInverse = True, doLoewdin = True ): 
        """Get the orthogonalizing transformation and, optionally, its inverse."""
        # . Ensure the overlap matrix is of the correct dimension for the grid evaluator.
        orbitalBases                 = target.qcState.orbitalBases
        overlap                      = Array.WithExtent ( orbitalBases.numberOfFunctions, storageType = StorageType.Symmetric )
        target.scratch.overlapMatrix = overlap
        self.gridPointEvaluator.f1Of1i ( target )
        super ( QCModelMNDO, self ).GetOrthogonalizer ( target, doInverse = doInverse, doLoewdin = doLoewdin )
        target.scratch.Pop ( "overlapMatrix" )

    def GetParameters ( self, target ):
        """Get the parameters for the model."""
        state = target.qcState
        state.mndoParameters = MNDOParametersContainer.FromParameterDirectory ( self.hamiltonian, state.atomicNumbers )
        state.energyBaseLine = state.mndoParameters.energyBaseLine
        state.nuclearCharges = state.mndoParameters.coreCharges
        state.orbitalBases   = state.mndoParameters.MakeOrbitalBasis ( state.atomicNumbers )

    @property
    def gridPointEvaluator ( self ):
        """The grid point evaluator for property calculations."""
        evaluator = self.__dict__.get ( "_gridPointEvaluator", None )
        if evaluator is None:
            evaluator = GaussianBasisIntegralEvaluator ( )
            self.__dict__["_gridPointEvaluator"] = evaluator
        return evaluator

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

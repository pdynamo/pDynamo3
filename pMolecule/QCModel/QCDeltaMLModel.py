"""Delta-ML (Machine Learning) correction to QC energy and gradients."""


from   pCore              import logFile               
from   pScientific        import Units , PeriodicTable
from  .QCModelError       import QCModelError
from ..EnergyModel        import EnergyClosurePriority , \
                                 EnergyModelState      , \
                                 EnergyModel

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

    def BuildModel ( self, target ):
        """Build the model."""
        if target.qcState is None:
            raise QCModelError ( "This QC Delta-ML model requires a pre-defined QC model." )
        elif target.qcModel.hamiltonian != 'pm6':     
            raise QCModelError ( "This QC Delta-ML correction is only valid for the PM6 hamiltonian." )
        else:
            state = super ( QCDeltaMLModel, self ).BuildModel ( target )

    def Energy ( self, target ):
        """Calculate ML and D3 corrections to the quantum chemical energy."""

        coordinates3  = target.scratch.qcCoordinates3AU
        atomicNumbers = target.qcState.atomicNumbers 
        # PeriodicTable.Symbol ( atomicNumbers[i] ) gives atomic symbol of atom rank i.
        # . Add energy corrections.
        eML = 0.0
        eD3 = 0.0
        #
        target.scratch.energyTerms["QC Delta-ML correction"]      = eML * Units.Energy_Hartrees_To_Kilojoules_Per_Mole
        target.scratch.energyTerms["QC Dispersion D3 correction"] = eD3 * Units.Energy_Hartrees_To_Kilojoules_Per_Mole

        doGradients  = target.scratch.doGradients
        if doGradients:
            n = len ( target.qcState.atomicNumbers )
            gradients3 = coordinates3.WithExtent ( n )
            gradients3.Set ( 0.0 )
            # . Add gradient corrections to gradients3.


            qcGradients3 = target.scratch.qcGradients3AU
            gradients3.ScatterAdd( 1.0 , qcGradients3 )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ):
            results = self.Energy ( target )
        return [ ( EnergyClosurePriority.QCEnergy, a, "QC Delta-ML and D3 correction" ) ] 

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCDeltaMLModel, self ).SummaryItems ( )
        return items

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

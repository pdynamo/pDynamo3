"""DFT grid integrator."""

from  pCore          import AttributableObject
from .DFTDefinitions import DFTGridAccuracy
from .DFTGrid        import DFTGrid
from .DFTIntegrator  import DFTIntegrator_Integrate

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Functions for a given order.
_DefaultFunctions  = 20
_NumberOfFunctions = { 0 : 1, 1 : 4, 2 : 10 }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DFTGridIntegrator ( AttributableObject ):
    """DFT grid integrator."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "accuracy" : DFTGridAccuracy.Medium ,
                             "inCore"   : False                  } )

    def BuildGrid ( self, target ):
        """Build the grid."""
        grid = DFTGrid.Construct ( self.accuracy, target.qcState.atomicNumbers, target.scratch.qcCoordinates3AU )
        target.scratch.dftGrid = grid
        target.scratch.qcEnergyReport["Quadrature Points"] = grid.numberOfPoints

    def EstimateMemory ( self, qcModel, qcState ):
        """Estimate grid memory (bytes)."""
        p = float ( DFTGrid.EstimatedPoints ( self.accuracy, qcState.atomicNumbers ) )
        m = 32.0 * p # 4 * R64 ( w, x, y, z ) for each point.
        if self.inCore:
            f  = p * float ( len ( qcState.orbitalBases ) ) # . Number of function values without sparsity.
            r  = _NumberOfFunctions.get ( qcModel.functionalModel.order, _DefaultFunctions ) # . Number of reals per value.
            m += f * ( 4.0 + 8.0 * r ) # . 1 I32 + r R64.
        return m

    def Fock ( self, target ):
        """Calculate the grid contributions to the Fock matrices."""
        scratch = target.scratch
        if hasattr ( scratch, "onePDMQ" ):
            dSpin = scratch.onePDMQ.density
            fSpin = scratch.onePDMQ.fock
        else:
            dSpin = None
            fSpin = None
        ( eQuad, rhoQuad ) = DFTIntegrator_Integrate ( target.qcModel.functionalModel              ,
                                                       scratch.dftGrid                             ,
                                                       target.qcState.orbitalBases                 ,
                                                       scratch.qcCoordinates3AU                    ,
                                                       scratch.onePDMP.density                     ,
                                                       dSpin                                       ,
                                                       self.inCore                                 ,
                                                       not target.electronicState.isSpinRestricted ,
                                                       scratch.onePDMP.fock                        ,
                                                       fSpin                                       ,
                                                       None                                        )
        scratch.qcEnergyReport["Quadrature Electron Count"] = ( rhoQuad, "{:.3f}" )
        scratch.qcEnergyReport["Quadrature Energy"        ] = eQuad
        if scratch.dftGrid.hasFunctionData:
            n        = scratch.dftGrid.numberOfFunctionValues
            o        = float ( len ( target.qcState.orbitalBases ) )
            ( s, m ) = scratch.dftGrid.functionByteSize
            scratch.qcEnergyReport["Quadrature BF Values"      ] = n
            scratch.qcEnergyReport["Quadrature BF Sparsity (%)"] = ( ( 1.0 - float ( n ) / ( float ( scratch.dftGrid.numberOfPoints ) * o ) ) * 100.0, "{:.1f}" )
            scratch.qcEnergyReport["Quadrature BF Storage ({:s}B)".format ( m.symbol )] = ( s, "{:.3f}" )
        return eQuad

    def Gradients ( self, target ):
        """Calculate the grid contributions to the gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            if hasattr ( scratch, "onePDMQ" ): dSpin = scratch.onePDMQ.density
            else:                              dSpin = None
            DFTIntegrator_Integrate ( target.qcModel.functionalModel              ,
                                      scratch.dftGrid                             ,
                                      target.qcState.orbitalBases                 ,
                                      scratch.qcCoordinates3AU                    ,
                                      scratch.onePDMP.density                     ,
                                      dSpin                                       ,
                                      self.inCore                                 ,
                                      not target.electronicState.isSpinRestricted ,
                                      None                                        ,
                                      None                                        ,
                                      scratch.qcGradients3AU                      )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Grid Accuracy"           , "{:s}".format ( self.accuracy.name   ) ) ,
                 ( "Save Grid Function Data" , "{:s}".format ( repr ( self.inCore ) ) ) ]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

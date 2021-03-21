"""Full QC/MM electrostatic density model for MNDO QC models."""

from  pCore                             import CrossPairList          , \
                                               logFile                , \
                                               LogFileActive
from  pMolecule.QCModel                 import MNDOQCMMEvaluator
from .NBDefaults                        import _NonUpdatablePairLists , \
                                               _PairListStatistics
from .QCMMElectrostaticModelDensityBase import QCMMElectrostaticModelDensityBase

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelDensityFullMNDO ( QCMMElectrostaticModelDensityBase ):
    """A full QC/MM electrostatic density model for MNDO models."""

    _attributable = dict ( QCMMElectrostaticModelDensityBase._attributable )
    _classLabel   = "MNDO Full Density QC/MM Electrostatic Model"
    _attributable.update ( { "evaluator" : MNDOQCMMEvaluator } )

    def _GetGradientDensity ( self, scratch ):
        """Get the appropriate gradient density."""
        if ( scratch.Get    ( "ci"     , None ) is not None ) and \
           ( scratch.ci.Get ( "dTotalZ", None ) is not None ): density = scratch.ci.dTotalZ
        else:                                                  density = scratch.onePDMP.density            
        return density

    def QCBPGradients ( self, target ):
        """QC/BP gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            state     = getattr ( target, self.__class__._stateName )
            bpCharges = getattr ( state, "bpCharges", None )
            if bpCharges is not None:
                self.evaluator.Gradients ( target.qcState.orbitalBases.functionCenters ,
                                           self._GetGradientDensity ( scratch )        ,
                                           scratch.bpQCMMDerivativeIntegrals           ,
                                           scratch.qcGradients3QCMM                    ,
                                           scratch.bpGradients3                        )

    def QCBPPotentials ( self, target ):
        """QC/BP potentials."""
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if bpCharges is not None:
            scratch        = target.scratch
            pNode          = scratch.GetSetNode ( _NonUpdatablePairLists )
            pairList       = pNode.Get ( "qcbpElectrostatic", None )
            qcCoordinates3 = scratch.qcCoordinates3QCMM
            if pairList is None:
                pairList                = CrossPairList.Full ( qcCoordinates3.rows, None, len ( bpCharges ), None )
                pNode.qcbpElectrostatic = pairList
                sNode                   = scratch.Get ( _PairListStatistics )
                sNode["QC/BP Electrostatic Pairs"] = float ( len ( pairList ) )
            if len ( pairList ) > 0:
                if scratch.doGradients:
                    qcGradients3 = scratch.qcGradients3QCMM
                    bpGradients3 = scratch.bpGradients3
                else:
                    qcGradients3 = None
                    bpGradients3 = None
                ( energy, integrals ) = self.evaluator.Integrals ( target.qcState.mndoParameters                      ,
                                                                   target.qcState.orbitalBases.centerFunctionPointers ,
                                                                   None                                               ,
                                                                   1.0 / self.dielectric                              ,
                                                                   qcCoordinates3                                     ,
                                                                   scratch.bpCoordinates3                             ,
                                                                   bpCharges                                          ,
                                                                   pairList                                           ,
                                                                   scratch.qcmmPotentials                             ,
                                                                   qcGradients3                                       ,
                                                                   bpGradients3                                       )
            scratch.energyTerms["QC/BP Core"] = energy
            scratch.bpQCMMDerivativeIntegrals = integrals

    def QCMMGradients ( self, target ):
        """QC/MM gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            self.evaluator.Gradients ( target.qcState.orbitalBases.functionCenters ,
                                       self._GetGradientDensity ( scratch )        ,
                                       scratch.mmQCMMDerivativeIntegrals           ,
                                       scratch.qcGradients3QCMM                    ,
                                       scratch.gradients3                          )

    def QCMMPotentials ( self, target ):
        """QC/MM potentials."""
        scratch        = target.scratch
        pNode          = scratch.GetSetNode ( _NonUpdatablePairLists )
        pairList       = pNode.Get ( "qcmmElectrostatic", None )
        qcCoordinates3 = scratch.qcCoordinates3QCMM
        if pairList is None:
            mmAtoms  = target.mmState.pureMMAtoms
            pairList = CrossPairList.Full ( qcCoordinates3.rows ,
                                            None                ,
                                            len ( mmAtoms )     ,
                                            mmAtoms             ,
                                            excluded = True     )
            pNode.qcmmElectrostatic = pairList
            sNode                   = scratch.Get ( _PairListStatistics )
            sNode["QC/MM Electrostatic Pairs"] = float ( len ( pairList ) )
        if len ( pairList ) > 0:
            if scratch.doGradients:
                qcGradients3 = scratch.qcGradients3QCMM
                mmGradients3 = scratch.gradients3
            else:
                qcGradients3 = None
                mmGradients3 = None
            ( energy, integrals ) = self.evaluator.Integrals ( target.qcState.mndoParameters                      ,
                                                               target.qcState.orbitalBases.centerFunctionPointers ,
                                                               None                                               ,
                                                               1.0 / self.dielectric                              ,
                                                               qcCoordinates3                                     ,
                                                               target.coordinates3                                ,
                                                               target.mmState.charges                             ,
                                                               pairList                                           ,
                                                               scratch.qcmmPotentials                             ,
                                                               qcGradients3                                       ,
                                                               mmGradients3                                       )
            scratch.energyTerms["QC/MM Core"] = energy
            scratch.mmQCMMDerivativeIntegrals = integrals

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

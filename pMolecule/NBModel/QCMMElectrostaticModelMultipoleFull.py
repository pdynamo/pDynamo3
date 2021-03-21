"""Defines a simple full QC/MM electrostatic multipole model."""

from   pCore                               import CrossPairList          , \
                                                  logFile                , \
                                                  LogFileActive
from   pScientific                         import Units
from  .NBDefaults                          import _NonUpdatablePairLists , \
                                                  _PairListStatistics
from  .PairwiseInteractionFull             import PairwiseInteractionFull
from  .QCMMElectrostaticModelMultipoleBase import QCMMElectrostaticModelMultipoleBase
from ..EnergyModel                         import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelMultipoleFull ( QCMMElectrostaticModelMultipoleBase ):
    """Define a full QC/MM electrostatic model."""

    _classLabel               = "Full Multipole QC/MM Electrostatic Model"
    _pairwiseInteractionClass = PairwiseInteractionFull

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def e ( ): self.QCBPPotentials ( target )
        def f ( ): self.QCBPGradients  ( target )
        closures = super ( QCMMElectrostaticModelMultipoleFull, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, e, "QC/BP Electrostatic Potentials" ) ,
                            ( EnergyClosurePriority.QCGradients, f, "QC/BP Electrostatic Gradients"  ) ] )
        return closures

    def Fock ( self, target ):
        """Energy and Fock matrix contributions."""
        scratch    = target.scratch
        multipoles = scratch.qcmmMultipoles
        potentials = scratch.qcmmPotentials
        target.qcModel.multipoleEvaluator.FockMultipoles ( target, self.multipoleOrder, multipoles )
        eQCMM      = multipoles.Dot ( potentials )
        scratch.energyTerms["QC/MM Electrostatic"] = ( eQCMM * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        target.qcModel.multipoleEvaluator.FockMultipoleDerivatives ( target, self.multipoleOrder, potentials, scratch.onePDMP.fock )
        return eQCMM

    def GetWeightedDensity ( self, target ): 
        """Get the weighted density."""
        scratch = target.scratch
        if scratch.doGradients:
            target.qcModel.multipoleEvaluator.WeightedDensity ( target                  ,
                                                                self.multipoleOrder     ,
                                                                scratch.qcmmPotentials  ,
                                                                scratch.weightedDensity )

    def QCBPGradients ( self, target ):
        """QC/BP gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            pNode    = scratch.Get  ( _NonUpdatablePairLists )
            pairList = pNode.Get ( "qcbpElectrostatic", None )
            if ( pairList is not None ) and ( len ( pairList ) > 0 ):
                state     = getattr ( target, self.__class__._stateName )
                bpCharges = getattr ( state, "bpCharges", None )
                self.pairwiseInteraction.QCMMGradients ( self.multipoleOrder        ,
                                                         scratch.qcmmMultipoles     ,
                                                         bpCharges                  ,
                                                         1.0 / self.dielectric      ,
                                                         scratch.qcCoordinates3QCMM ,
                                                         scratch.bpCoordinates3     ,
                                                         pairList                   ,
                                                         scratch.qcGradients3QCMM   ,
                                                         scratch.bpGradients3       )

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
                self.pairwiseInteraction.QCMMPotentials ( self.multipoleOrder    ,
                                                          bpCharges              ,
                                                          1.0 / self.dielectric  ,
                                                          qcCoordinates3         ,
                                                          scratch.bpCoordinates3 ,
                                                          pairList               ,
                                                          scratch.qcmmPotentials )

    def QCMMGradients ( self, target ):
        """QC/MM gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            pNode    = scratch.Get  ( _NonUpdatablePairLists )
            pairList = pNode.Get ( "qcmmElectrostatic", None )
            if len ( pairList ) > 0:
                self.pairwiseInteraction.QCMMGradients ( self.multipoleOrder        ,
                                                         scratch.qcmmMultipoles     ,
                                                         target.mmState.charges     ,
                                                         1.0 / self.dielectric      ,
                                                         scratch.qcCoordinates3QCMM ,
                                                         target.coordinates3        ,
                                                         pairList                   ,
                                                         scratch.qcGradients3QCMM   ,
                                                         scratch.gradients3         )

    def QCMMPotentials ( self, target ):
        """QC/MM potentials."""
        mmCharges      = target.mmState.charges
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
            self.pairwiseInteraction.QCMMPotentials ( self.multipoleOrder    ,
                                                      mmCharges              ,
                                                      1.0 / self.dielectric  ,
                                                      qcCoordinates3         ,
                                                      target.coordinates3    ,
                                                      pairList               ,
                                                      scratch.qcmmPotentials )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

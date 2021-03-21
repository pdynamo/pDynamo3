"""Defines a cut-off QC/MM electrostatic multipole model."""

from   pCore                               import CrossPairList                             , \
                                                  logFile                                   , \
                                                  LogFileActive
from   pScientific                         import Units
from   pScientific.Arrays                  import Array
from  .ImagePairListContainer              import ImagePairListContainer
from  .NBDefaults                          import _CheckCutOffs                             , \
                                                  _DefaultGeneratorCutOff                   , \
                                                  _DefaultPairwiseInteractionSplineABFSQCMM , \
                                                  _ImageScanData                            , \
                                                  _NonUpdatablePairLists                    , \
                                                  _QCGrid                                   , \
                                                  _QCOccupancy                              , \
                                                  _PairListStatistics                       , \
                                                  _UpdatablePairLists
from  .PairwiseInteractionABFS             import PairwiseInteractionABFS
from  .PairwiseInteractionSplineABFS       import PairwiseInteractionSplineABFS
from  .QCMMElectrostaticModelMultipoleBase import QCMMElectrostaticModelMultipoleBase
from ..EnergyModel                         import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelMultipoleCutOff ( QCMMElectrostaticModelMultipoleBase ):
    """Define a cut-off QC/MM electrostatic model."""

    _attributable             = dict ( QCMMElectrostaticModelMultipoleBase._attributable )
    _classLabel               = "CutOff Multipole QC/MM Electrostatic Model"
    _pairwiseInteractionClass = ( PairwiseInteractionABFS, PairwiseInteractionSplineABFS )
    _summarizable             = dict ( QCMMElectrostaticModelMultipoleBase._summarizable )
    _attributable.update ( { "generator" : None } )
    _summarizable.update ( { "generator" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        if self.generator is None:
            self.generator = _DefaultGeneratorCutOff ( )
        elif not isinstance ( self.generator, PairListGenerator ):
            raise TypeError ( "Invalid pairlist generator attribute." )
        if self.pairwiseInteraction is None:
            self.pairwiseInteraction = _DefaultPairwiseInteractionSplineABFSQCMM ( )
        elif not isinstance ( pairwiseInteraction, self.__class__._pairwiseInteractionClass ):
            raise TypeError ( "Invalid pairwise interaction attribute." )              
        _CheckCutOffs ( self )
        return self

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def e ( ): self.QCBPPotentials      ( target )
        def f ( ): self.QCBPPotentialsImage ( target )
        def g ( ): self.QCMMPotentialsImage ( target )
        def h ( ): self.QCBPGradients       ( target )
        def i ( ): self.QCBPGradientsImage  ( target )
        def j ( ): self.QCMMGradientsImage  ( target )
        closures = super ( QCMMElectrostaticModelMultipoleCutOff, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals , e, "QC/BP Electrostatic Potentials"       ) ,
                            ( EnergyClosurePriority.QCIntegrals , f, "QC/BP Image Electrostatic Potentials" ) ,
                            ( EnergyClosurePriority.QCIntegrals , g, "QC/MM Image Electrostatic Potentials" ) ,
                            ( EnergyClosurePriority.QCGradients , h, "QC/BP Electrostatic Gradients"        ) ,
                            ( EnergyClosurePriority.QCGradients , i, "QC/BP Image Electrostatic Gradients"  ) ,
                            ( EnergyClosurePriority.QCGradients , j, "QC/MM Image Electrostatic Gradients"  ) ] )
        return closures

    def EnergyInitialize ( self, target ):
        """Energy initialization"""
        super ( QCMMElectrostaticModelMultipoleCutOff, self ).EnergyInitialize ( target )
        self.ImageUpdateInitialize ( target )

    def Fock ( self, target ):
        """Energy and Fock matrix contributions for both QC/MM and QC/QC models."""
        scratch        = target.scratch
        multipoles     = scratch.qcmmMultipoles
        qcmmPotentials = scratch.Get ( "qcmmPotentials", None )
        qcqcPotentials = scratch.Get ( "qcqcPotentials", None )
        target.qcModel.multipoleEvaluator.FockMultipoles ( target, self.multipoleOrder, multipoles )
        dEdQ  = None
        eQCMM = 0.0
        eQCQC = 0.0
        if qcmmPotentials is not None:
            eQCMM = multipoles.Dot ( qcmmPotentials )
            scratch.energyTerms["QC/MM Electrostatic"] = ( eQCMM * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
            dEdQ = qcmmPotentials
        if qcqcPotentials is not None:
            dEdQ  = Array.WithExtent ( len ( multipoles ) )
            qcqcPotentials.VectorMultiply ( multipoles, dEdQ )
            eQCQC = multipoles.Dot ( dEdQ )
            scratch.energyTerms["QC/QC Electrostatic"] = ( eQCQC * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
            dEdQ.Scale ( 2.0 )
            if qcmmPotentials is not None: dEdQ.Add ( qcmmPotentials )
        if dEdQ is not None:
            target.qcModel.multipoleEvaluator.FockMultipoleDerivatives ( target, self.multipoleOrder, dEdQ, scratch.onePDMP.fock )
        return ( eQCMM + eQCQC )

    def GetWeightedDensity ( self, target ): 
        """Get the weighted density."""
        scratch = target.scratch
        if scratch.doGradients:
            qcmmPotentials = scratch.Get ( "qcmmPotentials", None )
            qcqcPotentials = scratch.Get ( "qcqcPotentials", None )
            dEdQ           = None
            if qcmmPotentials is not None:
                dEdQ = qcmmPotentials
            if qcqcPotentials is not None:
                multipoles = scratch.qcmmMultipoles
                dEdQ       = Array.WithExtent ( len ( multipoles ) )
                qcqcPotentials.VectorMultiply ( multipoles, dEdQ )
                dEdQ.Scale ( 2.0 )
                if qcmmPotentials is not None: dEdQ.Add ( qcmmPotentials )
            if dEdQ is not None:
                target.qcModel.multipoleEvaluator.WeightedDensity ( target, self.multipoleOrder, dEdQ, scratch.weightedDensity )

    def ImageUpdateInitialize ( self, target ):
        """Set up data for image calculation."""
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            scratch        = target.scratch
            coordinates3   = scratch.Get ( "coordinates3NB", target.coordinates3 )
            qcCoordinates3 = scratch.qcCoordinates3QCMM
            pNode          = scratch.GetSetNode ( _UpdatablePairLists  )
            grid           = pNode.Get          ( _QCGrid       , None )
            occupancy      = pNode.Get          ( _QCOccupancy  , None )
            if ( grid is None ) or ( occupancy is None ) and self.generator.UseGridSearch ( coordinates3 ):
                ( grid, occupancy ) = qcCoordinates3.MakeGridAndOccupancy ( self.generator.cellSize )
                pNode.Set ( _QCGrid     , grid      )
                pNode.Set ( _QCOccupancy, occupancy )

    def QCBPGradients ( self, target ):
        """QC/BP gradients."""
        scratch   = target.scratch
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if scratch.doGradients and ( bpCharges is not None ):
            pNode    = scratch.Get ( _NonUpdatablePairLists )
            pairList = pNode.Get ( "qcbpElectrostatic", None )
            if ( pairList is not None ) and ( len ( pairList ) > 0 ):
                self.pairwiseInteraction.QCMMGradients ( scratch.qcmmMultipoles     ,
                                                         bpCharges                  ,
                                                         1.0 / self.dielectric      ,
                                                         scratch.qcCoordinates3QCMM ,
                                                         scratch.bpCoordinates3     ,
                                                         pairList                   ,
                                                         scratch.qcGradients3QCMM   ,
                                                         scratch.bpGradients3       )

    def QCBPGradientsImage ( self, target ):
        """QC/BP image gradients."""
        scratch   = target.scratch
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        pNode     = scratch.Get ( _UpdatablePairLists )
        pairList  = pNode.Get ( "qcbpImageElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            self.pairwiseInteraction.QCMMGradientsImage ( scratch.qcmmMultipoles             ,
                                                          bpCharges                          ,
                                                          1.0 / self.dielectric              ,
                                                          scratch.qcCoordinates3QCMM         ,
                                                          scratch.bpCoordinates3             ,
                                                          target.symmetryParameters          ,
                                                          pairList                           ,
                                                          scratch.qcGradients3QCMM           ,
                                                          scratch.bpGradients3               ,
                                                          scratch.symmetryParameterGradients )

    def QCBPPotentials ( self, target ):
        """QC/BP potentials."""
        scratch   = target.scratch
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if bpCharges is not None:
            pNode          = scratch.GetSetNode ( _NonUpdatablePairLists )
            pairList       = pNode.Get ( "qcbpElectrostatic", None )
            coordinates3   = scratch.bpCoordinates3
            qcCoordinates3 = scratch.qcCoordinates3QCMM
            if pairList is None:
                pairList                = CrossPairList.Full ( qcCoordinates3.rows, None, coordinates3.rows, None )
                pNode.qcbpElectrostatic = pairList 
                sNode                   = scratch.Get ( _PairListStatistics )
                sNode["QC/BP Electrostatic Pairs"] = float ( len ( pairList ) )
            if len ( pairList ) > 0:
                self.pairwiseInteraction.QCMMPotentials ( bpCharges              ,
                                                          1.0 / self.dielectric  ,
                                                          qcCoordinates3         ,
                                                          coordinates3           ,
                                                          pairList               ,
                                                          scratch.qcmmPotentials )

    def QCBPPotentialsImage ( self, target ):
        """QC/BP image potentials."""
        scratch            = target.scratch
        state              = getattr ( target, self.__class__._stateName )
        bpCharges          = getattr ( state, "bpCharges", None )
        symmetryParameters = target.symmetryParameters
        if ( bpCharges is not None ) and ( symmetryParameters is not None ):
            pNode          = scratch.GetSetNode ( _UpdatablePairLists )
            pairList       = pNode.Get ( "qcbpImageElectrostatic", None )
            coordinates3   = scratch.bpCoordinates3
            qcCoordinates3 = scratch.qcCoordinates3QCMM
            if pairList is None:
                pairList = ImagePairListContainer.Constructor ( self.generator                      ,
                                                                None                                ,
                                                                None                                ,
                                                                None                                ,
                                                                qcCoordinates3                      ,
                                                                coordinates3                        ,
                                                                symmetryParameters                  ,
                                                                target.symmetry.transformations     ,
                                                                pNode.Get ( _ImageScanData , None ) ,
                                                                pNode.Get ( _QCGrid        , None ) ,
                                                                pNode.Get ( _QCOccupancy   , None ) ,
                                                                False                               )
                pNode.qcbpImageElectrostatic = pairList
                sNode = scratch.Get ( _PairListStatistics )
                nI    = float ( pairList.numberOfImages )
                nP    = float ( pairList.numberOfPairs  )
                sNode[ "QC/BP Image Electrostatic Images" ]  = nI
                sNode[ "QC/BP Image Electrostatic Pairs"  ]  = nP
                sNode["<QC/BP Image Electrostatic Images>"] += nI
                sNode["<QC/BP Image Electrostatic Pairs>" ] += nP
            if len ( pairList ) > 0:
                self.pairwiseInteraction.QCMMPotentialsImage ( bpCharges              ,
                                                               1.0 / self.dielectric  ,
                                                               qcCoordinates3         ,
                                                               coordinates3           ,
                                                               symmetryParameters     ,
                                                               pairList               ,
                                                               scratch.qcmmPotentials )

    def QCMMGradients ( self, target ):
        """QC/MM gradients."""
        scratch  = target.scratch
        pNode    = scratch.Get ( _UpdatablePairLists )
        pairList = pNode.Get ( "qcmmElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            self.pairwiseInteraction.QCMMGradients ( scratch.qcmmMultipoles     ,
                                                     target.mmState.charges     ,
                                                     1.0 / self.dielectric      ,
                                                     scratch.qcCoordinates3QCMM ,
                                                     coordinates3               ,
                                                     pairList                   ,
                                                     scratch.qcGradients3QCMM   ,
                                                     scratch.gradients3         )

    def QCMMGradientsImage ( self, target ):
        """QC/MM image gradients."""
        scratch  = target.scratch
        pNode    = scratch.Get ( _UpdatablePairLists )
        pairList = pNode.Get ( "qcmmImageElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            self.pairwiseInteraction.QCMMGradientsImage ( scratch.qcmmMultipoles             ,
                                                          target.mmState.charges             ,
                                                          1.0 / self.dielectric              ,
                                                          scratch.qcCoordinates3QCMM         ,
                                                          coordinates3                       ,
                                                          target.symmetryParameters          ,
                                                          pairList                           ,
                                                          scratch.qcGradients3QCMM           ,
                                                          scratch.gradients3                 ,
                                                          scratch.symmetryParameterGradients )

    def QCMMPotentials ( self, target ):
        """QC/MM potentials."""
        scratch        = target.scratch
        pNode          = scratch.GetSetNode ( _UpdatablePairLists )
        pairList       = pNode.Get ( "qcmmElectrostatic", None )
        coordinates3   = scratch.Get ( "coordinates3NB", target.coordinates3 )
        qcCoordinates3 = scratch.qcCoordinates3QCMM
        if pairList is None:
            pairList = self.generator.CrossPairListFromDoubleCoordinates3 ( qcCoordinates3                    ,
                                                                            coordinates3                      ,
                                                                            None                              ,
                                                                            None                              ,
                                                                            None                              ,
                                                                            target.mmState.pureMMAtoms        ,
                                                                            None                              ,
                                                                            None                              ,
                                                                            None                              ,
                                                                            pNode.Get ( _QCGrid      , None ) ,
                                                                            pNode.Get ( _QCOccupancy , None ) )
            pNode.qcmmElectrostatic = pairList
            sNode                   = scratch.Get ( _PairListStatistics )
            n                       = float ( pairList.numberOfPairs )
            sNode[ "QC/MM Electrostatic Pairs" ]  = n
            sNode["<QC/MM Electrostatic Pairs>"] += n
        if len ( pairList ) > 0:
            self.pairwiseInteraction.QCMMPotentials ( target.mmState.charges ,
                                                      1.0 / self.dielectric  ,
                                                      qcCoordinates3         ,
                                                      coordinates3           ,
                                                      pairList               ,
                                                      scratch.qcmmPotentials )

    def QCMMPotentialsImage ( self, target ):
        """QC/MM image potentials."""
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            scratch        = target.scratch
            pNode          = scratch.GetSetNode ( _UpdatablePairLists )
            pairList       = pNode.Get ( "qcmmImageElectrostatic", None )
            coordinates3   = scratch.Get ( "coordinates3NB", target.coordinates3 )
            qcCoordinates3 = scratch.qcCoordinates3QCMM
            if pairList is None:
                pairList = ImagePairListContainer.Constructor ( self.generator                      ,
                                                                None                                ,
                                                                target.mmState.pureMMAtoms          ,
                                                                None                                ,
                                                                qcCoordinates3                      ,
                                                                coordinates3                        ,
                                                                symmetryParameters                  ,
                                                                target.symmetry.transformations     ,
                                                                pNode.Get ( _ImageScanData , None ) ,
                                                                pNode.Get ( _QCGrid        , None ) ,
                                                                pNode.Get ( _QCOccupancy   , None ) ,
                                                                False                               )
                pNode.qcmmImageElectrostatic = pairList
                sNode = scratch.Get ( _PairListStatistics )
                nI    = float ( pairList.numberOfImages )
                nP    = float ( pairList.numberOfPairs  )
                sNode[ "QC/MM Image Electrostatic Images" ]  = nI
                sNode[ "QC/MM Image Electrostatic Pairs"  ]  = nP
                sNode["<QC/MM Image Electrostatic Images>"] += nI
                sNode["<QC/MM Image Electrostatic Pairs>" ] += nP
            if len ( pairList ) > 0:
                self.pairwiseInteraction.QCMMPotentialsImage ( target.mmState.charges ,
                                                               1.0 / self.dielectric  ,
                                                               qcCoordinates3         ,
                                                               coordinates3           ,
                                                               symmetryParameters     ,
                                                               pairList               ,
                                                               scratch.qcmmPotentials )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

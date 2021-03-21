"""Defines NB model classes that employ cut-offs."""

from   pCore                                 import logFile                         , \
                                                    LogFileActive                   , \
                                                    SelfPairList
from   pScientific.Geometry3                 import Coordinates3                    , \
                                                    PairListGenerator
from  .ImagePairListContainer                import ImagePairListContainer
from  .ImageScanContainer                    import ImageScanContainer
from  .NBDefaults                            import _CenteringTranslation3          , \
                                                    _CheckCutOffs                   , \
                                                    _DefaultGeneratorCutOff         , \
                                                    _DefaultPairwiseInteractionABFS , \
                                                    _ImageScanData                  , \
                                                    _MMGrid                         , \
                                                    _MMOccupancy                    , \
                                                    _NonUpdatablePairLists          , \
                                                    _PairListStatistics             , \
                                                    _UpdatablePairLists
from  .NBModel                               import NBModel
from  .PairwiseInteractionABFS               import PairwiseInteractionABFS
from  .PairwiseInteractionSplineABFS         import PairwiseInteractionSplineABFS
from  .QCMMElectrostaticModelMultipoleCutOff import QCMMElectrostaticModelMultipoleCutOff
from  .QCMMLennardJonesModelCutOff           import QCMMLennardJonesModelCutOff
from  .QCQCElectrostaticModelMultipoleCutOff import QCQCElectrostaticModelMultipoleCutOff
from  .QCQCLennardJonesModelCutOff           import QCQCLennardJonesModelCutOff
from  .UpdateChecker                         import UpdateChecker
from ..EnergyModel                           import EnergyClosurePriority

# . Currently changing the PLG or PWI options requires a manual clearance of the updatable pairlist node in scratch.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModelCutOff ( NBModel ):
    """The base class for NB models that employ cut-offs."""

    # . Check for inverses and image expand factor are for testing.
    # . Update checker need never be set explicitly.
    _attributable             = dict ( NBModel._attributable )
    _classLabel               = "CutOff NB Model"
    _pairwiseInteractionClass = ( PairwiseInteractionABFS, PairwiseInteractionSplineABFS )
    _summarizable             = dict ( NBModel._summarizable )
    _attributable.update ( { "checkForInverses"  : True ,
                             "generator"         : None ,
                             "imageExpandFactor" : 0    ,
                             "useCentering"      : True ,
                             "updateChecker"     : None } )
    _summarizable.update ( { "generator"         : None            ,
                             "useCentering"      : "Use Centering" } )

    def _CheckOptions ( self ):
        """Check options."""
        if self.generator is None:
            self.generator = _DefaultGeneratorCutOff ( )
        elif not isinstance ( self.generator, PairListGenerator ):
            raise TypeError ( "Invalid pairlist generator attribute." )
        if self.pairwiseInteraction is None:
            self.pairwiseInteraction = _DefaultPairwiseInteractionABFS ( )
        elif not isinstance ( self.pairwiseInteraction, self.__class__._pairwiseInteractionClass ):
            raise TypeError ( "Invalid pairwise interaction attribute." )              
        _CheckCutOffs ( self )
        return self

    def CenterCoordinates ( self, target ):
        """Center the coordinates."""
        symmetryParameters = target.symmetryParameters
        if self.useCentering and ( symmetryParameters is not None ):
            scratch         = target.scratch
            coordinates3    = target.coordinates3
            nbCoordinates3  = scratch.Get ( "coordinates3NB"      , None ) # . Centered coordinates.
            translation3    = scratch.Get ( _CenteringTranslation3, None )
            redoTranslation = ( nbCoordinates3 is None ) or ( translation3 is None ) or ( not hasattr ( scratch, _UpdatablePairLists ) )
            if nbCoordinates3 is None:
                nbCoordinates3         = Coordinates3.WithExtent ( coordinates3.rows )
                scratch.coordinates3NB = nbCoordinates3
            coordinates3.CopyTo ( nbCoordinates3 )
            if redoTranslation:
                if translation3 is None:
                    translation3 = Coordinates3.WithExtent ( coordinates3.rows )
                    scratch.Set ( _CenteringTranslation3, translation3 )
                if not hasattr ( scratch, "mmIsolates" ): target.mmState.MakeMMIsolates ( )
                freeIsolates = scratch.Get ( "freeMMIsolates", None )
                isolates     = scratch.mmIsolates
                symmetryParameters.CenterCoordinates3ByFreeIsolate ( isolates, freeIsolates, nbCoordinates3 )
                nbCoordinates3.CopyTo ( translation3 )
                translation3.Add      ( coordinates3 , scale = -1.0 )
            else:
                nbCoordinates3.Add ( translation3 )

    def CheckForUpdate ( self, target ):
        """Check for an update."""
        if self.updateChecker is None:
            buffer = self.generator.cutOff - self.pairwiseInteraction.range
            self.updateChecker = UpdateChecker.WithOptions ( buffer = buffer )
        self.updateChecker.Check ( target )

    def Energy ( self, target ):
        """Energy 1-5+."""
        energies     = {}
        scratch      = target.scratch
        coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
        pNode        = scratch.GetSetNode ( _UpdatablePairLists )
        pairList     = pNode.Get ( "mmmm", None )
        if pairList is None:
            pairList = self.generator.SelfPairListFromCoordinates3 ( coordinates3                       ,
                                                                     None                               ,
                                                                     target.mmState.mmAtoms             ,
                                                                     target.freeAtoms                   ,
                                                                     target.mmState.exclusions          ,
                                                                     pNode.Get ( _MMGrid       , None ) ,
                                                                     pNode.Get ( _MMOccupancy  , None ) )
            pNode.mmmm = pairList
            sNode      = scratch.Get ( _PairListStatistics )
            n          = float ( len ( pairList ) )
            sNode[ "MM/MM Pairs" ]  = n
            sNode["<MM/MM Pairs>"] += n
        if len ( pairList ) > 0:
            gradients3 = scratch.Get ( "gradients3", None )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( target.mmState.charges       ,
                                                                                      target.mmState.charges       ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljParameters  ,
                                                                                      ( 1.0 / self.dielectric )    ,
                                                                                        1.0                        ,
                                                                                      coordinates3                 ,
                                                                                      coordinates3                 ,
                                                                                      pairList                     ,
                                                                                      gradients3                   ,
                                                                                      gradients3                   )
            energies.update ( { "MM/MM Electrostatic" : eElectrostatic ,
                                "MM/MM Lennard-Jones" : eLennardJones  } )
        return energies

    def Energy14 ( self, target ):
        """Energy 1-4."""
        energies = {}
        scratch  = target.scratch
        pNode    = scratch.GetSetNode ( _NonUpdatablePairLists )
        pairList = pNode.Get ( "mmmm14", None )
        if pairList is None:
            pairList = SelfPairList.FromSelfPairList ( target.mmState.interactions14 ,
                                                       len ( target.atoms )          ,
                                                       target.mmState.mmAtoms        ,
                                                       target.freeAtoms              )
            pNode.mmmm14 = pairList
            sNode        = scratch.Get ( _PairListStatistics )
            sNode["MM/MM 1-4 Pairs"] = float ( len ( pairList ) )
        if pairList is not None:
            coordinates3 = scratch.Get ( "coordinates3NB" , target.coordinates3 )
            gradients3   = scratch.Get ( "gradients3"     , None                )
            scale        = ( target.mmModel.electrostaticScale14 / self.dielectric )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( target.mmState.charges        ,
                                                                                      target.mmState.charges        ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljParameters14 ,
                                                                                      scale                         ,
                                                                                      1.0                           ,
                                                                                      coordinates3                  ,
                                                                                      coordinates3                  ,
                                                                                      pairList                      ,
                                                                                      gradients3                    ,
                                                                                      gradients3                    )
            energies.update ( { "MM/MM 1-4 Electrostatic" : eElectrostatic ,
                                "MM/MM 1-4 Lennard-Jones" : eLennardJones  } )
        return energies

    # . Check whether target has symmetry and exclude EnergyImage if not?
    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): target.scratch.energyTerms.update ( self.Energy      ( target ) )
        def b ( ): target.scratch.energyTerms.update ( self.Energy14    ( target ) )
        def c ( ): target.scratch.energyTerms.update ( self.EnergyImage ( target ) )
        closures = super ( NBModelCutOff, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.IndependentEnergyTerm, a, "MM/MM NB Evaluation"       ) ,
                            ( EnergyClosurePriority.IndependentEnergyTerm, b, "MM/MM 1-4 NB Evaluation"   ) ,
                            ( EnergyClosurePriority.IndependentEnergyTerm, c, "MM/MM Image NB Evaluation" ) ] )
        return closures

    def EnergyImage ( self, target ):
        """Image energy."""
        energies           = {}
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            scratch      = target.scratch
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            pNode        = scratch.GetSetNode ( _UpdatablePairLists )
            pairList     = pNode.Get ( "mmmmImage", None )
            if pairList is None:
                pairList = ImagePairListContainer.Constructor ( self.generator                      ,
                                                                target.mmState.mmAtoms              ,
                                                                target.mmState.mmAtoms              ,
                                                                target.freeAtoms                    ,
                                                                coordinates3                        ,
                                                                coordinates3                        ,
                                                                symmetryParameters                  ,
                                                                target.symmetry.transformations     ,
                                                                pNode.Get ( _ImageScanData , None ) ,
                                                                pNode.Get ( _MMGrid        , None ) ,
                                                                pNode.Get ( _MMOccupancy   , None ) ,
                                                                self.checkForInverses               )
                pNode.mmmmImage = pairList
                sNode           = scratch.Get ( _PairListStatistics )
                nI              = float ( pairList.numberOfImages )
                nP              = float ( pairList.numberOfPairs  )
                sNode[ "MM/MM Image Images" ]  = nI
                sNode[ "MM/MM Image Pairs"  ]  = nP
                sNode["<MM/MM Image Images>"] += nI
                sNode["<MM/MM Image Pairs>" ] += nP
            if len ( pairList ) > 0:
                ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergyImage ( target.mmState.charges       ,
                                                                                               target.mmState.ljTypeIndices ,
                                                                                               target.mmState.ljParameters  ,
                                                                                               1.0 / self.dielectric        ,
                                                                                               coordinates3                 ,
                                                                                               symmetryParameters           ,
                                                                                               pairList                     ,
                                                                                               scratch.Get ( "gradients3"                , None ) ,
                                                                                               scratch.Get ( "symmetryParameterGradients", None ) )
                energies.update ( { "MM/MM Image Electrostatic" : eElectrostatic ,
                                    "MM/MM Image Lennard-Jones" : eLennardJones  } )
        return energies

    def EnergyInitialize ( self, target ):
        """Energy initialization"""
        super ( NBModelCutOff, self ).EnergyInitialize ( target )
        self.CheckForUpdate        ( target )
        self.CenterCoordinates     ( target )
        self.ImageUpdateInitialize ( target )

    def ImageUpdateInitialize ( self, target ):
        """Set up data for image calculation."""
        # . MM grids and scan data.
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            scratch      = target.scratch
            coordinates3 = scratch.Get        ( "coordinates3NB", target.coordinates3 )
            pNode        = scratch.GetSetNode ( _UpdatablePairLists )
            grid         = pNode.Get          ( _MMGrid       , None )
            occupancy    = pNode.Get          ( _MMOccupancy  , None )
            scanData     = pNode.Get          ( _ImageScanData, None )
            if ( ( grid is None ) or ( occupancy is None ) ) and self.generator.UseGridSearch ( coordinates3 ):
                ( grid, occupancy ) = coordinates3.MakeGridAndOccupancy ( self.generator.cellSize )
                pNode.Set ( _MMGrid     , grid      )
                pNode.Set ( _MMOccupancy, occupancy )
            if scanData is None:
                scanData = ImageScanContainer.Constructor ( coordinates3                    ,
                                                            symmetryParameters              ,
                                                            target.symmetry.transformations ,
                                                            self.generator.cutOff           ,
                                                            self.checkForInverses           ,
                                                            self.imageExpandFactor          )
                pNode.Set ( _ImageScanData, scanData )

    def QCMMModels ( self, qcModel = None, withSymmetry = False ):
        """Default companion QC/MM models for the model."""
        models = { "qcmmElectrostatic" : QCMMElectrostaticModelMultipoleCutOff ,
                   "qcmmLennardJones"  : QCMMLennardJonesModelCutOff           }
        if withSymmetry:
            models["qcqcElectrostatic"] = QCQCElectrostaticModelMultipoleCutOff
            models["qcqcLennardJones" ] = QCQCLennardJonesModelCutOff
        return models

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

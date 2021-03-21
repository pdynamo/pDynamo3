"""Cut-off QC/MM electrostatic density model for MNDO QC models."""

# . QC/QC interactions are handled using a zeroth-order multipole model which is simple
#   and symmetric. A QC(density)/QC(multipole) interaction could be implemented but is
#   more complicated and asymmetric as there are two types of QC atoms in the interaction
#   (although this might not be a problem). The full (and symmetric) QC/QC density model
#   could also be implemented but is way more complicated (need TEI contributions as
#   well) and probably not worth it (unless someone really wants to do it as they
#   want look at crystals, etc. ...).

from   pCore                             import CrossPairList           , \
                                                logFile                 , \
                                                LogFileActive
from   pMolecule.QCModel                 import _MNDOQCMMTermLabels     , \
                                                _MNDOQCMMTerms
from   pScientific                       import Units
from   pScientific.Arrays                import Array
from   pScientific.Splines               import CubicSpline             , \
                                                CubicSplineContainer  
from  .ABFSIntegrator                    import ABFSIntegrator
from  .ImagePairListContainer            import ImagePairListContainer
from  .MNDOQCMMImageEvaluator            import MNDOQCMMImageEvaluator
from  .NBDefaults                        import _DefaultGeneratorCutOff , \
                                                _ImageScanData          , \
                                                _NonUpdatablePairLists  , \
                                                _PairListStatistics     , \
                                                _QCGrid                 , \
                                                _QCOccupancy            , \
                                                _UpdatablePairLists
from  .NBModelError                      import NBModelError
from  .QCMMElectrostaticModelDensityBase import QCMMElectrostaticModelDensityBase
from ..EnergyModel                       import EnergyClosurePriority

# . Splines are in atomic units.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelDensityCutOffMNDO ( QCMMElectrostaticModelDensityBase ):
    """A cut-off QC/MM electrostatic density model for MNDO models."""

    _attributable = dict ( QCMMElectrostaticModelDensityBase._attributable )
    _classLabel   = "MNDO CutOff Density QC/MM Electrostatic Model"
    _attributable.update ( { "evaluator"  : MNDOQCMMImageEvaluator ,
                             "generator"  : None                   ,
                             "integrator" : None                   } )

    def _GetGradientDensity ( self, scratch ):
        """Get the appropriate gradient density."""
        if ( scratch.Get    ( "ci"     , None ) is not None ) and \
           ( scratch.ci.Get ( "dTotalZ", None ) is not None ): density = scratch.ci.dTotalZ
        else:                                                  density = scratch.onePDMP.density            
        return density

    def BuildModel ( self, target ):
        """Build the model."""
        state         = super ( QCMMElectrostaticModelDensityCutOffMNDO, self ).BuildModel ( target )
        state.splines = self.MakeSplines ( target )
        return state

    def _CheckOptions ( self ):
        """Check options."""
        if self.generator is None:
            self.generator = _DefaultGeneratorCutOff ( )
        elif not isinstance ( self.generator, PairListGenerator ):
            raise TypeError ( "Invalid pairlist generator attribute." )
        integrator = self.integrator
        if integrator is None:
            integrator = ABFSIntegrator.WithDefaults ( )
        elif not instance ( self.integrator, ABFSIntegrator ):
            raise TypeError ( "Invalid ABFS integrator attribute." )
        # . Check cutoffs.
        if ( self.generator.cutOff - integrator.outerCutOff ) < 0.0:
            raise NBModelError ( "Incompatible generator/ABFS integrator interaction cut-offs." )
        # . If default integrator convert to atomic units.
        if self.integrator is None:
            integrator.dampingCutOff *= Units.Length_Angstroms_To_Bohrs
            integrator.innerCutOff   *= Units.Length_Angstroms_To_Bohrs
            integrator.outerCutOff   *= Units.Length_Angstroms_To_Bohrs
            self.integrator = integrator
        return self

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def c ( ): self.QCBPPotentialsImage ( target )
        def d ( ): self.QCMMPotentialsImage ( target )
        def e ( ): self.QCBPGradientsImage  ( target )
        def f ( ): self.QCMMGradientsImage  ( target )
        closures = super ( QCMMElectrostaticModelDensityCutOffMNDO, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, c, "QC/BP Image Electrostatic Potentials" ) ,
                            ( EnergyClosurePriority.QCIntegrals, d, "QC/MM Image Electrostatic Potentials" ) ,
                            ( EnergyClosurePriority.QCGradients, e, "QC/BP Image Electrostatic Gradients"  ) ,
                            ( EnergyClosurePriority.QCGradients, f, "QC/MM Image Electrostatic Gradients"  ) ] )
        if hasattr ( target.qcModel.multipoleEvaluator, "WeightedDensity" ):
           def g ( ): self.GetWeightedDensity ( target )
           closures.append ( ( EnergyClosurePriority.QCPreGradients, g, "QC/MM Weighted Density" ) )
        return closures

    def EnergyInitialize ( self, target ):
        """Energy initialization"""
        super ( QCMMElectrostaticModelDensityCutOffMNDO, self ).EnergyInitialize ( target )
        self.ImageUpdateInitialize ( target ) # . Not the same as in other modules as includes multipole allocation.

    def Fock ( self, target ):
        """Energy and Fock matrix contributions for both QC/MM and QC/QC models."""
        eQCMM          = super ( QCMMElectrostaticModelDensityCutOffMNDO, self ).Fock ( target )
        scratch        = target.scratch
        qcqcPotentials = scratch.Get ( "qcqcPotentials", None )
        if qcqcPotentials is not None:
            multipoles = scratch.qcmmMultipoles
            dEdQ       = Array.WithExtent ( len ( multipoles ) )
            target.qcModel.multipoleEvaluator.FockMultipoles ( target, 0, multipoles )
            qcqcPotentials.VectorMultiply ( multipoles, dEdQ )
            eQCQC      = multipoles.Dot ( dEdQ )
            eQCMM     += eQCQC
            scratch.energyTerms["QC/QC Electrostatic"] = ( eQCQC * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
            dEdQ.Scale ( 2.0 )
            target.qcModel.multipoleEvaluator.FockMultipoleDerivatives ( target, 0, dEdQ, scratch.onePDMP.fock )
        return eQCMM

    def GetWeightedDensity ( self, target ): 
        """Get the weighted density."""
        scratch = target.scratch
        if scratch.doGradients:
            qcqcPotentials = scratch.Get ( "qcqcPotentials", None )
            if qcqcPotentials is not None:
                multipoles = scratch.qcmmMultipoles
                dEdQ       = Array.WithExtent ( len ( multipoles ) )
                qcqcPotentials.VectorMultiply ( multipoles, dEdQ )
                dEdQ.Scale ( 2.0 )
                target.qcModel.multipoleEvaluator.WeightedDensity ( target, 0, dEdQ, scratch.weightedDensity )

    def ImageUpdateInitialize ( self, target ):
        """Set up data for image calculation."""
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            # . Multipole data.
            multipoles = target.scratch.Get ( "qcmmMultipoles", None )
            if multipoles is None:
                multipoles = Array.WithExtent ( len ( target.qcState.qcAtoms ) )
                target.scratch.qcmmMultipoles = multipoles
            multipoles.Set ( 0.0 )
            # . Update data.
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

    def IntegratorClosure ( self, evaluator, qData, i, derivative ):
        """Get the function for integration."""
        def F ( r ):
            return evaluator.IntegralValue ( qData, r, i, derivative = derivative )
        return F

    def MakeSplines ( self, target, testing = False ):
        """Make the interaction splines."""
        # . Gather data.
        parameters = target.qcState.mndoParameters.uniqueEntries
        pInfo      = [ ( n, _MNDOQCMMTerms[parameters[n].numberOfOrbitals] ) for n in sorted ( parameters.keys ( ) ) ]
        # . Integration.
        x       = self.integrator.x
        results = {}
        for ( n, t ) in pInfo:
            qData = parameters[n]
            if testing: terms = []
            else:       y     = Array.WithExtents ( len ( x ), t )
            for i in range ( t ):
                F     = self.IntegratorClosure ( self.evaluator, qData, i, False )
                G     = self.IntegratorClosure ( self.evaluator, qData, i, True  )
                yABFS = self.integrator.Integrate ( F, G )
                if testing:
                    spline  = CubicSpline.FromArrays ( x, yT, lowerDerivative = 1, lowerValue = 0.0, upperDerivative = 1, upperValue = 0.0 )
                    yActual = [               F ( r )    for r in x ]
                    ySpline = [ spline.Evaluate ( r )[0] for r in x ]
                    terms.append ( ( _MNDOQCMMTermLabels[t], yActual, yABFS, ySpline ) ) # . Term label with actual, ABFS and spline values.
                else:
                    yABFS.CopyTo ( y[:,i] )
            if testing:
                results[n] = terms
            else:
                results[n] = CubicSpline.FromArrays ( x, y, lowerDerivative = 1, lowerValue = 0.0, upperDerivative = 1, upperValue = 0.0 )
        # . Finish up.
        if testing:
            return results
        else:
            return CubicSplineContainer.FromUniqueEntries ( results, target.qcState.atomicNumbers, label = "MNDO QC/MM Interaction Splines" )

    def QCBPGradients ( self, target ):
        """QC/BP gradients."""
        scratch   = target.scratch
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if scratch.doGradients and ( bpCharges is not None ):
            pNode    = scratch.Get ( _NonUpdatablePairLists )
            pairList = pNode.Get ( "qcbpElectrostatic", None )
            if ( pairList is not None ) and ( len ( pairList ) > 0 ):
                self.evaluator.Gradients ( target.qcState.orbitalBases.functionCenters ,
                                           self._GetGradientDensity ( scratch )        ,
                                           scratch.bpQCMMDerivativeIntegrals           ,
                                           scratch.qcGradients3QCMM                    ,
                                           scratch.bpGradients3                        )

    def QCBPGradientsImage ( self, target ):
        """QC/BP image gradients."""
        scratch   = target.scratch
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        pNode     = scratch.Get ( _UpdatablePairLists )
        pairList  = pNode.Get ( "qcbpImageElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            self.evaluator.GradientsImage ( target.qcState.orbitalBases.functionCenters ,
                                            self._GetGradientDensity ( scratch )        ,
                                            scratch.qcCoordinates3QCMM                  ,
                                            scratch.bpCoordinates3                      ,
                                            target.symmetryParameters                   ,
                                            pairList                                    ,
                                            scratch.bpQCMMImageDerivativeIntegrals      ,
                                            scratch.qcGradients3QCMM                    ,
                                            scratch.bpGradients3                        ,
                                            scratch.symmetryParameterGradients          )

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
                                                                   state.splines                                      ,
                                                                   1.0 / self.dielectric                              ,
                                                                   qcCoordinates3                                     ,
                                                                   scratch.bpCoordinates3                             ,
                                                                   bpCharges                                          ,
                                                                   pairList                                           ,
                                                                   scratch.qcmmPotentials                             ,
                                                                   qcGradients3                                       ,
                                                                   bpGradients3                                       ,
                                                                   cutOff = self.integrator.outerCutOff               )
                scratch.energyTerms["QC/BP Core"] = energy
                if integrals is None: scratch.Pop ( "bpQCMMDerivativeIntegrals" )
                else:                 scratch.bpQCMMDerivativeIntegrals = integrals

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
                if scratch.doGradients:
                    qcGradients3               = scratch.qcGradients3QCMM
                    bpGradients3               = scratch.bpGradients3
                    symmetryParameterGradients = scratch.symmetryParameterGradients
                else:
                    qcGradients3               = None
                    bpGradients3               = None
                    symmetryParameterGradients = None
                ( energy, integrals ) = self.evaluator.IntegralsImage ( target.qcState.mndoParameters                      ,
                                                                        target.qcState.orbitalBases.centerFunctionPointers ,
                                                                        state.splines                                      ,
                                                                        self.integrator.outerCutOff                        ,
                                                                        1.0 / self.dielectric                              ,
                                                                        qcCoordinates3                                     ,
                                                                        coordinates3                                       ,
                                                                        symmetryParameters                                 ,
                                                                        bpCharges                                          ,
                                                                        pairList                                           ,
                                                                        scratch.qcmmPotentials                             ,
                                                                        qcGradients3                                       ,
                                                                        bpGradients3                                       ,
                                                                        symmetryParameterGradients                         )
                scratch.energyTerms["QC/BP Image Core"] = energy
                if integrals is None: scratch.Pop ( "bpQCMMImageDerivativeIntegrals" )
                else:                 scratch.bpQCMMImageDerivativeIntegrals = integrals

    def QCMMGradients ( self, target ):
        """QC/MM gradients."""
        scratch  = target.scratch
        pNode    = scratch.Get ( _UpdatablePairLists )
        pairList = pNode.Get ( "qcmmElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            self.evaluator.Gradients ( target.qcState.orbitalBases.functionCenters ,
                                       self._GetGradientDensity ( scratch )        ,
                                       scratch.mmQCMMDerivativeIntegrals           ,
                                       scratch.qcGradients3QCMM                    ,
                                       scratch.gradients3                          )

    def QCMMGradientsImage ( self, target ):
        """QC/MM image gradients."""
        scratch  = target.scratch
        pNode    = scratch.Get ( _UpdatablePairLists )
        pairList = pNode.Get ( "qcmmImageElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            self.evaluator.GradientsImage ( target.qcState.orbitalBases.functionCenters ,
                                            self._GetGradientDensity ( scratch )        ,
                                            scratch.qcCoordinates3QCMM                  ,
                                            coordinates3                                ,
                                            target.symmetryParameters                   ,
                                            pairList                                    ,
                                            scratch.mmQCMMImageDerivativeIntegrals      ,
                                            scratch.qcGradients3QCMM                    ,
                                            scratch.gradients3                          ,
                                            scratch.symmetryParameterGradients          )

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
            if scratch.doGradients:
                qcGradients3 = scratch.qcGradients3QCMM
                mmGradients3 = scratch.gradients3
            else:
                qcGradients3 = None
                mmGradients3 = None
            state = getattr ( target, self.__class__._stateName )
            ( energy, integrals ) = self.evaluator.Integrals ( target.qcState.mndoParameters                      ,
                                                               target.qcState.orbitalBases.centerFunctionPointers ,
                                                               state.splines                                      ,
                                                               1.0 / self.dielectric                              ,
                                                               qcCoordinates3                                     ,
                                                               target.coordinates3                                ,
                                                               target.mmState.charges                             ,
                                                               pairList                                           ,
                                                               scratch.qcmmPotentials                             ,
                                                               qcGradients3                                       ,
                                                               mmGradients3                                       ,
                                                               cutOff = self.integrator.outerCutOff               )
            scratch.energyTerms["QC/MM Core"] = energy
            if integrals is None: scratch.Pop ( "mmQCMMDerivativeIntegrals" )
            else:                 scratch.mmQCMMDerivativeIntegrals = integrals

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
                if scratch.doGradients:
                    qcGradients3               = scratch.qcGradients3QCMM
                    mmGradients3               = scratch.gradients3
                    symmetryParameterGradients = scratch.symmetryParameterGradients
                else:
                    qcGradients3               = None
                    mmGradients3               = None
                    symmetryParameterGradients = None
                state = getattr ( target, self.__class__._stateName )
                ( energy, integrals ) = self.evaluator.IntegralsImage ( target.qcState.mndoParameters                      ,
                                                                        target.qcState.orbitalBases.centerFunctionPointers ,
                                                                        state.splines                                      ,
                                                                        self.integrator.outerCutOff                        ,
                                                                        1.0 / self.dielectric                              ,
                                                                        qcCoordinates3                                     ,
                                                                        coordinates3                                       ,
                                                                        symmetryParameters                                 ,
                                                                        target.mmState.charges                             ,
                                                                        pairList                                           ,
                                                                        scratch.qcmmPotentials                             ,
                                                                        qcGradients3                                       ,
                                                                        mmGradients3                                       ,
                                                                        symmetryParameterGradients                         )
                scratch.energyTerms["QC/MM Image Core"] = energy
                if integrals is None: scratch.Pop ( "mmQCMMImageDerivativeIntegrals" )
                else:                 scratch.mmQCMMImageDerivativeIntegrals = integrals

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

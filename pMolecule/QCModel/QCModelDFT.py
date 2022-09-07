"""pDynamo's in-built DFT QC model."""

from   pScientific.Arrays         import Array                                           , \
                                         StorageType
from  .DFTFunctionalModel         import DFTFunctionalModel
from  .DFTGridIntegrator          import DFTGridIntegrator
from  .FockConstruction           import FockConstruction_MakeFromFitIntegralsCoulomb    , \
                                         FockConstruction_MakeFromFitIntegralsNonCoulomb , \
                                         FockConstruction_MakeFromTEIsCoulomb            , \
                                         FockConstruction_MakeFromTEIsExchange
from  .GaussianBases              import GaussianBasisContainer                          , \
                                         GaussianBasisIntegralEvaluator                  , \
                                         GaussianBasisOperator
from  .LoewdinMultipoleEvaluator  import LoewdinMultipoleEvaluator
from  .MullikenMultipoleEvaluator import MullikenMultipoleEvaluator
from  .QCDefinitions              import ChargeModel                                     , \
                                         FockClosurePriority
from  .QCModelBase                import QCModelBase
from  .QCModelError               import QCModelError
from ..EnergyModel                import EnergyClosurePriority

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Fit operator.
_DefaultFitOperator = GaussianBasisOperator.Coulomb

# . Functional - defaults to regular HF.
_DefaultFunctional = "HF"

# . Maximum memory.
_DefaultMaximumMemory = 2.0 # GB.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelDFT ( QCModelBase ):
    """An DFT QC model."""

    _attributable = dict ( QCModelBase._attributable )
    _classLabel   = "DFT QC Model"
    _chargeModels = { ChargeModel.Loewdin  : LoewdinMultipoleEvaluator  , # . Note Loewdin multipoles are not rotationally invariant if Cartesian basis sets are used.
                      ChargeModel.Mulliken : MullikenMultipoleEvaluator } # . Older pDynamo versions corrected for this by orthogonalizing the basis but this introduced
                                                                          # . a host of other problems in the code and so has been removed. Use spherical harmonic basis
                                                                          # . sets instead!
    _summarizable = dict ( QCModelBase._summarizable )
    _attributable.update ( { "fitBasis"           : None                           , # . Can be None.
                             "fitOperator"        : _DefaultFitOperator            ,
                             "functional"         : _DefaultFunctional             ,
                             "functionalModel"    : None                           , # . Can be None.
                             "gridIntegrator"     : None                           , # . Need default if functional defined.
                             "integralEvaluator"  : GaussianBasisIntegralEvaluator ,
                             "multipoleEvaluator" : MullikenMultipoleEvaluator     ,
                             "maximumMemory"      : _DefaultMaximumMemory          ,
                             "orbitalBasis"       : "6-31g_st"                     } )
    _summarizable.update ( { "fitBasis"           : "Fit Basis"                    ,
                             "fitOperator"        : "Fit Operator"                 ,
                             "functional"         : "Functional"                   ,
                             "gridIntegrator"     : None                           ,
                             "maximumMemory"      : ( "Maximum Memory (GB)", "{:.3f}" ) ,
                             "orbitalBasis"       :   "Orbital Basis"              } )

    def _CheckOptions ( self ):
        """Check options."""
        DFTFunctionalModel.FindIDs ( self.functional ) # . Check functional option.

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelDFT, self ).BuildModel ( target, qcSelection = qcSelection )
        self.CheckMemory ( target )
        return state

    def CheckMemory ( self, target ):
        """A simple memory check."""
        # . No sparsity is assumed.
        # . Initialization.
        m = 0.0
        n = float ( len ( target.qcState.orbitalBases ) )
        # . Fit basis.
        if self.fitBasis is not None:
            f  = float ( len ( target.qcState.fitBases ) )
            m += 7.0 * ( f * n * ( n + 1.0 ) )   # . Electron-fit integrals with 1 Real64, 1 Integer32, 1 Integer16 ( 14 / 2 = 7 ).
            m += 4.0 * ( f + 1.0 ) * ( f + 2.0 ) # . Inverse fit matrix with Real64 ( 8 / 2 = 4 ).
        # . TEIs.
        if ( self.fitBasis is None ) or ( self.exchangeScaling != 0.0 ):
            p  = ( n * ( n + 1.0 ) ) / 2.0
            q  = ( p * ( p + 1.0 ) ) / 2.0
            m += 16.0 * q # . TEIs with 1 Real64, 4 Integer16.
        # . Grid.
        if self.functionalModel is not None: m += self.gridIntegrator.EstimateMemory ( self, target.qcState )
        # . Convert to GB and check.
        m /= 1.0e+09
        if m > self.maximumMemory:
            raise QCModelError ( "Estimated memory, {:.3f} GB, exceeds maximum memory, {:.3f} GB.".format ( m, self.maximumMemory ) )

    def EnergyClosureGradients ( self, target ):
        """Gradient energy closures."""
        # . One-electron terms.
        def a ( ): self.integralEvaluator.f2Cm1R1  ( target ) 
        def b ( ): self.integralEvaluator.f1KOf1R1 ( target ) 
        closures = [ ( EnergyClosurePriority.QCGradients, a, "QC Electron-Nuclear Gradients"    ) ,
                     ( EnergyClosurePriority.QCGradients, b, "QC Kinetic and Overlap Gradients" ) ]
        # . Fit basis.
        if self.fitBasis is not None:
            def c ( ): self.FitPreGradients ( target )
            def d ( ): self.integralEvaluator.f1Xg2R1 ( target, operator = self.fitOperator ) 
            def e ( ): self.integralEvaluator.f1Xf1R1 ( target, operator = self.fitOperator ) 
            closures.extend ( [ ( EnergyClosurePriority.QCPreGradients, c, "QC Fit Pregradients"       ) ,
                                ( EnergyClosurePriority.QCGradients   , d, "QC Electron-Fit Gradients" ) ,
                                ( EnergyClosurePriority.QCGradients   , e, "QC Fit-Fit Gradients"      ) ] )
            if self.fitOperator is not GaussianBasisOperator.Coulomb:
                def f ( ): self.integralEvaluator.f1Xf1R1 ( target, attributeX = "fitGradientVectorMc" ) # . Coulomb operator always.
                closures.append ( ( EnergyClosurePriority.QCGradients, f, "QC Coulomb Fit-Fit Integrals" ) )
        # . Functional.
        if self.functionalModel is not None:
            def g ( ): self.gridIntegrator.Gradients ( target )
            closures.append ( ( EnergyClosurePriority.QCGradients, g, "QC Grid Quadrature Gradients" ) )
        # . Two-electron integrals.
        if ( self.fitBasis is None ) or ( self.exchangeScaling != 0.0 ):
            def h ( ):
                self.integralEvaluator.f2Cf2R1 ( target                                      ,
                                                 doCoulomb       = ( self.fitBasis is None ) ,
                                                 exchangeScaling = self.exchangeScaling      )
            closures.append ( ( EnergyClosurePriority.QCGradients, h, "QC Two-Electron Gradients" ) )
        # . Weighted density.
        def i ( ): self.GetWeightedDensity ( target )
        closures.append ( ( EnergyClosurePriority.QCPreGradients, i, "QC Weighted Density" ) )
        return closures

    def EnergyClosures ( self, target ):
        """Energy closures."""
        # . Initialization.
        closures = super ( QCModelDFT, self ).EnergyClosures ( target )
        # . Nuclear and one-electron terms.
        def a ( ): self.integralEvaluator.m1Cm1ER1 ( target ) # . With derivatives if necessary.
        def b ( ): self.integralEvaluator.f2Cm1V   ( target )
        def c ( ): self.integralEvaluator.f1KOf1i  ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, a, "QC Nuclear"                       ) ,
                            ( EnergyClosurePriority.QCIntegrals, b, "QC Electron-Nuclear Integrals"    ) ,
                            ( EnergyClosurePriority.QCIntegrals, c, "QC Kinetic and Overlap Integrals" ) ] )
        # . Fit basis.
        if self.fitBasis is not None:
            def d ( ): self.integralEvaluator.f1Xg2i      ( target, operator = self.fitOperator )
            def e ( ): self.integralEvaluator.f1Xf1i_f1Oi ( target, operator = self.fitOperator )
            closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, d, "QC Electron-Fit Integrals" ) ,
                                ( EnergyClosurePriority.QCIntegrals, e, "QC Fit-Fit Integrals"      ) ] )
            if self.fitOperator is not GaussianBasisOperator.Coulomb:
                def f ( ): self.integralEvaluator.f1Xf1i_f1Oi ( target, attribute = "fitCoulombMatrix", withConstraints = False )
                closures.append ( ( EnergyClosurePriority.QCIntegrals, f, "QC Coulomb Fit-Fit Integrals" ) )
        # . Functional.
        if self.functionalModel is not None:
            def g ( ): self.gridIntegrator.BuildGrid ( target )
            closures.append ( ( EnergyClosurePriority.QCIntegrals, g, "QC Grid Quadrature Construction" ) )
        # . Two-electron integrals.
        if ( self.fitBasis is None ) or ( self.exchangeScaling != 0.0 ):
            def h ( ): self.integralEvaluator.f2Cf2i ( target ) 
            closures.append ( ( EnergyClosurePriority.QCIntegrals, h, "QC Two-Electron Integrals" ) )
        def i ( ): self.GetOrthogonalizer ( target )
        closures.append ( ( EnergyClosurePriority.QCOrthogonalizer, i, "QC Orthogonalizer" ) )
        # . Gradients.
        closures.extend ( self.EnergyClosureGradients ( target ) )
        # . Finish up.
        return closures

    def EnergyInitialize ( self, target ):
        """Energy initialization."""
        super ( QCModelDFT, self ).EnergyInitialize ( target )
        state   = getattr ( target, self.__class__._stateName )
        n       = len ( state.orbitalBases )
        scratch = target.scratch
        # . Fit arrays.
        if self.fitBasis is not None:
            f = len ( state.fitBases ) + 1 # . One less if no constraint.
            if scratch.Get ( "fitCoefficients", None ) is None:
                 scratch.fitCoefficients = Array.WithExtent ( f )
            if scratch.Get ( "fitGradientVectorB", None ) is None:
                scratch.fitGradientVectorB = Array.WithExtent ( f )
            scratch.fitGradientVectorB.Set ( 0.0 )
            if scratch.Get ( "fitGradientVectorM", None ) is None:
                scratch.fitGradientVectorM = Array.WithExtent ( f )
            scratch.fitGradientVectorM.Set ( 0.0 )
            if self.fitOperator is not GaussianBasisOperator.Coulomb:
                if scratch.Get ( "fitGradientVectorMc", None ) is None:
                    scratch.fitGradientVectorMc = Array.WithExtent ( f )
                scratch.fitGradientVectorMc.Set ( 0.0 )
                if scratch.Get ( "fitVectorD", None ) is None:
                    scratch.fitVectorD = Array.WithExtent ( f )
        # . Overlap.
        overlap = scratch.Get ( "overlapMatrix", None )
        if ( overlap is None ) or ( overlap.rows != n ):
            overlap = Array.WithExtent ( n, storageType = StorageType.Symmetric )
            scratch.overlapMatrix = overlap
        overlap.Set ( 0.0 )
        # . Weighted density.
        if scratch.doGradients:
            wDM = scratch.Get ( "weightedDensity", None )
            if wDM is None:
                wDM = Array.WithExtent ( n, storageType = StorageType.Symmetric )
                scratch.weightedDensity = wDM
            wDM.Set ( 0.0 )

    def FitPreGradients ( self, target ):
        """Fit pre-gradients."""
        scratch = target.scratch
        if self.fitOperator is GaussianBasisOperator.Coulomb:
            scratch.fitGradientVectorB.Add  ( scratch.fitCoefficients )
            scratch.fitGradientVectorM.Add  ( scratch.fitCoefficients )
        else:
            scratch.fitGradientVectorB.Add  ( scratch.fitVectorD                    )
            scratch.fitGradientVectorM.Add  ( scratch.fitVectorD     , scale =  2.0 )
            scratch.fitGradientVectorMc.Add ( scratch.fitCoefficients, scale = -1.0 )

    def FockClosures ( self, target ):
        """Fock closures."""
        # . OEI + finalization.
        def a ( ):
            return self.FockOne ( target )
        closures = [ ( FockClosurePriority.VeryLow , a ) ]
        # . Quadrature.
        if self.functionalModel is not None:
            def b ( ):
                return self.FockQuadrature ( target )
            closures.append ( ( FockClosurePriority.Medium, b ) )
        # . Coulomb via TEIs.
        if self.fitBasis is None:
            # . With exchange.
            if self.exchangeScaling != 0.0:
                def c ( ):
                    return self.FockTwo ( target )
                closures.append ( ( FockClosurePriority.VeryHigh, c ) )
            # . Without exchange.
            else:
                def d ( ):
                    return self.FockTwoCoulomb ( target )
                closures.append ( ( FockClosurePriority.VeryHigh, d ) )
        # . Coulomb via fit basis.
        else:
            # . With exchange.
            if self.exchangeScaling != 0.0:
                def e ( ):
                    return self.FockTwoExchange ( target )
                closures.append ( ( FockClosurePriority.VeryHigh, e ) )
            # . Without exchange.
            else:            
                def f ( ):
                    return self.FockInitialize ( target )
                closures.append ( ( FockClosurePriority.VeryHigh, f ) )
            # . Coulomb.
            def g ( ):
                return self.FockFit ( target )
            closures.append ( ( FockClosurePriority.Medium, g ) )
        return closures

    # . Beware - dTotal is temporarily modified here!
    def FockFit ( self, target ):
        """The fit Coulomb contribution to the Fock matrices."""
        scratch = target.scratch
        if self.fitOperator is GaussianBasisOperator.Coulomb:
            eFit = FockConstruction_MakeFromFitIntegralsCoulomb ( scratch.onePDMP.density      ,
                                                                  scratch.fitIntegrals         ,
                                                                  scratch.fitMatrix            ,
                                                                  scratch.onePDMP.totalCharge  ,
                                                                  scratch.fitCoefficients      ,
                                                                  scratch.onePDMP.fock         )
        else:
            eFit = FockConstruction_MakeFromFitIntegralsNonCoulomb ( scratch.onePDMP.density      ,
                                                                     scratch.fitIntegrals         ,
                                                                     scratch.fitMatrix            ,
                                                                     scratch.fitCoulombMatrix     ,
                                                                     scratch.onePDMP.totalCharge  ,
                                                                     scratch.fitCoefficients      ,
                                                                     scratch.fitVectorD           ,
                                                                     scratch.onePDMP.fock         )
        scratch.qcEnergyReport["Fit Coulomb Energy"       ]  = eFit
        scratch.qcEnergyReport["QC Electronic Accumulator"] += eFit
        return eFit

    # . This should always be used when the fit integral/TEI Fock methods are sorted out.
    def FockInitialize ( self, target ):
        """Initialization."""
        scratch = target.scratch
        scratch.onePDMP.fock.Set ( 0.0 )
        if hasattr ( scratch, "onePDMQ" ): scratch.onePDMQ.fock.Set ( 0.0 )
        scratch.qcEnergyReport["QC Electronic Accumulator"] = 0.0
        return 0.0

    def FockQuadrature ( self, target ):
        """The quadrature contribution to the Fock matrices."""
        scratch = target.scratch
        eQuad   = self.gridIntegrator.Fock ( target )
        scratch.qcEnergyReport["Quadrature Energy"        ]  = eQuad
        scratch.qcEnergyReport["QC Electronic Accumulator"] += eQuad
        return eQuad

    # . Accumulator and Fock matrices initialized in the next two methods.
    def FockTwoCoulomb ( self, target ):
        """The two-electron Coulomb contribution to the Fock matrices."""
        scratch = target.scratch
        dTotal  = scratch.onePDMP.density
        fTotal  = scratch.onePDMP.fock
        eTE     = FockConstruction_MakeFromTEIsCoulomb ( dTotal                       ,
                                                         scratch.twoElectronIntegrals ,
                                                         fTotal                       )
        scratch.qcEnergyReport["QC Electronic Accumulator"  ] = eTE
        scratch.qcEnergyReport["Two-Electron Coulomb Energy"] = eTE
        if hasattr ( scratch, "onePDMQ" ): scratch.onePDMQ.fock.Set ( 0.0 )
        return eTE

    def FockTwoExchange ( self, target ):
        """The two-electron exchange contribution to the Fock matrices."""
        scratch = target.scratch
        doSpin  = hasattr ( scratch, "onePDMQ" )
        dTotal  = scratch.onePDMP.density
        fTotal  = scratch.onePDMP.fock
        if doSpin:
            dSpin = scratch.onePDMQ.density
            fSpin = scratch.onePDMQ.fock
        else:
            dSpin = None
            fSpin = None
        eTE  = FockConstruction_MakeFromTEIsExchange ( dTotal                       ,
                                                       dSpin                        ,
                                                       scratch.twoElectronIntegrals ,
                                                       self.exchangeScaling         ,
                                                       fTotal                       ,
                                                       fSpin                        )
        scratch.qcEnergyReport["QC Electronic Accumulator"   ] = eTE
        scratch.qcEnergyReport["Two-Electron Exchange Energy"] = eTE
        return eTE

    def GetParameters ( self, target ):
        """Get the parameters for the model."""
        state = getattr ( target, self.__class__._stateName )
        if self.fitBasis is None:
            state.fitBases = None
        else:
            state.fitBases = GaussianBasisContainer.FromParameterDirectory ( self.fitBasis       ,
                                                                             state.atomicNumbers )
        state.orbitalBases   = GaussianBasisContainer.FromParameterDirectory ( self.orbitalBasis, state.atomicNumbers )
        state.nuclearCharges = state.orbitalBases.nuclearCharges
        self.functionalModel = DFTFunctionalModel.FromOptions ( self.functional, isSpinRestricted = target.electronicState.isSpinRestricted )
        if self.functionalModel is not None:
            self.exchangeScaling = self.functionalModel.exchangeScaling
            if self.gridIntegrator is None: self.gridIntegrator = DFTGridIntegrator ( )

    def GetWeightedDensity ( self, target ): 
        """Get the weighted density."""
        scratch = target.scratch
        if scratch.doGradients:
            wDM = scratch.weightedDensity
            for label in ( "orbitalsP", "orbitalsQ" ):
                orbitals = scratch.Get ( label, None )
                if orbitals is not None: orbitals.MakeWeightedDensity ( wDM )

    def MakeLabel ( self, target ):
        """Make a label."""
        if self.functionalModel is None: item = "HF"
        else:                            item = "DFT"
        if target.electronicState.isSpinRestricted: item = "R{:s}".format ( item )
        else:                                       item = "U{:s}".format ( item )
        items = [ item, self.orbitalBasis.upper ( ) ]
        if self.fitBasis is not None: items.append ( self.fitBasis.upper ( ) )
        return "/".join ( items )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

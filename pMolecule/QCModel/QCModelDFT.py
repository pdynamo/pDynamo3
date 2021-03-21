"""pDynamo's in-built DFT QC model."""

from   pScientific.Arrays             import Array                                 , \
                                             StorageType
from   pScientific.LinearAlgebra      import MatrixPower
from  .DFTFunctionalModel             import DFTFunctionalModel
from  .DFTGridIntegrator              import DFTGridIntegrator
from  .FockConstruction               import FockConstruction_MakeFromFitIntegrals , \
                                             FockConstruction_MakeFromTEIsCoulomb  , \
                                             FockConstruction_MakeFromTEIsExchange
from  .GaussianBasisContainer         import GaussianBasisContainer
from  .GaussianBasisIntegralEvaluator import GaussianBasisIntegralEvaluator
from  .LoewdinMultipoleEvaluator      import LoewdinMultipoleEvaluator
from  .MullikenMultipoleEvaluator     import MullikenMultipoleEvaluator
from  .QCDefinitions                  import ChargeModel                           , \
                                             FockClosurePriority
from  .QCModelBase                    import QCModelBase
from  .QCModelError                   import QCModelError
from ..EnergyModel                    import EnergyClosurePriority

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
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
    _chargeModels = { ChargeModel.Loewdin  : LoewdinMultipoleEvaluator  ,
                      ChargeModel.Mulliken : MullikenMultipoleEvaluator }
    _summarizable = dict ( QCModelBase._summarizable )
    _attributable.update ( { "fitBasis"           : None                           , # . Can be None.
                             "functional"         : _DefaultFunctional             ,
                             "functionalModel"    : None                           , # . Can be None.
                             "gridIntegrator"     : None                           , # . Need default if functional defined.
                             "integralEvaluator"  : GaussianBasisIntegralEvaluator ,
                             "multipoleEvaluator" : MullikenMultipoleEvaluator     ,
                             "maximumMemory"      : _DefaultMaximumMemory          ,
                             "orbitalBasis"       : "631gs"                        } )
    _summarizable.update ( { "fitBasis"           :   "Fit Basis"                       ,
                             "functional"         :   "Functional"                      ,
                             "gridIntegrator"     : None                                ,
                             "maximumMemory"      : ( "Maximum Memory (GB)", "{:.3f}" ) ,
                             "orbitalBasis"       :   "Orbital Basis"                   } )
                             
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
        if ( self.functional is not None ) and ( self.gridIntegrator is not None ): m += self.gridIntegrator.EstimateMemory ( self, target.qcState )
        # . Convert to GB and check.
        m /= 1.0e+09
        if m > self.maximumMemory:
            raise QCModelError ( "Estimated memory, {:.3f} GB, exceeds maximum memory, {:.3f} GB.".format ( m, self.maximumMemory ) )

    def EnergyClosureGradients ( self, target ):
        """Gradient energy closures."""
        # . One-electron terms.
        def a ( ): self.integralEvaluator.ElectronNuclearGradients ( target )
        def b ( ): self.integralEvaluator.KineticOverlapGradients  ( target )
        closures = [ ( EnergyClosurePriority.QCGradients, a, "QC Electron-Nuclear Gradients"    ) ,
                     ( EnergyClosurePriority.QCGradients, b, "QC Kinetic and Overlap Gradients" ) ]
        # . Fit basis.
        if self.fitBasis is not None:
            def c ( ): self.integralEvaluator.ElectronFitGradients ( target )
            def d ( ): self.integralEvaluator.FitFitGradients      ( target )
            closures.extend ( [ ( EnergyClosurePriority.QCGradients, c, "QC Electron-Fit Gradients" ) ,
                                ( EnergyClosurePriority.QCGradients, d, "QC Fit-Fit Gradients"      ) ] )
        # . Functional.
        if self.functionalModel is not None:
            def e ( ): self.gridIntegrator.Gradients ( target )
            closures.append ( ( EnergyClosurePriority.QCGradients, e, "QC Grid Quadrature Gradients" ) )
        # . Two-electron integrals.
        if ( self.fitBasis is None ) or ( self.exchangeScaling != 0.0 ):
            def f ( ):
                self.integralEvaluator.TwoElectronGradients ( target                                      ,
                                                              doCoulomb       = ( self.fitBasis is None ) ,
                                                              exchangeScaling = self.exchangeScaling      )
            closures.append ( ( EnergyClosurePriority.QCGradients, f, "QC Two-Electron Gradients" ) )
        # . Weighted density.
        def g ( ): self.GetWeightedDensity ( target )
        closures.append ( ( EnergyClosurePriority.QCPreGradients, g, "QC Weighted Density" ) )
        return closures

    def EnergyClosures ( self, target ):
        """Energy closures."""
        # . Initialization.
        closures = super ( QCModelDFT, self ).EnergyClosures ( target )
        # . Nuclear and one-electron terms.
        def a ( ): self.integralEvaluator.CoreCoreEnergy           ( target ) # . With derivatives if necessary.
        def b ( ): self.integralEvaluator.ElectronNuclearIntegrals ( target )
        def c ( ): self.integralEvaluator.KineticOverlapIntegrals  ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, a, "QC Nuclear"                       ) ,
                            ( EnergyClosurePriority.QCIntegrals, b, "QC Electron-Nuclear Integrals"    ) ,
                            ( EnergyClosurePriority.QCIntegrals, c, "QC Kinetic and Overlap Integrals" ) ] )
        # . Fit basis.
        if self.fitBasis is not None:
            def d ( ): self.integralEvaluator.ElectronFitIntegrals ( target )
            def e ( ): self.integralEvaluator.FitFitIntegrals      ( target )
            closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, d, "QC Electron-Fit Integrals" ) ,
                                ( EnergyClosurePriority.QCIntegrals, e, "QC Fit-Fit Integrals"      ) ] )
        # . Functional.
        if self.functionalModel is not None:
            def f ( ): self.gridIntegrator.BuildGrid ( target )
            closures.append ( ( EnergyClosurePriority.QCIntegrals, f, "QC Grid Quadrature Construction" ) )
        # . Two-electron integrals.
        if ( self.fitBasis is None ) or ( self.exchangeScaling != 0.0 ):
            def g ( ): self.integralEvaluator.TwoElectronIntegrals ( target )
            closures.append ( ( EnergyClosurePriority.QCIntegrals, g, "QC Two-Electron Integrals" ) )
        def h ( ): self.GetOrthogonalizer ( target )
        closures.append ( ( EnergyClosurePriority.QCPreEnergy, h, "QC Orthogonalizer" ) )
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
        # . Fit potential.
        if ( self.fitBasis is not None ) and ( scratch.Get ( "fitPotential", None ) is None ):
            scratch.fitPotential = Array.WithExtent ( len ( state.fitBases ) + 1 ) # . One less if no constraint.
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
        # . Coloumb via fit basis.
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
        eFit    = FockConstruction_MakeFromFitIntegrals ( scratch.onePDMP.density      ,
                                                          scratch.electronFitIntegrals ,
                                                          scratch.fitMatrix            ,
                                                          scratch.onePDMP.totalCharge  ,
                                                          scratch.fitPotential         ,
                                                          scratch.onePDMP.fock         )
        #print ( "\nFit Energy = {:.5f}".format ( eFit ) )
        #from pScientific.Arrays import ArrayPrint
        #ArrayPrint ( scratch.fitPotential, itemFormat = "{:.5f}", itemsPerRow = 10, title = "Fit Potential" )
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
        # . The fit and orbital bases are automatically recognized as Coulomb and Orbital from their definitions.
        state = getattr ( target, self.__class__._stateName )
        if self.fitBasis is None:
            state.fitBases = None
        else:
            state.fitBases = GaussianBasisContainer.FromParameterDirectory ( self.fitBasis, state.atomicNumbers )
            state.fitBases.basisRepresentation = self.integralEvaluator.basisRepresentation
        state.orbitalBases   = GaussianBasisContainer.FromParameterDirectory ( self.orbitalBasis, state.atomicNumbers )
        state.nuclearCharges = state.orbitalBases.nuclearCharges
        # . Functional model.
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

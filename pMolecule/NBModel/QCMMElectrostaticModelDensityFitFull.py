"""A full density-fitted QC/MM electrostatic model."""

from   pCore                                          import Clone
from   pMolecule.QCModel                              import FockConstruction_MakeCoefficientsFromFitIntegrals , \
                                                             FockConstruction_MakeFockFromFitIntegrals
from   pMolecule.QCModel.GaussianBases                import GaussianBasisContainer                            , \
                                                             GaussianBasisIntegralEvaluator                    , \
                                                             GaussianBasisOperator                             , \
                                                             GaussianBasisQCMMEvaluator
from   pScientific                                    import Units
from   pScientific.Arrays                             import Array                                             , \
                                                             StorageType
from   pScientific.Geometry3                          import Coordinates3
from   pScientific.LinearAlgebra                      import LinearEquations
from  .NBModelError                                   import NBModelError
from  .QCMMElectrostaticModelDensityFullGaussianBasis import QCMMElectrostaticModelDensityFullGaussianBasis
from ..EnergyModel                                    import EnergyClosurePriority

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Fit operator.
_DefaultFitOperator = GaussianBasisOperator.Coulomb

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelDensityFitFull ( QCMMElectrostaticModelDensityFullGaussianBasis ):
    """A full density-fitted QC/MM electrostatic model."""

    _attributable = dict ( QCMMElectrostaticModelDensityFullGaussianBasis._attributable )
    _summarizable = dict ( QCMMElectrostaticModelDensityFullGaussianBasis._summarizable )
    _classLabel   = "Gaussian Basis Full Density-Fitted QC/MM Electrostatic Model"
    _attributable.update ( { "fitBasis"     : None                           ,
                             "fitEvaluator" : GaussianBasisIntegralEvaluator ,
                             "fitOperator"  : None                           } )
    _summarizable.update ( { "fitBasis"     : "Fit Basis"                    ,
                             "fitOperator"  : "Fit Operator"                 } )

    def BuildModel ( self, target ):
        """Build the model."""
        state        = super ( QCMMElectrostaticModelDensityFullGaussianBasis, self ).BuildModel ( target )
        fitBasis0    = getattr ( target.qcModel, "fitBasis"   , None )
        fitOperator0 = getattr ( target.qcModel, "fitOperator", None )
        if self.fitBasis is None:
            if fitBasis0 is None: raise NBModelError ( "Fit basis not specified." )
            else:                 self.fitBasis = fitBasis0 # . Use the same as the QC model.
        if self.fitOperator is None:
            if fitOperator0 is None: raise NBModelError ( "Fit operator not specified." )
            else:                    self.fitOperator = fitOperator0
        if ( self.fitBasis == fitBasis0 ) and ( self.fitOperator is fitOperator0 ):
            state.fitBases    = target.qcState.fitBases
            state.isDuplicate = True
        else:
            state.fitBases    = GaussianBasisContainer.FromParameterDirectory ( self.fitBasis, target.qcState.atomicNumbers )
            state.isDuplicate = False
        return state

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        # . Initialization.
        closures = super ( QCMMElectrostaticModelDensityFitFull, self ).EnergyClosures ( target )
        state    = getattr ( target, self.__class__._stateName )
        fitBases = state.fitBases
        # . No need to do fit.
        if state.isDuplicate:
            # . Pre-gradients.
            def a ( ):
                scratch = target.scratch
                scratch.fitGradientVectorB.Add ( scratch.qcmmFitVectorW ) # . Will need to check these - need a factor (1/2, 2, ...)? AND not isDuplicate too!
                scratch.fitGradientVectorM.Add ( scratch.qcmmFitVectorW , scale = 2.0 )
            closures.append ( ( EnergyClosurePriority.QCPreGradients, a, "QC/MM Fit Pregradients" ) )
        # . Need to do fit.
        else:
            # . Integrals.
            def b ( ):
                self.fitEvaluator.f1Xg2i      ( target                             ,
                                                attribute  = "qcmmFitIntegrals"    ,
                                                fitBases   = fitBases              ,
                                                operator   = self.fitOperator      ,
                                                reportTag  = "QC/MM Fit"           )
            def c ( ):
                self.fitEvaluator.f1Xf1i_f1Oi ( target                             ,
                                                attribute  = "qcmmFitMatrix"       ,
                                                fitBases   = fitBases              ,
                                                operator   = self.fitOperator      )
            # . Gradients.
            def d ( ):
                self.fitEvaluator.f1Xg2R1     ( target                             ,
                                                attribute  = "qcmmFitVectorW"      ,
                                                fitBases   = fitBases              ,
                                                operator   = self.fitOperator      )
            def e ( ):
# . Fudge.
                target.scratch.qcmmFitCoefficients.Scale ( 2.0 )
                self.fitEvaluator.f1Xf1R1 ( target                             ,
                                            attributeA = "qcmmFitCoefficients" ,
                                            attributeX = "qcmmFitVectorW"      ,
                                            fitBases   = fitBases              ,
                                            operator   = self.fitOperator      )
                target.scratch.qcmmFitCoefficients.Scale ( 0.5 )
            closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, b, "QC/MM Electron/Fit Integrals" ) ,
                                ( EnergyClosurePriority.QCIntegrals, c, "QC/MM Fit/Fit Integrals"      ) ,
                                ( EnergyClosurePriority.QCGradients, d, "QC/MM Electron/Fit Gradients" ) ,
                                ( EnergyClosurePriority.QCGradients, e, "QC/MM Fit/Fit Gradients"      ) ] )
        # . Fock matrix.
        def f ( ): self.MakeFock ( target )
        closures.append ( ( EnergyClosurePriority.QCPreEnergy, f, "QC/MM Fit Fock" ) )
        # . Finish up.
        return closures

    def EnergyInitialize ( self, target ):
        """Energy initialization"""
        # . Initialization.
        super ( QCMMElectrostaticModelDensityFitFull, self ).EnergyInitialize ( target )
        scratch = target.scratch
        state   = getattr ( target, self.__class__._stateName )
        n       = len ( state.fitBases ) + 1 # . One less if no constraint.
        # . B vector.
        if scratch.Get ( "qcmmFitVectorB", None ) is None:
             scratch.qcmmFitVectorB = Array.WithExtent ( n )
        # . Coefficients.
        if not state.isDuplicate:
            if scratch.Get ( "qcmmFitCoefficients", None ) is None:
                 scratch.qcmmFitCoefficients = Array.WithExtent ( n )
        # . Fock matrix.
        if scratch.Get ( "qcmmFitFock", None ) is None:
            o = len ( target.qcState.orbitalBases )
            scratch.qcmmFitFock = Array.WithExtent ( o, storageType = StorageType.Symmetric )
        # . Potentials.
        potentials = target.scratch.Get ( "qcmmFitPotentials", None )
        if potentials is None:
            potentials = Array.WithExtent ( n )
            target.scratch.qcmmFitPotentials = potentials
        potentials.Set ( 0.0 )
        # . W vector.
        if scratch.Get ( "qcmmFitVectorW", None ) is None:
             scratch.qcmmFitVectorW = Array.WithExtent ( n )

    def Fock ( self, target ):
        """Energy and Fock matrix contributions."""
        scratch             = target.scratch
        state               = getattr ( target, self.__class__._stateName )
        qcmmFitCoefficients = scratch.Get ( "qcmmFitCoefficients",  scratch.Get ( "fitCoefficients", None ) )
        qcmmFitIntegrals    = scratch.Get ( "qcmmFitIntegrals"   ,  scratch.Get ( "fitIntegrals"   , None ) )
        qcmmFitMatrix       = scratch.Get ( "qcmmFitMatrix"      ,  scratch.Get ( "fitMatrix"      , None ) )
        # . Coefficients (unless already done by QC model).
        if not state.isDuplicate:
            FockConstruction_MakeCoefficientsFromFitIntegrals ( scratch.onePDMP.density     ,
                                                                qcmmFitIntegrals            ,
                                                                qcmmFitMatrix               ,
                                                                scratch.onePDMP.totalCharge ,
                                                                qcmmFitCoefficients         ,
                                                                scratch.qcmmFitVectorB      )
        # . Energy.
        eQCMM = qcmmFitCoefficients.Dot ( scratch.qcmmFitPotentials )
        # . Fock.
        scratch.onePDMP.fock.Add ( scratch.qcmmFitFock )
        # . Finish up.
        scratch.energyTerms["QC/MM Electrostatic"] = ( eQCMM * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        return eQCMM

    def MakeFock ( self, target ):
        """Make the Fock matrix - needs only to be done once."""
        scratch          = target.scratch
        qcmmFitIntegrals = scratch.Get ( "qcmmFitIntegrals", scratch.Get ( "fitIntegrals", None ) )
        qcmmFitMatrix    = scratch.Get ( "qcmmFitMatrix"   , scratch.Get ( "fitMatrix"   , None ) )
        # . Solve for the fit W-vector. */
        LinearEquations ( qcmmFitMatrix, scratch.qcmmFitPotentials, solution = scratch.qcmmFitVectorW )
        # . Fit Fock matrix.
        scratch.qcmmFitFock.Set ( 0.0 )
        FockConstruction_MakeFockFromFitIntegrals ( qcmmFitIntegrals       ,
                                                    scratch.qcmmFitVectorW ,
                                                    scratch.qcmmFitFock    )
        #from pScientific.Arrays import ArrayPrint, ArrayPrint2D
        #ArrayPrint   ( scratch.qcmmFitPotentials, itemFormat = "{:.3f}", title = "QC/MM Potentials" )
        #ArrayPrint   ( scratch.qcmmFitVectorW   , itemFormat = "{:.3f}", title = "QC/MM Vector W"   )
        #ArrayPrint2D ( scratch.qcmmFitFock      , itemFormat = "{:.3f}", title = "QC/MM Fock"       )
        #ArrayPrint2D ( qcmmFitMatrix            , itemFormat = "{:.3f}", title = "QC/MM Matrix M"   )
        #. Scale for gradients.
        #scratch.qcmmFitVectorW.Scale ( 2.0 )

    def QCBPGradients ( self, target ):
        """Calculate the QC/BP electrostatic gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            fitCoefficients = scratch.Get ( "qcmmFitCoefficients", scratch.Get ( "fitCoefficients", None ) )
            state           = getattr ( target, self.__class__._stateName )
            bpCharges       = getattr ( state, "bpCharges", None )
            if bpCharges is not None:
                self.evaluator.f1Cm1R1 ( state.fitBases               ,
                                         bpCharges                    ,
                                         scratch.qcCoordinates3QCMMAU ,
                                         scratch.bpCoordinates3AU     ,
                                         None                         ,
                                         fitCoefficients              ,
                                         scratch.qcGradients3QCMMAU   ,
                                         scratch.bpGradients3AU       )

    def QCBPPotentials ( self, target ):
        """Calculate the QC/BP electrostatic potentials."""
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if bpCharges is not None:
            scratch = target.scratch
            state   = getattr ( target, self.__class__._stateName )
            self.evaluator.f1Cm1V ( state.fitBases               , 
                                    bpCharges                    ,
                                    scratch.qcCoordinates3QCMMAU ,
                                    scratch.bpCoordinates3AU     ,
                                    None                         ,
                                    scratch.qcmmFitPotentials    )

    def QCMMGradients ( self, target ):
        """Calculate the QC/MM electrostatic gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            fitCoefficients = scratch.Get ( "qcmmFitCoefficients", scratch.Get ( "fitCoefficients", None ) )
            state           = getattr ( target, self.__class__._stateName )
            self.evaluator.f1Cm1R1 ( state.fitBases               ,
                                     target.mmState.charges       ,
                                     scratch.qcCoordinates3QCMMAU ,
                                     scratch.mmCoordinates3AU     ,
                                     target.mmState.pureMMAtoms   ,
                                     fitCoefficients              ,
                                     scratch.qcGradients3QCMMAU   ,
                                     scratch.mmGradients3AU       )

    def QCMMPotentials ( self, target ):
        """Calculate the QC/MM electrostatic potentials."""
        scratch = target.scratch
        state   = getattr ( target, self.__class__._stateName )
        self.evaluator.f1Cm1V ( state.fitBases               ,
                                target.mmState.charges       ,
                                scratch.qcCoordinates3QCMMAU ,
                                scratch.mmCoordinates3AU     ,
                                target.mmState.pureMMAtoms   ,
                                scratch.qcmmFitPotentials    )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

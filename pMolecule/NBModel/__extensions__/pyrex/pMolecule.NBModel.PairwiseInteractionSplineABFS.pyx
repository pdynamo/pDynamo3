"""Spline ABFS pairwise interactions."""

import math

from  enum               import Enum
from  pCore              import Clone          , \
                                logFile        , \
                                LogFileActive
from  pScientific.Arrays import Array
from .ABFSIntegrator     import ABFSIntegrator
from .NBModelError       import NBModelError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Spline models.
class SplineModel ( Enum ):
    """Spline models."""
    Delta_Delta       =  0
    Delta_Gaussian    = 10
    Gaussian_Gaussian = 20
    LennardJonesA     = 30
    LennardJonesB     = 40

# . Defaults.
_DefaultSplineModel = SplineModel.Delta_Delta

# . The minimum value of R allowed in integration functions (to prevent division by zero).
_MinimumR = 1.0e-2

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionSplineABFS ( PairwiseInteraction ):
    """Spline ABFS pairwise interactions."""

    def __copy__ ( self ):
        """Copying."""
        return self.__class__ ( **self.__getstate__ )

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairwiseInteractionSpline_DeassignSplines (  self.cObject )
            PairwiseInteractionSpline_Deallocate      ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        return { "dampingCutOff"      : self.dampingCutOff      ,
                 "electrostaticModel" : self.electrostaticModel ,
                 "innerCutOff"        : self.innerCutOff        ,
                 "outerCutOff"        : self.outerCutOff        ,
                 "pointDensity"       : self.pointDensity       ,
                 "width1"             : self.width1             ,
                 "width2"             : self.width2             }

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = PairwiseInteractionSpline_Allocate ( &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise NBModelError ( "Error allocating pairwise interaction." )

    def _CheckOptions ( self ):
        """Check the options."""
        if self.electrostaticModel is None: self.electrostaticModel = _DefaultSplineModel
        elif not isinstance ( self.electrostaticModel, SplineModel ):
            raise NBModelError ( "Invalid electrostatic model type: {:s}.".format ( repr ( self.electrostaticModel ) ) )
        if ( self.dampingCutOff < 0.0                ) or \
           ( self.innerCutOff   < self.dampingCutOff ) or \
           ( self.outerCutOff   < self.innerCutOff   ):
            raise NBModelError ( "Invalid cutOff values: damping - {:.3f}; inner - {:.3f}; outer - {:.3f}.".format ( self.dampingCutOff , \
                                                                                                                     self.innerCutOff   , \
                                                                                                                     self.outerCutOff   ) )
        self.integrator = ABFSIntegrator.WithOptions ( dampingCutOff = self.dampingCutOff ,
                                                       innerCutOff   = self.innerCutOff   ,
                                                       outerCutOff   = self.outerCutOff   ,
                                                       pointDensity  = self.pointDensity  )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject             = NULL
        self.isOwner             = False
        self.dampingCutOff       = 0.5
        self.electrostaticModel  = _DefaultSplineModel
        self.electrostaticSpline = None
        self.innerCutOff         = 8.0
        self.integrator          = None
        self.lennardJonesASpline = None
        self.lennardJonesBSpline = None
        self.outerCutOff         = 12.0
        self.pointDensity        = 50
        self.width1              = 0.0
        self.width2              = 0.0

    def _MakeSpline ( self, model ):
        """Make a spline."""
        ( F, G ) = self.IntegrationFunctions ( model )
        x = Clone ( self.integrator.x )
        y = self.integrator.Integrate ( F, G )
        x.Square ( ) # . Squaring X simplifies derivative calculation.
        spline = CubicSpline.FromArrays ( x, y, lowerDerivative = 1, lowerValue = 0.0, upperDerivative = 1, upperValue = 0.0 )
        return spline

    # . Gaussian models:
    # . f(r) = (a/pi)**(3/2) * exp ( - a * r**2 ).
    # . <r> = 2 / sqrt ( a * pi ) so a = 4 / ( pi * <r>**2 ).
    # . p = "sqrt ( a )"                        Delta_Gaussian
    # . p = "sqrt ( a * b / ( a + b ) )"        Gaussian_Gaussian
    # . Small x expansions up to 10th order. At switch over point differences are of the order of 10^-16 (F) and 10^-14/15 (G).
    def IntegrationFunctions ( self, model ):
        """Return the integration functions."""
        factor = 2.0 / math.sqrt ( math.pi )
        if model is SplineModel.Delta_Delta:            # . 1/r.
            def F ( r ):
                r = max ( r, _MinimumR )
                return (   1.0 / r     )
            def G ( r ):
                r = max ( r, _MinimumR )
                return ( - 1.0 / r**2  )
        elif model is SplineModel.LennardJonesA:        # . 1/r^12.
            def F ( r ):
                r = max ( r, _MinimumR )
                return (   1.0 / r**12 )
            def G ( r ):
                r = max ( r, _MinimumR )
                return ( -12.0 / r**13 )
        elif model is SplineModel.LennardJonesB:        # . -1/r^6.
            def F ( r ):
                r = max ( r, _MinimumR )
                return (  -1.0 / r**6  )
            def G ( r ):
                r = max ( r, _MinimumR )
                return (   6.0 / r**7  )
        else:                                           # . Gaussian models (Delta/Gaussian, Gaussian/Gaussian) - only differ in widths.
            if   model is SplineModel.Delta_Gaussian   : p = factor / self.width2
            elif model is SplineModel.Gaussian_Gaussian: p = factor / math.sqrt ( self.width1**2 + self.width2**2 )
            else: raise NBModelError ( "Invalid spline model: {:s}.".format ( model ) )
            def FGaussian ( p ):
                def F ( r ):
                    x  = p * r
                    x2 = x * x
                    if x > _MinimumR: f = math.erf ( x ) / x
                    else:             f = ( 2.0 - x2 * ( 2.0 / 3.0 - x2 * ( 0.2 - x2 * ( 1.0 / 21.0 - x2 * ( 1.0 / 108.0 - x2 / 660.0 ) ) ) ) ) / math.sqrt ( math.pi )
                    return ( f * p )
                return F
            def GGaussian ( p ):
                def G ( r ):
                    sPi = math.sqrt ( math.pi )
                    x   = p * r
                    x2  = x * x
                    if x > _MinimumR: g = 2.0 * math.exp ( - x2 ) / ( sPi * x ) - math.erf ( x ) / x2
                    else:             g = - x * ( 4.0 / 3.0 - x2 * ( 0.8 - x2 * ( 2.0 / 7.0 - x2 * ( 2.0 / 27.0 - x2 / 66.0 ) ) ) ) / sPi
                    return ( g * p * p )
                return G
            F = FGaussian ( p )
            G = GGaussian ( p )
        return ( F, G )

    def Interactions ( self, RealArray1D r not None ):
        """Calculate the interaction components for an array of distances."""
        cdef RealArray1D electrostatic, lennardJonesA, lennardJonesB
        extent        = len ( r )
        electrostatic = Array.WithExtent ( extent )
        lennardJonesA = Array.WithExtent ( extent )
        lennardJonesB = Array.WithExtent ( extent )
        PairwiseInteractionSpline_Interactions ( self.cObject          ,       
                                                 r.cObject             ,
                                                 electrostatic.cObject ,
                                                 lennardJonesA.cObject ,
                                                 lennardJonesB.cObject )
        return { "Electrostatic" : electrostatic ,
                 "LennardJonesA" : lennardJonesA ,
                 "LennardJonesB" : lennardJonesB }

    def MakeElectrostaticSpline ( self ):
        """Make the electrostatic spline."""
        cdef CubicSpline spline
        spline = self._MakeSpline ( self.electrostaticModel )
        PairwiseInteractionSpline_AssignElectrostaticSpline ( self.cObject, spline.cObject, self.outerCutOff )
        return spline

    def MakeLennardJonesASpline ( self ):
        """Make the SplineModel.LennardJonesA spline."""
        cdef CubicSpline spline
        spline = self._MakeSpline ( SplineModel.LennardJonesA )
        PairwiseInteractionSpline_AssignLennardJonesASpline ( self.cObject, spline.cObject, self.outerCutOff )
        return spline

    def MakeLennardJonesBSpline ( self ):
        """Make the SplineModel.LennardJonesB spline."""
        cdef CubicSpline spline
        spline = self._MakeSpline ( SplineModel.LennardJonesB )
        PairwiseInteractionSpline_AssignLennardJonesBSpline ( self.cObject, spline.cObject, self.outerCutOff )
        return spline

    def MMMMEnergy ( self, RealArray1D          chargesI               ,
                           RealArray1D          chargesJ               ,
                           IntegerArray1D       ljTypesI               ,
                           IntegerArray1D       ljTypesJ               ,
                           LJParameterContainer ljParameters           ,
                                                electrostaticScale     ,
                                                lennardJonesScale      ,
                           Coordinates3         coordinates3I not None ,
                           Coordinates3         coordinates3J not None ,
                           PairList             pairList      not None ,
                           Coordinates3         gradients3I            ,
                           Coordinates3         gradients3J            ):
        """MM/MM energy."""
        cdef CIntegerArray1D       *cLJTypesI      = NULL 
        cdef CIntegerArray1D       *cLJTypesJ      = NULL 
        cdef CRealArray2D          *cGradients3I   = NULL 
        cdef CRealArray2D          *cGradients3J   = NULL 
        cdef CLJParameterContainer *cLJParameters  = NULL
        cdef CReal                  eElectrostatic
        cdef CReal                  eLennardJones
        cdef CRealArray1D          *cChargesI      = NULL
        cdef CRealArray1D          *cChargesJ      = NULL
        cdef CStatus                cStatus        = CStatus_OK
        if chargesI     is not None: cChargesI     = chargesI.cObject
        if chargesJ     is not None: cChargesJ     = chargesJ.cObject
        if gradients3I  is not None: cGradients3I  = gradients3I.cObject
        if gradients3J  is not None: cGradients3J  = gradients3J.cObject
        if ljParameters is not None: cLJParameters = ljParameters.cObject
        if ljTypesI     is not None: cLJTypesI     = ljTypesI.cObject
        if ljTypesJ     is not None: cLJTypesJ     = ljTypesJ.cObject
        PairwiseInteractionSpline_MMMMEnergy ( self.cObject          ,
                                               cChargesI             ,
                                               cChargesJ             ,
                                               cLJTypesI             ,
                                               cLJTypesJ             ,
                                               cLJParameters         ,
                                               electrostaticScale    ,
                                               lennardJonesScale     ,
                                               coordinates3I.cObject ,
                                               coordinates3J.cObject ,
                                               pairList.cObject      ,
                                               &eElectrostatic       ,
                                               &eLennardJones        ,
                                               cGradients3I          ,
                                               cGradients3J          ,
                                               &cStatus               )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating MM energy." )
        return ( eElectrostatic, eLennardJones )

    def MMMMEnergyImage ( self, RealArray1D                charges                     ,
                                IntegerArray1D             ljTypes                     ,
                                LJParameterContainer       ljParameters                ,
                                                           electrostaticScale          ,
                                Coordinates3               coordinates3       not None ,
                                SymmetryParameters         symmetryParameters not None ,
                                ImagePairListContainer     imagePairLists     not None ,
                                Coordinates3               gradients3                  ,
                                SymmetryParameterGradients symmetryParameterGradients  ):
        """MM/MM image energy."""
        cdef CIntegerArray1D             *cLJTypes                    = NULL 
        cdef CRealArray2D                *cGradients3                 = NULL
        cdef CLJParameterContainer       *cLJParameters               = NULL
        cdef CReal                        eElectrostatic
        cdef CReal                        eLennardJones
        cdef CRealArray1D                *cCharges                    = NULL
        cdef CStatus                      cStatus                     = CStatus_OK
        cdef CSymmetryParameterGradients *cSymmetryParameterGradients = NULL
        if charges                    is not None: cCharges                    = charges.cObject
        if gradients3                 is not None: cGradients3                 = gradients3.cObject
        if ljParameters               is not None: cLJParameters               = ljParameters.cObject
        if ljTypes                    is not None: cLJTypes                    = ljTypes.cObject
        if symmetryParameterGradients is not None: cSymmetryParameterGradients = symmetryParameterGradients.cObject
        PairwiseInteractionSpline_MMMMEnergyImage ( self.cObject                ,
                                                    cCharges                    ,
                                                    cLJTypes                    ,
                                                    cLJParameters               ,
                                                    electrostaticScale          ,
                                                    coordinates3.cObject        ,
                                                    symmetryParameters.cObject  , 
                                                    imagePairLists.cObject      , 
                                                    &eElectrostatic             ,
                                                    &eLennardJones              ,
                                                    cGradients3                 ,
                                                    cSymmetryParameterGradients ,  
                                                    &cStatus                     )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating MM/MM image energy." )
        return ( eElectrostatic, eLennardJones )

    def MMMMEnergyMI ( self, RealArray1D                chargesI                    ,
                             RealArray1D                chargesJ                    ,
                             IntegerArray1D             ljTypesI                    ,
                             IntegerArray1D             ljTypesJ                    ,
                             LJParameterContainer       ljParameters                ,
                                                        electrostaticScale          ,
                                                        lennardJonesScale           ,
                             Coordinates3               coordinates3I      not None ,
                             Coordinates3               coordinates3J      not None ,
                             SymmetryParameters         symmetryParameters not None ,
                             PairList                   pairList           not None ,
                             Coordinates3               gradients3I                 ,
                             Coordinates3               gradients3J                 ,
                             SymmetryParameterGradients symmetryParameterGradients  ):
        """MM/MM energy within the minimum image convention."""
        cdef CIntegerArray1D             *cLJTypesI                   = NULL 
        cdef CIntegerArray1D             *cLJTypesJ                   = NULL 
        cdef CRealArray2D                *cGradients3I                = NULL 
        cdef CRealArray2D                *cGradients3J                = NULL 
        cdef CLJParameterContainer       *cLJParameters               = NULL
        cdef CReal                        eElectrostatic
        cdef CReal                        eLennardJones
        cdef CRealArray1D                *cChargesI                   = NULL
        cdef CRealArray1D                *cChargesJ                   = NULL
        cdef CStatus                      cStatus                     = CStatus_OK
        cdef CSymmetryParameterGradients *cSymmetryParameterGradients = NULL
        if chargesI                   is not None: cChargesI                   = chargesI.cObject
        if chargesJ                   is not None: cChargesJ                   = chargesJ.cObject
        if gradients3I                is not None: cGradients3I                = gradients3I.cObject
        if gradients3J                is not None: cGradients3J                = gradients3J.cObject
        if ljParameters               is not None: cLJParameters               = ljParameters.cObject
        if ljTypesI                   is not None: cLJTypesI                   = ljTypesI.cObject
        if ljTypesJ                   is not None: cLJTypesJ                   = ljTypesJ.cObject
        if symmetryParameterGradients is not None: cSymmetryParameterGradients = symmetryParameterGradients.cObject
        PairwiseInteractionSpline_MMMMEnergyMI ( self.cObject                ,
                                                 cChargesI                   ,
                                                 cChargesJ                   ,
                                                 cLJTypesI                   ,
                                                 cLJTypesJ                   ,
                                                 cLJParameters               ,
                                                 electrostaticScale          ,
                                                 lennardJonesScale           ,
                                                 coordinates3I.cObject       ,
                                                 coordinates3J.cObject       ,
                                                 symmetryParameters.cObject  , 
                                                 pairList.cObject            ,
                                                 &eElectrostatic             ,
                                                 &eLennardJones              ,
                                                 cGradients3I                ,
                                                 cGradients3J                ,
                                                 cSymmetryParameterGradients ,
                                                 &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating minimum image MM energy." )
        return ( eElectrostatic, eLennardJones )

    def QCMMGradients ( self, RealArray1D  chargesQ           not None ,
                              RealArray1D  chargesM           not None ,
                                           electrostaticScale          ,
                              Coordinates3 coordinates3Q      not None ,
                              Coordinates3 coordinates3M      not None ,
                              PairList     pairList           not None ,
                              Coordinates3 gradients3Q        not None ,
                              Coordinates3 gradients3M        not None ):
        """QC/MM gradients."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCMMGradients ( self.cObject          ,
                                                  chargesQ.cObject      ,
                                                  chargesM.cObject      ,
                                                  electrostaticScale    ,
                                                  coordinates3Q.cObject ,
                                                  coordinates3M.cObject ,
                                                  pairList.cObject      ,
                                                  gradients3Q.cObject   ,
                                                  gradients3M.cObject   ,
                                                  &cStatus               )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM gradients." )

    def QCMMGradientsImage ( self, RealArray1D                chargesA                   not None ,
                                   RealArray1D                chargesB                   not None ,
                                                              electrostaticScale                  ,
                                   Coordinates3               coordinates3A              not None ,
                                   Coordinates3               coordinates3B              not None ,
                                   SymmetryParameters         symmetryParameters         not None ,
                                   ImagePairListContainer     imagePairLists             not None ,
                                   Coordinates3               gradients3A                not None ,
                                   Coordinates3               gradients3B                not None ,
                                   SymmetryParameterGradients symmetryParameterGradients not None ):
        """QC/MM image gradients."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCMMGradientsImage ( self.cObject                       ,
                                                       chargesA.cObject                   ,
                                                       chargesB.cObject                   ,
                                                       electrostaticScale                 ,
                                                       coordinates3A.cObject              ,
                                                       coordinates3B.cObject              ,
                                                       symmetryParameters.cObject         , 
                                                       imagePairLists.cObject             ,
                                                       gradients3A.cObject                ,
                                                       gradients3B.cObject                ,
                                                       symmetryParameterGradients.cObject ,
                                                       &cStatus                            )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM image gradients." )

    def QCMMGradientsMI ( self, RealArray1D                chargesQ                   not None ,
                                RealArray1D                chargesM                   not None ,
                                                           electrostaticScale                  ,
                                Coordinates3               coordinates3Q              not None ,
                                Coordinates3               coordinates3M              not None ,
                                SymmetryParameters         symmetryParameters         not None ,
                                PairList                   pairList                   not None , 
                                Coordinates3               gradients3Q                not None , 
                                Coordinates3               gradients3M                not None ,
                                SymmetryParameterGradients symmetryParameterGradients not None ):
        """QC/MM gradients within the minimum image convention."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCMMGradientsMI ( self.cObject                       ,
                                                    chargesQ.cObject                   ,
                                                    chargesM.cObject                   ,
                                                    electrostaticScale                 ,
                                                    coordinates3Q.cObject              ,
                                                    coordinates3M.cObject              ,
                                                    symmetryParameters.cObject         ,
                                                    pairList.cObject                   ,
                                                    gradients3Q.cObject                ,
                                                    gradients3M.cObject                ,
                                                    symmetryParameterGradients.cObject ,
                                                    &cStatus                           )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating minimum image QC/MM gradients." )

    def QCMMPotentials ( self, RealArray1D  chargesM           not None ,
                                            electrostaticScale          ,
                               Coordinates3 coordinates3Q      not None ,
                               Coordinates3 coordinates3M      not None ,
                               PairList     pairList           not None ,
                               RealArray1D  potentials         not None ):
        """QC/MM potentials."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCMMPotentials ( self.cObject          ,
                                                   chargesM.cObject      ,
                                                   electrostaticScale    ,
                                                   coordinates3Q.cObject ,
                                                   coordinates3M.cObject ,
                                                   pairList.cObject      ,
                                                   potentials.cObject    ,
                                                   &cStatus               )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM potentials." )

    def QCMMPotentialsImage ( self, RealArray1D            charges            not None ,
                                                           electrostaticScale          ,
                                    Coordinates3           coordinates3A      not None ,
                                    Coordinates3           coordinates3B      not None ,
                                    SymmetryParameters     symmetryParameters not None ,
                                    ImagePairListContainer imagePairLists     not None ,
                                    RealArray1D            potentials         not None ):
        """QC/MM image potentials."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCMMPotentialsImage ( self.cObject               ,
                                                        charges.cObject            ,
                                                        electrostaticScale         ,
                                                        coordinates3A.cObject      ,
                                                        coordinates3B.cObject      ,
                                                        symmetryParameters.cObject , 
                                                        imagePairLists.cObject     ,
                                                        potentials.cObject         ,
                                                        &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM image potentials." )

    def QCQCGradientsImage ( self, RealArray1D                charges                    not None ,
                                                              electrostaticScale                  ,
                                   Coordinates3               coordinates3               not None ,
                                   SymmetryParameters         symmetryParameters         not None ,
                                   ImagePairListContainer     imagePairLists             not None ,
                                   Coordinates3               gradients3                 not None ,
                                   SymmetryParameterGradients symmetryParameterGradients not None ):
        """QC/QC image gradients."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCQCGradientsImage ( self.cObject                       ,
                                                       charges.cObject                    ,
                                                       electrostaticScale                 ,
                                                       coordinates3.cObject               ,
                                                       symmetryParameters.cObject         , 
                                                       imagePairLists.cObject             ,
                                                       gradients3.cObject                 ,
                                                       symmetryParameterGradients.cObject ,
                                                       &cStatus                            )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/QC image gradients." )

    def QCQCPotentialsImage ( self,                        electrostaticScale          ,
                                    Coordinates3           coordinates3       not None ,
                                    SymmetryParameters     symmetryParameters not None ,
                                    ImagePairListContainer imagePairLists     not None ,
                                    SymmetricMatrix        potentials         not None ):
        """QC/QC image potentials."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCQCPotentialsImage ( self.cObject               ,
                                                        electrostaticScale         ,
                                                        coordinates3.cObject       ,
                                                        symmetryParameters.cObject , 
                                                        imagePairLists.cObject     ,
                                                        potentials.cObject         ,
                                                        &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/QC image potentials." )

    def QCMMPotentialsMI ( self, RealArray1D        chargesM           not None ,
                                                    electrostaticScale          ,
                                 Coordinates3       coordinates3Q      not None ,
                                 Coordinates3       coordinates3M      not None ,
                                 SymmetryParameters symmetryParameters not None ,
                                 PairList           pairList           not None ,
                                 RealArray1D        potentials         not None ):
        """QC/MM potentials within the minimum image convention."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionSpline_QCMMPotentialsMI ( self.cObject               ,
                                                     chargesM.cObject           ,
                                                     electrostaticScale         ,
                                                     coordinates3Q.cObject      ,
                                                     coordinates3M.cObject      ,
                                                     symmetryParameters.cObject , 
                                                     pairList.cObject           ,
                                                     potentials.cObject         ,
                                                     &cStatus                   )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating minimum image QC/MM potentials." )

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "dampingCutOff"      in options: self.dampingCutOff      = options.pop ( "dampingCutOff"      )
        if "electrostaticModel" in options: self.electrostaticModel = options.pop ( "electrostaticModel" )
        if "innerCutOff"        in options: self.innerCutOff        = options.pop ( "innerCutOff"        )
        if "outerCutOff"        in options: self.outerCutOff        = options.pop ( "outerCutOff"        )
        if "pointDensity"       in options: self.pointDensity       = options.pop ( "pointDensity"       )
        if "width1"             in options: self.width1             = options.pop ( "width1"             )
        if "width2"             in options: self.width2             = options.pop ( "width2"             )
        if len ( options ) > 0: raise NBModelError ( "Invalid options: " + ", ".join ( sorted ( options.keys ( ) ) ) + "." )
        self._CheckOptions ( )
        self.electrostaticSpline = self.MakeElectrostaticSpline ( )
        self.lennardJonesASpline = self.MakeLennardJonesASpline ( )
        self.lennardJonesBSpline = self.MakeLennardJonesBSpline ( )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Spline ABFS Pairwise Interaction Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Damping CutOff"       , "{:.3f}".format ( self.dampingCutOff            ) ) ,
                 ( "Electrostatic Model"  , self.electrostaticModel.name.replace ( "_", "/" ) ) ,
                 ( "Inner CutOff"         , "{:.3f}".format ( self.innerCutOff              ) ) ,
                 ( "Outer CutOff"         , "{:.3f}".format ( self.outerCutOff              ) ) ,
                 ( "Spline Point Density" , "{:d}"  .format ( self.pointDensity             ) ) ,
                 ( "Width 1"              , "{:.3f}".format ( self.width1                   ) ) ,
                 ( "Width 2"              , "{:.3f}".format ( self.width2                   ) ) ]

    @property
    def range ( self ): return self.outerCutOff

"""Monte Carlo pairwise interactions."""

from pCore              import logFile              , \
                               LogFileActive        , \
                               RawObjectConstructor
from pMolecule.NBModel  import NBModelError
from pScientific.Arrays import Array

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionMonteCarlo ( PairwiseInteraction ):
    """Monte Carlo pairwise interactions."""

    def __copy__ ( self ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        cdef PairwiseInteractionMonteCarlo new
        new = self.__class__.Raw ( )
        new.cObject = PairwiseInteractionMonteCarlo_Clone ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise NBModelError ( "Error cloning pairwise interaction." )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairwiseInteractionMonteCarlo_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        return { "isolateScale" : self.cObject.isolateScale ,
                 "buffer"       : self.cObject.buffer       ,
                 "chargeScale"  : self.cObject.chargeScale  ,
                 "cutOff"       : self.cObject.cutOff       ,
                 "epsilonScale" : self.cObject.epsilonScale ,
                 "sigmaScale"   : self.cObject.sigmaScale   ,
                 "underFlowL"   : self.cObject.underFlowL   ,
                 "underFlowQ"   : self.cObject.underFlowQ   }

    def __str__ ( self ): return "Monte Carlo Pairwise Interaction"

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = PairwiseInteractionMonteCarlo_Allocate ( &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise NBModelError ( "Error allocating pairwise interaction." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def Interactions ( self, RealArray1D r not None ):
        """Calculate the interaction components for an array of distances."""
        cdef RealArray1D electrostatic, lennardJonesA, lennardJonesB
        extent        = len ( r )
        electrostatic = Array.WithExtent ( extent )
        lennardJonesA = Array.WithExtent ( extent )
        lennardJonesB = Array.WithExtent ( extent )
        PairwiseInteractionMonteCarlo_Interactions ( self.cObject          ,       
                                                     r.cObject             ,
                                                     electrostatic.cObject ,
                                                     lennardJonesA.cObject ,
                                                     lennardJonesB.cObject )
        return { "Electrostatic"   : electrostatic ,
                 "Lennard-Jones A" : lennardJonesA ,
                 "Lennard-Jones B" : lennardJonesB }

    def MMMMEnergy ( self, RealArray1D          charges            not None ,
                           IntegerArray1D       ljTypes            not None ,
                           LJParameterContainer ljParameters       not None ,
                                                electrostaticScale          ,
                           SelectionContainer   isolates           not None ,
                           Coordinates3         coordinates3       not None ,
                           SymmetryParameters   symmetryParameters not None ):
        """MM/MM energy."""
        cdef CReal   eElectrostatic
        cdef CReal   eLennardJones
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionMonteCarlo_MMMMEnergy ( self.cObject               ,
                                                   charges.cObject            ,
                                                   ljTypes.cObject            ,
                                                   ljParameters.cObject       ,
                                                   electrostaticScale         ,
                                                   1.0e+00                    ,
                                                   isolates.cObject           ,
                                                   NULL                       ,
                                                   coordinates3.cObject       ,
                                                   symmetryParameters.cObject ,
                                                   &eElectrostatic            ,
                                                   &eLennardJones             ,
                                                   &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating MM energy." )
        return ( eElectrostatic, eLennardJones )

    def MMMMIsolateEnergy ( self,                      isolate                     ,
                                  RealArray1D          charges            not None ,
                                  IntegerArray1D       ljTypes            not None ,
                                  LJParameterContainer ljParameters       not None ,
                                                       electrostaticScale          ,
                                  SelectionContainer   isolates           not None ,
                                  Coordinates3         coordinates3       not None ,
                                  SymmetryParameters   symmetryParameters not None ):
        """MM/MM energy of a single isolate."""
        cdef CReal   eElectrostatic
        cdef CReal   eLennardJones
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionMonteCarlo_MMMMIsolateEnergy ( self.cObject               ,
                                                          isolate                    ,
                                                          charges.cObject            ,
                                                          ljTypes.cObject            ,
                                                          ljParameters.cObject       ,
                                                          electrostaticScale         ,
                                                          1.0e+00                    ,
                                                          isolates.cObject           ,
                                                          NULL                       ,
                                                          coordinates3.cObject       ,
                                                          symmetryParameters.cObject ,
                                                          &eElectrostatic            ,
                                                          &eLennardJones             ,
                                                          &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating MM isolate energy." )
        return ( eElectrostatic, eLennardJones )

    def OptionRecords ( self ):
        """Option records and subobjects that also have options."""
        return ( [ ( "cutOff"      , "Cut-Off"                 , "float", "{:.3f}".format ( self.cObject.cutOff       ) ) ,
                   ( "buffer"      , "Buffer"                  , "float", "{:.3f}".format ( self.cObject.buffer       ) ) ,
                   ( "underFlowQ"  , "Electrostatic Underflow" , "float", "{:.3f}".format ( self.cObject.underFlowQ   ) ) ,
                   ( "underFlowL"  , "Lennard-Jones Underflow" , "float", "{:.3f}".format ( self.cObject.underFlowL   ) ) ,
                   ( "isolateScale", "Isolate for Scaling"     , "int"  , "{:d}".format   ( self.cObject.isolateScale ) ) ,
                   ( "chargeScale" , "Charge Scaling"          , "float", "{:.3f}".format ( self.cObject.chargeScale  ) ) ,
                   ( "epsilonScale", "Epsilon Scaling"         , "float", "{:.3f}".format ( self.cObject.epsilonScale ) ) ,
                   ( "sigmaScale"  , "Sigma Scaling"           , "float", "{:.3f}".format ( self.cObject.sigmaScale   ) ) ], [] )

    def ScaleIsolateInteractionParameters ( self, isolate, chargeScale, epsilonScale, sigmaScale, log = logFile ):
        """Scale the isolate interaction parameters."""
        if self.cObject != NULL:
            PairwiseInteractionMonteCarlo_ScaleIsolateInteractionParameters ( self.cObject           ,
                                                                              int   ( isolate      ) ,
                                                                              float ( chargeScale  ) ,
                                                                              float ( epsilonScale ) ,
                                                                              float ( sigmaScale   ) )
            if LogFileActive ( log ): log.Paragraph ( "Interaction parameters for isolate {:d} scaled.".format ( isolate ) )

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "isolateScale" in options: self.cObject.isolateScale = options.pop ( "isolateScale" )
        if "buffer"       in options: self.cObject.buffer       = options.pop ( "buffer"       )
        if "chargeScale"  in options: self.cObject.chargeScale  = options.pop ( "chargeScale"  )
        if "cutOff"       in options: self.cObject.cutOff       = options.pop ( "cutOff"       )
        if "epsilonScale" in options: self.cObject.epsilonScale = options.pop ( "epsilonScale" )
        if "sigmaScale"   in options: self.cObject.sigmaScale   = options.pop ( "sigmaScale"   )
        if "underFlowL"   in options: self.cObject.underFlowL   = options.pop ( "underFlowL"   )
        if "underFlowQ"   in options: self.cObject.underFlowQ   = options.pop ( "underFlowQ"   )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Monte Carlo Pairwise Interaction Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Cut-Off"                 , "{:.3f}".format ( self.cObject.cutOff       ) ) ,
                 ( "Buffer"                  , "{:.3f}".format ( self.cObject.buffer       ) ) ,
                 ( "Electrostatic Underflow" , "{:.3f}".format ( self.cObject.underFlowQ   ) ) ,
                 ( "Lennard-Jones Underflow" , "{:.3f}".format ( self.cObject.underFlowL   ) ) ,
                 ( "Isolate for Scaling"     , "{:d}".format   ( self.cObject.isolateScale ) ) ,
                 ( "Charge Scaling"          , "{:.3f}".format ( self.cObject.chargeScale  ) ) ,
                 ( "Epsilon Scaling"         , "{:.3f}".format ( self.cObject.epsilonScale ) ) ,
                 ( "Sigma Scaling"           , "{:.3f}".format ( self.cObject.sigmaScale   ) ) ]

    @property
    def  range ( self ):
        if self.cObject == NULL: return 0.0
        else:                    return self.cObject.cutOff

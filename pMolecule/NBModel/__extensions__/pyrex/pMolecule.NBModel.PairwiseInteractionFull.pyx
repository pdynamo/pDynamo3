"""Full pairwise interactions."""

from pCore              import logFile       , \
                               LogFileActive , \
                               RawObjectConstructor
from pMolecule.NBModel  import NBModelError
from pScientific.Arrays import Array

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionFull ( PairwiseInteraction ):
    """Base class for pairwise interactions."""

    def __copy__ ( self ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        cdef PairwiseInteractionFull new
        new = self.__class__.Raw ( )
        new.cObject = PairwiseInteractionFull_Clone ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise NBModelError ( "Error cloning pairwise interaction." )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairwiseInteractionFull_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        return { "dampingCutOff" : self.cObject.dampingCutOff }

    def __str__ ( self ): return "Full Pairwise Interaction"

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = PairwiseInteractionFull_Allocate ( &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise NBModelError ( "Error allocating pairwise interaction." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def Interactions ( self, RealArray1D r not None ):
        """Calculate the interaction components for an array of distances."""
        cdef RealArray1D lennardJonesA
        cdef RealArray1D lennardJonesB
        cdef RealArray1D multipole0
        cdef RealArray1D multipole1
        cdef RealArray1D multipole2
        extent        = len ( r )
        lennardJonesA = Array.WithExtent ( extent )
        lennardJonesB = Array.WithExtent ( extent )
        multipole0    = Array.WithExtent ( extent )
        multipole1    = Array.WithExtent ( extent )
        multipole2    = Array.WithExtent ( extent )
        PairwiseInteractionFull_Interactions ( self.cObject          ,       
                                               r.cObject             ,
                                               lennardJonesA.cObject ,
                                               lennardJonesB.cObject ,
                                               multipole0.cObject    ,
                                               multipole1.cObject    ,
                                               multipole2.cObject    )
        return { "Lennard-Jones A" : lennardJonesA ,
                 "Lennard-Jones B" : lennardJonesB ,
                 "Multipole 0"     : multipole0    ,
                 "Multipole 1"     : multipole1    ,
                 "Multipole 2"     : multipole2    }

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
        cdef CIntegerArray1D       *cLJTypesI     = NULL
        cdef CIntegerArray1D       *cLJTypesJ     = NULL
        cdef CRealArray2D          *cGradients3I  = NULL
        cdef CRealArray2D          *cGradients3J  = NULL
        cdef CLJParameterContainer *cLJParameters = NULL
        cdef CReal                  eElectrostatic
        cdef CReal                  eLennardJones
        cdef CRealArray1D          *cChargesI     = NULL
        cdef CRealArray1D          *cChargesJ     = NULL
        cdef CStatus                cStatus        = CStatus_OK
        if chargesI     is not None: cChargesI     = chargesI.cObject
        if chargesJ     is not None: cChargesJ     = chargesJ.cObject
        if gradients3I  is not None: cGradients3I  = gradients3I.cObject
        if gradients3J  is not None: cGradients3J  = gradients3J.cObject
        if ljParameters is not None: cLJParameters = ljParameters.cObject
        if ljTypesI     is not None: cLJTypesI     = ljTypesI.cObject
        if ljTypesJ     is not None: cLJTypesJ     = ljTypesJ.cObject
        PairwiseInteractionFull_MMMMEnergy ( self.cObject          ,
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

    def OptionRecords ( self ):
        """Option records and subobjects that also have options."""
        return ( [ ( "dampingCutOff", "Damping Cut-Off", "float", "{:.3f}".format ( self.cObject.dampingCutOff ) ) ], [] )

    def QCMMGradients ( self,              multipoleOrder              ,
                              RealArray1D  multipolesQ        not None ,
                              RealArray1D  chargesM           not None ,
                                           electrostaticScale          ,
                              Coordinates3 coordinates3Q      not None ,
                              Coordinates3 coordinates3M      not None ,
                              PairList     pairList           not None ,
                              Coordinates3 gradients3Q        not None ,
                              Coordinates3 gradients3M        not None ):
        """QC/MM gradients."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionFull_QCMMGradients ( self.cObject          ,
                                                multipoleOrder        ,
                                                multipolesQ.cObject   ,
                                                chargesM.cObject      ,
                                                electrostaticScale    ,
                                                coordinates3Q.cObject ,
                                                coordinates3M.cObject ,
                                                pairList.cObject      ,
                                                gradients3Q.cObject   ,
                                                gradients3M.cObject   ,
                                                &cStatus               )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM gradients." )

    def QCMMPotentials ( self,              multipoleOrder              ,
                               RealArray1D  chargesM           not None ,
                                            electrostaticScale          ,
                               Coordinates3 coordinates3Q      not None ,
                               Coordinates3 coordinates3M      not None ,
                               PairList     pairList           not None ,
                               RealArray1D  potentials         not None ):
        """QC/MM potentials."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionFull_QCMMPotentials ( self.cObject          ,
                                                 multipoleOrder        ,
                                                 chargesM.cObject      ,
                                                 electrostaticScale    ,
                                                 coordinates3Q.cObject ,
                                                 coordinates3M.cObject ,
                                                 pairList.cObject      ,
                                                 potentials.cObject    ,
                                                 &cStatus               )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM potentials." )

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "dampingCutOff" in options: self.cObject.dampingCutOff = options.pop ( "dampingCutOff" )
        PairwiseInteractionFull_InitializeDependent ( self.cObject )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Full Pairwise Interaction Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Damping Cut-Off", "{:.3f}".format ( self.cObject.dampingCutOff  ) ) ]

    @property
    def range ( self ): return 0.0


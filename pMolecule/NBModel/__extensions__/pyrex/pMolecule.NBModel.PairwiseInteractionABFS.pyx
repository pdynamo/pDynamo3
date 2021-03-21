"""Analytic ABFS pairwise interactions."""

from pCore              import logFile       , \
                               LogFileActive , \
                               RawObjectConstructor
from pMolecule.NBModel  import NBModelError
from pScientific.Arrays import Array

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionABFS ( PairwiseInteraction ):

    def __copy__ ( self ):
        """Copying."""
        cdef CStatus cStatus = CStatus_OK
        cdef PairwiseInteractionABFS new
        new = self.__class__.Raw ( )
        new.cObject = PairwiseInteractionABFS_Clone ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise NBModelError ( "Error cloning pairwise interaction." )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairwiseInteractionABFS_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        return { "dampingCutOff" : self.cObject.dampingCutOff ,
                 "innerCutOff"   : self.cObject.innerCutOff   ,
                 "outerCutOff"   : self.cObject.outerCutOff   }

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = PairwiseInteractionABFS_Allocate ( &cStatus )
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
        PairwiseInteractionABFS_Interactions ( self.cObject          ,       
                                               r.cObject             ,
                                               electrostatic.cObject ,
                                               lennardJonesA.cObject ,
                                               lennardJonesB.cObject )
        return { "Electrostatic"   : electrostatic ,
                 "Lennard-Jones A" : lennardJonesA ,
                 "Lennard-Jones B" : lennardJonesB }

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
        cdef CStatus                cStatus         = CStatus_OK
        if chargesI     is not None: cChargesI     = chargesI.cObject
        if chargesJ     is not None: cChargesJ     = chargesJ.cObject
        if gradients3I  is not None: cGradients3I  = gradients3I.cObject
        if gradients3J  is not None: cGradients3J  = gradients3J.cObject
        if ljParameters is not None: cLJParameters = ljParameters.cObject
        if ljTypesI     is not None: cLJTypesI     = ljTypesI.cObject
        if ljTypesJ     is not None: cLJTypesJ     = ljTypesJ.cObject
        PairwiseInteractionABFS_MMMMEnergy ( self.cObject          ,
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
                                             &cStatus              )
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
        PairwiseInteractionABFS_MMMMEnergyImage ( self.cObject                ,
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
                                                  &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating MM image energy." )
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
        PairwiseInteractionABFS_MMMMEnergyMI ( self.cObject                ,
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
        PairwiseInteractionABFS_QCMMGradients ( self.cObject          ,
                                                chargesQ.cObject      ,
                                                chargesM.cObject      ,
                                                electrostaticScale    ,
                                                coordinates3Q.cObject ,
                                                coordinates3M.cObject ,
                                                pairList.cObject      ,
                                                gradients3Q.cObject   ,
                                                gradients3M.cObject   ,
                                                &cStatus              )
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
        PairwiseInteractionABFS_QCMMGradientsImage ( self.cObject                       ,
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
                                                     &cStatus                           )
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
        PairwiseInteractionABFS_QCMMGradientsMI ( self.cObject                       ,
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
        PairwiseInteractionABFS_QCMMPotentials ( self.cObject          ,
                                                 chargesM.cObject      ,
                                                 electrostaticScale    ,
                                                 coordinates3Q.cObject ,
                                                 coordinates3M.cObject ,
                                                 pairList.cObject      ,
                                                 potentials.cObject    ,
                                                 &cStatus              )
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
        PairwiseInteractionABFS_QCMMPotentialsImage ( self.cObject               ,
                                                      charges.cObject            ,
                                                      electrostaticScale         ,
                                                      coordinates3A.cObject      ,
                                                      coordinates3B.cObject      ,
                                                      symmetryParameters.cObject , 
                                                      imagePairLists.cObject     ,
                                                      potentials.cObject         ,
                                                      &cStatus                   )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/MM image potentials." )

    def QCMMPotentialsMI ( self, RealArray1D        chargesM           not None ,
                                                    electrostaticScale          ,
                                 Coordinates3       coordinates3Q      not None ,
                                 Coordinates3       coordinates3M      not None ,
                                 SymmetryParameters symmetryParameters not None ,
                                 PairList           pairList           not None ,
                                 RealArray1D        potentials         not None ):
        """QC/MM potentials within the minimum image convention."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionABFS_QCMMPotentialsMI ( self.cObject               ,
                                                   chargesM.cObject           ,
                                                   electrostaticScale         ,
                                                   coordinates3Q.cObject      ,
                                                   coordinates3M.cObject      ,
                                                   symmetryParameters.cObject , 
                                                   pairList.cObject           ,
                                                   potentials.cObject         ,
                                                   &cStatus                   )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating minimum image QC/MM potentials." )

    def QCQCGradientsImage ( self, RealArray1D                charges                    not None ,
                                                              electrostaticScale                  ,
                                   Coordinates3               coordinates3               not None ,
                                   SymmetryParameters         symmetryParameters         not None ,
                                   ImagePairListContainer     imagePairLists             not None ,
                                   Coordinates3               gradients3                 not None ,
                                   SymmetryParameterGradients symmetryParameterGradients not None ):
        """QC/QC image gradients."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionABFS_QCQCGradientsImage ( self.cObject                       ,
                                                     charges.cObject                    ,
                                                     electrostaticScale                 ,
                                                     coordinates3.cObject               ,
                                                     symmetryParameters.cObject         , 
                                                     imagePairLists.cObject             ,
                                                     gradients3.cObject                 ,
                                                     symmetryParameterGradients.cObject ,
                                                     &cStatus                           )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/QC image gradients." )

    def QCQCPotentialsImage ( self,                        electrostaticScale          ,
                                    Coordinates3           coordinates3       not None ,
                                    SymmetryParameters     symmetryParameters not None ,
                                    ImagePairListContainer imagePairLists     not None ,
                                    SymmetricMatrix        potentials         not None ):
        """QC/QC image potentials."""
        cdef CStatus cStatus = CStatus_OK
        PairwiseInteractionABFS_QCQCPotentialsImage ( self.cObject               ,
                                                      electrostaticScale         ,
                                                      coordinates3.cObject       ,
                                                      symmetryParameters.cObject , 
                                                      imagePairLists.cObject     ,
                                                      potentials.cObject         ,
                                                      &cStatus                   )
        if cStatus != CStatus_OK: raise NBModelError ( "Error calculating QC/QC image potentials." )

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "dampingCutOff" in options: self.cObject.dampingCutOff = options.pop ( "dampingCutOff" )
        if "innerCutOff"   in options: self.cObject.innerCutOff   = options.pop ( "innerCutOff"   )
        if "outerCutOff"   in options: self.cObject.outerCutOff   = options.pop ( "outerCutOff"   )
        if ( self.cObject.dampingCutOff < 0.0                        ) or \
           ( self.cObject.innerCutOff   < self.cObject.dampingCutOff ) or \
           ( self.cObject.outerCutOff   < self.cObject.innerCutOff   ):
            raise NBModelError ( "Invalid cutOff values: damping - {:.3f}; inner - {:.3f}; outer - {:.3f}.".format ( self.cObject.dampingCutOff ,
                                                                                                                     self.cObject.innerCutOff   ,
                                                                                                                     self.cObject.outerCutOff ) )
 
    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "ABFS Pairwise Interaction Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Damping Cut-Off", "{:.3f}".format ( self.cObject.dampingCutOff ) ) ,
                 ( "Inner Cut-Off"  , "{:.3f}".format ( self.cObject.innerCutOff   ) ) ,
                 ( "Outer Cut-Off"  , "{:.3f}".format ( self.cObject.outerCutOff   ) ) ]

    @property
    def  range ( self ):
        if self.cObject == NULL: return 0.0
        else:                    return self.cObject.outerCutOff

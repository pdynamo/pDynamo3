"""Image pairlist container."""

from pCore             import logFile       , \
                              LogFileActive
from pMolecule.NBModel import NBModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ImagePairListContainer:
    """Image pairlist container."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            ImagePairListContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

    def __len__ ( self ): return self.numberOfPairs

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Constructor ( selfClass, PairListGenerator        generator          not None ,
                                 Selection                atomsA                      ,
                                 Selection                atomsB                      ,
                                 Selection                freeAtoms                   ,
                                 Coordinates3             coordinates3A      not None ,
                                 Coordinates3             coordinates3B      not None ,
                                 SymmetryParameters       symmetryParameters not None ,
                                 Transformation3Container transformations    not None ,
                                 ImageScanContainer       scanData           not None ,
                                 RegularGrid              gridA                       ,
                                 RegularGridOccupancy     occupancyA                  ,
                                                          checkForInverses            ):
        """Constructor given coordinates, symmetry parameters and transformations."""
        cdef ImagePairListContainer self
        cdef CBoolean                 cCheckForInverses
        cdef CImagePairListContainer *cObject           = NULL
        cdef CRegularGrid            *cGridA            = NULL
        cdef CRegularGridOccupancy   *cOccupancyA       = NULL
        cdef CSelection              *cAtomsA           = NULL
        cdef CSelection              *cAtomsB           = NULL
        cdef CSelection              *cFreeAtoms        = NULL
        cdef CStatus                  cStatus            = CStatus_OK
        if checkForInverses      : cCheckForInverses = CTrue
        else:                      cCheckForInverses = CFalse
        if atomsA     is not None: cAtomsA           = atomsA.cObject
        if atomsB     is not None: cAtomsB           = atomsB.cObject
        if freeAtoms  is not None: cFreeAtoms        = freeAtoms.cObject
        if gridA      is not None: cGridA            = gridA.cObject
        if occupancyA is not None: cOccupancyA       = occupancyA.cObject
        cObject = ImagePairListContainer_Constructor ( generator.cObject          ,
                                                       cAtomsA                    ,
                                                       cAtomsB                    ,
                                                       cFreeAtoms                 ,
                                                       coordinates3A.cObject      ,
                                                       coordinates3B.cObject      ,
                                                       symmetryParameters.cObject ,
                                                       transformations.cObject    ,
                                                       scanData.cObject           ,
                                                       cGridA                     ,
                                                       cOccupancyA                ,
                                                       cCheckForInverses          ,
                                                       &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error creating image pairlist container." )
        if cObject != NULL:
            self         = selfClass.Raw ( )
            self.cObject = cObject
            self.isOwner = True
        return self

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @property
    def numberOfImages ( self ): return ImagePairListContainer_NumberOfImages ( self.cObject )
    @property
    def numberOfPairs  ( self ): return ImagePairListContainer_NumberOfPairs  ( self.cObject )

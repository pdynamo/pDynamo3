"""Image pairlist container."""

from pCore             import logFile       , \
                              LogFileActive
from pMolecule.NBModel import NBModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ImageScanContainer:
    """Image pairlist container."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            ImageScanContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

    def _GetImages ( self ):
        """Return a list of images."""
        cdef CImageScan *record
        cdef CInteger    i
        images = []
        if ( self.cObject != NULL ) and ( self.cObject.records != NULL ):
            for i from 0 <= i < self.cObject.count:
                record = &(self.cObject.records[i])
                images.append ( ( record.doSkip == CTrue ,
                                  record.t               ,
                                  record.a               ,
                                  record.b               ,
                                  record.c               ,
                                  record.scale           ) )
        return images

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Constructor ( selfClass, Coordinates3             coordinates3       not None ,
                                 SymmetryParameters       symmetryParameters not None ,
                                 Transformation3Container transformations    not None ,
                                                          cutOff                      ,
                                                          checkForInverses            ,
                                                          expandFactor                ):
        """Constructor given coordinates, symmetry parameters and transformations."""
        cdef ImageScanContainer   self
        cdef CBoolean             cCheckForInverses
        cdef CImageScanContainer *cObject           = NULL
        cdef CInteger             cExpandFactor     = int ( expandFactor )
        cdef CReal                cCutOff           = float ( cutOff )
        cdef CStatus              cStatus           = CStatus_OK
        if checkForInverses: cCheckForInverses = CTrue
        else:                cCheckForInverses = CFalse
        cObject = ImageScanContainer_Constructor ( cCutOff                    ,
                                                   coordinates3.cObject       ,
                                                   symmetryParameters.cObject ,
                                                   transformations.cObject    ,
                                                   cCheckForInverses          ,
                                                   cExpandFactor              ,
                                                   &cStatus                   )
        if cStatus != CStatus_OK: raise NBModelError ( "Error creating image scan container." )
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

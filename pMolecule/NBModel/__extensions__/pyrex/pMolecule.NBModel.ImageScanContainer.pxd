from pCore.CPrimitiveTypes                          cimport CBoolean                  , \
                                                            CFalse                    , \
                                                            CInteger                  , \
                                                            CReal                     , \
                                                            CTrue
from pCore.Status                                   cimport CStatus                   , \
                                                            CStatus_OK
from pScientific.Arrays.RealArray2D                 cimport CRealArray2D
from pScientific.Geometry3.Coordinates3             cimport Coordinates3
from pScientific.Geometry3.Transformation3Container cimport CTransformation3Container , \
                                                            Transformation3Container
from pScientific.Symmetry.SymmetryParameters        cimport CSymmetryParameters       , \
                                                            SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ImageScanContainer.h":

    ctypedef struct CImageScan "ImageScan":
        CBoolean  doSkip
        CInteger  t
        CInteger  a
        CInteger  b
        CInteger  c
        CReal     scale

    ctypedef struct CImageScanContainer "ImageScanContainer":
        CInteger    capacity
        CInteger    count
        CReal       cutOff
        CImageScan *records

    cdef CImageScanContainer *ImageScanContainer_Allocate    ( CStatus                   *status             )
    cdef void                 ImageScanContainer_Deallocate  ( CImageScanContainer      **self               )
    cdef CImageScanContainer *ImageScanContainer_Constructor ( CReal                      cutOff             ,
                                                               CRealArray2D              *coordinates3       ,
                                                               CSymmetryParameters       *symmetryParameters ,
                                                               CTransformation3Container *transformations    ,
                                                               CBoolean                   checkForInverses   ,
                                                               CInteger                   expandFactor       ,
                                                               CStatus                   *status             )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ImageScanContainer:

    cdef CImageScanContainer *cObject
    cdef public object        isOwner

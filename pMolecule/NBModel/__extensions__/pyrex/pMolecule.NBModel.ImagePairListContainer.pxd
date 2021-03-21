from pCore.CPrimitiveTypes                          cimport CBoolean                  , \
                                                            CFalse                    , \
                                                            CInteger                  , \
                                                            CReal                     , \
                                                            CTrue
from pCore.Selection                                cimport CSelection                , \
                                                            Selection
from pCore.Status                                   cimport CStatus                   , \
                                                            CStatus_OK
from pMolecule.NBModel.ImageScanContainer           cimport CImageScanContainer       , \
                                                            ImageScanContainer
from pScientific.Arrays.RealArray2D                 cimport CRealArray2D
from pScientific.Geometry3.Coordinates3             cimport Coordinates3
from pScientific.Geometry3.PairListGenerator        cimport CPairListGenerator        , \
                                                            PairListGenerator
from pScientific.Geometry3.RegularGrid              cimport CRegularGrid              , \
                                                            RegularGrid
from pScientific.Geometry3.RegularGridOccupancy     cimport CRegularGridOccupancy     , \
                                                            RegularGridOccupancy
from pScientific.Geometry3.Transformation3Container cimport CTransformation3Container , \
                                                            Transformation3Container
from pScientific.Symmetry.SymmetryParameters        cimport CSymmetryParameters       , \
                                                            SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ImagePairListContainer.h":

    ctypedef struct CImagePairListContainer "ImagePairListContainer":
        pass

    cdef CImagePairListContainer *ImagePairListContainer_Allocate       ( CStatus                   *status             )
    cdef void                     ImagePairListContainer_Deallocate     ( CImagePairListContainer  **self               )
    cdef CInteger                 ImagePairListContainer_NumberOfImages ( CImagePairListContainer   *self               )
    cdef CInteger                 ImagePairListContainer_NumberOfPairs  ( CImagePairListContainer   *self               )
    cdef CImagePairListContainer *ImagePairListContainer_Constructor    ( CPairListGenerator        *generator          ,
                                                                          CSelection                *atomsA             ,
                                                                          CSelection                *atomsB             ,
                                                                          CSelection                *freeAtoms          ,
                                                                          CRealArray2D              *coordinates3A      ,
                                                                          CRealArray2D              *coordinates3B      ,
                                                                          CSymmetryParameters       *symmetryParameters ,
                                                                          CTransformation3Container *transformations    ,
                                                                          CImageScanContainer       *scanData           ,
                                                                          CRegularGrid              *gridA              ,
                                                                          CRegularGridOccupancy     *occupancyA         ,
                                                                          CBoolean                   checkForInverses   ,
                                                                          CStatus                   *status             )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class ImagePairListContainer:

    cdef CImagePairListContainer *cObject
    cdef public object            isOwner

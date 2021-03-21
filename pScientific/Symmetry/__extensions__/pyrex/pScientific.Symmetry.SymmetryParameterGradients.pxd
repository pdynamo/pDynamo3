from pCore.CPrimitiveTypes                   cimport CBoolean            , \
                                                     CFalse              , \
                                                     CInteger            , \
                                                     CReal               , \
                                                     CTrue
from pCore.Status                            cimport CStatus             , \
                                                     CStatus_OK
from pScientific.Arrays.RealArray2D          cimport CRealArray2D
from pScientific.Geometry3.Coordinates3      cimport Coordinates3
from pScientific.Geometry3.Matrix33          cimport Matrix33
from pScientific.Symmetry.SymmetryParameters cimport CSymmetryParameters , \
                                                     SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SymmetryParameterGradients.h":

    ctypedef struct CSymmetryParameterGradients "SymmetryParameterGradients":
        CReal         dEda
        CReal         dEdb
        CReal         dEdc
        CReal         dEdalpha
        CReal         dEdbeta
        CReal         dEdgamma
        CRealArray2D *dEdH

    cdef CSymmetryParameterGradients *SymmetryParameterGradients_AllocateWithMatrix    ( CRealArray2D                 *dEdH               ,
                                                                                         CStatus                      *status             )
    cdef void                         SymmetryParameterGradients_Deallocate            ( CSymmetryParameterGradients **self               )
    cdef void                         SymmetryParameterGradients_CrystalDerivatives    ( CSymmetryParameterGradients  *self               ,
                                                                                         CSymmetryParameters          *symmetryParameters )
    cdef void                         SymmetryParameterGradients_FractionalDerivatives ( CSymmetryParameterGradients  *self               ,
                                                                                         CSymmetryParameters          *symmetryParameters ,
                                                                                         CRealArray2D                 *coordinates3       ,
                                                                                         CRealArray2D                 *gradients3         )
    cdef void                         SymmetryParameterGradients_Initialize            ( CSymmetryParameterGradients  *self               )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetryParameterGradients:

    cdef CSymmetryParameterGradients *cObject
    cdef public object                isOwner
    cdef public Matrix33              dEdH

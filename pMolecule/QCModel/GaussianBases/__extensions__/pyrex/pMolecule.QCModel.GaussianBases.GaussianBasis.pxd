from pCore.CPrimitiveTypes                                   cimport CBoolean                 , \
                                                                     CFalse                   , \
                                                                     CInteger                 , \
                                                                     CReal                    , \
                                                                     CTrue
from pCore.Status                                            cimport CStatus                  , \
                                                                     CStatus_OK
from pScientific.Arrays.BaseArray2D                          cimport CView2D                  , \
                                                                     View2D_GetColumns        , \
                                                                     View2D_GetRows
from pScientific.Arrays.RealArray2D                          cimport CRealArray2D             , \
                                                                     RealArray2D              , \
                                                                     RealArray2D_CopyTo       , \
                                                                     RealArray2D_Deallocate
from pScientific.LinearAlgebra.OrthogonalizingTransformation cimport COrthogonalizationMethod , \
                                                                     OrthogonalizationMethod_Symmetric

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasis.h":

    ctypedef enum CGaussianBasisOperator "GaussianBasisOperator":
        GaussianBasisOperator_AntiCoulomb = 1 ,
        GaussianBasisOperator_Coulomb     = 2 ,
        GaussianBasisOperator_Dipole      = 3 ,
        GaussianBasisOperator_Kinetic     = 4 ,
        GaussianBasisOperator_Overlap     = 5 ,
        GaussianBasisOperator_Poisson     = 6 ,
        GaussianBasisOperator_Quadrupole  = 7

    ctypedef enum CGaussianBasisType "GaussianBasisType":
        GaussianBasisType_Density = 1 ,
        GaussianBasisType_Orbital = 2

    ctypedef struct CPrimitive "Primitive":
        CReal *cCBF
        CReal *coefficients
        CReal *coefficients0
        CReal  exponent
        CReal  exponent0

    ctypedef struct CShell "Shell":
        CInteger      lHigh
        CInteger      lLow
        CInteger      nBasis
        CInteger      nCBF
        CInteger      nPrimitives
        CInteger      nStart
        CInteger      nStartC
        CInteger     *cbfPowX
        CInteger     *cbfPowY
        CInteger     *cbfPowZ
        CRealArray2D *c2s
        CRealArray2D *s2c
        CPrimitive   *primitives

    ctypedef struct CGaussianBasis "GaussianBasis":
        CBoolean            pNormalized
        CBoolean            isSpherical
        CGaussianBasisType  basisType
        CInteger            atomicNumber
        CInteger            lHigh
        CInteger            nBasis
        CInteger            nCBF
        CInteger            nShells
        CInteger           *cbfPowX
        CInteger           *cbfPowY
        CInteger           *cbfPowZ
        CShell             *shells

    cdef CGaussianBasis *GaussianBasis_Allocate                           ( CInteger                 nShells         )
    cdef void            GaussianBasis_AllocateShell                      ( CGaussianBasis          *self            ,
                                                                            CInteger                 iShell          ,
                                                                            CInteger                 lHigh           ,
                                                                            CInteger                 lLow            ,
                                                                            CInteger                 nPrimitives     )
    cdef CGaussianBasis *GaussianBasis_Clone                              ( CGaussianBasis          *self            ,
                                                                            CStatus                 *status          )
    cdef void            GaussianBasis_Deallocate                         ( CGaussianBasis         **self            )
    cdef void            GaussianBasis_Finalize                           ( CGaussianBasis          *self            ,
                                                                            CStatus                 *status          )
    cdef void            GaussianBasis_ScaleShellExponents                ( CGaussianBasis          *self            ,
                                                                            CInteger                 index           ,
                                                                            CReal                    zeta            ,
                                                                            CStatus                 *status          )
    cdef CRealArray2D   *GaussianBasis_TransformationCartesianToSpherical ( CInteger                 lLow            ,
                                                                            CInteger                 lHigh           ,
                                                                            CStatus                 *status          )
    cdef CRealArray2D   *GaussianBasis_TransformationSphericalToCartesian ( CInteger                 lLow            ,
                                                                            CInteger                 lHigh           ,
                                                                            CRealArray2D            *c2s             ,
                                                                            CStatus                 *status          )
    cdef void            GaussianBasis_UnnormalizePrimitives              ( CGaussianBasis          *self            )   

cdef extern from "GaussianBasisOrthonormalize.h":

    cdef void            GaussianBasis_Orthonormalize                     ( CGaussianBasis          *self            ,
                                                                            CGaussianBasisOperator   operator        ,
                                                                            COrthogonalizationMethod method          ,
                                                                            CInteger                *nIndependent    ,
                                                                            CReal                   *deviation       ,
                                                                            CRealArray2D            *MOut            ,
                                                                            CRealArray2D            *XOut            ,
                                                                            CRealArray2D            *YOut            ,
                                                                            CStatus                 *status          )

cdef extern from "GaussianBasisRotate.h":

    cdef void GaussianBasis_MakeLRotations                                ( CInteger         L                       ,
                                                                            CRealArray2D    *R                       ,
                                                                            CRealArray2D    *Tc                      ,
                                                                            CStatus         *status                  )
    cdef void GaussianBasis_MakeRotationMatrix                            ( CGaussianBasis  *self                    ,
                                                                            CRealArray2D    *Tc                      ,
                                                                            CRealArray2D    *T                       ,
                                                                            CStatus         *status                  )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasis:

    cdef CGaussianBasis *cObject
    cdef public object   isOwner

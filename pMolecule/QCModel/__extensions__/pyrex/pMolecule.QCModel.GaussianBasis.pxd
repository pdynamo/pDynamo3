from pCore.CPrimitiveTypes          cimport CBoolean     , \
                                            CFalse       , \
                                            CInteger     , \
                                            CReal        , \
                                            CTrue
from pCore.Status                   cimport CStatus, CStatus_OK
from pScientific.Arrays.RealArray2D cimport CRealArray2D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasis.h":

    ctypedef enum CGaussianBasisType "GaussianBasisType":
        GaussianBasisType_Coulomb = 1 ,
        GaussianBasisType_Orbital = 2 ,
        GaussianBasisType_Poisson = 3

    ctypedef enum CNormalizationType "NormalizationType":
        NormalizationType_Canonical = 1 ,
        NormalizationType_Diagonal  = 2 ,
        NormalizationType_Symmetric = 3

    ctypedef struct CPrimitive "Primitive":
        CReal *ccbf
        CReal *coefficients
        CReal *coefficients0
        CReal  exponent
        CReal  exponent0

    ctypedef struct CShellDefinition "ShellDefinition":
        CInteger   angularmomentum_low
        CInteger   angularmomentum_high
        CInteger   cbfindex
        CInteger   nbasis
        CInteger   ncbf

    ctypedef struct CShell "Shell":
        CInteger          nbasisw
        CInteger          nprimitives
        CInteger          nstart
        CInteger          nstartw
        CRealArray2D     *c2s
        CRealArray2D     *s2c
        CPrimitive       *primitives
        CShellDefinition *type

    ctypedef struct CGaussianBasis "GaussianBasis":
        CBoolean           QNORMALIZEDPRIMITIVES
        CBoolean           QSPHERICAL
        CBoolean           QTOSPHERICAL
        CInteger           atomicNumber
        CInteger           maximum_angularmomentum
        CInteger           nbasis
        CInteger           nbasisw
        CInteger           nshells
        CGaussianBasisType basisType
        CNormalizationType normalizationType
        CRealArray2D      *c2o
        CRealArray2D      *o2c
        CShell            *shells

    cdef CGaussianBasis *GaussianBasis_Allocate              ( CInteger         nshells            )
    cdef void            GaussianBasis_AllocateShell         ( CGaussianBasis  *self               ,
                                                               CInteger         ishell             ,
                                                               CInteger         nprimitives        ,
                                                               CInteger         typeindex          )
    cdef CGaussianBasis *GaussianBasis_Clone                 ( CGaussianBasis  *self               )
    cdef void            GaussianBasis_Deallocate            ( CGaussianBasis **self               )
    cdef void            GaussianBasis_ScaleShellExponents   ( CGaussianBasis  *self               ,
                                                               CInteger         index              ,
                                                               CReal            zeta               ,
                                                               CStatus         *status             )
    cdef void            GaussianBasis_UnnormalizePrimitives ( CGaussianBasis  *self               )

cdef extern from "GaussianBasisNormalize.h":

    cdef CReal           GaussianBasis_Normalize             ( CGaussianBasis  *self               ,
                                                               CBoolean         checkNormalization ,
                                                               CStatus         *status             )

cdef extern from "GaussianBasisRotate.h":

    cdef void GaussianBasis_MakeLRotations                   ( CInteger         L                  ,
                                                               CRealArray2D    *R                  ,
                                                               CRealArray2D    *Tc                 ,
                                                               CStatus         *status             )
    cdef void GaussianBasis_MakeRotationMatrix               ( CGaussianBasis  *self               ,
                                                               CRealArray2D    *Tc                 ,
                                                               CBoolean         doC2O              ,
                                                               CRealArray2D    *T                  ,
                                                               CStatus         *status             )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasis:

    cdef CGaussianBasis *cObject
    cdef public object   isOwner

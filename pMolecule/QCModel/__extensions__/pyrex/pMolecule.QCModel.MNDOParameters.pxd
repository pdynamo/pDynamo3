from pCore.CPrimitiveTypes                         cimport CBoolean       , \
                                                           CFalse         , \
                                                           CInteger       , \
                                                           CReal          , \
                                                           CTrue
from pCore.Status                                  cimport CStatus        , \
                                                           CStatus_OK
from pMolecule.QCModel.GaussianBases.GaussianBasis cimport CGaussianBasis , \
                                                           GaussianBasis

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MNDOParameters.h":

    ctypedef struct CMNDOParameters "MNDOParameters":
        CBoolean     QDIATOMIC
        CBoolean    *QDIATOMICFLAGS
        CInteger     atomicNumber
        CInteger     iii
        CInteger     iiid
        CInteger     ir016
        CInteger     ir066
        CInteger     ir244
        CInteger     ir266
        CInteger     ir466
        CInteger     nam1pm3g
        CInteger     ndiatomic
        CInteger     norbitals
        CInteger     npddg
        CInteger     qnd
        CInteger     qnp
        CInteger     qns

        CInteger     nocteis
#        CCardinal16 *octeiindices
        CReal        hpp
        CReal       *octeivalues

        CReal        ad0
        CReal        alp0
        CReal        am0
        CReal        aq0
        CReal        betad0
        CReal        betap0
        CReal        betas0
        CReal        dd0
        CReal        eheat0
        CReal        eisol0
        CReal        f0sd0
        CReal        gphot0
        CReal        gpp0
        CReal        gp20
        CReal        gsp0
        CReal        gss0
        CReal        g2sd0
        CReal        hsp0
        CReal        pcore0
        CReal        qq0
        CReal        udd0
        CReal        upp0
        CReal        uss0
        CReal        zcore0
        CReal        zetad0
        CReal        zetap0
        CReal        zetas0
        CReal        zdn0
        CReal        zpn0
        CReal        zsn0

        CReal       *beta0
        CReal       *diatomica0
        CReal       *diatomicx0
        CReal       *fn10
        CReal       *fn20
        CReal       *fn30
        CReal       *pddgc0
        CReal       *pddge0
        CReal       *uspd0

        CReal        ad
        CReal        alp
        CReal        am
        CReal        aq
        CReal        betad
        CReal        betap
        CReal        betas
        CReal        dd
        CReal        eheat
        CReal        eisol
        CReal        f0sd
        CReal        gphot
        CReal        gpp
        CReal        gp2
        CReal        gsp
        CReal        gss
        CReal        g2sd
        CReal        hsp
        CReal        pcore
        CReal        qq
        CReal        udd
        CReal        upp
        CReal        uss
        CReal        zcore
        CReal        zetad
        CReal        zetap
        CReal        zetas
        CReal        zdn
        CReal        zpn
        CReal        zsn

        CReal       *beta
        CReal       *diatomica
        CReal       *diatomicx
        CReal       *fn1
        CReal       *fn2
        CReal       *fn3
        CReal       *pddgc
        CReal       *pddge
        CReal       *uspd

        CReal        ddp[6]
        CReal        po[9]

    cdef CMNDOParameters *MNDOParameters_Allocate               ( )
    cdef void             MNDOParameters_AllocateDiatomic       ( CMNDOParameters  *self , CInteger  n )
    cdef void             MNDOParameters_AllocateFN123          ( CMNDOParameters  *self , CInteger  n )
    cdef void             MNDOParameters_AllocatePDDG           ( CMNDOParameters  *self , CInteger  n )
    cdef void             MNDOParameters_CalculateOneCenterTEIs ( CMNDOParameters  *self )
    cdef CMNDOParameters *MNDOParameters_Clone                  ( CMNDOParameters  *self )
    cdef void             MNDOParameters_Deallocate             ( CMNDOParameters **self )
    cdef void             MNDOParameters_DetermineNormalization ( CMNDOParameters  *self , CGaussianBasis *basis, CStatus *status )
    cdef void             MNDOParameters_FillBetaUspd           ( CMNDOParameters  *self )
    cdef void             MNDOParameters_ToAtomicUnits          ( CMNDOParameters  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOParameters:

    cdef CMNDOParameters *cObject
    cdef public object    isOwner

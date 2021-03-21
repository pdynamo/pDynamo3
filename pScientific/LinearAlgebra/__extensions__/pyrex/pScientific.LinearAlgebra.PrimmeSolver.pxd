
#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "cprimme.h":

    ctypedef enum CPrimmePresetMethod "primme_preset_method":
        DYNAMIC                  ,
        DEFAULT_MIN_TIME         ,
        DEFAULT_MIN_MATVECS      ,
        Arnoldi                  ,
        GD                       ,
        GD_plusK                 ,
        GD_Olsen_plusK           ,
        JD_Olsen_plusK           ,
        RQI                      ,
        JDQR                     ,
        JDQMR                    ,
        JDQMR_ETol               ,
        SUBSPACE_ITERATION       ,
        LOBPCG_OrthoBasis        ,
        LOBPCG_OrthoBasis_Window

    ctypedef struct CPrimmeCorrectionParams "correction_params":
        int precondition

    ctypedef struct CPrimmeStats "primme_stats":
        int    numOuterIterations
        int    numRestarts
        int    numMatvecs
        int    numPreconds
        double elapsedTime

    ctypedef struct CPrimmeParams "primme_params":
        void ( *matrixMatvec        ) ( void *x,  void *y, int *blockSize, CPrimmeParams *primme )
        void ( *applyPreconditioner ) ( void *x,  void *y, int *blockSize, CPrimmeParams *primme )
        void ( *massMatrixMatvec    ) ( void *x,  void *y, int *blockSize, CPrimmeParams *primme )
        int                     n
        int                     numEvals
        int                     dynamicMethodSwitch
        int                     locking
        int                     initSize
        int                     numOrthoConst
        int                     maxBasisSize
        int                     minRestartSize
        int                     maxBlockSize
        int                     maxMatvecs
        int                     maxOuterIterations
        int                     intWorkSize
        long int                realWorkSize
        int                    *intWork
        void                   *realWork
        double                  aNorm
        double                  eps
        int                     printLevel
        void                   *matrix
        void                   *preconditioner
        CPrimmeCorrectionParams correctionParams
        CPrimmeStats            stats

    cdef int  dprimme               ( double *evals, double *evecs, double *resNorms, CPrimmeParams *primme )
    cdef void primme_display_params ( CPrimmeParams  primme )
    cdef void primme_initialize     ( CPrimmeParams *primme )
    cdef int  primme_set_method     ( CPrimmePresetMethod method, CPrimmeParams *params )

from pCore.CPrimitiveTypes          cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Status                   cimport CStatus, CStatus_OK
from pScientific.Arrays.RealArray1D cimport CRealArray1D, RealArray1D
from pScientific.Arrays.RealArray2D cimport CRealArray2D, RealArray2D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Correlation.h":

    cdef void Correlation_MakeDotProduct ( CRealArray2D *x            ,
                                           CRealArray2D *y            ,
                                           CBoolean      useFFT       ,
                                           CBoolean      normalize    ,
                                           CBoolean      removeMean   ,
                                           CInteger      tCorrelation ,
                                           CReal        *tolerance    ,
                                           CRealArray1D *f            ,
                                           CStatus      *status       )
    cdef void Correlation_MakeSimple     ( CRealArray1D *x            ,
                                           CRealArray1D *y            ,
                                           CBoolean      useFFT       ,
                                           CBoolean      normalize    ,
                                           CBoolean      removeMean   ,
                                           CInteger      tCorrelation ,
                                           CReal        *tolerance    ,
                                           CRealArray1D *f            ,
                                           CStatus      *status       )

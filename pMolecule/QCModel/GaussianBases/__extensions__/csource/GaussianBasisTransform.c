/*==================================================================================================================================
! . Functions for the transformation of Gaussian integrals.
! . Straightforward implementations with little optimization and no checking!
!=================================================================================================================================*/

# include <stdio.h>

# include "Array_Macros.h"
# include "GaussianBasisTransform.h"

/* . Options - by default use unoptimized versions. */
/*# define _UseCBLAS_*/         /* . Use CBLAS directly. */
/*# define _UseRealArray_*/     /* . Use pDynamo real arrays. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Static procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _UseRealArray_
static void Transform1 ( const Integer  dC       ,
                         const Integer  dS       ,
                         const Real    *T        ,
                         const Real    *In       ,
                               Real    *Out      ) ;
# endif
static void Transform2 ( const Integer  dI       ,
                         const Integer  dC       ,
                         const Integer  dS       ,
                         const Real    *T        ,
                         const Real    *In       ,
                         const Integer  iStrideI ,
                         const Integer  iStrideJ ,
                               Real    *Out      ,
                         const Integer  oStrideI ,
                         const Integer  oStrideJ ) ;
static void Transform3 ( const Integer  dI       ,
                         const Integer  dJ       ,
                         const Integer  dC       ,
                         const Integer  dS       ,
                         const Real    *T        ,
                         const Real    *In       ,
                         const Integer  iStrideI ,
                         const Integer  iStrideJ ,
                         const Integer  iStrideK ,
                               Real    *Out      ,
                         const Integer  oStrideI ,
                         const Integer  oStrideJ ,
                         const Integer  oStrideK ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . One transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisTransform1 ( const RealArray2D *tI     ,
                                     Real       **values ,
                                     Real       **work   )
{
    /* . Initialization. */
    if ( tI != NULL )
    {
        auto Real *In, *Out ;
# ifdef _UseRealArray_
        auto RealArray1D viewIn, viewOut ;
# else
        auto Integer dC, dS ;
        auto Real   *T ;
# endif
        In  = (*values) ;
        Out = (*work  ) ;
# ifdef _UseRealArray_
        RealArray1D_ViewOfRaw ( &viewIn , 0, tI->extent0, 1, In  ) ;
        RealArray1D_ViewOfRaw ( &viewOut, 0, tI->extent1, 1, Out ) ;
        RealArray2D_VectorMultiply ( True, 1.0e+00, tI, &viewIn, 0.0e+00, &viewOut, NULL ) ;
# else
        dC  = tI->extent0 ; dS = tI->extent1 ; T = Array_DataPointer ( tI ) ;
        Transform1 ( dC, dS, T, In, Out ) ;
# endif
        (*values) = Out ;
        (*work  ) = In  ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . One transformation for a matrix.
!   Assumed i x m unless transpose is true.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisTransform1M (       Integer      dI        ,
                                      Integer      dM        ,
                                const RealArray2D *tI        ,
                                const Boolean      transpose ,
                                      Real       **values    ,
                                      Real       **work      )
{
    if ( tI != NULL )
    {
        auto Real *In, *Out ;
# ifdef _UseRealArray_
        auto RealArray2D viewIn, viewOut ;
# else
        auto Integer dC, dS, iStrideI, iStrideM, oStrideI, oStrideM ;
        auto Real   *T ;
# endif
        In  = (*values) ;
        Out = (*work  ) ;
# ifdef _UseRealArray_
        if ( transpose )
        {
            RealArray2D_ViewOfRaw ( &viewIn , 0, dM, dI         , dI         , 1, In  ) ;
            RealArray2D_ViewOfRaw ( &viewOut, 0, dM, tI->extent1, tI->extent1, 1, Out ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, &viewIn, tI, 0.0e+00, &viewOut, NULL ) ;
        }
        else
        {
            RealArray2D_ViewOfRaw ( &viewIn , 0, dI         , dM, dM, 1, In  ) ;
            RealArray2D_ViewOfRaw ( &viewOut, 0, tI->extent1, dM, dM, 1, Out ) ;
            RealArray2D_MatrixMultiply ( True , False, 1.0e+00, tI, &viewIn, 0.0e+00, &viewOut, NULL ) ;
        }
# else
        dC  = tI->extent0 ; dS = tI->extent1 ; T = Array_DataPointer ( tI ) ;
        if ( transpose ) { iStrideI = 1 ; iStrideM = dI * iStrideI ; oStrideI = 1 ; oStrideM = dS * oStrideI ; }
        else             { iStrideM = 1 ; iStrideI = dM * iStrideM ; oStrideM = 1 ; oStrideI = dM * oStrideM ; }
        Transform2 ( dM, dC, dS, T, In, iStrideM, iStrideI, Out, oStrideM, oStrideI ) ;
# endif
        (*values) = Out ;
        (*work  ) = In  ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Two transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _strideJ_ 1
void GaussianBasisTransform2 (       Integer      dI     ,
                                     Integer      dJ     ,
                               const RealArray2D *tI     ,
                               const RealArray2D *tJ     ,
                                     Real       **values ,
                                     Real       **work   )
{
    Integer  dC, dS ;
    Real    *In, *Swap, *Out, *T ;
    /* . Initialization. */
    In  = (*values) ;
    Out = (*work  ) ;
    /* . I transformation. */
    if ( tI != NULL )
    {
        dC = tI->extent0 ; dS = tI->extent1 ; T = Array_DataPointer ( tI ) ;
        Transform2 ( dJ, dC, dS, T, In, _strideJ_, dJ * _strideJ_, Out, _strideJ_, dJ * _strideJ_ ) ;
        Swap = In ; In = Out ; Out = Swap ; dI = dS ;
    }
    /* . J transformation. */
    if ( tJ != NULL )
    {
        dC = tJ->extent0 ; dS = tJ->extent1 ; T = Array_DataPointer ( tJ ) ;
        Transform2 ( dI, dC, dS, T, In, dJ * _strideJ_, _strideJ_, Out, dS * _strideJ_, _strideJ_ ) ;
        Swap = In ; In = Out ; Out = Swap ;
    }
    /* . Finish up. */
    (*values) = In  ;
    (*work  ) = Out ;
}
# undef _strideJ_

/*----------------------------------------------------------------------------------------------------------------------------------
! . Three transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _strideK_ 1
void GaussianBasisTransform3 (       Integer      dI     ,
                                     Integer      dJ     ,
                                     Integer      dK     ,
                               const RealArray2D *tI     ,
                               const RealArray2D *tJ     ,
                               const RealArray2D *tK     ,
                                     Real       **values ,
                                     Real       **work   )
{
    Integer  d2, dC, dS, iStrideI, iStrideJ, oStrideI, oStrideJ ;
    Real    *In, *Swap, *Out, *T ;
    /* . Initialization. */
    In  = (*values) ;
    Out = (*work  ) ;
    iStrideJ = dK * _strideK_ ; iStrideI = dJ * iStrideJ ;
/*printf ( "T3>  Starting: dI %d, dJ %d, dK %d, iStrideJ %d, iStrideI %d\n", dI, dJ, dK, iStrideJ, iStrideI ) ; fflush ( stdout ) ;*/
    /* . I transformation. */
    if ( tI != NULL )
    {
        dC = tI->extent0 ; dS = tI->extent1 ; d2 = dJ * dK ; T = Array_DataPointer ( tI ) ;
        Transform2 ( d2, dC, dS, T, In, _strideK_, d2 * _strideK_, Out, _strideK_, d2 * _strideK_ ) ;
        Swap = In ; In = Out ; Out = Swap ; dI = dS ;
/*printf ( "T3> Done I: dI %d, dJ %d, dK %d, dC %d, dS %d, d2 %d, iStrideJ %d, iStrideI %d\n", dI, dJ, dK, dC, dS, d2, iStrideJ, iStrideI ) ; fflush ( stdout ) ;*/
/*printf ( "TI:\n" ) ; RealArray2D_Print ( tI ) ; fflush ( stdout ) ;*/
    }
    /* . J transformation. */
    if ( tJ != NULL )
    {
        dC       = tJ->extent0    ; dS       = tJ->extent1   ; T = Array_DataPointer ( tJ ) ;
        oStrideJ = dK * _strideK_ ; oStrideI = dS * oStrideJ ;
        Transform3 ( dI, dK, dC, dS, T, In , iStrideI, _strideK_, iStrideJ ,
                                        Out, oStrideI, _strideK_, oStrideJ ) ;
        Swap     = In       ; In       = Out      ; Out      = Swap     ; dJ = dS ;
        iStrideI = oStrideI ; iStrideJ = oStrideJ ;
/*printf ( "T3> Done J: dI %d, dJ %d, dK %d, dC %d, dS %d, iStrideJ %d, iStrideI %d\n", dI, dJ, dK, dC, dS, iStrideJ, iStrideI ) ; fflush ( stdout ) ;*/
/*printf ( "TJ:\n" ) ; RealArray2D_Print ( tJ ) ; fflush ( stdout ) ;*/
    }
    /* . K transformation. */
    if ( tK != NULL )
    {
/*
auto Real maxAbs ;
auto RealArray2D *test = NULL, viewIn, viewOut ;
*/
        dC = tK->extent0 ; dS = tK->extent1 ; d2 = dI * dJ ; T = Array_DataPointer ( tK ) ;
/*
if ( dC != dK ) printf ( "T3> Array mismatch T and In = %d  %d\n", dC, dK ) ;
test = RealArray2D_AllocateWithExtents ( d2, dS, NULL ) ;
RealArray2D_ViewOfRaw ( &viewIn , 0, d2, dC, dC, 1, In  ) ;
RealArray2D_ViewOfRaw ( &viewOut, 0, d2, dS, dS, 1, Out ) ;
RealArray2D_MatrixMultiply ( False, False, 1.0e+00, &viewIn, tK, 0.0e+00, test, NULL ) ;
*/
        Transform2 ( d2, dC, dS, T, In, iStrideJ, _strideK_, Out, dS * _strideK_, _strideK_ ) ;
/*
RealArray2D_Add ( test, -1.0e+00, &viewOut, NULL ) ;
maxAbs = RealArray2D_AbsoluteMaximum ( test ) ;
*/
        Swap = In ; In = Out ; Out = Swap ;
/*printf ( "T3> Done K: dI %d, dJ %d, dK %d, dC %d, dS %d, d2 %d, iStrideJ %d, iStrideI %d\n", dI, dJ, dK, dC, dS, d2, iStrideJ, iStrideI ) ; fflush ( stdout ) ;*/
/*printf ( "TK:\n" ) ; RealArray2D_Print ( tK ) ; fflush ( stdout ) ;*/
/*
printf ( "T3> Absolute maximum = %20.10f\n", maxAbs ) ;
if ( maxAbs > 1.0e-10 )
{
printf ( "Starting:\n"    ) ; RealArray2D_Print ( &viewIn  ) ; fflush ( stdout ) ;
printf ( "Actual:\n"      ) ; RealArray2D_Print ( &viewOut ) ; fflush ( stdout ) ;
printf ( "Differences:\n" ) ; RealArray2D_Print ( test     ) ; fflush ( stdout ) ;
RealArray2D_Add ( test, 1.0e+00, &viewOut, NULL ) ;
printf ( "Reference:\n"   ) ; RealArray2D_Print ( test     ) ; fflush ( stdout ) ;
}
RealArray2D_Deallocate ( &test ) ;*/
    }
    /* . Finish up. */
    (*values) = In  ;
    (*work  ) = Out ;
}
# undef _strideK_

/*----------------------------------------------------------------------------------------------------------------------------------
! . Four transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _strideL_ 1
void GaussianBasisTransform4 (       Integer      dI     ,
                                     Integer      dJ     ,
                                     Integer      dK     ,
                                     Integer      dL     ,
                               const RealArray2D *tI     ,
                               const RealArray2D *tJ     ,
                               const RealArray2D *tK     ,
                               const RealArray2D *tL     ,
                                     Real       **values ,
                                     Real       **work   )
{
    Integer  d2, d3, dC, dS, iStrideI, iStrideJ, iStrideK, oStrideI, oStrideJ, oStrideK ;
    Real    *In, *Swap, *Out, *T ;
    /* . Initialization. */
    In  = (*values) ;
    Out = (*work  ) ;
    iStrideK = dL * _strideL_ ; iStrideJ = dK * iStrideK ; iStrideI = dJ * iStrideJ ;
    /* . I transformation. */
    if ( tI != NULL )
    {
        dC = tI->extent0 ; dS = tI->extent1 ; d3 = dJ * dK * dL ; T = Array_DataPointer ( tI ) ;
        Transform2 ( d3, dC, dS, T, In, _strideL_, d3 * _strideL_, Out, _strideL_, d3 * _strideL_ ) ;
        Swap = In ; In = Out ; Out = Swap ; dI = dS ;
    }
    /* . J transformation. */
    if ( tJ != NULL )
    {
        dC       = tJ->extent0    ; dS       = tJ->extent1   ; d2 = dK * dL ; T = Array_DataPointer ( tJ ) ;
        oStrideK = dL * _strideL_ ; oStrideJ = dK * oStrideK ; oStrideI = dS * oStrideJ ;
        Transform3 ( dI, d2, dC, dS, T, In , iStrideI, _strideL_, iStrideJ ,
                                        Out, oStrideI, _strideL_, oStrideJ ) ;
        Swap     = In       ; In       = Out      ; Out      = Swap     ; dJ = dS ;
        iStrideI = oStrideI ; iStrideJ = oStrideJ ; iStrideK = oStrideK ;
    }
    /* . K transformation. */
    if ( tK != NULL )
    {
        dC       = tK->extent0    ; dS       = tK->extent1   ; d2 = dI * dJ ; T = Array_DataPointer ( tK ) ;
        oStrideK = dL * _strideL_ ; oStrideJ = dS * oStrideK ; oStrideI = dJ * oStrideJ ;
        Transform3 ( d2, dL, dC, dS, T, In , iStrideJ, _strideL_, iStrideK ,
                                        Out, oStrideJ, _strideL_, oStrideK ) ;
        Swap     = In       ; In       = Out      ; Out      = Swap     ; dK = dS ;
        iStrideI = oStrideI ; iStrideJ = oStrideJ ; iStrideK = oStrideK ;
    }
    /* . L transformation. */
    if ( tL != NULL )
    {
        dC       = tL->extent0    ; dS       = tL->extent1   ; d3 = dI * dJ * dK ; T = Array_DataPointer ( tL ) ;
        oStrideK = dS * _strideL_ ; oStrideJ = dK * oStrideK ; oStrideI = dJ * oStrideJ ;
        Transform2 ( d3, dC, dS, T, In, iStrideK, _strideL_, Out, oStrideK, _strideL_ ) ;
        Swap = In ; In = Out ; Out = Swap ;
    }
    /* . Finish up. */
    (*values) = In  ;
    (*work  ) = Out ;
}
# undef _strideL_

/*----------------------------------------------------------------------------------------------------------------------------------
! . Subsidiary transformations.
! . All arrays are compact including T (strideC = dS, strideS = 1).
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _UseRealArray_
/* . 1-transform. */
static void Transform1 ( const Integer  dC  ,
                         const Integer  dS  ,
                         const Real    *T   ,
                         const Real    *In  ,
                               Real    *Out )
{
    Integer c, s ;
    Real    v ;
    for ( s = 0 ; s < dS ; s++ )
    {
        for ( c = 0, v = 0.0e+00 ; c < dC ; c++ ) { v += ( T[c*dS+s] * In[c] ) ; }
        Out[s] = v ;
    }
}
# endif

/* . 2-transform. */
static void Transform2 ( const Integer  dI       ,
                         const Integer  dC       ,
                         const Integer  dS       ,
                         const Real    *T        ,
                         const Real    *In       ,
                         const Integer  iStrideI ,
                         const Integer  iStrideJ ,
                               Real    *Out      ,
                         const Integer  oStrideI ,
                         const Integer  oStrideJ )
{
    Integer c, i, iI, oI, s ;
    Real    v ;
    for ( i = 0 ; i < dI ; i++ )
    {
        iI = i * iStrideI ;
        oI = i * oStrideI ;
        for ( s = 0 ; s < dS ; s++ )
        {
            for ( c = 0, v = 0.0e+00 ; c < dC ; c++ ) { v += ( T[c*dS+s] * In[iI+c*iStrideJ] ) ; }
            Out[oI+s*oStrideJ] = v ;
        }
    }
}

/* . 3-transform. */
static void Transform3 ( const Integer  dI       ,
                         const Integer  dJ       ,
                         const Integer  dC       ,
                         const Integer  dS       ,
                         const Real    *T        ,
                         const Real    *In       ,
                         const Integer  iStrideI ,
                         const Integer  iStrideJ ,
                         const Integer  iStrideK ,
                               Real    *Out      ,
                         const Integer  oStrideI ,
                         const Integer  oStrideJ ,
                         const Integer  oStrideK )
{
    Integer c, i, iI, iIJ, j, oI, oIJ, s ;
    Real    v ;
    for ( i = 0 ; i < dI ; i++ )
    {
        iI = i * iStrideI ;
        oI = i * oStrideI ;
        for ( j = 0 ; j < dJ ; j++ )
        {
            iIJ = iI + j * iStrideJ ;
            oIJ = oI + j * oStrideJ ;
            for ( s = 0 ; s < dS ; s++ )
            {
                for ( c = 0, v = 0.0e+00 ; c < dC ; c++ ) { v += ( T[c*dS+s] * In[iIJ+c*iStrideK] ) ; }
                Out[oIJ+s*oStrideK] = v ;
            }
        }
    }
}

/*==================================================================================================================================
! . Functions for Fock construction.
!=================================================================================================================================*/

# include "stdio.h"

# include "Boolean.h"
# include "DenseLinearEquationSolvers.h"
# include "FockConstruction.h"
# include "Integer.h"

/* . Formulae:

     Pa = ( Pt + Ps ) / 2 ; Pb = ( Pt - Ps ) / 2
     Fa =   Ft + Fs       ; Fb =   Ft - Fs
     Pt =   Pa + Pb       ; Ps =   Pa - Pb
     Ft = ( Fa + Fb ) / 2 ; Fs = ( Fa - Fb ) / 2

*/

/* . Macros. */
# define BFINDEX(i) ( i * ( i + 1 ) ) / 2

/* . Options. */
/*# define _NOFITCONSTRAINTS*/
/*# define _USEINVERSEFITMATRIX*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the coefficients from the fit matrix and the fit integrals and density.
! . The b-vector is also returned.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Fock_MakeCoefficientsFromFitIntegrals (       SymmetricMatrix *fitMatrix       ,
                                                   BlockStorage    *fitIntegrals    ,
                                                   SymmetricMatrix *dTotal          ,
                                             const Real             totalCharge     ,
                                                   RealArray1D     *fitCoefficients ,
                                                   RealArray1D     *bVector         ,
                                                   Status          *status          )
{
    if ( ( fitIntegrals    != NULL ) &&
         ( fitMatrix       != NULL ) &&
         ( dTotal          != NULL ) &&
         ( fitCoefficients != NULL ) &&
         ( bVector         != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block   *block ;
        auto Integer  i ;
        /* . Scale diagonal elements of the density by 1/2. */
        SymmetricMatrix_ScaleDiagonal ( dTotal, 0.5e+00 ) ;
        /* . Determine b. */
        RealArray1D_Set ( bVector, 0.0e+00 ) ;
        List_Iterate_Initialize ( fitIntegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( fitIntegrals ) ) != NULL )
        {
            for ( i = 0 ; i < block->count ; i++ ) bVector->data[block->indices16[i]] += dTotal->data[block->indices32[i]] * block->data[i] ;
        }
        RealArray1D_Scale ( bVector, 2.0e+00 ) ;
# ifndef _NOFITCONSTRAINTS
        Array1D_Item ( bVector, View1D_Extent ( fitCoefficients ) - 1 ) = totalCharge ;
# endif
        /* . Find the fit coefficients. */
# ifdef _USEINVERSEFITMATRIX
        SymmetricMatrix_VectorMultiply ( fitMatrix, bVector, fitCoefficients, status ) ;
# else
        SymmetricMatrix_LinearEquationsSolve ( fitMatrix, bVector, fitCoefficients, status ) ;
# endif
        /* . Scale diagonal elements of the density by 2. */
        SymmetricMatrix_ScaleDiagonal ( dTotal, 2.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the Fock matrix from the fit integrals and the fit vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Fock_MakeFockFromFitIntegrals (       BlockStorage    *fitIntegrals ,
                                     const RealArray1D     *fitVector    ,
                                           SymmetricMatrix *fTotal       ,
                                           Status          *status       )
{
    if ( ( fitIntegrals != NULL ) &&
         ( fitVector    != NULL ) &&
         ( fTotal       != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block   *block ;
        auto Integer  i ;
        /* . Construct Fock matrix. */
        List_Iterate_Initialize ( fitIntegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( fitIntegrals ) ) != NULL )
        {
            for ( i = 0 ; i < block->count ; i++ )
            {
                fTotal->data[block->indices32[i]] += fitVector->data[block->indices16[i]] * block->data[i] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the fit-integral parts of the Fock matrices.
! . The fit energy and coefficients are also computed.
! . Fit-integrals computed using the Coulomb operator.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Fock_MakeFromFitIntegralsCoulomb (       BlockStorage    *fitIntegrals    ,
                                              SymmetricMatrix *fitMatrix       ,
                                        const Real             totalCharge     ,
                                              RealArray1D     *fitCoefficients ,
                                              SymmetricMatrix *dTotal          ,
                                              SymmetricMatrix *fTotal          ,
                                              Status          *status          )
{
    Real eFit = 0.0e+00 ;
    if ( ( fitIntegrals    != NULL ) &&
         ( fitMatrix       != NULL ) &&
         ( fitCoefficients != NULL ) &&
         ( dTotal          != NULL ) &&
         ( fTotal          != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n    = View1D_Extent ( fitCoefficients ) ;
        auto RealArray1D *work = RealArray1D_AllocateWithExtent ( n, status ) ;
        if ( work != NULL )
        {
            /* . Fit coefficients. */
            Fock_MakeCoefficientsFromFitIntegrals ( fitMatrix       ,
                                                    fitIntegrals    ,
                                                    dTotal          ,
                                                    totalCharge     ,
                                                    fitCoefficients ,
                                                    work            ,
                                                    status          ) ;
            /* . Fit energy. */
            eFit = 0.5e+00 * RealArray1D_Dot ( fitCoefficients, work, status ) ;
            /* . Fit Fock matrix. */
            Fock_MakeFockFromFitIntegrals ( fitIntegrals    ,
                                            fitCoefficients ,
                                            fTotal          ,
                                            status          ) ;
            /* . Finish up. */
            RealArray1D_Deallocate ( &work ) ;
        }
    }
    return eFit ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the fit-integral parts of the Fock matrices.
! . The fit energy and coefficients are also computed.
! . Fit-integrals computed using a non-Coulomb operator.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Fock_MakeFromFitIntegralsNonCoulomb (       BlockStorage    *fitIntegrals     ,
                                                 SymmetricMatrix *fitMatrix        ,
                                                 SymmetricMatrix *fitCoulombMatrix ,
                                           const Real             totalCharge      ,
                                                 RealArray1D     *fitCoefficients  ,
                                                 RealArray1D     *fitVectorD       ,
                                                 SymmetricMatrix *dTotal           ,
                                                 SymmetricMatrix *fTotal           ,
                                                 Status          *status           )
{
    Real eFit = 0.0e+00 ;
    if ( ( fitIntegrals     != NULL ) &&
         ( fitMatrix        != NULL ) &&
         ( fitCoulombMatrix != NULL ) &&
         ( fitCoefficients  != NULL ) &&
         ( fitVectorD       != NULL ) &&
         ( dTotal           != NULL ) &&
         ( fTotal           != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n    = View1D_Extent ( fitCoefficients ) ;
        auto RealArray1D *work = RealArray1D_AllocateWithExtent ( n, status ) ;
        if ( work != NULL )
        {
            /* . Fit coefficients. */
            Fock_MakeCoefficientsFromFitIntegrals ( fitMatrix       ,
                                                    fitIntegrals    ,
                                                    dTotal          ,
                                                    totalCharge     ,
                                                    fitCoefficients ,
                                                    work            ,
                                                    status          ) ;
            /* . Compute T = Mc * A. */
            SymmetricMatrix_VectorMultiply ( fitCoulombMatrix, fitCoefficients, work, status ) ;
            /* . Fit energy. */
            eFit = 0.5e+00 * RealArray1D_Dot ( fitCoefficients, work, status ) ;
            /* . Solve for the fit D-vector. */
# ifdef _USEINVERSEFITMATRIX
            SymmetricMatrix_VectorMultiply ( fitMatrix, work, fitVectorD, status ) ;
# else
            SymmetricMatrix_LinearEquationsSolve ( fitMatrix, work, fitVectorD, status ) ;
# endif
            /* . Fit Fock matrix. */
            Fock_MakeFockFromFitIntegrals ( fitIntegrals ,
                                            fitVectorD   ,
                                            fTotal       ,
                                            status       ) ;
            /* . Finish up. */
            RealArray1D_Deallocate ( &work ) ;
        }
    }
    return eFit ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the two-electron Coulomb and exchange parts of the Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Fock_MakeFromTEIs (       BlockStorage    *twoElectronIntegrals ,
                         const SymmetricMatrix *dTotal               ,
                         const SymmetricMatrix *dSpin                ,
                         const Real             exchangeScaling      ,
                               SymmetricMatrix *fTotal               ,
                               SymmetricMatrix *fSpin                )
{
    Real eTEI = 0.0e+00 ;
    if  ( ( twoElectronIntegrals != NULL ) && ( dTotal != NULL ) && ( fTotal != NULL ) )
    {
        auto Block   *block ;
        auto Boolean  doSpin = ( ( dSpin != NULL ) && ( fSpin != NULL ) ) ;
        auto Integer  i, i1, i2, i3, i4, n, nIJ, nIK, nIL, nJK, nJL, nKL, t ;
        auto Real     value ;
        SymmetricMatrix_Set ( fTotal, 0.0e+00 ) ;
        SymmetricMatrix_Set ( fSpin , 0.0e+00 ) ;
        List_Iterate_Initialize ( twoElectronIntegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( twoElectronIntegrals ) ) != NULL )
        {
            for ( i = 0, n = 0 ; i < block->count ; i++, n += 4 )
            {
                i1    = block->indices16[n  ] ;
                i2    = block->indices16[n+1] ;
                i3    = block->indices16[n+2] ;
                i4    = block->indices16[n+3] ;
                value = block->data[i] ;
	        if ( i1 < i2 ) { t = i1 ; i1 = i2 ; i2 = t ; }
                if ( i3 < i4 ) { t = i3 ; i3 = i4 ; i4 = t ; }
                if ( ( i1 < i3 ) || ( ( i1 == i3 ) && ( i2 < i4  ) ) ) { t = i1 ; i1 = i3 ; i3 = t ; t = i2 ; i2 = i4 ; i4 = t ; }
	        if ( i1 == i2 ) value *= 0.5e+00 ;
	        if ( i3 == i4 ) value *= 0.5e+00 ;
                if ( ( i1 == i3 ) && ( i2 == i4 ) ) value *= 0.5e+00 ;
                nIJ = BFINDEX ( i1 ) + i2 ;
                nKL = BFINDEX ( i3 ) + i4 ;
                nIK = BFINDEX ( i1 ) + i3 ;
                nIL = BFINDEX ( i1 ) + i4 ;
                if ( i2 > i3 ) nJK = BFINDEX ( i2 ) + i3 ;
                else           nJK = BFINDEX ( i3 ) + i2 ;
                if ( i2 > i4 ) nJL = BFINDEX ( i2 ) + i4 ;
                else           nJL = BFINDEX ( i4 ) + i2 ;
                /* . Coulomb. */
                fTotal->data[nIJ] += 4.0e+00 * value * dTotal->data[nKL] ;
                fTotal->data[nKL] += 4.0e+00 * value * dTotal->data[nIJ] ;
                /* . Exchange. */
                value *= exchangeScaling ;
                fTotal->data[nIK] -= value * dTotal->data[nJL] ;
                fTotal->data[nIL] -= value * dTotal->data[nJK] ;
                fTotal->data[nJK] -= value * dTotal->data[nIL] ;
                fTotal->data[nJL] -= value * dTotal->data[nIK] ;
                if ( doSpin )
                {
                    fSpin->data[nIK] -= value * dSpin->data[nJL] ;
                    fSpin->data[nIL] -= value * dSpin->data[nJK] ;
                    fSpin->data[nJK] -= value * dSpin->data[nIL] ;
                    fSpin->data[nJL] -= value * dSpin->data[nIK] ;
                }
            }
        }
        SymmetricMatrix_ScaleOffDiagonal ( fTotal, 0.5e+00 ) ;
        SymmetricMatrix_ScaleOffDiagonal ( fSpin , 0.5e+00 ) ;
        eTEI = 0.5e+00 * SymmetricMatrix_TraceOfProduct ( dTotal, fTotal, NULL ) ;
        if ( doSpin ) eTEI += 0.5 * SymmetricMatrix_TraceOfProduct ( dSpin, fSpin, NULL ) ;
    }
    return eTEI ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the Coulomb two-electron part of the Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Fock_MakeFromTEIsCoulomb (       BlockStorage    *twoElectronIntegrals ,
                                const SymmetricMatrix *dTotal               ,
                                      SymmetricMatrix *fTotal               )
{
    Real eTEI = 0.0e+00 ;
    if  ( ( twoElectronIntegrals != NULL ) && ( dTotal != NULL ) && ( fTotal != NULL ) )
    {
        auto Block   *block ;
        auto Integer  i, i1, i2, i3, i4, n, nIJ, nKL, t ;
        auto Real     value ;
        SymmetricMatrix_Set ( fTotal, 0.0e+00 ) ;
        List_Iterate_Initialize ( twoElectronIntegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( twoElectronIntegrals ) ) != NULL )
        {
            for ( i = 0, n = 0 ; i < block->count ; i++, n += 4 )
            {
                i1    = block->indices16[n  ] ;
                i2    = block->indices16[n+1] ;
                i3    = block->indices16[n+2] ;
                i4    = block->indices16[n+3] ;
                value = 4.0e+00 * block->data[i] ;
	        if ( i1 < i2 ) { t = i1 ; i1 = i2 ; i2 = t ; }
                if ( i3 < i4 ) { t = i3 ; i3 = i4 ; i4 = t ; }
                if ( ( i1 < i3 ) || ( ( i1 == i3 ) && ( i2 < i4  ) ) ) { t = i1 ; i1 = i3 ; i3 = t ; t = i2 ; i2 = i4 ; i4 = t ; }
	        if ( i1 == i2 ) value *= 0.5e+00 ;
	        if ( i3 == i4 ) value *= 0.5e+00 ;
                if ( ( i1 == i3 ) && ( i2 == i4 ) ) value *= 0.5e+00 ;
                nIJ = BFINDEX ( i1 ) + i2 ;
                nKL = BFINDEX ( i3 ) + i4 ;
                fTotal->data[nIJ] += value * dTotal->data[nKL] ;
                fTotal->data[nKL] += value * dTotal->data[nIJ] ;
            }
        }
        SymmetricMatrix_ScaleOffDiagonal ( fTotal, 0.5e+00 ) ;
        eTEI = 0.5e+00 * SymmetricMatrix_TraceOfProduct ( dTotal, fTotal, NULL ) ;
    }
    return eTEI ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the two-electron exchange part of the Fock matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Fock_MakeFromTEIsExchange (       BlockStorage    *twoElectronIntegrals ,
                                 const SymmetricMatrix *dTotal               ,
                                 const SymmetricMatrix *dSpin                ,
                                 const Real             exchangeScaling      ,
                                       SymmetricMatrix *fTotal               ,
                                       SymmetricMatrix *fSpin                )
{
    Real eTEI = 0.0e+00 ;
    if  ( ( twoElectronIntegrals != NULL ) && ( dTotal != NULL ) && ( fTotal != NULL ) )
    {
        auto Block   *block ;
        auto Boolean  doSpin = ( ( dSpin != NULL ) && ( fSpin != NULL ) ) ;
        auto Integer  i, i1, i2, i3, i4, n, nIK, nIL, nJK, nJL, t ;
        auto Real     value ;
        SymmetricMatrix_Set ( fTotal, 0.0e+00 ) ;
        SymmetricMatrix_Set ( fSpin , 0.0e+00 ) ;
        List_Iterate_Initialize ( twoElectronIntegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( twoElectronIntegrals ) ) != NULL )
        {
            for ( i = 0, n = 0 ; i < block->count ; i++, n += 4 )
            {
                i1    = block->indices16[n  ] ;
                i2    = block->indices16[n+1] ;
                i3    = block->indices16[n+2] ;
                i4    = block->indices16[n+3] ;
                value = exchangeScaling * block->data[i] ;
	        if ( i1 < i2 ) { t = i1 ; i1 = i2 ; i2 = t ; }
                if ( i3 < i4 ) { t = i3 ; i3 = i4 ; i4 = t ; }
                if ( ( i1 < i3 ) || ( ( i1 == i3 ) && ( i2 < i4  ) ) ) { t = i1 ; i1 = i3 ; i3 = t ; t = i2 ; i2 = i4 ; i4 = t ; }
	        if ( i1 == i2 ) value *= 0.5e+00 ;
	        if ( i3 == i4 ) value *= 0.5e+00 ;
                if ( ( i1 == i3 ) && ( i2 == i4 ) ) value *= 0.5e+00 ;
                nIK = BFINDEX ( i1 ) + i3 ;
                nIL = BFINDEX ( i1 ) + i4 ;
                if ( i2 > i3 ) nJK = BFINDEX ( i2 ) + i3 ;
                else           nJK = BFINDEX ( i3 ) + i2 ;
                if ( i2 > i4 ) nJL = BFINDEX ( i2 ) + i4 ;
                else           nJL = BFINDEX ( i4 ) + i2 ;
                fTotal->data[nIK] -= value * dTotal->data[nJL] ;
                fTotal->data[nIL] -= value * dTotal->data[nJK] ;
                fTotal->data[nJK] -= value * dTotal->data[nIL] ;
                fTotal->data[nJL] -= value * dTotal->data[nIK] ;
                if ( doSpin )
                {
                    fSpin->data[nIK] -= value * dSpin->data[nJL] ;
                    fSpin->data[nIL] -= value * dSpin->data[nJK] ;
                    fSpin->data[nJK] -= value * dSpin->data[nIL] ;
                    fSpin->data[nJL] -= value * dSpin->data[nIK] ;
                }
            }
        }
        SymmetricMatrix_ScaleOffDiagonal ( fTotal, 0.5e+00 ) ;
        SymmetricMatrix_ScaleOffDiagonal ( fSpin , 0.5e+00 ) ;
        eTEI = 0.5e+00 * SymmetricMatrix_TraceOfProduct ( dTotal, fTotal, NULL ) ;
        if ( doSpin ) eTEI += 0.5 * SymmetricMatrix_TraceOfProduct ( dSpin, fSpin, NULL ) ;
    }
    return eTEI ;
}

# undef BFINDEX

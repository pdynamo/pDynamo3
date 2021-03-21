/*==================================================================================================================================
! . Functions for Fock construction.
!=================================================================================================================================*/

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

/* . Inverse fit matrix is less accurate especially for derivatives. */
/*# define _USEINVERSEFITMATRIX*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Form the fit-integral parts of the Fock matrices.
! . The fit energy and potential are also computed.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Fock_MakeFromFitIntegrals (       BlockStorage    *fitIntegrals ,
                                       SymmetricMatrix *fitMatrix    ,
                                 const Real             totalCharge  ,
                                       RealArray1D     *fitPotential ,
                                       SymmetricMatrix *dTotal       ,
                                       SymmetricMatrix *fTotal       ,
                                       Status          *status       )
{
    Real eTEI = 0.0e+00 ;
    if ( ( fitIntegrals != NULL ) &&
         ( fitPotential != NULL ) &&
         ( fitMatrix    != NULL ) &&
         ( dTotal       != NULL ) &&
         ( fTotal       != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n = View1D_Extent ( fitPotential ) ;
        auto RealArray1D *b = RealArray1D_AllocateWithExtent ( n, status ) ;
        if ( b != NULL )
        {
            auto Block   *block ;
            auto Integer  i ;
            /* . Initialization -  diagonal elements of density need scaling by 1/2. */
            /* . SymmetricMatrix_Set           ( fTotal, 0.0e+00 ) ; */
            SymmetricMatrix_ScaleDiagonal ( dTotal, 0.5e+00 ) ;
            /* . Determine b. */
            RealArray1D_Set ( b, 0.0e+00 ) ;
            List_Iterate_Initialize ( fitIntegrals->blocks ) ;
            while ( ( block = BlockStorage_Iterate ( fitIntegrals ) ) != NULL )
            {
                for ( i = 0 ; i < block->count ; i++ ) b->data[block->indices16[i]] += dTotal->data[block->indices32[i]] * block->data[i] ;
            }
            RealArray1D_Scale ( b, 2.0e+00 ) ;
/* . Remove following line if no constraint. */
            Array1D_Item ( b, n-1 ) = totalCharge ;
            /* . Find the fit potential. */
# ifdef _USEINVERSEFITMATRIX
            SymmetricMatrix_VectorMultiply ( fitMatrix, b, fitPotential, NULL ) ;
# else
            SymmetricMatrix_LinearEquationsSolve ( fitMatrix, b, fitPotential, status ) ;
            if ( Status_IsOK ( status ) )
            {
# endif
            /* . Construct Fock matrix. */
            List_Iterate_Initialize ( fitIntegrals->blocks ) ;
            while ( ( block = BlockStorage_Iterate ( fitIntegrals ) ) != NULL )
            {
                for ( i = 0 ; i < block->count ; i++ ) fTotal->data[block->indices32[i]] += fitPotential->data[block->indices16[i]] * block->data[i] ;
            }
            /* . Find the two-electron energy. */
            eTEI = 0.5e+00 * RealArray1D_Dot ( b, fitPotential, NULL ) ;
# ifndef _USEINVERSEFITMATRIX
            }
# endif
            /* . Finish up. */
            SymmetricMatrix_ScaleDiagonal ( dTotal, 2.0e+00 ) ;
            RealArray1D_Deallocate ( &b ) ;
        }
    }
    return eTEI ;
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

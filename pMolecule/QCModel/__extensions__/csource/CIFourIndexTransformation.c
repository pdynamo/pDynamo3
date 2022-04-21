/*==================================================================================================================================
! . A simple two-electron integral four index transformation module.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "CIFourIndexTransformation.h"
# include "Integer.h"
# include "Real.h"
# include "SymmetricMatrix.h"

/*==================================================================================================================================
! . Local procedures declarations.
!=================================================================================================================================*/
static void CIFIT_TransformIndex1    ( const RealArray2D *mos, const RealArrayND *tei234, DoubleSymmetricMatrix *moTEIs ) ;
static void CIFIT_TransformIndex2    ( const RealArray2D *mos, const RealArray2D *tei34, RealArrayND *tei234 ) ;
static void CIFIT_TransformIndices34 ( const RealArray2D *mos, BlockStorage *twoElectronIntegrals, RealArray2D *tei34 ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . A no-nonsense four index transformation for small numbers of MOs only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CIFourIndexTransformation ( const RealArray2D           *activeMOs            ,
                                       BlockStorage          *twoElectronIntegrals ,
                                       RealArray2D           *moTEI34              ,
                                       RealArrayND           *moTEI234             ,
                                       DoubleSymmetricMatrix *moTEIs               )
{
    if ( ( activeMOs            != NULL ) &&
         ( twoElectronIntegrals != NULL ) &&
         ( moTEI34              != NULL ) &&
         ( moTEI234             != NULL ) &&
         ( moTEIs               != NULL ) )
    {
        CIFIT_TransformIndices34 ( activeMOs, twoElectronIntegrals, moTEI34  ) ;
        CIFIT_TransformIndex2    ( activeMOs, moTEI34             , moTEI234 ) ;
        CIFIT_TransformIndex1    ( activeMOs, moTEI234            , moTEIs   ) ;
    }
}

/*==================================================================================================================================
! . CI four index transformation procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform index 1 by reading the hybrid integrals already with indices 2, 3 and 4 transformed.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIFIT_TransformIndex1 ( const RealArray2D *mos, const RealArrayND *tei234, DoubleSymmetricMatrix *moTEIs )
{
    if  ( ( mos != NULL ) && ( tei234 != NULL ) && ( moTEIs != NULL ) )
    {
        auto Integer  i, nActive, nBasis, p, pq, q, r, rs, s, upper ;
        auto Real     sum ;

        /* . Initialization. */
        nActive = View2D_Columns ( mos ) ;
        nBasis  = View2D_Rows    ( mos ) ;
        DoubleSymmetricMatrix_Set ( moTEIs, 0.0e+00 ) ;

        /* . Loop over MOs. */
        for ( p = pq = 0 ; p < nActive ; p++ )
        {
            for ( q = 0 ; q <= p ; q++, pq++ )
            {
                for ( r = rs = 0 ; r <= p ; r++ )
                {
                    if ( r == p ) upper = q ;
                    else          upper = r ;
                    for ( s = 0 ; s <= upper ; s++, rs++ )
                    {
                        /* . This is a dot-product so ultimately should use slices. */
                        for ( i = 0, sum = 0.0e+00 ; i < nBasis ; i++ ) sum += Array2D_Item ( mos, i, p ) * ArrayND_Item3D ( tei234, i, q, rs ) ;
                        DoubleSymmetricMatrix_SetItem ( moTEIs, p, q, r, s, sum, NULL ) ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform index 2 by reading the hybrid integrals already with indices 3 and 4 transformed.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIFIT_TransformIndex2 ( const RealArray2D *mos, const RealArray2D *tei34, RealArrayND *tei234 )
{
    if  ( ( mos != NULL ) && ( tei34 != NULL ) && ( tei234 != NULL ) )
    {
        auto Integer  i, ij, j, nActive, nBasis,q, r, rs, s ;
        auto Real     t ;

        /* . Initialization. */
        nActive = View2D_Columns ( mos ) ;
        nBasis  = View2D_Rows    ( mos ) ;
        RealArrayND_Set ( tei234, 0.0e+00, NULL ) ;

        /* . Loop over MO pairs. */
        for ( r = rs = 0 ; r < nActive ; r++ )
        {
            for ( s = 0 ; s <= r ; s++, rs++ )
            {
                /* . Loop over AOs. */
                for ( i = ij = 0 ; i < nBasis ; i++ )
                {
                    for ( j = 0 ; j <= i ; ij++, j++ )
                    {
                        t  = Array2D_Item ( tei34, ij, rs ) ;
                        if ( i == j ) t *= 0.5e+00 ;
                        for ( q = 0 ; q < nActive ; q++ )
                        {
                            ArrayND_Item3D ( tei234, i, q, rs ) += t * Array2D_Item ( mos, j, q ) ;
                            ArrayND_Item3D ( tei234, j, q, rs ) += t * Array2D_Item ( mos, i, q ) ;
                        }
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform indices 3 and 4 together by reading the A.O. integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIFIT_TransformIndices34 ( const RealArray2D *mos, BlockStorage *twoElectronIntegrals, RealArray2D *tei34 )
{
    if  ( ( mos != NULL ) && ( twoElectronIntegrals != NULL ) && ( tei34 != NULL ) )
    {
        auto Integer  i, ii, ij, j, k, kl, l, m, n, nActive, r, rs, s ;
        auto Real     t, wij, wkl ;
        auto Block   *block ;

        /* . Initialization. */
        nActive = View2D_Columns ( mos ) ;
        RealArray2D_Set ( tei34, 0.0e+00 ) ;

        /* . Loop over the integral blocks. */
        List_Iterate_Initialize ( twoElectronIntegrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( twoElectronIntegrals ) ) != NULL )
        {
            /* . Loop over the integrals. */
            for ( ii = 0, n = 0 ; ii < block->count ; ii++, n += 4 )
            {
                /* . Get the data. */
                i = block->indices16[n  ] ;
                j = block->indices16[n+1] ;
                k = block->indices16[n+2] ;
                l = block->indices16[n+3] ;
                t = block->data[ii] ;

                /* . Shuffle the index pairs. */
	        if ( i < j ) { m = i ; i = j ; j = m ; }
                if ( k < l ) { m = k ; k = l ; l = m ; }

                /* . Shuffle the indices of both pairs. */
                if ( ( i < k ) || ( ( i == k ) && ( j < l  ) ) ) { m = i ; i = k ; k = m ; m = j ; j = l ; l = m ; }

                /* . RealArray2D indices. */
                ij = SymmetricMatrix_ItemIndex ( i, j ) ;
                kl = SymmetricMatrix_ItemIndex ( k, l ) ;

                /* . Get the appropriate scaling factors. */
                wij = 1.0e+00 ;
                wkl = 1.0e+00 ;
	        if ( i == j ) wij *= 0.5e+00 ;
	        if ( k == l ) wkl *= 0.5e+00 ;
                if ( ij == kl ) { wij *= 0.5e+00 ; wkl *= 0.5e+00 ; }

                /* . Loop over MO pairs. */
                for ( r = rs = 0 ; r < nActive ; r++ )
                {
                    for ( s = 0 ; s <= r ; rs++, s++ )
                    {
                        Array2D_Item ( tei34, ij, rs ) += t * ( Array2D_Item ( mos, k, r ) * Array2D_Item ( mos, l, s ) + Array2D_Item ( mos, l, r ) * Array2D_Item ( mos, k, s ) ) * wkl ;
                        Array2D_Item ( tei34, kl, rs ) += t * ( Array2D_Item ( mos, i, r ) * Array2D_Item ( mos, j, s ) + Array2D_Item ( mos, j, r ) * Array2D_Item ( mos, i, s ) ) * wij ;
                    }
                }
            }
        }
    }
}

/*==================================================================================================================================
! . Integrals - 1 basis, 1 electron, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_f1Cp1.h"
# include "GaussianBasisSubsidiary.h"
# include "GaussianBasisTransform.h"
# include "GaussianNucleus.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Selected( selection, i ) ( (selection) == NULL ? True : Block_Item ( selection->flags, i ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fit-nuclear/point derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s1 and real 5 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Cm1R1 ( const GaussianBasis *fBasis ,
                                      const Real          *rF         ,
                                      const RealArray1D   *charges    ,
                                      const RealArray1D   *widthsE    ,
                                      const RealArray1D   *widthsN    ,
                                      const Coordinates3  *rNP        ,
                                      const Selection     *selectionN ,
                                      const RealArray1D   *dOneF      ,
                                      const Integer        s1         ,
                                            Integer       *iWork      ,
                                            Real          *rWork      ,
                                            Real          *gF         ,
                                            Coordinates3  *gN         )
{
    Integer       f, fammax, ff, fp, fShell, k, m, n, nCFuncF, nRoots ;
    Real          aa, aandb, ab, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dGx, dGy, dGz, dnuc,
                  expfac, expN, fac, facN, fac2, f00, qN, rho, u2,
                  xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
    Real         *pGx, *pGy, *pGz, *values = NULL, *work = NULL ;     
    Integer      *Ix, *Iy, *Iz ;
    Real         *Ci, *gT, *gX, *gY, *gZ,
                  xidG[MAXAMP1], yidG[MAXAMP1], zidG[MAXAMP1] ,
                  xint[MAXAMP2], yint[MAXAMP2], zint[MAXAMP2] ;
    Real         *rN ;
    RealArray2D  *fC2S ;
    RysQuadrature roots ;
    /* . Initialization. */
    for ( m = 0 ; m < 3 ; m++ ) { gF[m] = 0.0e+00 ; }
    /* . Set pointers. */
    Ci = &rWork[   0] ; gT = &rWork[  s1] ;
    gX = &rWork[2*s1] ; gY = &rWork[3*s1] ; gZ = &rWork[4*s1] ;
    Ix = &iWork[   0] ; Iy = &iWork[  s1] ; Iz = &iWork[2*s1] ;
    /* . Loop over the nuclear densities. */
    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
    {
        if ( _Selected ( selectionN, k ) )
        {
            expN = _GetWidthE ( widthsE, k ) ;
            facN = _GetWidthN ( widthsN, k ) ;
            qN   = - Array1D_Item ( charges, k ) ; /* . Negative as electrons. */
            rN   = Coordinates3_RowPointer ( rNP, k ) ;
            /* . Initialize some accumulators. */
            dGx = dGy = dGz = 0.0e+00 ;
            /* . Loop over shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Initialization. */
                fammax  = fBasis->shells[fShell].lHigh ;
                fC2S    = fBasis->shells[fShell].c2s ;
                nCFuncF = fBasis->shells[fShell].nCBF     ;
                nRoots  = ( fammax + 1 ) / 2 + 1 ;
                /* . Index arrays. */
                for ( f = 0 ; f < nCFuncF ; f++ )
                {
   	            Ix[f] = fBasis->shells[fShell].cbfPowX[f] ;
	            Iy[f] = fBasis->shells[fShell].cbfPowY[f] ;
	            Iz[f] = fBasis->shells[fShell].cbfPowZ[f] ;
                }
                /* . Initialize the integral blocks. */
                for ( n = 0 ; n < nCFuncF ; n++ ) gX[n] = gY[n] = gZ[n] = 0.0e+00 ;
                /* . Loop over primitives. */
                for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                {
                    /* . Get some information for the primitive. */
	            aa     = fBasis->shells[fShell].primitives[fp].exponent ;
                    expfac = PI252 / aa ;
                    /* . Calculate some factors. */
                    ab    = aa * expN ;
                    aandb = aa + expN ;
                    rho   = ab / aandb ;
                    dnuc  = expfac * ( facN * qN ) / ( expN * sqrt ( aandb ) ) ;
                    /* . Start of point-specific code. */
                    /* . Calculate the rys polynomial roots. */
                    c1x = ( rF[0] - rN[0] ) ;
                    c1y = ( rF[1] - rN[1] ) ;
                    c1z = ( rF[2] - rN[2] ) ;
                    RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                    c1x *= aa ;
                    c1y *= aa ;
                    c1z *= aa ;
                    c3x  = expN * ( rN[0] - rF[0] ) ;
                    c3y  = expN * ( rN[1] - rF[1] ) ;
                    c3z  = expN * ( rN[2] - rF[2] ) ;
                    /* . Coefficient array. */
                    for ( n = 0 ; n < nCFuncF ; n++ ) Ci[n] = dnuc * fBasis->shells[fShell].primitives[fp].cCBF[n] ;
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( aa   + u2 ) * fac2 ;
                        b00   =          u2   * fac2 ;
                        b10   = ( expN + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        xc00  = u2 * c3x * fac ;
                        yc00  = u2 * c3y * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1  ( fammax+1 ,
                                                         0        ,
                                                         b00      ,
                                                         b10      ,
                                                         bp01     ,
                                                         f00      ,
                                                         xc00     ,
                                                         xcp00    ,
                                                         yc00     ,
                                                         ycp00    ,
                                                         zc00     ,
                                                         zcp00    ,
                                                         1        ,
                                                         xint     ,
                                                         yint     ,
                                                         zint     );
                        GaussianBasisSubsidiary_f1Xg1r ( xint     ,
                                                         yint     ,
                                                         zint     ,
                                                         aa       ,
                                                         fammax   ,
                                                         0        ,
                                                         1        ,
                                                         1        ,
                                                         xidG     ,
                                                         yidG     ,
                                                         zidG     );
                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFuncF ; n++ )
                        {
                            gX[n] += Ci[n] * ( xidG[Ix[n]] * yint[Iy[n]] * zint[Iz[n]] ) ;
                            gY[n] += Ci[n] * ( xint[Ix[n]] * yidG[Iy[n]] * zint[Iz[n]] ) ;
                            gZ[n] += Ci[n] * ( xint[Ix[n]] * yint[Iy[n]] * zidG[Iz[n]] ) ;
                        }
                    } /* . nRoots. */
                } /* . fP. */
                /* . Transform the integrals. */
                work   = gT ;
                values = gX ; GaussianBasisTransform1 ( fC2S, &values, &work ) ; pGx = values ;
                values = gY ; GaussianBasisTransform1 ( fC2S, &values, &work ) ; pGy = values ;
                values = gZ ; GaussianBasisTransform1 ( fC2S, &values, &work ) ; pGz = values ;
                /* .  Add in the blocks of integrals to the derivatives. */
                for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++ )
                {
                    ff   = fBasis->shells[fShell].nStart + f ;
                    fac  = Array1D_Item ( dOneF, ff ) ;
                    dGx += fac * pGx[f] ;
                    dGy += fac * pGy[f] ;
                    dGz += fac * pGz[f] ;
                }
            } /* . fShell. */
            /* . Sum in the contributions to the gradients. */
            gF[0] += dGx ; gF[1] += dGy ; gF[2] += dGz ;
            Coordinates3_DecrementRow ( gN, k, dGx, dGy, dGz ) ;
        } /* . Selected. */
    } /* . k. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fit-nuclear/point integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s1 and real 3 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Cm1V ( const GaussianBasis *fBasis ,
                                     const Real          *rF         ,
                                     const RealArray1D   *charges    ,
                                     const RealArray1D   *widthsE    ,
                                     const RealArray1D   *widthsN    ,
                                     const Coordinates3  *rNP        ,
                                     const Selection     *selectionN ,
                                     const Integer        s1         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           RealArray1D   *integrals  )
{
    Integer       f, fammax, ff, fp, fShell, k, m, n, nCFuncF, nRoots ;
    Real          aa, aandb, ab, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dnuc, expfac, expN,
                  fac, facN, fac2, f00, qN, rho, u2, xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
    Real         *pG, *values = NULL, *work = NULL ;     
    Integer      *Ix, *Iy, *Iz ;
    Real         *Ci, *g , *gT ,
                  xint[MAXAMP1], yint[MAXAMP1], zint[MAXAMP1] ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Set pointers. */
    Ci = &rWork[0] ; g  = &rWork[s1] ; gT = &rWork[2*s1] ;
    Ix = &iWork[0] ; Iy = &iWork[s1] ; Iz = &iWork[2*s1] ;
    /* . Loop over shells. */
    for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
    {
        /* . Initialization. */
        fammax  = fBasis->shells[fShell].lHigh ;
        nCFuncF = fBasis->shells[fShell].nCBF     ;
        nRoots  = fammax / 2 + 1 ;
        /* . Index arrays. */
        for ( f = 0 ; f < nCFuncF ; f++ )
        {
   	    Ix[f] = fBasis->shells[fShell].cbfPowX[f] ;
	    Iy[f] = fBasis->shells[fShell].cbfPowY[f] ;
	    Iz[f] = fBasis->shells[fShell].cbfPowZ[f] ;
        }
        /* . Initialize the integral blocks. */
        for ( f = 0 ; f < nCFuncF ; f++ ) g[f] = 0.0e+00 ;
        /* . Loop over primitives. */
        for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
        {
            /* . Get some information for the primitive. */
	    aa     = fBasis->shells[fShell].primitives[fp].exponent ;
            expfac = PI252 / aa ;
            /* . Loop over the nuclear densities. */
            for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
            {
                if ( _Selected ( selectionN, k ) )
                {
                    expN  = _GetWidthE ( widthsE, k ) ;
                    facN  = _GetWidthN ( widthsN, k ) ;
                    qN    = - Array1D_Item ( charges, k ) ; /* . Negative as electrons. */
                    rN    = Coordinates3_RowPointer ( rNP, k ) ;
                    /* . Calculate some factors. */
                    ab    = aa * expN ;
                    aandb = aa + expN ;
                    rho   = ab / aandb ;
                    dnuc  = expfac * ( facN * qN ) / ( expN * sqrt ( aandb ) ) ;
                    /* . Calculate the rys polynomial roots. */
                    c1x   = ( rF[0] - rN[0] ) ;
                    c1y   = ( rF[1] - rN[1] ) ;
                    c1z   = ( rF[2] - rN[2] ) ;
                    RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                    /* . Calculate some displacements. */
                    c1x  *= aa ;
                    c1y  *= aa ;
                    c1z  *= aa ;
                    c3x   = expN * ( rN[0] - rF[0] ) ;
                    c3y   = expN * ( rN[1] - rF[1] ) ;
                    c3z   = expN * ( rN[2] - rF[2] ) ;
                    /* . Coefficient array. */
                    for ( n = 0 ; n < nCFuncF ; n++ ) Ci[n] = dnuc * fBasis->shells[fShell].primitives[fp].cCBF[n] ;
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( aa   + u2 ) * fac2 ;
                        b00   =          u2   * fac2 ;
                        b10   = ( expN + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        xc00  = u2 * c3x * fac ;
                        yc00  = u2 * c3y * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1 ( fammax ,
                                                        0      ,
                                                        b00    ,
                                                        b10    ,
                                                        bp01   ,
                                                        f00    ,
                                                        xc00   ,
                                                        xcp00  ,
                                                        yc00   ,
                                                        ycp00  ,
                                                        zc00   ,
                                                        zcp00  ,
                                                        1      ,
                                                        xint   ,
                                                        yint   ,
                                                        zint   ) ;
                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFuncF ; n++ ) g[n] += ( Ci[n] * xint[Ix[n]] * yint[Iy[n]] * zint[Iz[n]] ) ;
                    } /* . nRoots. */
                } /* . Selected. */
            } /* . k. */
        } /* . fp. */
        /* . Transform the integrals. */
        values = g  ;
        work   = gT ;
        GaussianBasisTransform1 ( fBasis->shells[fShell].c2s, &values, &work ) ;
        pG = values ;
        /* . Save the integrals. */
        for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++ )
        {
            ff = fBasis->shells[fShell].nStart + f ;
            Array1D_Item ( integrals, ff ) += pG[f] ;
        }
    } /* . fShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fit-nuclear/point potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Work space: integer 3 * s1 and real 3 * s1 where s1 = ( maximum shell size ). */
void GaussianBasisIntegrals_f1Cp1V ( const GaussianBasis *fBasis     ,
                                     const Real          *rF         ,
                                     const RealArray1D   *widthsE    ,
                                     const RealArray1D   *widthsN    ,
                                     const Coordinates3  *rNP        ,
                                     const Selection     *selectionN ,
                                     const RealArray1D   *dOneF      ,
                                     const Integer        s1         ,
                                           Integer       *iWork      ,
                                           Real          *rWork      ,
                                           RealArray1D   *potentials )
{
    Integer       f, fammax, ff, fp, fShell, k, m, n, nCFuncF, nRoots ;
    Real          aa, aandb, ab, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dnuc, expfac, expN,
                  fac, facN, fac2, f00, pot, rho, u2, xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
    Real         *pG, *values = NULL, *work = NULL ;     
    Integer      *Ix, *Iy, *Iz ;
    Real         *Ci, *g , *gT ,
                  xint[MAXAMP1], yint[MAXAMP1], zint[MAXAMP1] ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Set pointers. */
    Ci = &rWork[0] ; g  = &rWork[s1] ; gT = &rWork[2*s1] ;
    Ix = &iWork[0] ; Iy = &iWork[s1] ; Iz = &iWork[2*s1] ;
    /* . Loop over the points. */
    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
    {
        if ( _Selected ( selectionN, k ) )
        {
            expN = _GetWidthE ( widthsE, k ) ;
            facN = _GetWidthN ( widthsN, k ) ;
            rN   = Coordinates3_RowPointer ( rNP, k ) ;
            pot  = 0.0e+00 ;
            /* . Loop over shells. */
            for ( fShell = 0 ; fShell < fBasis->nShells ; fShell++ )
            {
                /* . Initialization. */
                fammax  = fBasis->shells[fShell].lHigh ;
                nCFuncF = fBasis->shells[fShell].nCBF  ;
                nRoots  = fammax / 2 + 1 ;
                /* . Index arrays. */
                for ( f = 0 ; f < nCFuncF ; f++ )
                {
   	            Ix[f] = fBasis->shells[fShell].cbfPowX[f] ;
	            Iy[f] = fBasis->shells[fShell].cbfPowY[f] ;
	            Iz[f] = fBasis->shells[fShell].cbfPowZ[f] ;
                }
                /* . Initialize the integral blocks. */
                for ( f = 0 ; f < nCFuncF ; f++ ) g[f] = 0.0e+00 ;
                /* . Loop over primitives. */
                for ( fp = 0 ; fp < fBasis->shells[fShell].nPrimitives ; fp++ )
                {
                    /* . Get some information for the primitive. */
	            aa     = fBasis->shells[fShell].primitives[fp].exponent ;
                    expfac = PI252 / aa ;
                    /* . Start of point-specific code. */
                    /* . Calculate some factors. */
                    ab     = aa * expN ;
                    aandb  = aa + expN ;
                    rho    = ab / aandb ;
                    dnuc   = expfac * facN / ( expN * sqrt ( aandb ) ) ;
                    /* . Calculate the rys polynomial roots. */
                    c1x = ( rF[0] - rN[0] ) ;
                    c1y = ( rF[1] - rN[1] ) ;
                    c1z = ( rF[2] - rN[2] ) ;
                    RysQuadrature_Roots ( &roots, nRoots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                    c1x *= aa ;
                    c1y *= aa ;
                    c1z *= aa ;
                    c3x  = expN * ( rN[0] - rF[0] ) ;
                    c3y  = expN * ( rN[1] - rF[1] ) ;
                    c3z  = expN * ( rN[2] - rF[2] ) ;
                    /* . Coefficient array. */
                    for ( n = 0 ; n < nCFuncF ; n++ ) Ci[n] = dnuc * fBasis->shells[fShell].primitives[fp].cCBF[n] ;
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nRoots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( aa   + u2 ) * fac2 ;
                        b00   =          u2   * fac2 ;
                        b10   = ( expN + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        xc00  = u2 * c3x * fac ;
                        yc00  = u2 * c3y * fac ;
                        zc00  = u2 * c3z * fac ;
                        GaussianBasisSubsidiary_f1Cg1 ( fammax ,
                                                        0      ,
                                                        b00    ,
                                                        b10    ,
                                                        bp01   ,
                                                        f00    ,
                                                        xc00   ,
                                                        xcp00  ,
                                                        yc00   ,
                                                        ycp00  ,
                                                        zc00   ,
                                                        zcp00  ,
                                                        1      ,
                                                        xint   ,
                                                        yint   ,
                                                        zint   ) ;

                        /* . Assemble the integrals. */
                        for ( n = 0 ; n < nCFuncF ; n++ ) g[n] += ( Ci[n] * xint[Ix[n]] * yint[Iy[n]] * zint[Iz[n]] ) ;
                    } /* . nRoots. */
                } /* . fP. */
                /* . Transform the integrals. */
                values = g  ;
                work   = gT ;
                GaussianBasisTransform1 ( fBasis->shells[fShell].c2s, &values, &work ) ;
                pG = values ;
                /* .  Add in the block of integrals to the potential. */
                for ( f = 0 ; f < fBasis->shells[fShell].nBasis ; f++ )
                {
                    ff   = fBasis->shells[fShell].nStart + f ;
                    pot += ( Array1D_Item ( dOneF, ff ) * pG[f] ) ;
                }
            } /* . fShell. */
            /* . Save the potential - negative as electrons. */
            Array1D_Item ( potentials, k ) -= pot ;
        } /* . Selected. */
    } /* . k. */
}

# undef _Selected

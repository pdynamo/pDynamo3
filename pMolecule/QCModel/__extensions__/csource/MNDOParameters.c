/*==================================================================================================================================
! . Procedures for dealing with parameters for MNDO methods.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "Memory.h"
# include "MNDODefinitions.h"
# include "MNDOParameters.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real AIJL                        ( const Real z1, const Real z2, const Integer  n1, const Integer  n2, const Integer  l ) ;
static Real Binomial                    ( const Integer  n, const Integer  k ) ;
static Real Factorial                   ( const Integer  n ) ;
static void MNDOParameters_FillDDPPO    ( MNDOParameters *self, const Real *octeis ) ;
static Real POIJ                        ( const Integer  l, const Real d, const Real fg ) ;
static Real RadialSlaterCondonParameter ( const Integer  k, const Integer  na, const Real ea, const Integer  nb, const Real eb, const Integer  nc, const Real ec, const Integer  nd, const Real ed ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOParameters *MNDOParameters_Allocate ( void )
{
    MNDOParameters *self = Memory_AllocateType ( MNDOParameters ) ;
    /* . Counters. */
    self->QDIATOMIC      = False ;
    self->QDIATOMICFLAGS = NULL  ;
    self->atomicNumber   =  0 ;
    self->iii            =  0 ;
    self->iiid           =  0 ;
    self->ir016          =  0 ;
    self->ir066          =  0 ;
    self->ir244          =  0 ;
    self->ir266          =  0 ;
    self->ir466          =  0 ;
    self->nam1pm3g       =  0 ;
    self->ndiatomic      =  0 ;
    self->norbitals      =  0 ;
    self->npddg          =  0 ;
    self->qnd            =  0 ;
    self->qnp            =  0 ;
    self->qns            =  0 ;
    /* . Derived quantities. */
    self->nocteis      =    0    ;
    self->octeiindices = NULL    ;
    self->hpp          = 0.0e+00 ;
    self->octeivalues  = NULL    ;
    /* . Input parameters. */
    self->ad0        = 0.0e+00 ;
    self->alp0       = 0.0e+00 ;
    self->am0        = 0.0e+00 ;
    self->aq0        = 0.0e+00 ;
    self->betad0     = 0.0e+00 ;
    self->betap0     = 0.0e+00 ;
    self->betas0     = 0.0e+00 ;
    self->dd0        = 0.0e+00 ;
    self->eheat0     = 0.0e+00 ;
    self->eisol0     = 0.0e+00 ;
    self->f0sd0      = 0.0e+00 ;
    self->gphot0     = 1.0e+00 ;
    self->gpp0       = 0.0e+00 ;
    self->gp20       = 0.0e+00 ;
    self->gsp0       = 0.0e+00 ;
    self->gss0       = 0.0e+00 ;
    self->g2sd0      = 0.0e+00 ;
    self->hsp0       = 0.0e+00 ;
    self->pcore0     = 0.0e+00 ;
    self->qq0        = 0.0e+00 ;
    self->udd0       = 0.0e+00 ;
    self->upp0       = 0.0e+00 ;
    self->uss0       = 0.0e+00 ;
    self->zcore0     = 0.0e+00 ;
    self->zetad0     = 0.0e+00 ;
    self->zetap0     = 0.0e+00 ;
    self->zetas0     = 0.0e+00 ;
    self->zdn0       = 0.0e+00 ;
    self->zpn0       = 0.0e+00 ;
    self->zsn0       = 0.0e+00 ;
    self->beta0      = NULL ;
    self->diatomica0 = NULL ;
    self->diatomicx0 = NULL ;
    self->fn10       = NULL ;
    self->fn20       = NULL ;
    self->fn30       = NULL ;
    self->pddgc0     = NULL ;
    self->pddge0     = NULL ;
    self->uspd0      = NULL ;
    /* . Internal parameters. */
    self->ad        = 0.0e+00 ;
    self->alp       = 0.0e+00 ;
    self->am        = 0.0e+00 ;
    self->aq        = 0.0e+00 ;
    self->betad     = 0.0e+00 ;
    self->betap     = 0.0e+00 ;
    self->betas     = 0.0e+00 ;
    self->dd        = 0.0e+00 ;
    self->eheat     = 0.0e+00 ;
    self->eisol     = 0.0e+00 ;
    self->f0sd      = 0.0e+00 ;
    self->gphot     = 1.0e+00 ;
    self->gpp       = 0.0e+00 ;
    self->gp2       = 0.0e+00 ;
    self->gsp       = 0.0e+00 ;
    self->gss       = 0.0e+00 ;
    self->g2sd      = 0.0e+00 ;
    self->hsp       = 0.0e+00 ;
    self->pcore     = 0.0e+00 ;
    self->qq        = 0.0e+00 ;
    self->udd       = 0.0e+00 ;
    self->upp       = 0.0e+00 ;
    self->uss       = 0.0e+00 ;
    self->zcore     = 0.0e+00 ;
    self->zetad     = 0.0e+00 ;
    self->zetap     = 0.0e+00 ;
    self->zetas     = 0.0e+00 ;
    self->zdn       = 0.0e+00 ;
    self->zpn       = 0.0e+00 ;
    self->zsn       = 0.0e+00 ;
    self->beta      = NULL ;
    self->diatomica = NULL ;
    self->diatomicx = NULL ;
    self->fn1       = NULL ;
    self->fn2       = NULL ;
    self->fn3       = NULL ;
    self->pddgc     = NULL ;
    self->pddge     = NULL ;
    self->uspd      = NULL ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation of individual sets of arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOParameters_AllocateDiatomic ( MNDOParameters *self, const Integer  n )
{
    if ( ( self != NULL ) && ( n > 0 ) )
    {
        auto Integer  i ;
        self->ndiatomic      = n    ;
        self->QDIATOMIC      = True ;
        self->QDIATOMICFLAGS = Memory_AllocateArrayOfTypes ( n, Boolean ) ;
        self->diatomica0     = Memory_AllocateArrayOfTypes ( n, Real    ) ;
        self->diatomicx0     = Memory_AllocateArrayOfTypes ( n, Real    ) ;
        for ( i = 0 ; i < n ; i++ )
        {
            self->QDIATOMICFLAGS[i] = False   ;
            self->diatomica0    [i] = 0.0e+00 ;
            self->diatomicx0    [i] = 0.0e+00 ;
        }
    }
}
void MNDOParameters_AllocateFN123 ( MNDOParameters *self, const Integer  n )
{
    if ( ( self != NULL ) && ( n > 0 ) )
    {
        self->nam1pm3g = n ;
        self->fn10     = Memory_AllocateArrayOfTypes ( n, Real ) ;
        self->fn20     = Memory_AllocateArrayOfTypes ( n, Real ) ;
        self->fn30     = Memory_AllocateArrayOfTypes ( n, Real ) ;
    }
}
void MNDOParameters_AllocatePDDG ( MNDOParameters *self, const Integer  n )
{
    if ( ( self != NULL ) && ( n > 0 ) )
    {
        self->npddg  = n ;
        self->pddgc0 = Memory_AllocateArrayOfTypes ( n, Real ) ;
        self->pddge0 = Memory_AllocateArrayOfTypes ( n, Real ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the one-center TEIs.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define NREPD 53 /* . For repd only the last 52 are needed as Fortran indexing is kept here. */
void MNDOParameters_CalculateOneCenterTEIs ( MNDOParameters *self )
{
    /* . Check that there are orbitals. */
    if ( self->norbitals > 0 )
    {
        auto Integer     i, n = 0, t = 0 ;
        auto Cardinal16 *indices ;
        auto Real        repd[NREPD], *values ;

        /* . Initialization. */
        for ( i = 0 ; i < NREPD ; i++ ) repd[i] = 0.0e+00 ;

        /* . Determine the number of integrals. */
        if      ( self->norbitals == 1 ) n =   1 ;
        else if ( self->norbitals == 4 ) n =  16 ;
        else if ( self->norbitals == 9 ) n = 155 ;

        /* . Allocate the internal arrays. */
        self->nocteis      = n ;
        self->octeiindices = Memory_AllocateArrayOfTypes ( 4 * n, Cardinal16  ) ;
        self->octeivalues  = Memory_AllocateArrayOfTypes (     n, Real       ) ;

        /* . Calculate the integrals. */
        n = 0 ;
        indices = self->octeiindices ;
        values  = self->octeivalues  ;

        /* . ( s s | s s ). */
        values[n] = self->gss ; indices[t] = 0 ; indices[t+1] = 0 ; indices[t+2] = 0 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;

# ifdef PRINTMOPACPARAMETERS
{
printf ( "\n\nMopac Parameters for Element %d:\n", self->atomicNumber ) ;
printf ( "MOPACPARAMETERS> %d   EHEAT   %.15f   ATOMIC\n", self->atomicNumber, self->eheat ) ;
}
# endif

        /* . Check for p-orbitals. */
        /* . Order is z, x, y. */
        if ( self->norbitals >= 4 )
        {
            /* . (pp'|pp') Integral. */
            self->hpp = 0.5e+00 * ( self->gpp - self->gp2 ) ;

            /* . ( x x | x x ), etc. */
            values[n] = self->gpp ; indices[t] = 1 ; indices[t+1] = 1 ; indices[t+2] = 1 ; indices[t+3] = 1 ; n += 1 ; t += 4 ;
            values[n] = self->gpp ; indices[t] = 2 ; indices[t+1] = 2 ; indices[t+2] = 2 ; indices[t+3] = 2 ; n += 1 ; t += 4 ;
            values[n] = self->gpp ; indices[t] = 3 ; indices[t+1] = 3 ; indices[t+2] = 3 ; indices[t+3] = 3 ; n += 1 ; t += 4 ;

            /* . ( x x | s s ), etc. */
            values[n] = self->gsp ; indices[t] = 1 ; indices[t+1] = 1 ; indices[t+2] = 0 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;
            values[n] = self->gsp ; indices[t] = 2 ; indices[t+1] = 2 ; indices[t+2] = 0 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;
            values[n] = self->gsp ; indices[t] = 3 ; indices[t+1] = 3 ; indices[t+2] = 0 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;

            /* . ( x s | x s ), etc. */
            values[n] = self->hsp ; indices[t] = 1 ; indices[t+1] = 0 ; indices[t+2] = 1 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;
            values[n] = self->hsp ; indices[t] = 2 ; indices[t+1] = 0 ; indices[t+2] = 2 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;
            values[n] = self->hsp ; indices[t] = 3 ; indices[t+1] = 0 ; indices[t+2] = 3 ; indices[t+3] = 0 ; n += 1 ; t += 4 ;

            /* . ( y y | x x ), etc. */
            values[n] = self->gp2 ; indices[t] = 2 ; indices[t+1] = 2 ; indices[t+2] = 1 ; indices[t+3] = 1 ; n += 1 ; t += 4 ;
            values[n] = self->gp2 ; indices[t] = 3 ; indices[t+1] = 3 ; indices[t+2] = 1 ; indices[t+3] = 1 ; n += 1 ; t += 4 ;
            values[n] = self->gp2 ; indices[t] = 3 ; indices[t+1] = 3 ; indices[t+2] = 2 ; indices[t+3] = 2 ; n += 1 ; t += 4 ;

            /* . ( y x | y x ), etc. */
            values[n] = self->hpp ; indices[t] = 2 ; indices[t+1] = 1 ; indices[t+2] = 2 ; indices[t+3] = 1 ; n += 1 ; t += 4 ;
            values[n] = self->hpp ; indices[t] = 3 ; indices[t+1] = 1 ; indices[t+2] = 3 ; indices[t+3] = 1 ; n += 1 ; t += 4 ;
            values[n] = self->hpp ; indices[t] = 3 ; indices[t+1] = 2 ; indices[t+2] = 3 ; indices[t+3] = 2 ; n += 1 ; t += 4 ;
        }

        /* . Check for d-orbitals. */
        /* . Order is z2, xz, yz,  x2y2, xy. */
        if ( self->norbitals >= 9 )
        {
            auto Integer  nd, ns ;
            auto Real     ed, ep, es, r016, r036, r066, r125, r155, r234, r236, r244, r246, r266, r355, r466,  s3, s5, s15 ;

            /* . Constants. */
            s3  = sqrt (  3.0e+00 ) ;
            s5  = sqrt (  5.0e+00 ) ;
            s15 = sqrt ( 15.0e+00 ) ;

            /* . Calculate Slater-Condon parameters - rLIJ. */
            /* . L -  L quantum number of Slater-Condon parameter.
            !  . I -  ss 1, sp 2, pp 3, sd 4, pd 5, dd 6 - electron 1.
            !  . J -  ss 1, sp 2, pp 3, sd 4, pd 5, dd 6 - electron 2. */
            ns   = self->iii  ;
            nd   = self->iiid ;
            es   = self->zsn  ;
            ep   = self->zpn  ;
            ed   = self->zdn  ;
            r016 = RadialSlaterCondonParameter ( 0, ns, es, ns, es, nd, ed, nd, ed ) ;
            r036 = RadialSlaterCondonParameter ( 0, ns, ep, ns, ep, nd, ed, nd, ed ) ;
            r066 = RadialSlaterCondonParameter ( 0, nd, ed, nd, ed, nd, ed, nd, ed ) ;
            r155 = RadialSlaterCondonParameter ( 1, ns, ep, nd, ed, ns, ep, nd, ed ) ;
            r125 = RadialSlaterCondonParameter ( 1, ns, es, ns, ep, ns, ep, nd, ed ) ;
            r244 = RadialSlaterCondonParameter ( 2, ns, es, nd, ed, ns, es, nd, ed ) ;
            r236 = RadialSlaterCondonParameter ( 2, ns, ep, ns, ep, nd, ed, nd, ed ) ;
            r266 = RadialSlaterCondonParameter ( 2, nd, ed, nd, ed, nd, ed, nd, ed ) ;
            r234 = RadialSlaterCondonParameter ( 2, ns, ep, ns, ep, ns, es, nd, ed ) ;
            r246 = RadialSlaterCondonParameter ( 2, ns, es, nd, ed, nd, ed, nd, ed ) ;
            r355 = RadialSlaterCondonParameter ( 3, ns, ep, nd, ed, ns, ep, nd, ed ) ;
            r466 = RadialSlaterCondonParameter ( 4, nd, ed, nd, ed, nd, ed, nd, ed ) ;

            /* . Override r016 and r244 if these are input. */
            if ( self->f0sd0 > 1.0e-3 ) r016 = self->f0sd ;
            if ( self->g2sd0 > 1.0e-3 ) r244 = self->g2sd ;

            /* . Modify the atomic energy for those elements with partially-filled d-orbitals. */
            /* . r016: <ss|dd>.
            !  . r066: <dd|dd> "0" term.
            !  . r244: <sd|sd>.
            !  . r266: <dd|dd> "2" term.
            !  . r466: <dd|dd> "4" term. */
            self->eisol += ( Real ) self->ir016 * r016 + ( Real ) self->ir066 * r066 - ( Real ) self->ir244 * r244 / 5.0e+00 - ( Real ) self->ir266 * r266 / 49.0e+00 - ( Real ) self->ir466 * r466 / 49.0e+00 ;

            /* . Determine the integrals. */
            repd[1 ] = r016 ;
            repd[2 ] = 2.0e+00/(3.0e+00*s5)*r125 ;
            repd[3 ] = 1.0e+00/s15*r125 ;
            repd[4 ] = 2.0e+00/(5.0e+00*s5)*r234 ;
            repd[5 ] = r036 + 4.0e+00/35.0e+00*r236 ;
            repd[6 ] = r036 + 2.0e+00/35.0e+00*r236 ;
            repd[7 ] = r036 - 4.0e+00/35.0e+00*r236 ;
            repd[8 ] = -1.0e+00/(3.0e+00*s5)*r125 ;
            repd[9 ] = sqrt(3.0e+00/125.0e+00)*r234 ;
            repd[10] = s3/35.0e+00*r236 ;
            repd[11] = 3.0e+00/35.0e+00*r236 ;
            repd[12] = -1.0e+00/(5.0e+00*s5)*r234 ;
            repd[13] = r036 - 2.0e+00/35.0e+00*r236 ;
            repd[14] = -2.0e+00*s3/35.0e+00*r236 ;
            repd[15] = -repd[3] ;
            repd[16] = -repd[11] ;
            repd[17] = -repd[9] ;
            repd[18] = -repd[14] ;
            repd[19] = 1.0e+00/5.0e+00*r244 ;
            repd[20] = 2.0e+00/(7.0e+00*s5)*r246 ;
            repd[21] = repd[20]/2.0e+00 ;
            repd[22] = -repd[20] ;
            repd[23] = 4.0e+00/15.0e+00*r155 + 27.0e+00/245.0e+00*r355 ;
            repd[24] = 2.0e+00*s3/15.0e+00*r155 - 9.0e+00*s3/245.0e+00*r355 ;
            repd[25] = 1.0e+00/15.0e+00*r155 + 18.0e+00/245.0e+00*r355 ;
            repd[26] = (-s3/15.0e+00*r155) + 12.0e+00*s3/245.0e+00*r355 ;
            repd[27] = (-s3/15.0e+00*r155) - 3.0e+00*s3/245.0e+00*r355 ;
            repd[28] = -repd[27] ;
            repd[29] = r066 + 4.0e+00/49.0e+00*r266 + 4.0e+00/49.0e+00*r466 ;
            repd[30] = r066 + 2.0e+00/49.0e+00*r266 - 24.0e+00/441.0e+00*r466 ;
            repd[31] = r066 - 4.0e+00/49.0e+00*r266 + 6.0e+00/441.0e+00*r466 ;
            repd[32] = sqrt(3.0e+00/245.0e+00)*r246 ;
            repd[33] = 1.0e+00/5.0e+00*r155 + 24.0e+00/245.0e+00*r355 ;
            repd[34] = 1.0e+00/5.0e+00*r155 - 6.0e+00/245.0e+00*r355 ;
            repd[35] = 3.0e+00/49.0e+00*r355 ;
            repd[36] = 1.0e+00/49.0e+00*r266 + 30.0e+00/441.0e+00*r466 ;
            repd[37] = s3/49.0e+00*r266 - 5.0e+00*s3/441.0e+00*r466 ;
            repd[38] = r066 - 2.0e+00/49.0e+00*r266 - 4.0e+00/441.0e+00*r466 ;
            repd[39] = (-2.0e+00*s3/49.0e+00*r266) + 10.0e+00*s3/441.0e+00*r466 ;
            repd[40] = -repd[32] ;
            repd[41] = -repd[34] ;
            repd[42] = -repd[35] ;
            repd[43] = -repd[37] ;
            repd[44] = 3.0e+00/49.0e+00*r266 + 20.0e+00/441.0e+00*r466 ;
            repd[45] = -repd[39] ;
            repd[46] = 1.0e+00/5.0e+00*r155 - 3.0e+00/35.0e+00*r355 ;
            repd[47] = -repd[46] ;
            repd[48] = 4.0e+00/49.0e+00*r266 + 15.0e+00/441.0e+00*r466 ;
            repd[49] = 3.0e+00/49.0e+00*r266 - 5.0e+00/147.0e+00*r466 ;
            repd[50] = -repd[49] ;
            repd[51] = r066 + 4.0e+00/49.0e+00*r266 - 34.0e+00/441.0e+00*r466 ;
            repd[52] = 35.0e+00/441.0e+00*r466 ;

            /* . Save the resulting values. */
            values[n] = repd[ 1] ; indices[t] = DX2Y2 ; indices[t+1] = DX2Y2 ; indices[t+2] =     S ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 1] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] =     S ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 1] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =     S ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 1] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =     S ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 1] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =     S ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 2] ; indices[t] =   DZ2 ; indices[t+1] =    PZ ; indices[t+2] =    PZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] = DX2Y2 ; indices[t+1] =    PX ; indices[t+2] =    PX ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] =   DXZ ; indices[t+1] =    PX ; indices[t+2] =    PZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] =   DXZ ; indices[t+1] =    PZ ; indices[t+2] =    PX ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] =   DYZ ; indices[t+1] =    PY ; indices[t+2] =    PZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] =   DYZ ; indices[t+1] =    PZ ; indices[t+2] =    PY ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] =   DXY ; indices[t+1] =    PX ; indices[t+2] =    PY ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 3] ; indices[t] =   DXY ; indices[t+1] =    PY ; indices[t+2] =    PX ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 4] ; indices[t] =   DZ2 ; indices[t+1] =     S ; indices[t+2] =    PZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[ 5] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =    PZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] = DX2Y2 ; indices[t+1] = DX2Y2 ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] = DX2Y2 ; indices[t+1] = DX2Y2 ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] =    PZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =    PZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[ 6] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[ 7] ; indices[t] = DX2Y2 ; indices[t+1] = DX2Y2 ; indices[t+2] =    PZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[ 7] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[ 7] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[ 7] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =    PZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[ 8] ; indices[t] =   DZ2 ; indices[t+1] =    PX ; indices[t+2] =    PX ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 8] ; indices[t] =   DZ2 ; indices[t+1] =    PY ; indices[t+2] =    PY ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[ 9] ; indices[t] = DX2Y2 ; indices[t+1] =     S ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[ 9] ; indices[t] =   DXZ ; indices[t+1] =     S ; indices[t+2] =    PZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[ 9] ; indices[t] =   DYZ ; indices[t+1] =     S ; indices[t+2] =    PZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[ 9] ; indices[t] =   DXY ; indices[t+1] =     S ; indices[t+2] =    PY ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[10] ; indices[t] =   DZ2 ; indices[t+1] =   DXZ ; indices[t+2] =    PZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[10] ; indices[t] =   DYZ ; indices[t+1] =   DZ2 ; indices[t+2] =    PZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[11] ; indices[t] =   DXZ ; indices[t+1] = DX2Y2 ; indices[t+2] =    PZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[11] ; indices[t] =   DYZ ; indices[t+1] =   DXZ ; indices[t+2] =    PY ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[11] ; indices[t] =   DXY ; indices[t+1] =   DXZ ; indices[t+2] =    PZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[11] ; indices[t] =   DXY ; indices[t+1] =   DYZ ; indices[t+2] =    PZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[12] ; indices[t] =   DZ2 ; indices[t+1] =     S ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[12] ; indices[t] =   DZ2 ; indices[t+1] =     S ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[13] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[13] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[14] ; indices[t] =   DZ2 ; indices[t+1] = DX2Y2 ; indices[t+2] =    PX ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[14] ; indices[t] =   DXY ; indices[t+1] =   DZ2 ; indices[t+2] =    PY ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[15] ; indices[t] = DX2Y2 ; indices[t+1] =    PY ; indices[t+2] =    PY ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[16] ; indices[t] =   DYZ ; indices[t+1] = DX2Y2 ; indices[t+2] =    PZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[17] ; indices[t] = DX2Y2 ; indices[t+1] =     S ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[18] ; indices[t] =   DZ2 ; indices[t+1] = DX2Y2 ; indices[t+2] =    PY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[19] ; indices[t] = DX2Y2 ; indices[t+1] =     S ; indices[t+2] = DX2Y2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[19] ; indices[t] =   DXZ ; indices[t+1] =     S ; indices[t+2] =   DXZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[19] ; indices[t] =   DZ2 ; indices[t+1] =     S ; indices[t+2] =   DZ2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[19] ; indices[t] =   DYZ ; indices[t+1] =     S ; indices[t+2] =   DYZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[19] ; indices[t] =   DXY ; indices[t+1] =     S ; indices[t+2] =   DXY ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[20] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =   DZ2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[21] ; indices[t] =   DZ2 ; indices[t+1] =     S ; indices[t+2] =   DXZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[21] ; indices[t] =   DZ2 ; indices[t+1] =   DXZ ; indices[t+2] =   DXZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[21] ; indices[t] =   DYZ ; indices[t+1] =   DZ2 ; indices[t+2] =   DYZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[21] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =   DZ2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[22] ; indices[t] =   DZ2 ; indices[t+1] =     S ; indices[t+2] = DX2Y2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[22] ; indices[t] =   DZ2 ; indices[t+1] = DX2Y2 ; indices[t+2] = DX2Y2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[22] ; indices[t] =   DXY ; indices[t+1] =   DZ2 ; indices[t+2] =   DXY ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[22] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =   DZ2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[23] ; indices[t] =   DZ2 ; indices[t+1] =    PZ ; indices[t+2] =   DZ2 ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[24] ; indices[t] =   DZ2 ; indices[t+1] =    PZ ; indices[t+2] =   DXZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[24] ; indices[t] =   DYZ ; indices[t+1] =    PY ; indices[t+2] =   DZ2 ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[25] ; indices[t] =   DZ2 ; indices[t+1] =    PX ; indices[t+2] =   DZ2 ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[25] ; indices[t] =   DZ2 ; indices[t+1] =    PY ; indices[t+2] =   DZ2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[26] ; indices[t] =   DZ2 ; indices[t+1] =    PX ; indices[t+2] =   DXZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[26] ; indices[t] =   DYZ ; indices[t+1] =    PZ ; indices[t+2] =   DZ2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[27] ; indices[t] =   DZ2 ; indices[t+1] =    PX ; indices[t+2] = DX2Y2 ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[27] ; indices[t] =   DXY ; indices[t+1] =    PX ; indices[t+2] =   DZ2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[27] ; indices[t] =   DXY ; indices[t+1] =    PY ; indices[t+2] =   DZ2 ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[28] ; indices[t] =   DZ2 ; indices[t+1] =    PY ; indices[t+2] = DX2Y2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[29] ; indices[t] = DX2Y2 ; indices[t+1] = DX2Y2 ; indices[t+2] = DX2Y2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[29] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] =   DXZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[29] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =   DZ2 ; indices[t+3] =   DZ2 ; n += 1 ; t += 4 ;
            values[n] = repd[29] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =   DYZ ; indices[t+3] =   DYZ ; n += 1 ; t += 4 ;
            values[n] = repd[29] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =   DXY ; indices[t+3] =   DXY ; n += 1 ; t += 4 ;
            values[n] = repd[30] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] =   DXZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[30] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =   DZ2 ; indices[t+3] =   DZ2 ; n += 1 ; t += 4 ;
            values[n] = repd[31] ; indices[t] =   DZ2 ; indices[t+1] =   DZ2 ; indices[t+2] = DX2Y2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[31] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =   DZ2 ; indices[t+3] =   DZ2 ; n += 1 ; t += 4 ;
            values[n] = repd[32] ; indices[t] =   DXZ ; indices[t+1] = DX2Y2 ; indices[t+2] =   DXZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[32] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] = DX2Y2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[32] ; indices[t] =   DXY ; indices[t+1] =     S ; indices[t+2] =   DYZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[32] ; indices[t] =   DXY ; indices[t+1] =   DXZ ; indices[t+2] =   DYZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[32] ; indices[t] =   DXY ; indices[t+1] =   DYZ ; indices[t+2] =   DXZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] = DX2Y2 ; indices[t+1] =    PX ; indices[t+2] = DX2Y2 ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] = DX2Y2 ; indices[t+1] =    PY ; indices[t+2] = DX2Y2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] =   DXZ ; indices[t+1] =    PX ; indices[t+2] =   DXZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] =   DXZ ; indices[t+1] =    PZ ; indices[t+2] =   DXZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] =   DYZ ; indices[t+1] =    PY ; indices[t+2] =   DYZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] =   DYZ ; indices[t+1] =    PZ ; indices[t+2] =   DYZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] =   DXY ; indices[t+1] =    PX ; indices[t+2] =   DXY ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[33] ; indices[t] =   DXY ; indices[t+1] =    PY ; indices[t+2] =   DXY ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[34] ; indices[t] =   DXZ ; indices[t+1] =    PZ ; indices[t+2] = DX2Y2 ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[34] ; indices[t] =   DYZ ; indices[t+1] =    PY ; indices[t+2] =   DXZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[34] ; indices[t] =   DXY ; indices[t+1] =    PX ; indices[t+2] =   DYZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[34] ; indices[t] =   DXY ; indices[t+1] =    PY ; indices[t+2] =   DXZ ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] = DX2Y2 ; indices[t+1] =    PZ ; indices[t+2] = DX2Y2 ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DXZ ; indices[t+1] =    PX ; indices[t+2] = DX2Y2 ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DXZ ; indices[t+1] =    PY ; indices[t+2] =   DXZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DYZ ; indices[t+1] =    PX ; indices[t+2] =   DXZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DYZ ; indices[t+1] =    PX ; indices[t+2] =   DYZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DXY ; indices[t+1] =    PZ ; indices[t+2] =   DXZ ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DXY ; indices[t+1] =    PZ ; indices[t+2] =   DYZ ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[35] ; indices[t] =   DXY ; indices[t+1] =    PZ ; indices[t+2] =   DXY ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[36] ; indices[t] =   DZ2 ; indices[t+1] =   DXZ ; indices[t+2] =   DZ2 ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[36] ; indices[t] =   DYZ ; indices[t+1] =   DZ2 ; indices[t+2] =   DYZ ; indices[t+3] =   DZ2 ; n += 1 ; t += 4 ;
            values[n] = repd[37] ; indices[t] =   DZ2 ; indices[t+1] =   DXZ ; indices[t+2] =   DXZ ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[37] ; indices[t] =   DXY ; indices[t+1] =   DXZ ; indices[t+2] =   DYZ ; indices[t+3] =   DZ2 ; n += 1 ; t += 4 ;
            values[n] = repd[37] ; indices[t] =   DXY ; indices[t+1] =   DYZ ; indices[t+2] =   DZ2 ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[38] ; indices[t] =   DXZ ; indices[t+1] =   DXZ ; indices[t+2] = DX2Y2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[38] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] = DX2Y2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[38] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =   DXZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[38] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =   DXZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[38] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] =   DYZ ; indices[t+3] =   DYZ ; n += 1 ; t += 4 ;
            values[n] = repd[39] ; indices[t] =   DZ2 ; indices[t+1] = DX2Y2 ; indices[t+2] =   DXZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[39] ; indices[t] =   DXY ; indices[t+1] =   DZ2 ; indices[t+2] =   DYZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[40] ; indices[t] =   DYZ ; indices[t+1] = DX2Y2 ; indices[t+2] =   DYZ ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[40] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] = DX2Y2 ; indices[t+3] =     S ; n += 1 ; t += 4 ;
            values[n] = repd[41] ; indices[t] =   DYZ ; indices[t+1] =    PZ ; indices[t+2] = DX2Y2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[42] ; indices[t] =   DYZ ; indices[t+1] =    PY ; indices[t+2] = DX2Y2 ; indices[t+3] =    PZ ; n += 1 ; t += 4 ;
            values[n] = repd[43] ; indices[t] =   DYZ ; indices[t+1] =   DZ2 ; indices[t+2] =   DYZ ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[44] ; indices[t] =   DXZ ; indices[t+1] = DX2Y2 ; indices[t+2] =   DXZ ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[44] ; indices[t] =   DYZ ; indices[t+1] = DX2Y2 ; indices[t+2] =   DYZ ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[44] ; indices[t] =   DYZ ; indices[t+1] =   DXZ ; indices[t+2] =   DYZ ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[44] ; indices[t] =   DXY ; indices[t+1] =   DXZ ; indices[t+2] =   DXY ; indices[t+3] =   DXZ ; n += 1 ; t += 4 ;
            values[n] = repd[44] ; indices[t] =   DXY ; indices[t+1] =   DYZ ; indices[t+2] =   DXY ; indices[t+3] =   DYZ ; n += 1 ; t += 4 ;
            values[n] = repd[45] ; indices[t] =   DYZ ; indices[t+1] =   DYZ ; indices[t+2] =   DZ2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[46] ; indices[t] =   DXY ; indices[t+1] =    PY ; indices[t+2] = DX2Y2 ; indices[t+3] =    PX ; n += 1 ; t += 4 ;
            values[n] = repd[47] ; indices[t] =   DXY ; indices[t+1] =    PX ; indices[t+2] = DX2Y2 ; indices[t+3] =    PY ; n += 1 ; t += 4 ;
            values[n] = repd[48] ; indices[t] =   DZ2 ; indices[t+1] = DX2Y2 ; indices[t+2] =   DZ2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[48] ; indices[t] =   DXY ; indices[t+1] =   DZ2 ; indices[t+2] =   DXY ; indices[t+3] =   DZ2 ; n += 1 ; t += 4 ;
            values[n] = repd[49] ; indices[t] =   DXY ; indices[t+1] =   DYZ ; indices[t+2] =   DXZ ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[50] ; indices[t] =   DXY ; indices[t+1] =   DXZ ; indices[t+2] =   DYZ ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[51] ; indices[t] =   DXY ; indices[t+1] =   DXY ; indices[t+2] = DX2Y2 ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;
            values[n] = repd[52] ; indices[t] =   DXY ; indices[t+1] = DX2Y2 ; indices[t+2] =   DXY ; indices[t+3] = DX2Y2 ; n += 1 ; t += 4 ;

            /* . Save the Slater-Condon parameters. */
/*
            self->f0dd = r066 ;
            self->f2dd = r266 ;
            self->f4dd = r466 ;
*/
            self->f0sd = r016 ;
            self->g2sd = r244 ;
/*
            self->f0pd = r036 ;
            self->f2pd = r236 ;
            self->g1pd = r155 ;
            self->g3pd = r355 ;
*/

# ifdef PRINTMOPACPARAMETERS
{
printf ( "MOPACPARAMETERS> %d   F0DD    %.15f   ATOMIC\n", self->atomicNumber, r066 ) ;
printf ( "MOPACPARAMETERS> %d   F2DD    %.15f   ATOMIC\n", self->atomicNumber, r266 ) ;
printf ( "MOPACPARAMETERS> %d   F4DD    %.15f   ATOMIC\n", self->atomicNumber, r466 ) ;
printf ( "MOPACPARAMETERS> %d   F0SD    %.15f   ATOMIC\n", self->atomicNumber, r016 ) ;
printf ( "MOPACPARAMETERS> %d   G2SD    %.15f   ATOMIC\n", self->atomicNumber, r244 ) ;
printf ( "MOPACPARAMETERS> %d   F0PD    %.15f   ATOMIC\n", self->atomicNumber, r036 ) ;
printf ( "MOPACPARAMETERS> %d   F2PD    %.15f   ATOMIC\n", self->atomicNumber, r236 ) ;
printf ( "MOPACPARAMETERS> %d   G1PD    %.15f   ATOMIC\n", self->atomicNumber, r155 ) ;
printf ( "MOPACPARAMETERS> %d   G3PD    %.15f   ATOMIC\n", self->atomicNumber, r355 ) ;
}
# endif

        }

        /* . Calculate the terms required for the calculation of the two-center two-electron integrals. */
        MNDOParameters_FillDDPPO ( self, repd ) ;

# ifdef PRINTMOPACPARAMETERS
{
auto Integer  i ;
for ( i = 0 ; i < 6 ; i++ ) printf ( "MOPACPARAMETERS> %d   DD%d     %.15f   ATOMIC\n", self->atomicNumber, i+1, self->ddp[i] ) ;
for ( i = 0 ; i < 9 ; i++ ) printf ( "MOPACPARAMETERS> %d   PO%d     %.15f   ATOMIC\n", self->atomicNumber, i+1, self->po[i] ) ;
printf ( "MOPACPARAMETERS> %d   EISOL   %.15f   ATOMIC\n", self->atomicNumber, self->eisol ) ;
}
# endif

# ifdef PRINTMOPACOEIS
{
/* . Printing. */
printf ( "\n\nONE-CENTER TEIS:\n" ) ;
printf ( "%d\n", self->atomicNumber ) ;
for ( t = 0 ; t < n ; t++ ) printf ( "%5d %5d %5d %5d %5d %20.10f\n", t, indices[4*t], indices[4*t+1], indices[4*t+2], indices[4*t+3], values[t] ) ;
}
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOParameters *MNDOParameters_Clone ( const MNDOParameters *self )
{
    MNDOParameters *new = NULL ;
    if ( self != NULL )
    {
        auto Integer  i ;
        new = MNDOParameters_Allocate ( ) ;
        /* . Counters. */
        new->QDIATOMIC    = self->QDIATOMIC    ;
        new->atomicNumber = self->atomicNumber ;
        new->iii          = self->iii          ;
        new->iiid         = self->iiid         ;
        new->ir016        = self->ir016        ;
        new->ir066        = self->ir066        ;
        new->ir244        = self->ir244        ;
        new->ir266        = self->ir266        ;
        new->ir466        = self->ir466        ;
        new->nam1pm3g     = self->nam1pm3g     ;
        new->ndiatomic    = self->ndiatomic    ;
        new->norbitals    = self->norbitals    ;
        new->npddg        = self->npddg        ;
        new->qnd          = self->qnd          ;
        new->qnp          = self->qnp          ;
        new->qns          = self->qns          ;
        /* . Derived quantities. */
        new->nocteis      = self->nocteis      ;
        new->hpp          = self->hpp          ;
        /* . Scalars. */
        /* . Input parameters. */
        new->ad0    = self->ad0    ;
        new->alp0   = self->alp0   ;
        new->am0    = self->am0    ;
        new->aq0    = self->aq0    ;
        new->betad0 = self->betad0 ;
        new->betap0 = self->betap0 ;
        new->betas0 = self->betas0 ;
        new->dd0    = self->dd0    ;
        new->eheat0 = self->eheat0 ;
        new->eisol0 = self->eisol0 ;
        new->f0sd0  = self->f0sd0  ;
	new->gphot0 = self->gphot0 ;
        new->gpp0   = self->gpp0   ;
        new->gp20   = self->gp20   ;
        new->gsp0   = self->gsp0   ;
        new->gss0   = self->gss0   ;
        new->g2sd0  = self->g2sd0  ;
        new->hsp0   = self->hsp0   ;
        new->pcore0 = self->pcore0 ;
        new->qq0    = self->qq0    ;
        new->udd0   = self->udd0   ;
        new->upp0   = self->upp0   ;
        new->uss0   = self->uss0   ;
        new->zcore0 = self->zcore0 ;
        new->zetad0 = self->zetad0 ;
        new->zetap0 = self->zetap0 ;
        new->zetas0 = self->zetas0 ;
        new->zdn0   = self->zdn0   ;
        new->zpn0   = self->zpn0   ;
        new->zsn0   = self->zsn0   ;
        /* . Internal parameters. */
        new->ad    = self->ad    ;
        new->alp   = self->alp   ;
        new->am    = self->am    ;
        new->aq    = self->aq    ;
        new->betad = self->betad ;
        new->betap = self->betap ;
        new->betas = self->betas ;
        new->dd    = self->dd    ;
        new->eheat = self->eheat ;
        new->eisol = self->eisol ;
        new->f0sd  = self->f0sd  ;
	new->gphot = self->gphot ;
        new->gpp   = self->gpp   ;
        new->gp2   = self->gp2   ;
        new->gsp   = self->gsp   ;
        new->gss   = self->gss   ;
        new->g2sd  = self->g2sd  ;
        new->hsp   = self->hsp   ;
        new->pcore = self->pcore ;
        new->qq    = self->qq    ;
        new->udd   = self->udd   ;
        new->upp   = self->upp   ;
        new->uss   = self->uss   ;
        new->zcore = self->zcore ;
        new->zetad = self->zetad ;
        new->zetap = self->zetap ;
        new->zetas = self->zetas ;
        new->zdn   = self->zdn   ;
        new->zpn   = self->zpn   ;
        new->zsn   = self->zsn   ;
        /* . Arrays. */
        for ( i = 0 ; i < 6 ; i++ ) new->ddp[i] = self->ddp[i] ;
        for ( i = 0 ; i < 9 ; i++ ) new->po [i] = self->po [i] ;
        if ( new->nocteis > 0 )
        {
            new->octeiindices = Memory_AllocateArrayOfTypes ( 4 * new->nocteis, Cardinal16  ) ;
            new->octeivalues  = Memory_AllocateArrayOfTypes (     new->nocteis, Real       ) ;
            for ( i = 0 ; i < 4 * new->nocteis ; i++ ) new->octeiindices[i] = self->octeiindices[i] ;
            for ( i = 0 ; i <     new->nocteis ; i++ ) new->octeivalues [i] = self->octeivalues [i] ;
        }
        if ( new->norbitals > 0 )
        {
            new->beta0 = Memory_AllocateArrayOfTypes ( new->norbitals, Real ) ;
            new->uspd0 = Memory_AllocateArrayOfTypes ( new->norbitals, Real ) ;
            new->beta  = Memory_AllocateArrayOfTypes ( new->norbitals, Real ) ;
            new->uspd  = Memory_AllocateArrayOfTypes ( new->norbitals, Real ) ;
            for ( i = 0 ; i < new->norbitals ; i++ )
            {
                new->beta0 [i] = self->beta0 [i] ;
                new->uspd0 [i] = self->uspd0 [i] ;
                new->beta  [i] = self->beta  [i] ;
                new->uspd  [i] = self->uspd  [i] ;
            }
        }
        if ( new->nam1pm3g > 0 )
        {
            new->fn10 = Memory_AllocateArrayOfTypes ( new->nam1pm3g, Real ) ;
            new->fn20 = Memory_AllocateArrayOfTypes ( new->nam1pm3g, Real ) ;
            new->fn30 = Memory_AllocateArrayOfTypes ( new->nam1pm3g, Real ) ;
            new->fn1  = Memory_AllocateArrayOfTypes ( new->nam1pm3g, Real ) ;
            new->fn2  = Memory_AllocateArrayOfTypes ( new->nam1pm3g, Real ) ;
            new->fn3  = Memory_AllocateArrayOfTypes ( new->nam1pm3g, Real ) ;
            for ( i = 0 ; i < new->nam1pm3g ; i++ )
            {
                new->fn10 [i] = self->fn10 [i] ;
                new->fn20 [i] = self->fn20 [i] ;
                new->fn30 [i] = self->fn30 [i] ;
                new->fn1  [i] = self->fn1  [i] ;
                new->fn2  [i] = self->fn2  [i] ;
                new->fn3  [i] = self->fn3  [i] ;
            }
        }
        if ( new->ndiatomic > 0 )
        {
            new->QDIATOMICFLAGS = Memory_AllocateArrayOfTypes ( new->ndiatomic, Boolean ) ;
            new->diatomica0     = Memory_AllocateArrayOfTypes ( new->ndiatomic, Real    ) ;
            new->diatomicx0     = Memory_AllocateArrayOfTypes ( new->ndiatomic, Real    ) ;
            new->diatomica      = Memory_AllocateArrayOfTypes ( new->ndiatomic, Real    ) ;
            new->diatomicx      = Memory_AllocateArrayOfTypes ( new->ndiatomic, Real    ) ;
            for ( i = 0 ; i < new->ndiatomic ; i++ )
            {
                new->QDIATOMICFLAGS [i] = self->QDIATOMICFLAGS [i] ;
                new->diatomica0     [i] = self->diatomica0     [i] ;
                new->diatomicx0     [i] = self->diatomicx0     [i] ;
                new->diatomica      [i] = self->diatomica      [i] ;
                new->diatomicx      [i] = self->diatomicx      [i] ;
            }
        }
        if ( new->npddg > 0 )
        {
            new->pddgc0 = Memory_AllocateArrayOfTypes ( new->npddg, Real ) ;
            new->pddge0 = Memory_AllocateArrayOfTypes ( new->npddg, Real ) ;
            new->pddgc  = Memory_AllocateArrayOfTypes ( new->npddg, Real ) ;
            new->pddge  = Memory_AllocateArrayOfTypes ( new->npddg, Real ) ;
            for ( i = 0 ; i < new->npddg ; i++ )
            {
                new->pddgc0 [i] = self->pddgc0 [i] ;
                new->pddge0 [i] = self->pddge0 [i] ;
                new->pddgc  [i] = self->pddgc  [i] ;
                new->pddge  [i] = self->pddge  [i] ;
            }
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOParameters_Deallocate ( MNDOParameters **self )
{
    if ( (*self) != NULL )
    {
        MNDOParameters_DeallocateAtomicUnitArrays ( *self ) ;
        Memory_Deallocate ( (*self)->QDIATOMICFLAGS ) ;
        Memory_Deallocate ( (*self)->octeiindices   ) ;
        Memory_Deallocate ( (*self)->octeivalues    ) ;
        Memory_Deallocate ( (*self)->beta0          ) ;
        Memory_Deallocate ( (*self)->diatomica0     ) ;
        Memory_Deallocate ( (*self)->diatomicx0     ) ;
        Memory_Deallocate ( (*self)->fn10           ) ;
        Memory_Deallocate ( (*self)->fn20           ) ;
        Memory_Deallocate ( (*self)->fn30           ) ;
        Memory_Deallocate ( (*self)->pddgc0         ) ;
        Memory_Deallocate ( (*self)->pddge0         ) ;
        Memory_Deallocate ( (*self)->uspd0          ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOParameters_DeallocateAtomicUnitArrays ( MNDOParameters *self )
{
    if ( self != NULL )
    {
        Memory_Deallocate ( self->beta      ) ;
        Memory_Deallocate ( self->diatomica ) ;
        Memory_Deallocate ( self->diatomicx ) ;
        Memory_Deallocate ( self->fn1       ) ;
        Memory_Deallocate ( self->fn2       ) ;
        Memory_Deallocate ( self->fn3       ) ;
        Memory_Deallocate ( self->pddgc     ) ;
        Memory_Deallocate ( self->pddge     ) ;
        Memory_Deallocate ( self->uspd      ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fill beta and uspd arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOParameters_FillBetaUspd ( MNDOParameters *self )
{
    if ( self != NULL )
    {
        if ( self->norbitals > 0 )
        {
            auto Integer  i ;
            self->beta0 = Memory_AllocateArrayOfTypes ( self->norbitals, Real ) ;
            self->uspd0 = Memory_AllocateArrayOfTypes ( self->norbitals, Real ) ;
            /* . s. */
            self->beta0[0] = self->betas0 ;
            self->uspd0[0] = self->uss0   ;
            /* . p. */
            if ( self->norbitals >= 4 )
            {
                for ( i = 1 ; i < 4 ; i++ )
                {
                    self->beta0[i] = self->betap0 ;
                    self->uspd0[i] = self->upp0   ;
                }
            }
            /* . d. */
            if ( self->norbitals >= 9 )
            {
                for ( i = 4 ; i < 9 ; i++ )
                {
                    self->beta0[i] = self->betad0 ;
                    self->uspd0[i] = self->udd0   ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert to atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _BohrsToAngstroms            ( Units_Length_Bohrs_To_Angstroms          )
# define _MNDOHartreesToElectronVolts ( 27.2113834e+00                           ) /* . New value. */
# define _MNDOHartreesToKCal          ( _MNDOHartreesToElectronVolts * 23.060529 ) /* . New value. */
void MNDOParameters_ToAtomicUnits ( MNDOParameters *self )
{
    if ( self != NULL )
    {
        auto Integer  i ;
        MNDOParameters_DeallocateAtomicUnitArrays ( self ) ;
        /* . ad, am, aq, dd, qq, zcore, zetap and zetas are already in atomic units. */
        self->ad    = self->ad0    ;
        self->alp   = self->alp0   * _BohrsToAngstroms            ; /* . A^-1. */
        self->am    = self->am0    ;
        self->aq    = self->aq0    ;
        self->betad = self->betad0 / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->betap = self->betap0 / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->betas = self->betas0 / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->dd    = self->dd0    ;
        self->eheat = self->eheat0 / _MNDOHartreesToKCal          ; /* . kcal mol^-1. */
        self->eisol = self->eisol0 / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->f0sd  = self->f0sd0  / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->gphot = self->gphot0 ;
        self->gpp   = self->gpp0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->gp2   = self->gp20   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->gsp   = self->gsp0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->gss   = self->gss0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->g2sd  = self->g2sd0  / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->hsp   = self->hsp0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->pcore = self->pcore0 ;
        self->qq    = self->qq0    ;
        self->udd   = self->udd0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->upp   = self->upp0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->uss   = self->uss0   / _MNDOHartreesToElectronVolts ; /* . eV. */
        self->zcore = self->zcore0 ;
        self->zetad = self->zetad0 ;
        self->zetap = self->zetap0 ;
        self->zetas = self->zetas0 ;
        self->zdn   = self->zdn0   ;
        self->zpn   = self->zpn0   ;
        self->zsn   = self->zsn0   ;
        if ( self->norbitals > 0 )
        {
            self->beta = Memory_AllocateArrayOfTypes ( self->norbitals, Real ) ;
            self->uspd = Memory_AllocateArrayOfTypes ( self->norbitals, Real ) ;
            for ( i = 0 ; i < self->norbitals ; i++ )
            {
                self->beta[i] = self->beta0[i] / _MNDOHartreesToElectronVolts ; /* . eV. */
                self->uspd[i] = self->uspd0[i] / _MNDOHartreesToElectronVolts ; /* . eV. */
            }
        }
        if ( self->nam1pm3g > 0 )
        {
            self->fn1 = Memory_AllocateArrayOfTypes ( self->nam1pm3g, Real ) ;
            self->fn2 = Memory_AllocateArrayOfTypes ( self->nam1pm3g, Real ) ;
            self->fn3 = Memory_AllocateArrayOfTypes ( self->nam1pm3g, Real ) ;
            for ( i = 0 ; i < self->nam1pm3g ; i++ )
            {
                self->fn1[i] = self->fn10[i] / ( _BohrsToAngstroms * _MNDOHartreesToElectronVolts ) ; /* . A eV. */
                self->fn2[i] = self->fn20[i] * _BohrsToAngstroms * _BohrsToAngstroms                ; /* . A^-2. */
                self->fn3[i] = self->fn30[i] / _BohrsToAngstroms                                    ; /* . A. */
            }
        }
        if ( self->ndiatomic > 0 )
        {
            self->diatomica = Memory_AllocateArrayOfTypes ( self->ndiatomic, Real ) ;
            self->diatomicx = Memory_AllocateArrayOfTypes ( self->ndiatomic, Real ) ;
            for ( i = 0 ; i < self->ndiatomic ; i++ )
            {
                self->diatomica[i] = self->diatomica0[i] * _BohrsToAngstroms ; /* . A^-1. */
                self->diatomicx[i] = self->diatomicx0[i]                     ; /* . Dimensionless. */
            }
        }
        if ( self->npddg > 0 )
        {
            self->pddgc = Memory_AllocateArrayOfTypes ( self->npddg, Real ) ;
            self->pddge = Memory_AllocateArrayOfTypes ( self->npddg, Real ) ;
            for ( i = 0 ; i < self->npddg ; i++ )
            {
                self->pddgc[i] = self->pddgc0[i] / _MNDOHartreesToElectronVolts ; /* . eV. */
                self->pddge[i] = self->pddge0[i] / _BohrsToAngstroms            ; /* . A. */
            }
        }
    }
}
# undef _BohrsToAngstroms
# undef _MNDOHartreesToElectronVolts
# undef _MNDOHartreesToKCal

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate AIJL.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . There is a disagreement here with the paper by addition of an extra factor of 2^l. */
static Real AIJL ( const Real z1, const Real z2, const Integer  n1, const Integer  n2, const Integer  l )
{
    auto Real value ;
    value = Factorial ( n1 + n2 + l ) / sqrt ( Factorial ( 2 * n1 ) * Factorial ( 2 * n2 ) ) *
            pow ( ( 2.0e+00 * z1 / ( z1 + z2 ) ), n1 ) * sqrt ( 2.0e+00 * z1 / ( z1 + z2 ) ) *
            pow ( ( 2.0e+00 * z2 / ( z1 + z2 ) ), n2 ) * sqrt ( 2.0e+00 * z2 / ( z1 + z2 ) ) * pow ( 2.0e+00, l ) / pow ( ( z1 + z2 ), l ) ;
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Binomial.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Binomial ( const Integer  n, const Integer  k )
{
    if ( n >= k ) return Factorial ( n ) / ( Factorial ( k ) * Factorial ( n - k ) ) ;
    else          return 0.0e+00 ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Factorial.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Factorial ( const Integer  n )
{
    Integer  i ;
    Real     value = 1.0e+00 ;
    for ( i = 1 ; i <= n ; i++ ) value *= ( Real ) i ;
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize ddp and po.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
  Compared to MNDOD, there is a disagreement in the definition of AIJL and in the ddp/DD values.
  The conversions are:

  ddp[1] = DD(2)
  ddp[2] = DD(3) * sqrt ( 2 )
  ddp[3] = DD(4)
  ddp[4] = DD(5)
  ddp[5] = DD(6)

  Originally ddp[3] was DD(4) / sqrt ( 2 ) but this was changed for integral evaluation. DD(4) was always needed for POIJ.

*/
static void MNDOParameters_FillDDPPO ( MNDOParameters *self, const Real *octeis )
{
    if ( self != NULL )
    {
        auto Boolean  hasDOrbitals ;
        auto Integer  i, ni ;
        auto Real     aij[6], d, fg ;
        /* . Initialization. */
        hasDOrbitals = ( self->norbitals == 9 ) ;
        ni           =   self->atomicNumber ;
        for ( i = 0 ; i < 6 ; i++ ) aij[i]       = 0.0e+00 ;
        for ( i = 0 ; i < 6 ; i++ ) self->ddp[i] = 0.0e+00 ;
        for ( i = 0 ; i < 9 ; i++ ) self->po [i] = 0.0e+00 ;
        /* . Calculate aij. */
        /* . 0 ss, 1 sp, 2 pp, 3 sd, 4 pd, 5 dd. */
        if ( ni > 2 )
        {
            auto Integer  nsp ;
            auto Real     z1, z2, z3, zz ;
            z1  = self->zetas ;
            z2  = self->zetap ;
            z3  = self->zetad ;
            nsp = self->iii   ;
            zz  = z1 * z2 ;
            if ( zz >= 1.0e-2 )
            {
                aij[1] = AIJL ( z1, z2, nsp, nsp, 1 ) ;
                aij[2] = AIJL ( z2, z2, nsp, nsp, 2 ) ;
                if ( hasDOrbitals )
                {
                    auto Integer  nd ;
                    nd = self->iiid ;
                    aij[3] = AIJL ( z1, z3, nsp, nd, 2 ) ;
                    aij[4] = AIJL ( z2, z3, nsp, nd, 1 ) ;
                    aij[5] = AIJL ( z3, z3, nd,  nd, 2 ) ;
                }
            }
        }
        /* . Calculate ddp and po. */
        /* . ss. */
        fg = self->gss ;
        self->po[0] = POIJ ( 0, 1.0e+00, fg ) ;
        self->po[8] = self->po[0] ;
        if ( ni > 2 )
        {
            /* . sp. */
            d  = aij[1] / sqrt ( 12.0e+00 ) ;
            fg = self->hsp ;
            self->ddp[1] = d ;
            self->po [1] = POIJ ( 1, d, fg ) ;
            /* . pp. */
            self->po[6] = self->po[0] ;
            d  = sqrt ( aij[2] / 10.0e+00 ) ;
            fg = self->hpp ;
            self->ddp[2] = d ; /* . DD(3) * sqrt ( 2 ) in MNDOD. */
            self->po [2] = POIJ ( 2, d, fg ) ;
            if ( hasDOrbitals )
            {
                /* . sd. */
                d  = sqrt ( aij[3] / 30.0e+00 ) ;
                fg = octeis[19] ;
                self->ddp[3] = d ;
                self->po [3] = POIJ ( 2, d, fg ) ;
                /* . pd. */
                d  = aij[4] / sqrt ( 20.0e+00 ) ;
                fg = octeis[23] - 1.8e+00 * octeis[35] ;
                self->ddp[4] = d ;
                self->po [4] = POIJ ( 1, d, fg ) ;
                /* . dd. */
                fg = 0.2e+00 * ( octeis[29] + 2.0e+00 * octeis[30] + 2.0e+00 * octeis[31] ) ;
                self->po[7] = POIJ ( 0, 1.0e+00, fg ) ;
                d  = sqrt ( aij[5] / 14.0e+00 ) ;
                fg = octeis[44] - ( 20.0e+00 / 35.0e+00 ) * octeis[52] ;
                self->ddp[5] = d ;
                self->po [5] = POIJ ( 2, d, fg ) ;
            }
        }

        /* . Non-d-orbital elements. */
        if ( ! hasDOrbitals )
        {
            if ( self->am < 1.0e-4 ) self->am = 1.0e+00 ;
            self->po[0] = 0.5e+00 / self->am ;
            if ( self->ad > 1.0e-5 ) self->po[1] = 0.5e+00 / self->ad ;
            if ( self->aq > 1.0e-5 ) self->po[2] = 0.5e+00 / self->aq ;
            self->po[6] = self->po[0] ;
            self->po[8] = self->po[0] ;
            self->ddp[1] = self->dd ;
            self->ddp[2] = self->qq * sqrt ( 2.0e+00 ) ;
        }
        /* . A core term has been specified so use it instead of po[0]. */
/*        if ( useDOrbitals && ( self->pcore > 1.0e-5 ) ) self->po[8] = self->pcore ; */
        if ( self->pcore > 1.0e-5 ) self->po[8] = self->pcore ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate POIJ.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EPSILON      1.0e-08
# define G1           0.382e+00
# define G2           0.618e+00
# define NITERATIONS  100
static Real POIJ ( const Integer  l, const Real d, const Real fg )
{
    Real value = 0.0e+00 ;
    if ( l == 0 )
    {
        value = 0.5e+00 / fg ;
    }
    else
    {
        auto Boolean  QOK = False ;
        auto Integer  i ;
        auto Real     a1, a2, delta, dsq, ev4, ev8, f1 = 0.0e+00, f2 = 0.0e+00, y1, y2 ;
        dsq = d * d ;
        ev4 = 1.0e+00 / 4.0e+00 ;
        ev8 = 1.0e+00 / 8.0e+00 ;
        a1  = 1.0e-01 ;
        a2  = 5.0e+00 ;
        for ( i = 0 ; i < NITERATIONS ; i++ )
        {
            delta = a2 - a1 ;
            if ( delta < EPSILON ) { QOK = True ; break ; }
            y1 = a1 + delta * G1 ;
            y2 = a1 + delta * G2 ;
            if ( l == 1 )
            {
                f1 = pow ( ( ev4 * ( 1.0e+00 / y1 - 1.0e+00 / sqrt ( y1*y1 + dsq ) ) - fg ), 2 ) ;
                f2 = pow ( ( ev4 * ( 1.0e+00 / y2 - 1.0e+00 / sqrt ( y2*y2 + dsq ) ) - fg ), 2 ) ;
            }
            else if ( l == 2 )
            {
                f1 = pow ( ( ev8 * ( 1.0e+00 / y1 - 2.0e+00 / sqrt ( y1*y1 + dsq*0.5e+00 ) + 1.0e+00 / sqrt ( y1*y1 + dsq ) ) - fg ), 2 ) ;
                f2 = pow ( ( ev8 * ( 1.0e+00 / y2 - 2.0e+00 / sqrt ( y2*y2 + dsq*0.5e+00 ) + 1.0e+00 / sqrt ( y2*y2 + dsq ) ) - fg ), 2 ) ;
            }
            if ( f1 < f2 ) a2 = y2 ;
            else           a1 = y1 ;
        }
        if ( ! QOK ) printf ( "\nNon-convergence in POIJ.\n" ) ;
        if ( f1 >= f2 ) value = a2 ;
        else            value = a1 ;
    }
    return value ;
}
# undef EPSILON
# undef G1
# undef G2
# undef NITERATIONS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the radial part of the Slater-Condon parameter.
! . k     - type of integral, can be equal to 0,1,2,3,4 in spd-basis.
! . na,nb - principle quantum number of ao,corresponding electron 1.
! . ea,eb - exponents of ao,corresponding electron 1.
! . nc,nd - principle quantum number of ao,corresponding electron 2.
! . ec,ed - exponents of ao,corresponding electron 2.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real RadialSlaterCondonParameter ( const Integer  k, const Integer  na, const Real ea, const Integer  nb, const Real eb, const Integer  nc, const Real ec, const Integer  nd, const Real ed )
{
    Integer  nab, ncd, n, m, i, m2 ;
    Real     aea, aeb, aec, aed, ecd, eab, e, ae, a2, acd, aab, ff, c, s0, s1, s2, s3 ;
    aea = log ( ea ) ;
    aeb = log ( eb ) ;
    aec = log ( ec ) ;
    aed = log ( ed ) ;
    nab = na + nb ;
    ncd = nc + nd ;
    ecd = ec + ed ;
    eab = ea + eb ;
    e   = ecd + eab ;
    n   = nab + ncd ;
    ae  = log ( e ) ;
    a2  = log ( 2.0e+00 ) ;
    acd = log ( ecd ) ;
    aab = log ( eab ) ;
    ff  = Factorial ( n-1 ) / sqrt ( Factorial ( 2*na ) * Factorial ( 2*nb ) * Factorial ( 2*nc ) * Factorial ( 2*nd ) )  ;
/*    c   = ev*ff*exp(na*aea + nb*aeb + nc*aec + nd*aed + 0.5e+00*(aea + aeb + aec + aed) + a2*(n + 2) - ae*n) ; */
    c   = ff*exp(na*aea + nb*aeb + nc*aec + nd*aed + 0.5e+00*(aea + aeb + aec + aed) + a2*(n + 2) - ae*n) ;
    s0  = 1.0e+00 / e ;
    s1  = 0.0e+00 ;
    s2  = 0.0e+00 ;
    m   = ncd - k ;
    for ( i = 0 ; i < m ; i++ )
    {
        s0 = s0 * e / ecd ;
        s1 = s1 + s0 * ( Binomial ( ncd-k-1, i ) - Binomial ( ncd+k, i ) ) / Binomial ( n-1, i ) ;
    }
    m2 = ncd + k ;
    for ( i = m ; i <= m2 ; i++ )
    {
        s0 = s0 * e / ecd ;
        s2 = s2 + s0 * Binomial ( m2, i ) / Binomial ( n-1, i ) ;
    }
    s3 = exp ( ae * n - acd * ( m2 + 1 ) - aab * ( nab - k ) ) / Binomial ( n-1, m2 ) ;
    return  c * ( s1 - s2 + s3 ) ;
}

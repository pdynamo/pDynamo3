/*==================================================================================================================================
! . This module handles the DFT grid.
!
! . This module uses a modified Mura-Knowles method for radial integration and Lebedev grids for angular integration.
!
! . All units are atomic.
!===================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Boolean.h"
# include "DFTGrid.h"
# include "Lebedev.h"
# include "Memory.h"
# include "RealUtilities.h"
# include "Units.h"

/* . Nothing complicated done here until decide what to do about the grid. */

# define _PRINTWARNINGS False

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The elemental bragg radii. */
/* . Units are Angstroms and must be converted to Bohrs before use. */
# define NELEMENTS 120
const Real BraggRadii[NELEMENTS] = { 0.75, 0.35, 0.35, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50 ,
                                     0.50, 1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.00, 2.20 ,
                                     1.80, 1.60, 1.40, 1.35, 1.40, 1.40, 1.40, 1.35, 1.35, 1.35 ,
                                     1.35, 1.30, 1.25, 1.15, 1.15, 1.15, 1.15, 2.35, 2.00, 1.80 ,
                                     1.55, 1.45, 1.45, 1.35, 1.30, 1.35, 1.40, 1.60, 1.55, 1.55 ,
                                     1.45, 1.45, 1.40, 1.40, 1.40, 2.60, 2.15, 1.95, 1.85, 1.85 ,
                                     1.85, 1.85, 1.85, 1.85, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75 ,
                                     1.75, 1.75, 1.55, 1.45, 1.35, 1.35, 1.30, 1.35, 1.35, 1.35 ,
                                     1.50, 1.90, 1.80, 1.60, 1.90, 1.90, 1.90, 2.60, 2.15, 1.95 ,
                                     1.80, 1.80, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75, 1.75 ,
                                     1.75, 1.75, 1.75, 1.75, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55 ,
                                     1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55, 1.55 } ;

/* . Tolerances for basis function calculation. */
const Real BFTolerances[NDFTGRID_ACCURACY] = { 1.0e-8, 1.0e-9, 1.0e-10, 1.0e-12, 1.0e-15 } ;

/* . Tolerances for density calculation. */
const Real RhoTolerances[NDFTGRID_ACCURACY] = { 1.0e-13, 1.0e-14, 1.0e-15, 1.0e-17, 1.0e-20 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void               DFTGrid_Atom_Parameters      ( const DFTGridAccuracy accuracy       ,
                                                         const Integer         ni             ,
                                                               Integer        *nR             ,
                                                               Integer        *lValue         ,
                                                               Real           *maximumRadius  ) ;
static Real               DFTGrid_Bragg_Radius         ( const Integer         atomicNumber   ) ;
static void               DFTGrid_Radial_Points        ( const Integer         nR             ,
                                                         const Real            range          ,
                                                               Real           *r              ,
                                                               Real           *w              ) ;

static DFTGridPointBlock *DFTGridPointBlock_Allocate   ( const Integer         numberOfPoints ,
                                                         const Integer         atom           ,
                                                               Coordinates3  **rG             ,
                                                               RealArray1D   **wG             ,
                                                               Status         *status         ) ;
static void               DFTGridPointBlock_Deallocate (       void           *vSelf          ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTGrid *DFTGrid_Allocate ( const DFTGridAccuracy  accuracy ,
                                  Status          *status   )
{
    DFTGrid *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( DFTGrid ) ;
        if ( self != NULL )
        {
            self->accuracy        = accuracy ;
            self->blockSize       = 128 ;
            self->numberOfPoints  =   0 ;
            self->numberOfRecords =  -1 ;
            self->bfTolerance     = BFTolerances [accuracy] ;
            self->rhoTolerance    = RhoTolerances[accuracy] ;
            self->points          = List_Allocate ( ) ;
            self->records         = NULL ;
            self->weights         = NULL ;
            if ( self->points != NULL ) self->points->Element_Deallocate = DFTGridPointBlock_Deallocate ;
            else DFTGrid_Deallocate ( &self, status ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Construct a grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MINIMUM_LVALUE       9
# define RADIAL_CUTOFF_FACTOR 0.2e+00
# define WEIGHT_TOLERANCE     1.0e-30

DFTGrid *DFTGrid_Construct ( const DFTGridAccuracy  accuracy       ,
                             const IntegerArray1D  *atomicNumbers  ,
                             const Coordinates3    *qcCoordinates3 ,
                                   Status          *status         )
{
    DFTGrid *self = NULL ;
    if ( ( atomicNumbers                   != NULL ) &&
         ( qcCoordinates3                  != NULL ) &&
         ( View1D_Extent ( atomicNumbers ) >  0    ) &&
         Status_IsOK ( status ) )
    {
        auto Integer  nAtoms = View1D_Extent ( atomicNumbers ) ;
        auto Real    *work1, *work2 ;
        auto Status  *localStatus = NULL ;
        Status_Set ( localStatus, Status_OK ) ;
        self  = DFTGrid_Allocate ( accuracy , localStatus ) ;
        work1 = Real_Allocate    ( nAtoms   , localStatus ) ;
        work2 = Real_Allocate    ( nAtoms   , localStatus ) ;
        if ( Status_IsOK ( localStatus ) )
        {
            auto Integer iqm ;
            /* . Temporarily fill work1 with the Bragg radii for weights allocation. */
            for ( iqm = 0 ; iqm < nAtoms ; iqm++ ) work1[iqm] = DFTGrid_Bragg_Radius ( Array1D_Item ( atomicNumbers, iqm ) ) ;
            self->weights = DFTGridWeights_Allocate ( qcCoordinates3, work1, localStatus ) ;
            if ( Status_IsOK ( localStatus ) )
            {
                auto Integer            ia, ipt, ir, lMax = -1, lOld, lStart, lStop, lVal, n, nAngMax, nAngMin, nApts = 0, nLocal, nR = 0 ;
                auto Real               maximumRadius = 0.0e+00, range, rCutoff, w, wfac, xqm, yqm, zqm ;
                auto Real               pg[3], *rr, *wa, *wr, *xa, *ya, *za ;
                auto Coordinates3      *rG, *rLocal, rView ;
                auto DFTGridPointBlock *block ;
                auto RealArray1D       *wG, *wLocal, wView ;
                /* . Loop over the atoms. */
                for ( iqm = 0 ; iqm < nAtoms ; iqm++ )
                {
                    /* . Get the data for the atom. */
                    Coordinates3_GetRow ( qcCoordinates3, iqm, xqm, yqm, zqm ) ;
                    DFTGrid_Atom_Parameters ( self->accuracy, Array1D_Item ( atomicNumbers, iqm ), &nR, &lMax, &maximumRadius ) ;
                    nAngMax = LebedevLaikov_Number_Of_Points ( lMax           ) ;
                    nAngMin = LebedevLaikov_Number_Of_Points ( MINIMUM_LVALUE ) ;
                    range   = DFTGrid_Bragg_Radius ( Array1D_Item ( atomicNumbers, iqm ) ) ;
                    /* . Allocate space. */
                    rr = Real_Allocate ( nR     , localStatus ) ;
                    wr = Real_Allocate ( nR     , localStatus ) ;
                    wa = Real_Allocate ( nAngMax, localStatus ) ;
                    xa = Real_Allocate ( nAngMax, localStatus ) ;
                    ya = Real_Allocate ( nAngMax, localStatus ) ;
                    za = Real_Allocate ( nAngMax, localStatus ) ;
                    rG = Coordinates3_Allocate          ( nAngMax * nR, localStatus ) ;
                    wG = RealArray1D_AllocateWithExtent ( nAngMax * nR, localStatus ) ;
                    if ( Status_IsOK ( localStatus ) )
                    {
                        /* . Get the radial grid points. */
                        rCutoff = RADIAL_CUTOFF_FACTOR * range ;
                        DFTGrid_Radial_Points ( nR, range, rr, wr ) ;
                        /* . Loop  over the radial grid points. */
                        for ( ir = 0, ipt = 0, lOld = -1 ; ir < nR ; ir++ )
                        {
                            /* . Check for the maximum value of r. */
                            if ( rr[ir] > maximumRadius ) break ;
                            /* . Get the angular points. */
                            if ( rr[ir] > rCutoff )
                            {
                                lVal = lMax ;
                            }
                            else
                            {
                                n = ( Integer ) ceil ( ( ( Real ) nAngMax ) * rr[ir] / rCutoff ) ;
                                if ( n < nAngMin ) lVal = MINIMUM_LVALUE ;
                                else               lVal = LebedevLaikov_Angular_Momentum_Value ( n ) ;
                            }
                            /* . Get the angular points. */
                            if ( lVal != lOld )
                            {
                                nApts = LebedevLaikov_Number_Of_Points ( lVal ) ;
                                nApts = LebedevLaikov_Points ( nApts, xa, ya, za, wa ) ;
                                lOld  = lVal ;
                            }
                            /* . Construct the integration points. */
                            {
                                auto Real sum = 0.0e+00, xdev ;
                                wfac = 4.0e+00 * M_PI * wr[ir] ;
                                for ( ia = 0 ; ia < nApts ; ia++ )
                                {
                                    /* . Get the point. */
                                    pg[0] = xa[ia] * rr[ir] + xqm ;
                                    pg[1] = ya[ia] * rr[ir] + yqm ;
                                    pg[2] = za[ia] * rr[ir] + zqm ;
                                    w = wfac * wa[ia] * DFTGridWeights_Weight ( self->weights, iqm, pg, work1, work2 ) ;
                                    if ( fabs ( w ) > WEIGHT_TOLERANCE )
                                    {
                                        Coordinates3_SetRow ( rG, ipt, pg[0], pg[1], pg[2] ) ;
                                        Array1D_Item ( wG, ipt ) = w ;
                                        ipt++ ;
                                    }
                                    xdev = fabs ( 1.0e+00 - xa[ia]*xa[ia] - ya[ia]*ya[ia] - za[ia]*za[ia] ) ;
# ifdef _PRINTWARNINGS
                                    if ( xdev > 1.0e-8 ) printf ( "Node Inaccuracy = %5d %5d %5d %25.15f\n", ir, nApts, lVal, xdev ) ;
# endif
                                    sum += wa[ia] ;
                                }
                                xdev = fabs ( sum - 1.0e+00 ) ;
# ifdef _PRINTWARNINGS
                                if ( xdev > 1.0e-9 ) printf ( "Weight Inaccuracy = %5d %5d %5d %25.15f\n", ir, nApts, lVal, xdev ) ;
# endif
                            }
                        }
                        /* . Save the grid points in blocks of the appropriate size. */
                        lStart = 0 ;
                        do
                        {
                            lStop                 = Minimum ( lStart + self->blockSize, ipt ) ;
                            nLocal                = lStop - lStart ;
                            self->numberOfPoints += nLocal ;
                            rLocal                = Coordinates3_Allocate          ( nLocal, localStatus ) ;
                            wLocal                = RealArray1D_AllocateWithExtent ( nLocal, localStatus ) ;
                            if ( Status_IsOK ( localStatus ) )
                            {
                                Coordinates3_View2D ( rG, lStart, 0, nLocal, 3, 1, 1, False, &rView, NULL ) ;
                                RealArray1D_View    ( wG, lStart,    nLocal,    1,    False, &wView, NULL ) ;
                                Coordinates3_CopyTo ( &rView, rLocal, NULL ) ;
                                RealArray1D_CopyTo  ( &wView, wLocal, NULL ) ;
                                block = DFTGridPointBlock_Allocate ( nLocal, iqm, &rLocal, &wLocal, localStatus ) ;
                                List_Element_Append ( self->points, ( void * ) block ) ;
                                lStart += self->blockSize ;
                            }
                            else
                            {
                                Coordinates3_Deallocate ( &rLocal ) ;
                                RealArray1D_Deallocate  ( &wLocal ) ;
                                goto EndOfLoop ;
                            }
                        }
                        while ( lStart < ipt ) ;
                    }
                    /* . Deallocation. */
                EndOfLoop:
                    Coordinates3_Deallocate ( &rG ) ;
                    RealArray1D_Deallocate  ( &wG ) ;
                    Real_Deallocate ( &rr ) ;
                    Real_Deallocate ( &wr ) ;
                    Real_Deallocate ( &wa ) ;
                    Real_Deallocate ( &xa ) ;
                    Real_Deallocate ( &ya ) ;
                    Real_Deallocate ( &za ) ;
                }
            }
        }
        /* . Deallocation. */
        Real_Deallocate ( &work1 ) ;
        Real_Deallocate ( &work2 ) ;
        if ( ! Status_IsOK ( localStatus ) )
        {
            DFTGrid_Deallocate ( &self, NULL ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

# undef MINIMUM_LVALUE
# undef RADIAL_CUTOFF_FACTOR
# undef WEIGHT_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTGrid_Deallocate ( DFTGrid **self, Status *status )
{
    if ( (*self) != NULL )
    {
        DFTGrid_DeallocateFunctionData (  (*self), status    ) ;
        DFTGridWeights_Deallocate      ( &((*self)->weights) ) ;
        List_Deallocate                ( &((*self)->points)  ) ;
        Memory_Deallocate              (  (*self)->records   ) ;
        Memory_Deallocate              (  (*self)            ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation of function data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTGrid_DeallocateFunctionData ( DFTGrid *self, Status *status )
{
    if ( ( self != NULL ) && DFTGrid_HasFunctionData ( self, status ) && Status_IsOK ( status ) )
    {
        auto Integer r ;
        for ( r = 0 ; r < self->numberOfRecords ; r++ ) GridFunctionDataBlock_Deallocate ( &(self->records[r]->functionData) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Estimate the number of points in the grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MINIMUM_LVALUE       9
# define RADIAL_CUTOFF_FACTOR 0.2e+00
Integer DFTGrid_EstimatedPoints ( const DFTGridAccuracy  accuracy      ,
                                  const IntegerArray1D  *atomicNumbers ,
                                        Status          *status        )
{
    Integer p = 0 ;
    if ( ( atomicNumbers                   != NULL ) &&
         ( View1D_Extent ( atomicNumbers ) >  0    ) &&
         Status_IsOK ( status )                    )
    {
        auto Integer iqm, ir, lMax = -1, lOld, lVal, n, nAngMax, nAngMin, nApts = 0, nR = 0 ;
        auto Real    maximumRadius = 0.0e+00, range, rCutoff, *rr, *wr ;
        /* . Loop over the atoms. */
        for ( iqm = 0 ; iqm < View1D_Extent ( atomicNumbers ) ; iqm++ )
        {
            DFTGrid_Atom_Parameters ( accuracy, Array1D_Item ( atomicNumbers, iqm ), &nR, &lMax, &maximumRadius ) ;
            nAngMax = LebedevLaikov_Number_Of_Points ( lMax           ) ;
            nAngMin = LebedevLaikov_Number_Of_Points ( MINIMUM_LVALUE ) ;
            range   = DFTGrid_Bragg_Radius ( Array1D_Item ( atomicNumbers, iqm ) ) ;
            rr      = Real_Allocate ( nR, status ) ;
            wr      = Real_Allocate ( nR, status ) ;
            if Status_IsOK ( status )
            {
                /* . Get the radial grid points. */
                rCutoff = RADIAL_CUTOFF_FACTOR * range ;
                DFTGrid_Radial_Points ( nR, range, rr, wr ) ;
                /* . Loop  over the radial grid points. */
                for ( ir = 0, lOld = -1 ; ir < nR ; ir++ )
                {
                    /* . Check for the maximum value of r. */
                    if ( rr[ir] > maximumRadius ) break ;
                    /* . Get the angular points. */
                    if ( rr[ir] > rCutoff )
                    {
                        lVal = lMax ;
                    }
                    else
                    {
                        n = ( Integer ) ceil ( ( ( Real ) nAngMax ) * rr[ir] / rCutoff ) ;
                        if ( n < nAngMin ) lVal = MINIMUM_LVALUE ;
                        else               lVal = LebedevLaikov_Angular_Momentum_Value ( n ) ;
                    }
                    /* . Get the angular points. */
                    if ( lVal != lOld )
                    {
                        nApts = LebedevLaikov_Number_Of_Points ( lVal ) ;
                        lOld  = lVal ;
                    }
                    p += nApts ;
                }
            }
            Real_Deallocate ( &rr ) ;
            Real_Deallocate ( &wr ) ;
        }
    }
    return p ;
}
# undef MINIMUM_LVALUE
# undef RADIAL_CUTOFF_FACTOR

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size of the block storage in bytes.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real DFTGrid_FunctionByteSize ( DFTGrid *self, Status *status )
{
    Real size = 0.0e+00 ;
    if ( ( self != NULL ) && ( self->numberOfPoints > 0 ) && DFTGrid_HasFunctionData ( self, status ) && Status_IsOK ( status ) )
    {
        auto Integer r ;
        for ( r = 0 ; r < self->numberOfRecords ; r++ ) size += GridFunctionDataBlock_ByteSize ( self->records[r]->functionData ) ;
    }
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Has the grid function data?
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean DFTGrid_HasFunctionData ( DFTGrid *self, Status *status )
{
    Boolean hasData = False ;
    if ( ( self != NULL ) && ( self->numberOfPoints > 0 ) && Status_IsOK ( status ) )
    {
        DFTGrid_MakeRecords ( self, status ) ;
        if ( ( self->records != NULL ) && ( self->numberOfRecords > 0 ) ) hasData = ( self->records[0]->functionData != NULL ) ;
    }
    return hasData ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTGridPointBlock *DFTGrid_Iterate ( DFTGrid *self )
{
   if ( self == NULL ) return NULL ;
   else                return ( DFTGridPointBlock * ) List_Iterate ( self->points ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the records representation of the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTGrid_MakeRecords ( DFTGrid *self, Status *status )
{
    if ( ( self != NULL ) && ( self->records == NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = DFTGrid_NumberOfRecords ( self ) ;
        self->records = Memory_AllocateArrayOfReferences ( n, DFTGridPointBlock ) ;
        if ( self->records == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        else
        {
            auto DFTGridPointBlock *record ;
            n = 0 ;
            List_Iterate_Initialize ( self->points ) ;
            while ( ( record = DFTGrid_Iterate ( self ) ) != NULL ) { self->records[n] = record ; n += 1 ; }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the number of stored function values.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DFTGrid_NumberOfFunctionValues ( DFTGrid *self, Status *status )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( self->numberOfPoints > 0 ) && Status_IsOK ( status ) )
    {
        auto Integer r ;
        DFTGrid_MakeRecords ( self, status ) ;
        if ( Status_IsOK ( status ) )
        {
            for ( r = 0 ; r < self->numberOfRecords ; r++ )
            {
                if ( self->records[r]->functionData != NULL ) n += ( self->records[r]->functionData->numberOfFunctions * self->records[r]->numberOfPoints ) ;
            }
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of points in the grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DFTGrid_NumberOfPoints ( DFTGrid *self )
{
    Integer n = 0 ;
    if ( self != NULL ) n = self->numberOfPoints ;
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of records in the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DFTGrid_NumberOfRecords ( DFTGrid *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        if ( self->numberOfRecords >= 0 ) n = self->numberOfRecords ;
        else
        {
            auto DFTGridPointBlock *record ;
            List_Iterate_Initialize ( self->points ) ;
            while ( ( record = DFTGrid_Iterate ( self ) ) != NULL ) n += 1 ;
            self->numberOfRecords = n ;
        }
    }
    return n ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the grid parameters for an element.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTGrid_Atom_Parameters ( const DFTGridAccuracy accuracy, const Integer ni, Integer *nR, Integer *lValue, Real *maximumRadius )
{
    switch ( accuracy )
    {
        case DFTGridAccuracy_VeryLow:
                 if ( ni < 11 ) (*nR) = 21 ;
            else if ( ni < 19 ) (*nR) = 42 ;
            else if ( ni < 37 ) (*nR) = 75 ;
            else                (*nR) = 84 ;
            (*maximumRadius) = 20.0e+00 ;
            (*lValue)        = 23 ;
            break ;
        case DFTGridAccuracy_Low:
                 if ( ni < 11 ) (*nR) =  35 ;
            else if ( ni < 19 ) (*nR) =  70 ;
            else if ( ni < 37 ) (*nR) =  95 ;
            else                (*nR) = 104 ;
            (*maximumRadius) = 25.0e+00 ;
            (*lValue)        = 35 ;
            break ;
        case DFTGridAccuracy_Medium:
                 if ( ni < 11 ) (*nR) =  49 ;
            else if ( ni < 19 ) (*nR) =  88 ;
            else if ( ni < 37 ) (*nR) = 112 ;
            else                (*nR) = 123 ;
            (*maximumRadius) = 30.0e+00 ;
            if ( ni < 11 ) (*lValue) = 35 ;
            else           (*lValue) = 41 ;
            break ;
        case DFTGridAccuracy_High:
                 if ( ni < 11 ) (*nR) =  70 ;
            else if ( ni < 19 ) (*nR) = 123 ;
            else if ( ni < 37 ) (*nR) = 130 ;
            else                (*nR) = 155 ;
            (*maximumRadius) = 35.0e+00 ;
                 if ( ni < 11 ) (*lValue) = 41 ;
            else if ( ni < 19 ) (*lValue) = 47 ;
            else if ( ni < 89 ) (*lValue) = 53 ;
            else                (*lValue) = 59 ;
            break ;
        case DFTGridAccuracy_VeryHigh:
                 if ( ni < 11 ) (*nR) = 100 ;
            else if ( ni < 19 ) (*nR) = 125 ;
            else if ( ni < 37 ) (*nR) = 160 ;
            else                (*nR) = 205 ;
            (*maximumRadius) = 35.0e+00 ;
            (*lValue)        = 65 ;
            break ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the bragg radius for an atom.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real DFTGrid_Bragg_Radius ( const Integer atomicNumber )
{
    Integer ni ;
    ni = atomicNumber ;
    if      ( ni <  0         ) ni = 0 ;
    else if ( ni >= NELEMENTS ) ni = NELEMENTS - 1 ;
    return Units_Length_Angstroms_To_Bohrs * BraggRadii[ni] ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the radial integration points.
! . Follows the Mura and Knowles scheme (JCP 104, 9848, 1996).
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef DFTGRID_MURAKNOWLES
# define EXPONENT 3.0e+00
# define SCALING  3.3e+00
static void DFTGrid_Radial_Points ( const Integer nR, const Real range, Real *r, Real *w )
{
    Integer i ;
    Real    alpha, fmn, qi, ri, wi ;
    alpha = SCALING * range ;
    fmn   = EXPONENT / ( 1.0e+00 + ( Real ) nR ) ;
    for ( i = 0 ; i < nR ; i++ )
    {
       qi = ( Real ) i / ( ( Real ) nR + 1.0e+00 ) ;
       ri = - alpha * log ( 1.0e+00 - pow ( qi, EXPONENT ) ) ;
       wi = fmn * alpha * ( ri * ri ) / ( 1.0e+00 - pow ( qi, EXPONENT ) ) * pow ( qi, EXPONENT - 1.0e+00 ) ;
       r[i] = ri ;
       w[i] = wi ;
    }
}
# undef EXPONENT
# undef SCALING
# else
static void DFTGrid_Radial_Points ( const Integer nR, const Real range, Real *r, Real *w )
{
    Integer i ;
    Real    dfac, ri, wi, xnode ;
    dfac = ( 1.0e+00 + ( Real ) nR ) ;
    for ( i = 0 ; i < nR ; i++ )
    {
       xnode = ( Real ) ( i + 1 ) / dfac ;
       ri    = range * pow ( ( xnode / ( 1.0e+00 - xnode ) ), 2.0e+00 ) ;
       wi    = 2.0e+00 * pow ( range, 3.0e+00 ) * pow ( xnode, 5.0e+00 ) / pow ( ( 1.0e+00 - xnode ), 7.0e+00 ) ;
       r[i] = ri ;
       w[i] = wi / dfac ;
    }
}
# endif

/*==================================================================================================================================
! . Private grid point block procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
static DFTGridPointBlock *DFTGridPointBlock_Allocate ( const Integer        numberOfPoints ,
                                                       const Integer        atom           ,
                                                             Coordinates3 **rG             ,
                                                             RealArray1D  **wG             ,
                                                             Status        *status         )
{
   DFTGridPointBlock *self = NULL ;
   if ( Status_IsOK ( status ) && ( numberOfPoints > 0 ) )
   {
       self = Memory_AllocateType ( DFTGridPointBlock ) ;
       if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
       else
       {
           self->atom           = atom  ;
           self->coordinates3   = (*rG) ;
           self->functionData   = NULL  ;
           self->numberOfPoints = numberOfPoints ;
           self->weights        = (*wG) ;
           /* . Take ownership of the arrays. */
           (*rG) = NULL ;
           (*wG) = NULL ;
       }
   }
   return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocate a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DFTGridPointBlock_Deallocate ( void *vSelf )
{
   DFTGridPointBlock *self ;
   self = ( DFTGridPointBlock * ) vSelf ;
   Coordinates3_Deallocate          ( &(self->coordinates3) ) ;
   GridFunctionDataBlock_Deallocate ( &(self->functionData) ) ;
   RealArray1D_Deallocate           ( &(self->weights     ) ) ;
   Memory_Deallocate ( self ) ;
}

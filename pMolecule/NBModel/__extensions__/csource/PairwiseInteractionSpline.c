/*==================================================================================================================================
! . Pairwise interactions of general form and limited range.
!=================================================================================================================================*/

# include "Boolean.h"
# include "Memory.h"
# include "MinimumImageUtilities.h"
# include "NumericalMacros.h"
# include "PairwiseInteractionSpline.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DefaultCutOff 12.0e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionSpline_Initialize ( PairwiseInteractionSpline *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionSpline *PairwiseInteractionSpline_Allocate ( Status *status )
{
    PairwiseInteractionSpline *self = NULL ;
    self = Memory_AllocateType ( PairwiseInteractionSpline ) ;
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    else PairwiseInteractionSpline_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Assign splines.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_AssignElectrostaticSpline ( PairwiseInteractionSpline *self, CubicSpline *spline, Real cutOff )
{
    if ( self != NULL )
    {
        CubicSpline_Deallocate ( &(self->electrostaticSpline) ) ;
        self->electrostaticSpline = spline ;
        self->cutOff              = Maximum ( cutOff, self->cutOff ) ;
        self->cutOff2             = self->cutOff * self->cutOff      ;
    }
}
void PairwiseInteractionSpline_AssignLennardJonesASpline ( PairwiseInteractionSpline *self, CubicSpline *spline, Real cutOff )
{
    if ( self != NULL )
    {
        CubicSpline_Deallocate ( &(self->lennardJonesASpline) ) ;
        self->lennardJonesASpline = spline ;
        self->cutOff              = Maximum ( cutOff, self->cutOff ) ;
        self->cutOff2             = self->cutOff * self->cutOff      ;
    }
}
void PairwiseInteractionSpline_AssignLennardJonesBSpline ( PairwiseInteractionSpline *self, CubicSpline *spline, Real cutOff )
{
    if ( self != NULL )
    {
        CubicSpline_Deallocate ( &(self->lennardJonesBSpline) ) ;
        self->lennardJonesBSpline = spline ;
        self->cutOff              = Maximum ( cutOff, self->cutOff ) ;
        self->cutOff2             = self->cutOff * self->cutOff      ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionSpline *PairwiseInteractionSpline_Clone ( PairwiseInteractionSpline *self, Status *status )
{
    PairwiseInteractionSpline *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = PairwiseInteractionSpline_Allocate ( status ) ;
        if ( clone != NULL )
        {
            clone->cutOff              = self->cutOff  ;
            clone->cutOff2             = self->cutOff2 ;
            clone->electrostaticSpline = CubicSpline_Clone ( self->electrostaticSpline, status ) ;
            clone->lennardJonesASpline = CubicSpline_Clone ( self->lennardJonesASpline, status ) ;
            clone->lennardJonesBSpline = CubicSpline_Clone ( self->lennardJonesBSpline, status ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_Deallocate ( PairwiseInteractionSpline **self )
{
    if ( (*self) != NULL )
    {
        CubicSpline_Deallocate ( &((*self)->electrostaticSpline) ) ;
        CubicSpline_Deallocate ( &((*self)->lennardJonesASpline) ) ;
        CubicSpline_Deallocate ( &((*self)->lennardJonesBSpline) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deassign splines.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_DeassignSplines ( PairwiseInteractionSpline *self )
{
    if ( self != NULL )
    {
        self->electrostaticSpline = NULL ;
        self->lennardJonesASpline = NULL ;
        self->lennardJonesBSpline = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionSpline_Initialize ( PairwiseInteractionSpline *self )
{
    if ( self != NULL )
    {
        self->cutOff              = _DefaultCutOff ;
        self->cutOff2             = self->cutOff * self->cutOff ;
        self->electrostaticSpline = NULL ;
        self->lennardJonesASpline = NULL ;
        self->lennardJonesBSpline = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_Interactions ( const PairwiseInteractionSpline *self ,
                                              const RealArray1D *r                  ,
                                                    RealArray1D *electrostatic      ,
                                                    RealArray1D *lennardJonesA      ,
                                                    RealArray1D *lennardJonesB      )
{
    Boolean doElectrostatic = ( electrostatic != NULL ) ,
            doLennardJonesA = ( lennardJonesA != NULL ) ,
            doLennardJonesB = ( lennardJonesB != NULL ) ;
    if ( ( self != NULL ) && ( r != NULL ) && ( doElectrostatic|| doLennardJonesA || doLennardJonesB ) )
    {
        auto Integer  i ;
        auto Integer  l, u ;
        auto Real     d, fA, fB, fQ, g, s, t, x, x2 ;
        RealArray1D_Set ( electrostatic, 0.0e+00 ) ;
        RealArray1D_Set ( lennardJonesA, 0.0e+00 ) ;
        RealArray1D_Set ( lennardJonesB, 0.0e+00 ) ;
        for ( i = 0 ; i < View1D_Extent ( r ) ; i++ )
        {
            x  = Array1D_Item ( r, i ) ;
            x2 = x * x ;
            if ( x2 > self->cutOff2 ) continue ;
            CubicSpline_EvaluateLUDST   ( self->electrostaticSpline, x2, &l, &u, &d, &s, &t        ) ;
            CubicSpline_FastEvaluateFGN ( self->electrostaticSpline,  0,  l,  u,  d,  s,  t, fQ, g ) ;
            CubicSpline_EvaluateLUDST   ( self->lennardJonesASpline, x2, &l, &u, &d, &s, &t        ) ;
            CubicSpline_FastEvaluateFGN ( self->lennardJonesASpline,  0,  l,  u,  d,  s,  t, fA, g ) ;
            CubicSpline_EvaluateLUDST   ( self->lennardJonesBSpline, x2, &l, &u, &d, &s, &t        ) ;
            CubicSpline_FastEvaluateFGN ( self->lennardJonesBSpline,  0,  l,  u,  d,  s,  t, fB, g ) ;
            if ( doElectrostatic ) Array1D_Item ( electrostatic, i ) = fQ ;
            if ( doLennardJonesA ) Array1D_Item ( lennardJonesA, i ) = fA ;
            if ( doLennardJonesB ) Array1D_Item ( lennardJonesB, i ) = fB ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM energy and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_MMMMEnergy ( const PairwiseInteractionSpline *self               ,
                                            const RealArray1D               *chargesI           ,
                                            const RealArray1D               *chargesJ           ,
                                                  IntegerArray1D            *ljTypesI           ,
                                                  IntegerArray1D            *ljTypesJ           ,
                                            const LJParameterContainer      *ljParameters       ,
                                            const Real                       electrostaticScale ,
                                            const Real                       lennardJonesScale  ,
                                            const Coordinates3              *coordinates3I      ,
                                            const Coordinates3              *coordinates3J      ,
                                                  PairList                  *pairList           ,
                                                  Real                      *eElectrostatic     ,
                                                  Real                      *eLennardJones      ,
                                                  Coordinates3              *gradients3I        ,
                                                  Coordinates3              *gradients3J        ,
                                                  Status                    *status             )
{
    if ( eElectrostatic != NULL ) (*eElectrostatic ) = 0.0e+00 ;
    if ( eLennardJones  != NULL ) (*eLennardJones ) = 0.0e+00 ;
    if ( ( self          != NULL ) &&
         ( coordinates3I != NULL ) &&
         ( coordinates3J != NULL ) &&
         ( pairList      != NULL ) &&
           Status_IsOK ( status ) )
    {
        auto Boolean doElectrostatic, doGradients, doLennardJones ;
        doElectrostatic = ( chargesI                  != NULL    ) &&
                          ( chargesJ                  != NULL    ) &&
                          ( eElectrostatic            != NULL    ) &&
                          ( electrostaticScale        != 0.0e+00 ) &&
                          ( self->electrostaticSpline != NULL    ) ;
        doGradients     = ( gradients3I               != NULL    ) &&
                          ( gradients3J               != NULL    ) ;
        doLennardJones  = ( eLennardJones             != NULL    ) &&
                          ( ljTypesI                  != NULL    ) &&
                          ( ljTypesJ                  != NULL    ) &&
                          ( ljParameters              != NULL    ) &&
                          ( lennardJonesScale         != 0.0e+00 ) &&
                          ( self->lennardJonesASpline != NULL    ) &&
                          ( self->lennardJonesBSpline != NULL    ) ;
        if ( doElectrostatic|| doLennardJones )
        {
            auto Integer           i, j, n, numberOfLJTypes = 0, tI = 0, tIJ ;
            auto CubicSpline      *referenceSpline ;
            auto Integer           l, u  ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              aIJ, bIJ, cutOff2 = self->cutOff2, d, eScale, f, g, gL, qI = 0.0e+00, qIJ,
                                   r2, s, t, xI, xIJ, xJ, yI, yIJ, yJ, zI, zIJ, zJ ;
            auto Real              eLJ = 0.0e+00, eQQ = 0.0e+00 ;
            /* . Initialization. */
            eScale = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            if ( doElectrostatic ) referenceSpline = self->electrostaticSpline ;
            else                   referenceSpline = self->lennardJonesASpline ;
            if ( doLennardJones  ) numberOfLJTypes = ljParameters->ntypes ;
            /* . Loop over record. */
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                /* . First atom. */
                i  = record->index ;
	        if ( doElectrostatic ) qI = eScale          * Array1D_Item   ( chargesI, i ) ;
                if ( doLennardJones  ) tI = numberOfLJTypes * Array1D_Item ( ljTypesI, i ) ;
	        Coordinates3_GetRow ( coordinates3I, i, xI, yI, zI ) ;
                /* . Second atom. */
       	        for ( n = 0 ; n < record->capacity ; n++ )
	        {
	            j   = record->indices[n] ;
	            Coordinates3_GetRow ( coordinates3J , j, xJ, yJ, zJ ) ;
                    xIJ = xI - xJ ;
                    yIJ = yI - yJ ;
                    zIJ = zI - zJ ;
                    r2  = ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_EvaluateLUDST ( referenceSpline, r2, &l, &u, &d, &s, &t ) ; /* . Assume all abscissae are the same. */
                    g   = 0.0e+00 ;
                    if ( doElectrostatic )
                    {
                        qIJ  = qI * Array1D_Item ( chargesJ, j ) ;
                        CubicSpline_FastEvaluateFGN ( self->electrostaticSpline, 0, l, u, d, s, t, f, gL ) ;
                        eQQ += ( qIJ * f ) ;
                        g   += ( 2.0e+00 * qIJ * gL ) ;
                    }
                    if ( doLennardJones )
                    {
	                tIJ  = ljParameters->tableindex[tI+Array1D_Item ( ljTypesJ, j )] ;
                        aIJ  = ljParameters->tableA[tIJ] * lennardJonesScale ;
                        bIJ  = ljParameters->tableB[tIJ] * lennardJonesScale ;
                        CubicSpline_FastEvaluateFGN ( self->lennardJonesASpline, 0, l, u, d, s, t, f, gL ) ;
                        eLJ += aIJ * f  ;
                        g   += ( 2.0e+00 * aIJ * gL ) ;
                        CubicSpline_FastEvaluateFGN ( self->lennardJonesBSpline, 0, l, u, d, s, t, f, gL ) ;
                        eLJ += bIJ * f  ;
                        g   += ( 2.0e+00 * bIJ * gL ) ;
                    }
                    if ( doGradients )
                    {
                        xIJ *= g ;
                        yIJ *= g ;
                        zIJ *= g ;
                        Coordinates3_IncrementRow ( gradients3I, i, xIJ, yIJ, zIJ ) ;
                        Coordinates3_DecrementRow ( gradients3J, j, xIJ, yIJ, zIJ ) ;
	            }
                }
            }
            if ( doElectrostatic ) (*eElectrostatic) = eQQ ;
            if ( doLennardJones  ) (*eLennardJones ) = eLJ ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image MM/MM energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_MMMMEnergyImage ( const PairwiseInteractionSpline  *self                       ,
                                                 const RealArray1D                *charges                    ,
                                                       IntegerArray1D             *ljTypes                    ,
                                                 const LJParameterContainer       *ljParameters               ,
                                                 const Real                        electrostaticScale         ,
                                                       Coordinates3               *coordinates3               ,
                                                       SymmetryParameters         *symmetryParameters         ,
                                                       ImagePairListContainer     *imagePairLists             ,
                                                       Real                       *eElectrostatic             ,
                                                       Real                       *eLennardJones              ,
                                                       Coordinates3               *gradients3                 ,
                                                       SymmetryParameterGradients *symmetryParameterGradients ,
                                                       Status                     *status                     )
{
    Boolean doElectrostatics, doLennardJones ;
    doElectrostatics = ( charges != NULL ) && ( electrostaticScale != 0.0e+00 ) && ( eElectrostatic != NULL ) ;
    doLennardJones   = ( ljTypes != NULL ) && ( ljParameters       != NULL    ) && ( eLennardJones  != NULL ) ;
    if ( eElectrostatic != NULL ) (*eElectrostatic) = 0.0e+00 ;
    if ( eLennardJones  != NULL ) (*eLennardJones ) = 0.0e+00 ;
    if ( ( self               != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( imagePairLists     != NULL ) &&
         ( doElectrostatics || doLennardJones ) &&
         ( Status_IsOK ( status ) ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator                  ,
                                           imagePairLists             ,
                                           coordinates3               ,
                                           symmetryParameters         ,
                                           gradients3                 ,
                                           symmetryParameterGradients ,
                                           status                     ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Real  eLJ     = 0.0e+00 ,  eLJ0    = 0.0e+00 ,
                       eQQ     = 0.0e+00 ,  eQQ0    = 0.0e+00 ,
                      *eTempLJ = NULL    , *eTempQQ = NULL    ;
            if ( doElectrostatics ) eTempQQ = &eQQ0 ;
            if ( doLennardJones   ) eTempLJ = &eLJ0 ;
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                PairwiseInteractionSpline_MMMMEnergy ( self                   ,
                                                       charges                ,
                                                       charges                ,
                                                       ljTypes                ,
                                                       ljTypes                ,
                                                       ljParameters           ,
                                                       electrostaticScale * iterator.scale ,
                                                       iterator.scale         ,
                                                       coordinates3           ,
                                                       iterator.iCoordinates3 ,
                                                       iterator.pairList      ,
                                                       eTempQQ                ,
                                                       eTempLJ                ,
                                                       gradients3             ,
                                                       iterator.iGradients3   ,
                                                       status                 ) ;
                if ( doElectrostatics ) eQQ += (*eTempQQ) ;
                if ( doLennardJones   ) eLJ += (*eTempLJ) ;
                ImagePairListIterator_Gradients ( &iterator ) ;
            }
            if ( doElectrostatics ) (*eElectrostatic) = eQQ ;
            if ( doLennardJones   ) (*eLennardJones ) = eLJ ;
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM energy and gradients within the minimum image convention.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_MMMMEnergyMI ( const PairwiseInteractionSpline  *self                       ,
                                              const RealArray1D                *chargesI                   ,
                                              const RealArray1D                *chargesJ                   ,
                                                    IntegerArray1D             *ljTypesI                   ,
                                                    IntegerArray1D             *ljTypesJ                   ,
                                              const LJParameterContainer       *ljParameters               ,
                                              const Real                        electrostaticScale         ,
                                              const Real                        lennardJonesScale          ,
                                              const Coordinates3               *coordinates3I              ,
                                              const Coordinates3               *coordinates3J              ,
                                                    SymmetryParameters         *symmetryParameters         ,
                                                    PairList                   *pairList                   ,
                                                    Real                       *eElectrostatic             ,
                                                    Real                       *eLennardJones              ,
                                                    Coordinates3               *gradients3I                ,
                                                    Coordinates3               *gradients3J                ,
                                                    SymmetryParameterGradients *symmetryParameterGradients ,
                                                    Status                     *status                     )
{
    if ( eElectrostatic != NULL ) (*eElectrostatic ) = 0.0e+00 ;
    if ( eLennardJones  != NULL ) (*eLennardJones ) = 0.0e+00 ;
    if ( ( self               != NULL ) &&
         ( coordinates3I      != NULL ) &&
         ( coordinates3J      != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( pairList           != NULL ) &&
           Status_IsOK ( status ) )
    {
        auto Boolean doElectrostatic, doGradients, doLennardJones ;
        doElectrostatic = ( chargesI                   != NULL    ) &&
                          ( chargesJ                   != NULL    ) &&
                          ( eElectrostatic             != NULL    ) &&
                          ( electrostaticScale         != 0.0e+00 ) &&
                          ( self->electrostaticSpline  != NULL    ) ;
        doGradients     = ( gradients3I                != NULL    ) &&
                          ( gradients3J                != NULL    ) &&
                          ( symmetryParameterGradients != NULL    ) ;
        doLennardJones  = ( eLennardJones              != NULL    ) &&
                          ( ljTypesI                   != NULL    ) &&
                          ( ljTypesJ                   != NULL    ) &&
                          ( ljParameters               != NULL    ) &&
                          ( lennardJonesScale          != 0.0e+00 ) &&
                          ( self->lennardJonesASpline  != NULL    ) &&
                          ( self->lennardJonesBSpline  != NULL    ) ;
        if ( doElectrostatic|| doLennardJones )
        {
            auto Integer           i, j, l, n, numberOfLJTypes = 0, tI = 0, tIJ, u ;
            auto CubicSpline      *referenceSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              aIJ, bIJ, cutOff2 = self->cutOff2, eScale, f, g, gL, qI = 0.0e+00, qIJ, r2, s, xIJ, yIJ, zIJ ;
            auto Real              eLJ = 0.0e+00, eQQ = 0.0e+00 ;
            /* . Minimum image setup. */
            MinimumImage_Allocate ( coordinates3I, coordinates3J ) ;
            /* . Initialization. */
            eScale = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            if ( doElectrostatic ) referenceSpline = self->electrostaticSpline ;
            else                   referenceSpline = self->lennardJonesASpline ;
            if ( doLennardJones  ) numberOfLJTypes = ljParameters->ntypes ;
            /* . Loop over record. */
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                /* . First atom. */
                i  = record->index ;
	        if ( doElectrostatic ) qI = eScale          * Array1D_Item ( chargesI, i ) ;
                if ( doLennardJones  ) tI = numberOfLJTypes * Array1D_Item ( ljTypesI, i ) ;
                /* . Displacements. */
                MinimumImage_Displacements ( i, j ) ;
                /* . Second atom. */
       	        for ( n = 0 ; n < record->capacity ; n++ )
	        {
	            j = record->indices[n] ;
	            Coordinates3_GetRow ( rDisplacements, n, xIJ, yIJ, zIJ ) ;
                    r2 = ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_EvaluateLUDST ( referenceSpline, r2, &l, &u, &d, &s, &t ) ; /* . Assume all abscissae are the same. */
                    g   = 0.0e+00 ;
                    if ( doElectrostatic )
                    {
                        qIJ  = qI * Array1D_Item ( chargesJ, j ) ;
                        CubicSpline_FastEvaluateFGN ( self->electrostaticSpline, 0, l, u, d, s, t, f, gL ) ;
                        eQQ += ( qIJ * f ) ;
                        g   += ( 2.0e+00 * qIJ * gL ) ;
                    }
                    if ( doLennardJones )
                    {
	                tIJ  = ljParameters->tableindex[tI+Array1D_Item ( ljTypesJ, j )] ;
                        aIJ  = ljParameters->tableA[tIJ] * lennardJonesScale ;
                        bIJ  = ljParameters->tableB[tIJ] * lennardJonesScale ;
                        CubicSpline_FastEvaluateFGN ( self->lennardJonesASpline, 0, l, u, d, s, t, f, gL ) ;
                        eLJ += aIJ * f  ;
                        g   += ( 2.0e+00 * aIJ * gL ) ;
                        CubicSpline_FastEvaluateFGN ( self->lennardJonesBSpline, 0, l, u, d, s, t, f, gL ) ;
                        eLJ += bIJ * f  ;
                        g   += ( 2.0e+00 * bIJ * gL ) ;
                    }
                    if ( doGradients )
                    {
                        xIJ *= g ;
                        yIJ *= g ;
                        zIJ *= g ;
                        Coordinates3_IncrementRow ( gradients3I   , i, xIJ, yIJ, zIJ ) ;
                        Coordinates3_DecrementRow ( gradients3J   , j, xIJ, yIJ, zIJ ) ;
                        Coordinates3_SetRow       ( fDisplacements, n, xIJ, yIJ, zIJ ) ;
	            }
                }
                MinimumImage_Gradients ;
            }
            if ( doElectrostatic ) (*eElectrostatic) = eQQ ;
            if ( doLennardJones  ) (*eLennardJones ) = eLJ ;
            MinimumImage_Deallocate ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCMMGradients ( const PairwiseInteractionSpline *self               ,
                                               const RealArray1D               *chargesQ           ,
                                               const RealArray1D               *chargesM           ,
                                               const Real                       electrostaticScale ,
                                               const Coordinates3              *coordinates3Q      ,
                                               const Coordinates3              *coordinates3M      ,
                                                     PairList                  *pairList           ,
                                               const Coordinates3              *gradients3Q        ,
                                               const Coordinates3              *gradients3M        ,
                                                     Status                    *status             )
{
    if ( ( self               != NULL    ) &&
         ( chargesM           != NULL    ) &&
         ( chargesQ           != NULL    ) &&
         ( coordinates3M      != NULL    ) &&
         ( coordinates3Q      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( gradients3M        != NULL    ) &&
         ( gradients3Q        != NULL    ) &&
         ( pairList           != NULL    ) &&
           Status_IsOK ( status ) )
    {
        if ( self->electrostaticSpline != NULL )
        {
            auto Integer           m, n, q ;
            auto CubicSpline      *spline = self->electrostaticSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              cutOff2 = self->cutOff2, eScale, g, gX, gY, gZ, qQ, qQM,
                                   r2, xM, xQ, xQM, yM, yQ, yQM, zM, zQ, zQM ;
            eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                q  = record->index ;
                qQ = eScale * Array1D_Item ( chargesQ, q ) ;
	        Coordinates3_GetRow ( coordinates3Q, q, xQ, yQ, zQ ) ;
       	        for ( n = 0, gX = gY = gZ = 0.0e+00 ; n < record->capacity ; n++ )
	        {
	            m   = record->indices[n] ;
                    qQM = 2.0e+00 * qQ * Array1D_Item ( chargesM, m ) ; /*. Note the extra factor of 2 here rather than later. */
	            Coordinates3_GetRow ( coordinates3M, m, xM, yM, zM ) ;
                    xQM = xQ - xM ;
                    yQM = yQ - yM ;
                    zQM = zQ - zM ;
                    r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_Evaluate ( spline, 0, r2, NULL, &g, NULL, NULL ) ;
                    g   *= qQM ;
                    xQM *= g   ;
                    yQM *= g   ;
                    zQM *= g   ;
                    Coordinates3_DecrementRow ( gradients3M, m, xQM, yQM, zQM ) ;
                    gX  += xQM ;
                    gY  += yQM ;
                    gZ  += zQM ;
	        }
                Coordinates3_IncrementRow ( gradients3Q, q, gX, gY, gZ ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/MM gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCMMGradientsImage ( const PairwiseInteractionSpline  *self                       ,
                                                    const RealArray1D                *chargesA                   ,
                                                    const RealArray1D                *chargesB                   ,
                                                    const Real                        electrostaticScale         ,
                                                          Coordinates3               *coordinates3A              ,
                                                          Coordinates3               *coordinates3B              ,
                                                          SymmetryParameters         *symmetryParameters         ,
                                                          ImagePairListContainer     *imagePairLists             ,
                                                          Coordinates3               *gradients3A                ,
                                                          Coordinates3               *gradients3B                ,
                                                          SymmetryParameterGradients *symmetryParameterGradients ,
                                                          Status                     *status                     )
{
    if ( ( self                       != NULL ) &&
         ( chargesA                   != NULL ) &&
         ( chargesB                   != NULL ) &&
         ( coordinates3A              != NULL ) &&
         ( coordinates3B              != NULL ) &&
         ( symmetryParameters         != NULL ) &&
         ( imagePairLists             != NULL ) &&
         ( gradients3A                != NULL ) &&
         ( gradients3B                != NULL ) &&
         ( symmetryParameterGradients != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator                  ,
                                           imagePairLists             ,
                                           coordinates3B              ,
                                           symmetryParameters         ,
                                           gradients3B                ,
                                           symmetryParameterGradients ,
                                           status                     ) ;
        if ( Status_IsOK ( status ) )
        {
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                PairwiseInteractionSpline_QCMMGradients ( self                   ,
                                                          chargesA               ,
                                                          chargesB               ,
                                                          electrostaticScale * iterator.scale ,
                                                          coordinates3A          ,
                                                          iterator.iCoordinates3 ,
                                                          iterator.pairList      ,
                                                          gradients3A            ,
                                                          iterator.iGradients3   ,
                                                          status                 ) ;
                ImagePairListIterator_Gradients ( &iterator ) ;
            }
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM gradients within the minimum image convention.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCMMGradientsMI ( const PairwiseInteractionSpline  *self                       ,
                                                 const RealArray1D                *chargesQ                   ,
                                                 const RealArray1D                *chargesM                   ,
                                                 const Real                        electrostaticScale         ,
                                                 const Coordinates3               *coordinates3Q              ,
                                                 const Coordinates3               *coordinates3M              ,
                                                       SymmetryParameters         *symmetryParameters         ,
                                                       PairList                   *pairList                   ,
                                                 const Coordinates3               *gradients3Q                ,
                                                 const Coordinates3               *gradients3M                ,
                                                       SymmetryParameterGradients *symmetryParameterGradients ,
                                                       Status                     *status                     )
{
    if ( ( self                       != NULL    ) &&
         ( chargesM                   != NULL    ) &&
         ( chargesQ                   != NULL    ) &&
         ( coordinates3M              != NULL    ) &&
         ( coordinates3Q              != NULL    ) &&
         ( symmetryParameters         != NULL    ) &&
         ( electrostaticScale         != 0.0e+00 ) &&
         ( gradients3M                != NULL    ) &&
         ( gradients3Q                != NULL    ) &&
         ( symmetryParameterGradients != NULL    ) &&
         ( pairList                   != NULL    ) &&
           Status_IsOK ( status ) )
    {
        if ( self->electrostaticSpline != NULL )
        {
            auto Integer           m, n, q ;
            auto CubicSpline      *spline = self->electrostaticSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              cutOff2 = self->cutOff2, eScale, g, gX, gY, gZ, qQ, qQM, r2, xQM, yQM, zQM ;
            /* . Minimum image setup. */
            MinimumImage_Allocate ( coordinates3Q, coordinates3M ) ;
            /* . Initialization. */
            eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                q  = record->index ;
                qQ = eScale * Array1D_Item ( chargesQ, q ) ;
                /* . Displacements. */
                MinimumImage_Displacements ( q, m ) ;
                /* . MM atom. */
       	        for ( n = 0, gX = gY = gZ = 0.0e+00 ; n < record->capacity ; n++ )
	        {
	            m   = record->indices[n] ;
                    qQM = 2.0e+00 * qQ * Array1D_Item ( chargesM, m ) ; /*. Note the extra factor of 2 here rather than later. */
	            Coordinates3_GetRow ( rDisplacements, m, xQM, yQM, zQM ) ;
                    r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_Evaluate ( spline, 0, r2, NULL, &g, NULL, NULL ) ;
                    g   *= qQM ;
                    xQM *= g   ;
                    yQM *= g   ;
                    zQM *= g   ;
                    Coordinates3_DecrementRow ( gradients3M   , m, xQM, yQM, zQM ) ;
                    Coordinates3_SetRow       ( fDisplacements, m, xQM, yQM, zQM ) ;
                    gX  += xQM ;
                    gY  += yQM ;
                    gZ  += zQM ;
	        }
                Coordinates3_IncrementRow ( gradients3Q, q, gX, gY, gZ ) ;
                MinimumImage_Gradients ;
            }
            MinimumImage_Deallocate ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCMMPotentials ( const PairwiseInteractionSpline *self               ,
                                                const RealArray1D               *chargesM           ,
                                                const Real                       electrostaticScale ,
                                                const Coordinates3              *coordinates3Q      ,
                                                const Coordinates3              *coordinates3M      ,
                                                      PairList                  *pairList           ,
                                                      RealArray1D               *potentials         ,
                                                      Status                    *status             )
{
    if ( ( self               != NULL    ) &&
         ( chargesM           != NULL    ) &&
         ( coordinates3M      != NULL    ) &&
         ( coordinates3Q      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( pairList           != NULL    ) &&
         ( potentials         != NULL    ) &&
           Status_IsOK ( status ) )
    {
        if ( self->electrostaticSpline != NULL )
        {
            auto Integer           m, n, q ;
            auto CubicSpline      *spline = self->electrostaticSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              cutOff2 = self->cutOff2, eScale, f, p, qM, r2, xM, xQ, xQM, yM, yQ, yQM, zM, zQ, zQM ;
            eScale  = electrostaticScale / Units_Length_Angstroms_To_Bohrs ;
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                q  = record->index ;
	        Coordinates3_GetRow ( coordinates3Q, q, xQ, yQ, zQ ) ;
       	        for ( n = 0, p = 0.0e+00 ; n < record->capacity ; n++ )
	        {
	            m   = record->indices[n] ;
                    qM  = Array1D_Item ( chargesM, m ) ;
	            Coordinates3_GetRow ( coordinates3M, m, xM, yM, zM ) ;
                    xQM = xQ - xM ;
                    yQM = yQ - yM ;
                    zQM = zQ - zM ;
                    r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_Evaluate ( spline, 0, r2, &f, NULL, NULL, NULL ) ;
                    p  += qM * f ;
	        }
                Array1D_Item ( potentials, q ) += ( p * eScale ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/MM potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCMMPotentialsImage ( const PairwiseInteractionSpline *self               ,
                                                     const RealArray1D               *charges            ,
                                                     const Real                       electrostaticScale ,
                                                           Coordinates3              *coordinates3A      ,
                                                           Coordinates3              *coordinates3B      ,
                                                           SymmetryParameters        *symmetryParameters ,
                                                           ImagePairListContainer    *imagePairLists     ,
                                                           RealArray1D               *potentials         ,
                                                           Status                    *status             )
{
    if ( ( self                != NULL ) &&
         ( charges             != NULL ) &&
         ( coordinates3A       != NULL ) &&
         ( coordinates3B       != NULL ) &&
         ( symmetryParameters  != NULL ) &&
         ( imagePairLists      != NULL ) &&
         ( potentials          != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator          ,
                                           imagePairLists     ,
                                           coordinates3B      ,
                                           symmetryParameters ,
                                           NULL               ,
                                           NULL               ,
                                           status             ) ;
        if ( Status_IsOK ( status ) )
        {
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                PairwiseInteractionSpline_QCMMPotentials ( self                   ,
                                                           charges                ,
                                                           electrostaticScale * iterator.scale ,
                                                           coordinates3A          ,
                                                           iterator.iCoordinates3 ,
                                                           iterator.pairList      ,
                                                           potentials             ,
                                                           status                 ) ;
            }
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM potentials in atomic units within the minimum image convention.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCMMPotentialsMI ( const PairwiseInteractionSpline *self               ,
                                                  const RealArray1D               *chargesM           ,
                                                  const Real                       electrostaticScale ,
                                                  const Coordinates3              *coordinates3Q      ,
                                                  const Coordinates3              *coordinates3M      ,
                                                        SymmetryParameters        *symmetryParameters ,
                                                        PairList                  *pairList           ,
                                                        RealArray1D               *potentials         ,
                                                        Status                    *status             )
{
    if ( ( self               != NULL    ) &&
         ( chargesM           != NULL    ) &&
         ( coordinates3M      != NULL    ) &&
         ( coordinates3Q      != NULL    ) &&
         ( symmetryParameters != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( pairList           != NULL    ) &&
         ( potentials         != NULL    ) &&
           Status_IsOK ( status ) )
    {
        if ( self->electrostaticSpline != NULL )
        {
            auto Integer           m, n, q ;
            auto CubicSpline      *spline = self->electrostaticSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              cutOff2 = self->cutOff2, eScale, f, p, qM, r2, xQM, yQM, zQM ;
            /* . Minimum image setup. */
            MinimumImage_Allocate ( coordinates3Q, coordinates3M ) ;
            /* . Initialization. */
            eScale  = electrostaticScale / Units_Length_Angstroms_To_Bohrs ;
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                q  = record->index ;
                /* . Displacements. */
                MinimumImage_Displacements ( q, m ) ;
                /* . MM atom. */
       	        for ( n = 0, p = 0.0e+00 ; n < record->capacity ; n++ )
	        {
	            m   = record->indices[n] ;
                    qM  = Array1D_Item ( chargesM, m ) ;
	            Coordinates3_GetRow ( rDisplacements, m, xQM, yQM, zQM ) ;
                    r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_Evaluate ( spline, 0, r2, &f, NULL, NULL, NULL ) ;
                    p  += qM * f ;
	        }
                Array1D_Item ( potentials, q ) += ( p * eScale ) ;
            }
            MinimumImage_Deallocate ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCQCGradients ( const PairwiseInteractionSpline *self               ,
                                               const RealArray1D               *charges            ,
                                               const Real                       electrostaticScale ,
                                               const Coordinates3              *coordinates3I      ,
                                               const Coordinates3              *coordinates3J      ,
                                                     PairList                  *pairList           ,
                                               const Coordinates3              *gradients3I        ,
                                               const Coordinates3              *gradients3J        ,
                                                     Status                    *status             )
{
    if ( ( self               != NULL    ) &&
         ( charges            != NULL    ) &&
         ( coordinates3I      != NULL    ) &&
         ( coordinates3J      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( gradients3I        != NULL    ) &&
         ( gradients3J        != NULL    ) &&
         ( pairList           != NULL    ) &&
         ( Status_IsOK ( status )        ) )
    {
        if ( self->electrostaticSpline != NULL )
        {
            auto Integer           i, j, n ;
            auto CubicSpline      *spline = self->electrostaticSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              cutOff2 = self->cutOff2, eScale, g, gX, gY, gZ, qI, qIJ, r2, xI, xJ, xIJ, yI, yJ, yIJ, zI, zJ, zIJ ;
            eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                i  = record->index ;
                qI = eScale * Array1D_Item ( charges, i ) ;
	        Coordinates3_GetRow ( coordinates3I, i, xI, yI, zI ) ;
       	        for ( n = 0, gX = gY = gZ = 0.0e+00 ; n < record->capacity ; n++ )
	        {
	            j   = record->indices[n] ;
                    qIJ = 2.0e+00 * qI * Array1D_Item ( charges, j ) ; /*. Note the extra factor of 2 here rather than later. */
	            Coordinates3_GetRow ( coordinates3J, j, xJ, yJ, zJ ) ;
                    xIJ = xI - xJ ;
                    yIJ = yI - yJ ;
                    zIJ = zI - zJ ;
                    r2  = ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_Evaluate ( spline, 0, r2, NULL, &g, NULL, NULL ) ;
                    g   *= qIJ ;
                    xIJ *= g   ;
                    yIJ *= g   ;
                    zIJ *= g   ;
                    Coordinates3_DecrementRow ( gradients3J, j, xIJ, yIJ, zIJ ) ;
                    gX  += xIJ ;
                    gY  += yIJ ;
                    gZ  += zIJ ;
	        }
                Coordinates3_IncrementRow ( gradients3I, i, gX, gY, gZ ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/QC gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCQCGradientsImage ( const PairwiseInteractionSpline  *self                       ,
                                                    const RealArray1D                *charges                    ,
                                                    const Real                        electrostaticScale         ,
                                                          Coordinates3               *coordinates3               ,
                                                          SymmetryParameters         *symmetryParameters         ,
                                                          ImagePairListContainer     *imagePairLists             ,
                                                          Coordinates3               *gradients3                 ,
                                                          SymmetryParameterGradients *symmetryParameterGradients ,
                                                          Status                     *status                     )
{
    if ( ( self                       != NULL ) &&
         ( charges                    != NULL ) &&
         ( coordinates3               != NULL ) &&
         ( symmetryParameters         != NULL ) &&
         ( imagePairLists             != NULL ) &&
         ( gradients3                 != NULL ) &&
         ( symmetryParameterGradients != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator                  ,
                                           imagePairLists             ,
                                           coordinates3               ,
                                           symmetryParameters         ,
                                           gradients3                 ,
                                           symmetryParameterGradients ,
                                           status                     ) ;
        if ( Status_IsOK ( status ) )
        {
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                PairwiseInteractionSpline_QCQCGradients ( self                   ,
                                                          charges                ,
                                                          electrostaticScale * iterator.scale ,
                                                          coordinates3           ,
                                                          iterator.iCoordinates3 ,
                                                          iterator.pairList      ,
                                                          gradients3             ,
                                                          iterator.iGradients3   ,
                                                          status                 ) ;
                ImagePairListIterator_Gradients ( &iterator ) ;
            }
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCQCPotentials ( const PairwiseInteractionSpline *self               ,
                                                const Real                       electrostaticScale ,
                                                const Coordinates3              *coordinates3I      ,
                                                const Coordinates3              *coordinates3J      ,
                                                      PairList                  *pairList           ,
                                                      SymmetricMatrix           *potentials         ,
                                                      Status                    *status             )
{
    if ( ( self               != NULL    ) &&
         ( coordinates3I      != NULL    ) &&
         ( coordinates3J      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( pairList           != NULL    ) &&
         ( potentials         != NULL    ) &&
         ( Status_IsOK ( status )        ) )
    {
        if ( self->electrostaticSpline != NULL )
        {
            auto Integer           i, j, n ;
            auto CubicSpline      *spline = self->electrostaticSpline ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              cutOff2 = self->cutOff2, eScale, f, r2, xI, xJ, xIJ, yI, yJ, yIJ, zI, zJ, zIJ ;
            eScale  = electrostaticScale / Units_Length_Angstroms_To_Bohrs ;
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                i = record->index ;
	        Coordinates3_GetRow ( coordinates3I, i, xI, yI, zI ) ;
       	        for ( n = 0 ; n < record->capacity ; n++ )
	        {
	            j   = record->indices[n] ;
	            Coordinates3_GetRow ( coordinates3J, j, xJ, yJ, zJ ) ;
                    xIJ = xI - xJ ;
                    yIJ = yI - yJ ;
                    zIJ = zI - zJ ;
                    r2  = ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                    if ( r2 > cutOff2 ) continue ;
                    CubicSpline_Evaluate ( spline, 0, r2, &f, NULL, NULL, NULL ) ;
                    if ( i != j ) f *= 0.5e+00 ;
                    f *= eScale ;
                    f += SymmetricMatrix_GetItem ( potentials, i, j,    NULL ) ;
                    SymmetricMatrix_SetItem      ( potentials, i, j, f, NULL ) ;
	        }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/QC potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionSpline_QCQCPotentialsImage ( const PairwiseInteractionSpline *self               ,
                                                     const Real                       electrostaticScale ,
                                                           Coordinates3              *coordinates3       ,
                                                           SymmetryParameters        *symmetryParameters ,
                                                           ImagePairListContainer    *imagePairLists     ,
                                                           SymmetricMatrix           *potentials         ,
                                                           Status                    *status             )
{
    if ( ( self               != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( imagePairLists     != NULL ) &&
         ( potentials         != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator          ,
                                           imagePairLists     ,
                                           coordinates3       ,
                                           symmetryParameters ,
                                           NULL               ,
                                           NULL               ,
                                           status             ) ;
        if ( Status_IsOK ( status ) )
        {
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                PairwiseInteractionSpline_QCQCPotentials ( self                   ,
                                                           electrostaticScale * iterator.scale ,
                                                           coordinates3           ,
                                                           iterator.iCoordinates3 ,
                                                           iterator.pairList      ,
                                                           potentials             ,
                                                           status                 ) ;
            }
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
}

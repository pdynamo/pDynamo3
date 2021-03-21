/*==================================================================================================================================
! . Pairwise interactions with ABFS smoothing.
!=================================================================================================================================*/

# include "Boolean.h"
# include "Memory.h"
# include "MinimumImageUtilities.h"
# include "NumericalMacros.h"
# include "PairwiseInteractionABFS.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DefaultDampingCutOff  0.5e+00
# define _DefaultInnerCutOff    8.0e+00
# define _DefaultOuterCutOff   12.0e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    /* . CutOff factors. */
    Real r2Damp   ;
    Real r2Off    ;
    Real r2On     ;
    /* . Electrostatic factors. */
    Real a        ;
    Real b        ;
    Real c        ;
    Real d        ;
    Real qShift1  ;
    Real qShift2  ;
    Real qF0      ;
    Real qAlpha   ;
    /* . Lennard-Jones A factors. */
    Real aF6      ;
    Real aK12     ;
    Real aShift12 ;
    Real aF0      ;
    Real aAlpha   ;
    /* . Lennard-Jones B factors. */
    Real bF3      ;
    Real bK6      ;
    Real bShift6  ;
    Real bF0      ;
    Real bAlpha   ;
} ABFSFactors ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Check distances. */
# define CheckDistances( self, r2, s, s2 ) \
    { \
        if ( r2 > self.r2Off  ) continue ; \
        else if ( r2 < self.r2Damp ) { s2 = 0.0e+00 ; s = 0.0e+00 ; } \
        else { s2 = 1.0e+00 / r2 ; s = sqrt ( s2 ) ; } \
    }

/* . Electrostaticinteraction. */
# define ElectrostaticTerm( self, r2, s, qIJ, f, g ) \
    { \
        if ( r2 > self.r2On ) \
        { \
            f  = qIJ * ( s * ( self.a - r2 * ( self.b + r2 * ( self.c + self.d * r2 ) ) ) + self.qShift2 ) ; \
            g += - qIJ * 0.5e+00 * s * ( self.a + r2 * ( self.b + r2 * ( 3.0e+00 * self.c + 5.0e+00 * self.d * r2 ) ) ) / r2 ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f  = qIJ * ( s + self.qShift1 ) ; \
            g += - qIJ * 0.5e+00 * s / r2 ; \
        } \
        else \
        { \
            f  = qIJ * ( self.qF0 - self.qAlpha * r2 ) ; \
            g += - qIJ * self.qAlpha ; \
        } \
    }

/* . Lennard-Jones interactions - A and B. */
# define LennardJonesTerm( self, r2, s, s2, aIJ, bIJ, f, g ) \
    { \
        auto Real s6 = s2 * s2 * s2 ; \
        if ( r2 > self.r2On ) \
        { \
            auto Real l1 = s6 - self.aF6, l2 = ( s / r2 ) - self.bF3 ; \
            f  = aIJ * self.aK12 * pow ( l1, 2 ) - bIJ * self.bK6 * pow ( l2, 2 ) ; \
            g += - 3.0e+00 * s6 * ( 2.0e+00 * aIJ * self.aK12 * l1 / r2 - bIJ * self.bK6 * l2 / s ) ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f  = aIJ * ( s6 * s6 - self.aShift12 ) - bIJ * ( s6 - self.bShift6 ) ; \
            g += - 3.0e+00 * s6 * ( 2.0e+00 * aIJ * s6 - bIJ ) / r2 ; \
        } \
        else \
        { \
            f  = aIJ * ( self.aF0 - self.aAlpha * r2 ) - bIJ * ( self.bF0 - self.bAlpha * r2 ) ; \
            g += - aIJ * self.aAlpha + bIJ * self.bAlpha ; \
        } \
    }

/* . Lennard-Jones A interaction. */
# define LennardJonesATerm( self, r2, s, s2, aIJ, f, g ) \
    { \
        auto Real s6 = s2 * s2 * s2 ; \
        if ( r2 > self.r2On ) \
        { \
            auto Real l1 = s6 - self.aF6 ; \
            f  =   aIJ * self.aK12 * pow ( l1, 2 ) ; \
            g += - aIJ * self.aK12 * 6.0e+00 * s6 * l1 / r2 ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f  = aIJ * ( s6 * s6 - self.aShift12 ) ; \
            g += - 6.0e+00 * s6 * aIJ * s6 / r2 ; \
        } \
        else \
        { \
            f  = aIJ * ( self.aF0 - self.aAlpha * r2 ) ; \
            g += - aIJ * self.aAlpha ; \
        } \
    }

/* . Lennard-Jones B interaction. */
# define LennardJonesBTerm( self, r2, s, s2, bIJ, f, g ) \
    { \
        auto Real s6 = s2 * s2 * s2 ; \
        if ( r2 > self.r2On ) \
        { \
            auto Real l2 = ( s / r2 ) - self.bF3 ; \
            f  = - bIJ * self.bK6 * pow ( l2, 2 ) ; \
            g +=   bIJ * self.bK6 * 3.0e+00 * s6 * l2 / s ; \
        } \
        else if ( r2 > self.r2Damp ) \
        { \
            f  = - bIJ * ( s6 - self.bShift6 ) ; \
            g +=   bIJ * 3.0e+00 * s6 / r2 ; \
        } \
        else \
        { \
            f  = - bIJ * ( self.bF0 - self.bAlpha * r2 ) ; \
            g +=   bIJ * self.bAlpha ; \
        } \
    }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionABFS_CopyTo            ( const PairwiseInteractionABFS *self, PairwiseInteractionABFS *other ) ;
static void PairwiseInteractionABFS_Initialize        (       PairwiseInteractionABFS *self ) ;
static void PairwiseInteractionABFS_InitializeFactors ( const PairwiseInteractionABFS *self, ABFSFactors *factors ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionABFS *PairwiseInteractionABFS_Allocate ( Status *status )
{
    PairwiseInteractionABFS *self = NULL ;
    self = Memory_AllocateType ( PairwiseInteractionABFS ) ;
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    else PairwiseInteractionABFS_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionABFS *PairwiseInteractionABFS_Clone ( PairwiseInteractionABFS *self, Status *status )
{
    PairwiseInteractionABFS *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = PairwiseInteractionABFS_Allocate ( status ) ;
        PairwiseInteractionABFS_CopyTo ( self, clone ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionABFS_CopyTo ( const PairwiseInteractionABFS *self, PairwiseInteractionABFS *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->dampingCutOff = self->dampingCutOff ;
        other->innerCutOff   = self->innerCutOff   ;
        other->outerCutOff   = self->outerCutOff   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_Deallocate ( PairwiseInteractionABFS **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionABFS_Initialize ( PairwiseInteractionABFS *self )
{
    if ( self != NULL )
    {
        self->dampingCutOff = _DefaultDampingCutOff ;
        self->innerCutOff   = _DefaultInnerCutOff   ;
        self->outerCutOff   = _DefaultOuterCutOff   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize factors.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionABFS_InitializeFactors ( const PairwiseInteractionABFS *self, ABFSFactors *factors )
{
    if ( ( self != NULL ) && ( factors != NULL ) )
    {
        auto Real f, g, gamma, s12, s6 ;
        /* . CutOff factors. */
        factors->r2Damp   = self->dampingCutOff * self->dampingCutOff ;
        factors->r2Off    = self->outerCutOff   * self->outerCutOff   ;
        factors->r2On     = self->innerCutOff   * self->innerCutOff   ;
        /* . Electrostaticfactors. */
        /* . Interaction. */
        gamma             = pow ( factors->r2Off - factors->r2On, 3 ) ;
        factors->a        = factors->r2Off * factors->r2Off * ( factors->r2Off - 3.0e+00 * factors->r2On ) / gamma ;
        factors->b        = 6.0e+00 * factors->r2Off * factors->r2On / gamma ;
        factors->c        = - ( factors->r2Off + factors->r2On )     / gamma ;
        factors->d        = 0.4e+00                                  / gamma ;
        factors->qShift1  = 8.0e+00 * ( factors->r2Off * factors->r2On * ( self->outerCutOff - self->innerCutOff ) -
                            0.2e+00 * ( self->outerCutOff * factors->r2Off * factors->r2Off - self->innerCutOff * factors->r2On * factors->r2On ) ) / gamma ;
        factors->qShift2  = - ( factors->a / self->outerCutOff ) + factors->b * self->outerCutOff + factors->c * self->outerCutOff * factors->r2Off +
                                                                                           factors->d * self->outerCutOff * factors->r2Off * factors->r2Off ;
        /* . Damping. */
        f                 =   1.0e+00 / self->dampingCutOff + factors->qShift1 ;
        g                 = - 1.0e+00 / factors->r2Damp ;
        factors->qF0      = f - 0.5e+00 * self->dampingCutOff * g ;
        factors->qAlpha   =   - 0.5e+00 * g / self->dampingCutOff ;
        /* . Lennard-Jones A factors. */
        /* . Interaction. */
        factors->aF6      = 1.0e+00 / ( factors->r2Off * factors->r2Off * factors->r2Off ) ;
        factors->aK12     = pow ( factors->r2Off, 3 ) / ( pow ( factors->r2Off, 3 ) - pow ( factors->r2On, 3 ) ) ;
        factors->aShift12 = 1.0e+00 / pow ( self->innerCutOff * self->outerCutOff, 6 ) ;
        /* . Damping. */
        s12               = 1.0e+00 / pow ( self->dampingCutOff, 12 ) ;
        f                 = s12 - factors->aShift12 ;
        g                 = - 12.0e+00 * s12 / self->dampingCutOff ;
        factors->aF0      = f - 0.5e+00 * self->dampingCutOff * g ;
        factors->aAlpha   =   - 0.5e+00 * g / self->dampingCutOff ;
        /* . Lennard-Jones B factors. */
        /* . Interaction. */
        factors->bF3      = 1.0e+00 / ( self->outerCutOff * factors->r2Off ) ;
        factors->bK6      = ( self->outerCutOff * factors->r2Off ) / ( self->outerCutOff * factors->r2Off - self->innerCutOff * factors->r2On ) ;
        factors->bShift6  = 1.0e+00 / pow ( self->innerCutOff * self->outerCutOff, 3 ) ;
        /* . Damping. */
        s6                = 1.0e+00 / pow ( self->dampingCutOff, 6 ) ;
        f                 = s6 - factors->bShift6 ;
        g                 = - 6.0e+00 * s6 / self->dampingCutOff ;
        factors->bF0      = f - 0.5e+00 * self->dampingCutOff * g ;
        factors->bAlpha   =   - 0.5e+00 * g / self->dampingCutOff ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_Interactions ( const PairwiseInteractionABFS *self , 
                                            const RealArray1D *r                ,
                                                  RealArray1D *electrostatic    ,
                                                  RealArray1D *lennardJonesA    ,
                                                  RealArray1D *lennardJonesB    )
{
    Boolean doElectrostatic = ( electrostatic != NULL ) ,
            doLennardJonesA = ( lennardJonesA != NULL ) ,
            doLennardJonesB = ( lennardJonesB != NULL ) ;
    if ( ( self != NULL ) && ( r != NULL ) && ( doElectrostatic|| doLennardJonesA || doLennardJonesB ) )
    {
        auto ABFSFactors factors ;
        auto Integer     i ;
        auto Real        fA, fB, fQ, g, s, s2, x, x2 ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
        RealArray1D_Set ( electrostatic, 0.0e+00 ) ;
        RealArray1D_Set ( lennardJonesA, 0.0e+00 ) ;
        RealArray1D_Set ( lennardJonesB, 0.0e+00 ) ;
        for ( i = 0 ; i < View1D_Extent ( r ) ; i++ )
        {
            x  = Array1D_Item ( r, i ) ;
            x2 = x * x ;
            CheckDistances ( factors, x2, s, s2 ) ;
            ElectrostaticTerm ( factors, x2, s,     1.0e+00, fQ, g ) ;
            LennardJonesATerm ( factors, x2, s, s2, 1.0e+00, fA, g ) ;
            LennardJonesBTerm ( factors, x2, s, s2, 1.0e+00, fB, g ) ;
            if ( doElectrostatic ) Array1D_Item ( electrostatic, i ) = fQ ;
            if ( doLennardJonesA ) Array1D_Item ( lennardJonesA, i ) = fA ;
            if ( doLennardJonesB ) Array1D_Item ( lennardJonesB, i ) = fB ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM energy and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_MMMMEnergy ( const PairwiseInteractionABFS *self               ,
                                          const RealArray1D             *chargesI           ,
                                          const RealArray1D             *chargesJ           ,
                                                IntegerArray1D          *ljTypesI           ,
                                                IntegerArray1D          *ljTypesJ           ,
                                          const LJParameterContainer    *ljParameters       ,
                                          const Real                     electrostaticScale ,
                                          const Real                     lennardJonesScale  ,
                                          const Coordinates3            *coordinates3I      ,
                                          const Coordinates3            *coordinates3J      ,
                                                PairList                *pairList           ,
                                                Real                    *eElectrostatic     ,
                                                Real                    *eLennardJones      ,
                                                Coordinates3            *gradients3I        ,
                                                Coordinates3            *gradients3J        ,
                                                Status                  *status             )
{
    if ( eElectrostatic!= NULL ) (*eElectrostatic) = 0.0e+00 ;
    if ( eLennardJones  != NULL ) (*eLennardJones ) = 0.0e+00 ;
    if ( ( self          != NULL ) &&
         ( coordinates3I != NULL ) &&
         ( coordinates3J != NULL ) &&
         ( pairList      != NULL ) &&
           Status_IsOK ( status ) )
    {
        auto Boolean doElectrostatic, doGradients, doLennardJones ;
        doElectrostatic = ( chargesI           != NULL    ) &&
                          ( chargesJ           != NULL    ) &&
                          ( eElectrostatic     != NULL    ) &&
                          ( electrostaticScale != 0.0e+00 ) ;
        doGradients     = ( gradients3I        != NULL    ) && 
                          ( gradients3J        != NULL    ) ;
        doLennardJones  = ( eLennardJones      != NULL    ) &&
                          ( ljTypesI           != NULL    ) &&
                          ( ljTypesJ           != NULL    ) &&
                          ( ljParameters       != NULL    ) &&
                          ( lennardJonesScale  != 0.0e+00 ) ;
        if ( doElectrostatic|| doLennardJones )
        {
            auto ABFSFactors       factors ;
            auto Integer           i, j, n, numberOfLJTypes = 0, tI = 0, tIJ ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              aIJ, bIJ, eScale, f, g, qI = 0.0e+00, qIJ, r2, s, s2, xI, xIJ, xJ, yI, yIJ, yJ, zI, zIJ, zJ ;
            auto Real              eLJ = 0.0e+00, eQQ = 0.0e+00 ;
            /* . Initialization. */
            eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            if ( doLennardJones  ) numberOfLJTypes = ljParameters->ntypes ;
            PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
            /* . Loop over record. */
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                /* . First atom. */
                i  = record->index ;
	        if ( doElectrostatic ) qI = eScale          * Array1D_Item ( chargesI, i ) ;
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
                    CheckDistances ( factors, r2, s, s2 ) ;
                    f = 0.0e+00 ;
                    g = 0.0e+00 ;
                    if ( doElectrostatic)
                    {
                        qIJ = qI * Array1D_Item ( chargesJ, j ) ;
                        ElectrostaticTerm ( factors, r2, s, qIJ, f, g ) ;
                        eQQ += f ;
                    }
                    if ( doLennardJones )
                    {
	                tIJ = ljParameters->tableindex[tI+Array1D_Item ( ljTypesJ, j )] ;
                        aIJ = ljParameters->tableA[tIJ] * lennardJonesScale ;
                        bIJ = ljParameters->tableB[tIJ] * lennardJonesScale ;
                        LennardJonesTerm ( factors, r2, s, s2, aIJ, bIJ, f, g ) ;
                        eLJ += f ;
                    }
                    if ( doGradients )
                    {
                        g   *= 2.0e+00 ;
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
void PairwiseInteractionABFS_MMMMEnergyImage ( const PairwiseInteractionABFS    *self                       ,
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
                PairwiseInteractionABFS_MMMMEnergy ( self                   ,
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
void PairwiseInteractionABFS_MMMMEnergyMI ( const PairwiseInteractionABFS    *self                       ,
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
    if ( eElectrostatic!= NULL ) (*eElectrostatic) = 0.0e+00 ;
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
                          ( electrostaticScale         != 0.0e+00 ) ;
        doGradients     = ( gradients3I                != NULL    ) && 
                          ( gradients3J                != NULL    ) &&
                          ( symmetryParameterGradients != NULL    ) ;
        doLennardJones  = ( eLennardJones              != NULL    ) &&
                          ( ljTypesI                   != NULL    ) &&
                          ( ljTypesJ                   != NULL    ) &&
                          ( ljParameters               != NULL    ) &&
                          ( lennardJonesScale          != 0.0e+00 ) ;
        if ( doElectrostatic|| doLennardJones )
        {
            auto ABFSFactors       factors ;
            auto Integer           i, j, n, numberOfLJTypes = 0, tI = 0, tIJ ;
            auto PairListIterator  iterator ;
            auto PairRecord       *record ;
            auto Real              aIJ, bIJ, eScale, f, g, qI = 0.0e+00, qIJ, r2, s, s2, xIJ, yIJ, zIJ ;
            auto Real              eLJ = 0.0e+00, eQQ = 0.0e+00 ;
            /* . Minimum image setup. */
            MinimumImage_Allocate ( coordinates3I, coordinates3J ) ;
            /* . Initialization. */
            eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            if ( doLennardJones  ) numberOfLJTypes = ljParameters->ntypes ;
            PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
            /* . Loop over record. */
            PairListIterator_Initialize ( &iterator, pairList ) ;
            while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
            {
                /* . First atom. */
                i = record->index ;
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
                    CheckDistances ( factors, r2, s, s2 ) ;
                    f = 0.0e+00 ;
                    g = 0.0e+00 ;
                    if ( doElectrostatic )
                    {
                        qIJ = qI * Array1D_Item ( chargesJ, j ) ;
                        ElectrostaticTerm ( factors, r2, s, qIJ, f, g ) ;
                        eQQ += f ;
                    }
                    if ( doLennardJones )
                    {
	                tIJ = ljParameters->tableindex[tI+Array1D_Item ( ljTypesJ, j )] ;
                        aIJ = ljParameters->tableA[tIJ] * lennardJonesScale ;
                        bIJ = ljParameters->tableB[tIJ] * lennardJonesScale ;
                        LennardJonesTerm ( factors, r2, s, s2, aIJ, bIJ, f, g ) ;
                        eLJ += f ;
                    }
                    if ( doGradients )
                    {
                        g   *= 2.0e+00 ;
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
void PairwiseInteractionABFS_QCMMGradients ( const PairwiseInteractionABFS *self               ,
                                             const RealArray1D             *chargesQ           ,
                                             const RealArray1D             *chargesM           ,
                                             const Real                     electrostaticScale ,
                                             const Coordinates3            *coordinates3Q      ,
                                             const Coordinates3            *coordinates3M      ,
                                                   PairList                *pairList           ,
                                             const Coordinates3            *gradients3Q        ,
                                             const Coordinates3            *gradients3M        ,
                                                   Status                  *status             )
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
        auto ABFSFactors              factors ;
        auto Integer                  m, n, q ;
        auto PairListIterator         iterator ;
        auto PairRecord              *record ;
        auto Real                     eScale, f, g, gX, gY, gZ, qQ, qQM, r2, s, s2, xM, xQ, xQM, yM, yQ, yQM, zM, zQ, zQM ;
        eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
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
                f   = 0.0e+00 ;
                g   = 0.0e+00 ;
                CheckDistances    ( factors, r2, s, s2        ) ;
                ElectrostaticTerm ( factors, r2, s, qQM, f, g ) ;
                xQM *= g ;
                yQM *= g ;
                zQM *= g ;
                Coordinates3_DecrementRow ( gradients3M, m, xQM, yQM, zQM ) ;
                gX  += xQM ;
                gY  += yQM ;
                gZ  += zQM ;
	    }
            Coordinates3_IncrementRow ( gradients3Q, q, gX, gY, gZ ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/MM gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCMMGradientsImage ( const PairwiseInteractionABFS    *self                       ,
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
                PairwiseInteractionABFS_QCMMGradients ( self                   ,
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
void PairwiseInteractionABFS_QCMMGradientsMI ( const PairwiseInteractionABFS    *self                       ,
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
        auto ABFSFactors              factors ;
        auto Integer                  m, n, q ;
        auto PairListIterator         iterator ;
        auto PairRecord              *record ;
        auto Real                     eScale, f, g, gX, gY, gZ, qQ, qQM, r2, s, s2, xQM, yQM, zQM ;
        /* . Minimum image setup. */
        MinimumImage_Allocate ( coordinates3Q, coordinates3M ) ;
        /* . Initialization. */
        eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
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
                f   = 0.0e+00 ;
                g   = 0.0e+00 ;
                CheckDistances    ( factors, r2, s, s2        ) ;
                ElectrostaticTerm ( factors, r2, s, qQM, f, g ) ;
                xQM *= g ;
                yQM *= g ;
                zQM *= g ;
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

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM potentials in atomic units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCMMPotentials ( const PairwiseInteractionABFS *self               ,
                                              const RealArray1D             *chargesM           ,
                                              const Real                     electrostaticScale ,
                                              const Coordinates3            *coordinates3Q      ,
                                              const Coordinates3            *coordinates3M      ,
                                                    PairList                *pairList           ,
                                                    RealArray1D             *potentials         ,
                                                    Status                  *status             )
{
    if ( ( self               != NULL    ) &&
         ( chargesM           != NULL    ) &&
         ( coordinates3M      != NULL    ) &&
         ( coordinates3Q      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( pairList           != NULL    ) &&
         ( potentials         != NULL    ) &&
           Status_IsOK ( status )        )
    {
        auto ABFSFactors              factors ;
        auto Integer                  m, n, q ;
        auto PairListIterator         iterator ;
        auto PairRecord              *record ;
        auto Real                     eScale, f, g, p, qM, r2, s, s2, xM, xQ, xQM, yM, yQ, yQM, zM, zQ, zQM ;
        eScale  = electrostaticScale / Units_Length_Angstroms_To_Bohrs ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
        PairListIterator_Initialize ( &iterator, pairList ) ;
        while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
        {
            q  = record->index ;
	    Coordinates3_GetRow ( coordinates3Q, q, xQ, yQ, zQ ) ;
       	    for ( n = 0, g = p = 0.0e+00 ; n < record->capacity ; n++ )
	    {
	        m   = record->indices[n] ;
                qM  = Array1D_Item ( chargesM, m ) ;
	        Coordinates3_GetRow ( coordinates3M, m, xM, yM, zM ) ;
                xQM = xQ - xM ;
                yQM = yQ - yM ;
                zQM = zQ - zM ;
                r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                CheckDistances    ( factors, r2, s, s2       ) ;
                ElectrostaticTerm ( factors, r2, s, qM, f, g ) ;
                p += f ;
	    }
            Array1D_Item ( potentials, q ) += ( p * eScale ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/MM potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCMMPotentialsImage ( const PairwiseInteractionABFS *self               ,
                                                   const RealArray1D             *charges            ,
                                                   const Real                     electrostaticScale ,
                                                         Coordinates3            *coordinates3A      ,
                                                         Coordinates3            *coordinates3B      ,
                                                         SymmetryParameters      *symmetryParameters ,
                                                         ImagePairListContainer  *imagePairLists     ,
                                                         RealArray1D             *potentials         ,
                                                         Status                  *status             )
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
                PairwiseInteractionABFS_QCMMPotentials ( self                   ,
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
void PairwiseInteractionABFS_QCMMPotentialsMI ( const PairwiseInteractionABFS *self               ,
                                                const RealArray1D             *chargesM           ,
                                                const Real                     electrostaticScale ,
                                                const Coordinates3            *coordinates3Q      ,
                                                const Coordinates3            *coordinates3M      ,
                                                      SymmetryParameters      *symmetryParameters ,
                                                      PairList                *pairList           ,
                                                      RealArray1D             *potentials         ,
                                                      Status                  *status             )
{
    if ( ( self               != NULL    ) &&
         ( chargesM           != NULL    ) &&
         ( coordinates3M      != NULL    ) &&
         ( coordinates3Q      != NULL    ) &&
         ( symmetryParameters != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( pairList           != NULL    ) &&
         ( potentials         != NULL    ) &&
           Status_IsOK ( status )        )
    {
        auto ABFSFactors              factors ;
        auto Integer                  m, n, q ;
        auto PairListIterator         iterator ;
        auto PairRecord              *record ;
        auto Real                     eScale, f, g, p, qM, r2, s, s2, xQM, yQM, zQM ;
        /* . Minimum image setup. */
        MinimumImage_Allocate ( coordinates3Q, coordinates3M ) ;
        /* . Initialization. */
        eScale  = electrostaticScale / Units_Length_Angstroms_To_Bohrs ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
        PairListIterator_Initialize ( &iterator, pairList ) ;
        while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
        {
            q  = record->index ;
            /* . Displacements. */
            MinimumImage_Displacements ( q, m ) ;
            /* . MM atom. */
       	    for ( n = 0, g = p = 0.0e+00 ; n < record->capacity ; n++ )
	    {
	        m   = record->indices[n] ;
                qM  = Array1D_Item ( chargesM, m ) ;
	        Coordinates3_GetRow ( rDisplacements, m, xQM, yQM, zQM ) ;
                r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                CheckDistances    ( factors, r2, s, s2       ) ;
                ElectrostaticTerm ( factors, r2, s, qM, f, g ) ;
                p += f ;
	    }
            Array1D_Item ( potentials, q ) += ( p * eScale ) ;
        }
        MinimumImage_Deallocate ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/QC gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCQCGradients ( const PairwiseInteractionABFS *self               ,
                                             const RealArray1D             *charges            ,
                                             const Real                     electrostaticScale ,
                                             const Coordinates3            *coordinates3I      ,
                                             const Coordinates3            *coordinates3J      ,
                                                   PairList                *pairList           ,
                                             const Coordinates3            *gradients3I        ,
                                             const Coordinates3            *gradients3J        ,
                                                   Status                  *status             )
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
        auto ABFSFactors              factors ;
        auto Integer                  i, j, n ;
        auto PairListIterator         iterator ;
        auto PairRecord              *record ;
        auto Real                     eScale, f, g, gX, gY, gZ, qI, qIJ, r2, s, s2, xI, xJ, xIJ, yI, yJ, yIJ, zI, zJ, zIJ ;
        eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
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
                f   = 0.0e+00 ;
                g   = 0.0e+00 ;
                CheckDistances    ( factors, r2, s, s2        ) ;
                ElectrostaticTerm ( factors, r2, s, qIJ, f, g ) ;
                xIJ *= g ;
                yIJ *= g ;
                zIJ *= g ;
                Coordinates3_DecrementRow ( gradients3J, j, xIJ, yIJ, zIJ ) ;
                gX  += xIJ ;
                gY  += yIJ ;
                gZ  += zIJ ;
	    }
            Coordinates3_IncrementRow ( gradients3I, i, gX, gY, gZ ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/QC gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCQCGradientsImage ( const PairwiseInteractionABFS    *self                       ,
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
                PairwiseInteractionABFS_QCQCGradients ( self                   ,
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
void PairwiseInteractionABFS_QCQCPotentials ( const PairwiseInteractionABFS *self               ,
                                              const Real                     electrostaticScale ,
                                              const Coordinates3            *coordinates3I      ,
                                              const Coordinates3            *coordinates3J      ,
                                                    PairList                *pairList           ,
                                                    SymmetricMatrix         *potentials         ,
                                                    Status                  *status             )
{
    if ( ( self               != NULL    ) &&
         ( coordinates3I      != NULL    ) &&
         ( coordinates3J      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( pairList           != NULL    ) &&
         ( potentials         != NULL    ) &&
         ( Status_IsOK ( status )        ) )
    {
        auto ABFSFactors              factors ;
        auto Integer                  i, j, n ;
        auto PairListIterator         iterator ;
        auto PairRecord              *record ;
        auto Real                     eScale, f, g, r2, s, s2, xI, xJ, xIJ, yI, yJ, yIJ, zI, zJ, zIJ ;
        eScale  = electrostaticScale / Units_Length_Angstroms_To_Bohrs ;
        PairwiseInteractionABFS_InitializeFactors ( self, &factors ) ;
        PairListIterator_Initialize ( &iterator, pairList ) ;
        while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
        {
            i = record->index ;
	    Coordinates3_GetRow ( coordinates3I, i, xI, yI, zI ) ;
       	    for ( n = 0, g = 0.0e+00 ; n < record->capacity ; n++ )
	    {
	        j   = record->indices[n] ;
	        Coordinates3_GetRow ( coordinates3J, j, xJ, yJ, zJ ) ;
                xIJ = xI - xJ ;
                yIJ = yI - yJ ;
                zIJ = zI - zJ ;
                r2  = ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                CheckDistances    ( factors, r2, s, s2           ) ;
                ElectrostaticTerm ( factors, r2, s, eScale, f, g ) ;
                if ( i != j ) f *= 0.5e+00 ;
                f += SymmetricMatrix_GetItem ( potentials, i, j,    NULL ) ;
                SymmetricMatrix_SetItem      ( potentials, i, j, f, NULL ) ;
	    }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/QC potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionABFS_QCQCPotentialsImage ( const PairwiseInteractionABFS *self               ,
                                                   const Real                     electrostaticScale ,
                                                         Coordinates3            *coordinates3       ,
                                                         SymmetryParameters      *symmetryParameters ,
                                                         ImagePairListContainer  *imagePairLists     ,
                                                         SymmetricMatrix         *potentials         ,
                                                         Status                  *status             )
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
                PairwiseInteractionABFS_QCQCPotentials ( self                   ,
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

/*==================================================================================================================================
! . Pairwise interactions with no energy modification except at short range.
!=================================================================================================================================*/

/*# include <stdio.h>*/

# include "Boolean.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "PairwiseInteractionFull.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DefaultDampingCutOff 0.5e+00 /* . Always assumed to be greater than zero. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionFull_CopyTo     ( const PairwiseInteractionFull *self , PairwiseInteractionFull *other ) ;
static void PairwiseInteractionFull_Initialize (       PairwiseInteractionFull *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionFull *PairwiseInteractionFull_Allocate ( Status *status )
{
    PairwiseInteractionFull *self = NULL ;
    self = Memory_AllocateType ( PairwiseInteractionFull ) ;
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    else PairwiseInteractionFull_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionFull *PairwiseInteractionFull_Clone ( PairwiseInteractionFull *self, Status *status )
{
    PairwiseInteractionFull *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = PairwiseInteractionFull_Allocate ( status ) ;
        PairwiseInteractionFull_CopyTo ( self, clone ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionFull_CopyTo ( const PairwiseInteractionFull *self, PairwiseInteractionFull *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->dampingCutOff = self->dampingCutOff ;
        PairwiseInteractionFull_InitializeDependent ( other ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionFull_Deallocate ( PairwiseInteractionFull **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PairwiseInteractionFull_Initialize ( PairwiseInteractionFull *self )
{
    if ( self != NULL )
    {
        self->dampingCutOff = _DefaultDampingCutOff  ;
        PairwiseInteractionFull_InitializeDependent ( self ) ;
    }
}

void PairwiseInteractionFull_InitializeDependent ( PairwiseInteractionFull *self )
{
    if ( self != NULL )
    {
        if ( self->dampingCutOff > 0.0 )
        {
            auto Real c = self->dampingCutOff, cF, f, g ;
            /* . 1. */
            cF = c ;
            f             =   1.0e+00 / cF ;
            g             = - f / c ;
            self->alpha1  =   0.5e+00 * g / c ;
            self->beta1   = - 0.5e+00 * g * c + f ;
            /* . 2. */
            cF *= c ;
            f             =   1.0e+00 / cF ;
            g             = - 2.0e+00 * f / c ;
            self->alpha2  =   0.5e+00 * g / c ;
            self->beta2   = - 0.5e+00 * g * c + f ;
            /* . 3. */
            cF *= c ;
            f             =   1.0e+00 / cF ;
            g             = - 3.0e+00 * f / c ;
            self->alpha3  =   0.5e+00 * g / c ;
            self->beta3   = - 0.5e+00 * g * c + f ;
            /* . 6 * -1. */
            cF = cF * cF ;
            f             = - 1.0e+00 / cF ;
            g             = - 6.0e+00  * f / c ;
            self->alpha6  =   0.5e+00  * g / c ;
            self->beta6   = - 0.5e+00  * g * c + f ;
            /* . 12. */
            cF = cF * cF ;
            f             =   1.0e+00 / cF ;
            g             = - 12.0e+00 * f / c ;
            self->alpha12 =   0.5e+00  * g / c ;
            self->beta12  = - 0.5e+00  * g * c + f ;
        }
        else
        {
            self->alpha1  = 0.0e+00 ; self->beta1  = 0.0e+00 ;
            self->alpha2  = 0.0e+00 ; self->beta2  = 0.0e+00 ;
            self->alpha3  = 0.0e+00 ; self->beta3  = 0.0e+00 ;
            self->alpha6  = 0.0e+00 ; self->beta6  = 0.0e+00 ;
            self->alpha12 = 0.0e+00 ; self->beta12 = 0.0e+00 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionFull_Interactions ( const PairwiseInteractionFull *self          , 
                                            const RealArray1D             *r             ,
                                                  RealArray1D             *lennardJonesA ,
                                                  RealArray1D             *lennardJonesB ,
                                                  RealArray1D             *multipole0    ,
                                                  RealArray1D             *multipole1    ,
                                                  RealArray1D             *multipole2    )
{
    if (   ( self          != NULL ) &&
           ( r             != NULL ) &&
         ( ( lennardJonesA != NULL ) ||
           ( lennardJonesB != NULL ) ||
           ( multipole0    != NULL ) ||
           ( multipole1    != NULL ) ||
           ( multipole2    != NULL ) ) )
    {
        auto Integer  i ;
        auto Real     f1, f2, f3, f6, f12, s, s2, s6, x, x2 ;
        for ( i = 0 ; i < View1D_Extent ( r ) ; i++ )
        {
            x  = Array1D_Item ( r, i ) ;
            x2 = x * x ;
            if ( x <= self->dampingCutOff )
            {
                f1  = self->alpha1  * x2 + self->beta1  ;
                f2  = self->alpha2  * x2 + self->beta2  ;
                f3  = self->alpha3  * x2 + self->beta3  ;
                f6  = self->alpha6  * x2 + self->beta6  ;
                f12 = self->alpha12 * x2 + self->beta12 ;
            }
            else
            {
                s2  = 1.0e+00 / x2 ;
                s   = sqrt ( s2 ) ;
                s6  = s2 * s2 * s2 ;
                f1  = s  ;
                f2  = s2 ;
                f3  = s * s2 ;
                f6  = - s6 ;
                f12 = s6 * s6 ;
            }
            if ( lennardJonesA != NULL ) Array1D_Item ( lennardJonesA, i ) = f12 ;
            if ( lennardJonesB != NULL ) Array1D_Item ( lennardJonesB, i ) = f6  ;
            if ( multipole0    != NULL ) Array1D_Item ( multipole0   , i ) = f1  ;
            if ( multipole1    != NULL ) Array1D_Item ( multipole1   , i ) = f2  ;
            if ( multipole2    != NULL ) Array1D_Item ( multipole2   , i ) = f3  ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM energy and gradients with single sets of parameters, coordinates and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionFull_MMMMEnergy ( const PairwiseInteractionFull *self               , 
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
    if ( eElectrostatic != NULL ) (*eElectrostatic) = 0.0e+00 ;
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
            auto Integer     i, j, n, numberOfLJTypes = 0, r, tI = 0, tIJ ;
            auto PairRecord *record ;
            auto Real        aIJ, alpha, beta, bIJ, cutOff2, eScale, f, g, qI = 0.0e+00, qIJ, r2, s, s2, s6, xI, xIJ, xJ, yI, yIJ, yJ, zI, zIJ, zJ ;
            auto Real        eLJ = 0.0e+00, eQQ = 0.0e+00 ;
            /* . Initialization. */
            cutOff2 = Maximum ( 0.0e+00 , self->dampingCutOff * self->dampingCutOff ) ;
            eScale  = electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ;
            if ( doLennardJones  ) numberOfLJTypes = ljParameters->ntypes ;
            /* . Loop over records. */
            for ( r = 0 ; r < pairList->count ; r++ )
            {
                record = PairList_GetRecord ( pairList, r ) ;
                /* . First atom. */
                i  = record->index ;
	        if ( doElectrostatic) qI = eScale          * Array1D_Item   ( chargesI, i ) ;
                if ( doLennardJones  ) tI = numberOfLJTypes * Array1D_Item ( ljTypesI, i ) ;
	        Coordinates3_GetRow ( coordinates3I, i, xI, yI, zI ) ;
                /* . Second atom. */
       	        for ( n = 0 ; n < record->capacity ; n++ )
	        {
	            j   = record->indices[n] ;
	            Coordinates3_GetRow ( coordinates3J, j, xJ, yJ, zJ ) ;
                    xIJ = xI - xJ ;
                    yIJ = yI - yJ ;
                    zIJ = zI - zJ ;
                    r2  = ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                    g   = 0.0e+00 ;
                    /* . Damping. */
                    if ( r2 <= cutOff2 )
                    {
                        if ( doElectrostatic)
                        {
                            qIJ   = qI  * Array1D_Item ( chargesJ, j ) ;
                            alpha = qIJ * self->alpha1 ;
                            beta  = qIJ * self->beta1  ;
                            eQQ  += ( alpha * r2 + beta ) ;
                            g    +=   alpha ;
                        }
                        if ( doLennardJones )
                        {
	                    tIJ   = ljParameters->tableindex[tI+Array1D_Item ( ljTypesJ, j )] ;
                            aIJ   = ljParameters->tableA[tIJ] * lennardJonesScale ;
                            bIJ   = ljParameters->tableB[tIJ] * lennardJonesScale ;
                            alpha = aIJ * self->alpha12 + bIJ * self->alpha6 ;
                            beta  = aIJ * self->beta12  + bIJ * self->beta6  ;
                            eLJ  += ( alpha * r2 + beta ) ;
                            g    +=   alpha ;
                        }
                    }
                    /* . Full. */
                    else
                    {
                        s2 = 1.0e+00 / r2 ;
                        s  = sqrt ( s2 ) ;
                        f  = 0.0e+00 ;
                        if ( doElectrostatic)
                        {
                            qIJ  = qI * Array1D_Item ( chargesJ, j ) ;
                            f    = qIJ * s  ;
                            g   -= 0.5e+00 * f * s2 ;
                            eQQ += f ;
                        }
                        if ( doLennardJones )
                        {

	                    tIJ  = ljParameters->tableindex[tI+Array1D_Item ( ljTypesJ, j )] ;
                            aIJ  = ljParameters->tableA[tIJ] * lennardJonesScale ;
                            bIJ  = ljParameters->tableB[tIJ] * lennardJonesScale ;
                            s6   = s2 * s2 * s2 ;
                            f    = ( aIJ * s6 - bIJ ) * s6 ;
                            g   -= 3.0e+00 * s2 * ( aIJ * s6 * s6 + f ) ;
                            eLJ += f ;
                        }
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
! . QC/MM gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionFull_QCMMGradients ( const PairwiseInteractionFull *self               ,
                                             const Integer                  multipoleOrder     ,
                                             const RealArray1D             *multipolesQ        ,
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
         ( multipolesQ        != NULL    ) &&
         ( chargesM           != NULL    ) &&
         ( coordinates3M      != NULL    ) &&
         ( coordinates3Q      != NULL    ) &&
         ( electrostaticScale != 0.0e+00 ) &&
         ( gradients3M        != NULL    ) &&
         ( gradients3Q        != NULL    ) &&
         ( pairList           != NULL    ) &&
         ( Status_IsOK ( status )        ) )
    {
        auto Integer     d = Coordinates3_Rows ( coordinates3Q ), m, n, order, q, r ;
        auto PairRecord *record ;

        auto Real        a, b, cutOff, eScale, dX, dY, dZ, f2, f3, g1, g2, g3, gX, gXt, gY, gYt, gZ, gZt,
                         q0, qM, qX  = 0.0e+00, qY  = 0.0e+00, qZ  = 0.0e+00,
                                 qXX = 0.0e+00, qXY = 0.0e+00, qXZ = 0.0e+00,
                                 qYY = 0.0e+00, qYZ = 0.0e+00, qZZ = 0.0e+00,
                         r1, r2, scale1, scale2, scale3, t, xM, xQ, yM, yQ, zM, zQ ;

        cutOff = Maximum ( 0.0e+00, self->dampingCutOff ) ;
        eScale = electrostaticScale * Units_Energy_Hartrees_To_Kilojoules_Per_Mole ;
        order  = Maximum ( Minimum ( multipoleOrder, 2 ), 0 ) ;
        scale1 =       Units_Length_Angstroms_To_Bohrs      ;
        scale2 = pow ( Units_Length_Angstroms_To_Bohrs, 2 ) ;
        scale3 = pow ( Units_Length_Angstroms_To_Bohrs, 3 ) ;
        for ( r = 0 ; r < pairList->count ; r++ )
        {
            record = PairList_GetRecord ( pairList, r ) ;
            q      = record->index ;
            q0     = eScale * Array1D_Item ( multipolesQ, q ) ;
            if ( order > 0 )
            {
                qX = eScale * Array1D_Item ( multipolesQ, q +   d ) ;
                qY = eScale * Array1D_Item ( multipolesQ, q + 2*d ) ;
                qZ = eScale * Array1D_Item ( multipolesQ, q + 3*d ) ;
                if ( order > 1 )
                {
                    qXX = eScale * Array1D_Item ( multipolesQ, q + 4*d ) ;
                    qXY = eScale * Array1D_Item ( multipolesQ, q + 5*d ) ;
                    qXZ = eScale * Array1D_Item ( multipolesQ, q + 6*d ) ;
                    qYY = eScale * Array1D_Item ( multipolesQ, q + 7*d ) ;
                    qYZ = eScale * Array1D_Item ( multipolesQ, q + 8*d ) ;
                    qZZ = eScale * Array1D_Item ( multipolesQ, q + 9*d ) ;
                }
            }
	    Coordinates3_GetRow ( coordinates3Q, q, xQ, yQ, zQ ) ;
       	    for ( n = 0, gXt = gYt = gZt = 0.0e+00 ; n < record->capacity ; n++ )
	    {
	        m  = record->indices[n] ;
                qM = Array1D_Item ( chargesM, m ) ;
	        Coordinates3_GetRow ( coordinates3M, m, xM, yM, zM ) ;
                dX = xQ - xM ;
                dY = yQ - yM ;
                dZ = zQ - zM ;
                r2  = ( dX * dX + dY * dY + dZ * dZ ) ;
                r1  = sqrt ( r2 ) ;
                if ( r1 <= cutOff )
                {
                    f2 = qM * ( self->alpha2 * r2 + self->beta2 ) / scale2 ; ;
                    f3 = qM * ( self->alpha3 * r2 + self->beta3 ) / scale3 ; ;
                    g1 = 2.0e+00 * qM * self->alpha1 * r1 / scale1 ;
                    g2 = 2.0e+00 * qM * self->alpha2 * r1 / scale2 ;
                    g3 = 2.0e+00 * qM * self->alpha3 * r1 / scale3 ;
                    if ( r1 == 0.0e+00 ) { dX = dY = dZ = 0.0e+00 ; r1 = 1.0e+00 ; }
                }
                else
                {
                    f2 = qM / (      r2 * scale2 ) ;
                    f3 = qM / ( r1 * r2 * scale3 ) ;
                    g1 = -           qM / (      r2 * scale1 ) ;
                    g2 = - 2.0e+00 * qM / ( r1 * r2 * scale2 ) ;
                    g3 = - 3.0e+00 * qM / ( r2 * r2 * scale3 ) ;
                }
                dX /= r1 ;
                dY /= r1 ;
                dZ /= r1 ;
                g1 *= q0 ;
                gX  = g1 * dX ;
                gY  = g1 * dY ;
                gZ  = g1 * dZ ;
                if ( order > 0 )
                {
                    a   = f2 / r1 ;
                    b   = ( g2 - a ) * ( dX * qX + dY * qY + dZ * qZ ) ;
                    gX -= ( qX * a + dX * b ) ;
                    gY -= ( qY * a + dY * b ) ;
                    gZ -= ( qZ * a + dZ * b ) ;
                    if ( order > 1 )
                    {
                        t   = qXX + qYY + qZZ ;
                        a   = 2.0e+00 * f3 / r1 ;
                        b   = 0.5e+00 * ( 3.0e+00 * ( g3 - a ) * ( qXX * dX * dX +
                                                                   qYY * dY * dY +
                                                                   qZZ * dZ * dZ +
                                                       2.0e+00 * ( qXY * dX * dY +
                                                                   qXZ * dX * dZ +
                                                                   qYZ * dY * dZ ) ) - g3 * t ) ;
                        gX += ( 1.5e+00 * a * ( qXX * dX + qXY * dY + qXZ * dZ ) + dX * b ) ;
                        gY += ( 1.5e+00 * a * ( qXY * dX + qYY * dY + qYZ * dZ ) + dY * b ) ;
                        gZ += ( 1.5e+00 * a * ( qXZ * dX + qYZ * dY + qZZ * dZ ) + dZ * b ) ;
                    }
                }
                gXt += gX ; gYt += gY ; gZt += gZ ;
                Coordinates3_DecrementRow ( gradients3M, m, gX, gY, gZ ) ;
	    }
            Coordinates3_IncrementRow ( gradients3Q, q, gXt, gYt, gZt ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . QC/MM potentials in atomic units.
! . These are not the true potentials but are formulated so that the energy can be written as C^T * P where C is the vector of
! . unique Cartesian multipoles. This means, for example, that the potentials for symmetric multipoles (e.g. qXY = qYX) are doubled.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionFull_QCMMPotentials ( const PairwiseInteractionFull *self               ,
                                              const Integer                  multipoleOrder     ,
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
         ( Status_IsOK ( status )        ) )
    {
        auto Integer     d = Coordinates3_Rows ( coordinates3Q ), m, n, order, q, r ;
        auto PairRecord *record ;
        auto Real        cutOff, dX, dY, dZ, eScale, f1, f2, f3,
                         p0, pX, pY, pZ, pXX, pXY, pXZ, pYY, pYZ, pZZ,
                         qM, r1, r2, scale1, scale2, scale3,
                         xM, xQ, yM, yQ, zM, zQ ;
        cutOff = Maximum ( 0.0e+00, self->dampingCutOff ) ;
        eScale = electrostaticScale ;
        order  = Maximum ( Minimum ( multipoleOrder, 2 ), 0 ) ;
        scale1 =       Units_Length_Angstroms_To_Bohrs      ;
        scale2 = pow ( Units_Length_Angstroms_To_Bohrs, 2 ) ;
        scale3 = pow ( Units_Length_Angstroms_To_Bohrs, 3 ) ;
        for ( r = 0 ; r < pairList->count ; r++ )
        {
            record = PairList_GetRecord ( pairList, r ) ;
            q      = record->index ;
	    Coordinates3_GetRow ( coordinates3Q, q, xQ, yQ, zQ ) ;
            p0 = pX = pY = pZ = pXX = pXY = pXZ = pYY = pYZ = pZZ = 0.0e+00 ;
       	    for ( n = 0 ; n < record->capacity ; n++ )
	    {
	        m   = record->indices[n] ;
                qM  = Array1D_Item ( chargesM, m ) ;
	        Coordinates3_GetRow ( coordinates3M, m, xM, yM, zM ) ;
                dX = xQ - xM ;
                dY = yQ - yM ;
                dZ = zQ - zM ;
                r2  = ( dX * dX + dY * dY + dZ * dZ ) ;
                r1  = sqrt ( r2 ) ;
                if ( r1 <= cutOff )
                {
                    f1 = qM * ( self->alpha1 * r2 + self->beta1 ) / scale1 ;
                    f2 = qM * ( self->alpha2 * r2 + self->beta2 ) / scale2 ;
                    f3 = qM * ( self->alpha3 * r2 + self->beta3 ) / scale3 ;
                    if ( r1 == 0.0e+00 ) { dX = dY = dZ = 0.0e+00 ; r1 = 1.0e+00 ; }
                }
                else
                {
                    f1 = qM / ( r1      * scale1 ) ;
                    f2 = qM / (      r2 * scale2 ) ;
                    f3 = qM / ( r1 * r2 * scale3 ) ;
                }
                p0 += f1 ;
                if ( order > 0 )
                {
                    dX /= r1 ;
                    dY /= r1 ;
                    dZ /= r1 ;
                    pX -= dX * f2 ;
                    pY -= dY * f2 ;
                    pZ -= dZ * f2 ;
                    if ( order > 1 )
                    {
                        pXX += 0.5e+00 * ( 3.0e+00 * dX * dX - 1.0e+00 ) * f3 ;
                        pXY +=           ( 3.0e+00 * dX * dY           ) * f3 ; /* times 2 */
                        pXZ +=           ( 3.0e+00 * dX * dZ           ) * f3 ; /* times 2 */
                        pYY += 0.5e+00 * ( 3.0e+00 * dY * dY - 1.0e+00 ) * f3 ;
                        pYZ +=           ( 3.0e+00 * dY * dZ           ) * f3 ; /* times 2 */
                        pZZ += 0.5e+00 * ( 3.0e+00 * dZ * dZ - 1.0e+00 ) * f3 ;
                    }
                }
	    }
            Array1D_Item ( potentials, q ) += ( eScale * p0 ) ;
            if ( order > 0 )
            {
                Array1D_Item ( potentials, q +   d ) += ( eScale * pX ) ;
                Array1D_Item ( potentials, q + 2*d ) += ( eScale * pY ) ;
                Array1D_Item ( potentials, q + 3*d ) += ( eScale * pZ ) ;
                if ( order > 1 )
                {
                    Array1D_Item ( potentials, q + 4*d ) += ( eScale * pXX ) ;
                    Array1D_Item ( potentials, q + 5*d ) += ( eScale * pXY ) ;
                    Array1D_Item ( potentials, q + 6*d ) += ( eScale * pXZ ) ;
                    Array1D_Item ( potentials, q + 7*d ) += ( eScale * pYY ) ;
                    Array1D_Item ( potentials, q + 8*d ) += ( eScale * pYZ ) ;
                    Array1D_Item ( potentials, q + 9*d ) += ( eScale * pZZ ) ;
                }
            }
        }
    }
}

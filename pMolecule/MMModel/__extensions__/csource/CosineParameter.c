/*==================================================================================================================================
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "CosineParameter.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Binomial  ( const Integer n, const Integer k ) ;
static Real Factorial ( const Integer n ) ;

/*==================================================================================================================================
! . Procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineParameter_Allocate ( CosineParameter *self, const Integer nTerms )
{
    if ( ( self != NULL ) && ( nTerms > 0 ) )
    {
        self->nTerms            = nTerms ;
        self->termCoefficients  = Memory_AllocateArrayOfTypes ( nTerms, Real    ) ;
        self->periods           = Memory_AllocateArrayOfTypes ( nTerms, Integer ) ;
        self->nPowers           = -1 ;
        self->powerCoefficients = NULL ;
        Memory_Set ( self->periods         , nTerms, 0       ) ;
        Memory_Set ( self->termCoefficients, nTerms, 0.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal clone.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineParameter_Clone ( CosineParameter *self, const CosineParameter *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        CosineParameter_Allocate ( self, other->nTerms ) ;
        for ( i = 0 ; i < other->nTerms ; i++ )
        {
            self->periods         [i] = other->periods         [i] ;
            self->termCoefficients[i] = other->termCoefficients[i] ;
        }
        CosineParameter_MakePowers ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineParameter_Deallocate ( CosineParameter *self )
{
    if ( self != NULL )
    {
        self->nPowers = -1 ;
        self->nTerms  =  0 ;
        Memory_Deallocate ( self->periods           ) ;
        Memory_Deallocate ( self->powerCoefficients ) ;
        Memory_Deallocate ( self->termCoefficients  ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineParameter_Initialize ( CosineParameter *self )
{
    if ( self != NULL )
    {
        self->nPowers           = -1   ;
        self->nTerms            =  0   ;
        self->periods           = NULL ;
        self->powerCoefficients = NULL ;
        self->termCoefficients  = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameter internal make powers.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real Binomial ( const Integer n, const Integer k )
{
    return Factorial ( n ) / ( Factorial ( n - k ) * Factorial ( k ) ) ;
}

static Real Factorial ( const Integer n )
{
    auto Integer i = n ;
    auto Real    f = 1.0e+00 ;
    while ( i > 0 )
    {
        f *= ( Real ) i ;
        i -= 1 ;
    }
    return f ;
}

void CosineParameter_MakePowers ( CosineParameter *self )
{
    if ( ( self != NULL ) && ( self->nPowers < 0 ) )
    {
        auto Integer i, nPowers ;
        /* . Find the largest power. */
        nPowers = -1 ;
        for ( i = 0 ; i < self->nTerms ; i++ ) nPowers = Maximum ( nPowers, self->periods[i] ) ;
        /* . Determine the power coefficients. */
        if ( nPowers >= 0 )
        {
            auto Integer l, m, n, nby2 ;
            auto Real    c, f, phase ;
            /* . Initialization. */
            self->nPowers           = nPowers ;
            self->powerCoefficients = Memory_AllocateArrayOfTypes ( nPowers + 1, Real ) ;
            Memory_Set ( self->powerCoefficients, nPowers + 1, 0.0e+00 ) ;
            /* . Loop over terms. */
            for ( i = 0 ; i < self->nTerms ; i++ )
            {
                c     = self->termCoefficients[i] ;
                n     = self->periods[i]          ;
                nby2  = n / 2 ;
                phase = 1.0e+00 ;
                for ( l = 0 ; l <= nby2 ; l++ )
                {
                    f = 0.0e+00 ;
                    for ( m = l ; m <= nby2 ; m++ ) f += Binomial ( n, 2 * m ) * Binomial ( m, l ) ;
                    self->powerCoefficients[n-2*l] += phase * f * c ;
                    /* printf ( "\nCoefficient for n = %d and power = %d = %10.3f\n", n, n-2*l, phase * f ) ; */
                    phase *= -1.0e+00 ;
                }
            }
        }
    }
}

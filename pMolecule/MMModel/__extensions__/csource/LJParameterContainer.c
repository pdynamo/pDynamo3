/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean.h"
# include "LJParameterContainer.h"
# include "Memory.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The tolerance for determining whether an LJ interaction is zero. */
# define NULL_INTERACTION 1.0e-20

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static void LJParameterContainer_MakeEpsilonSigma ( LJParameterContainer *self, const Boolean QSIGMA ) ;
static void LJParameterContainer_MakeTable        ( LJParameterContainer *self, const Boolean QGEOMETRIC, const Boolean QSIGMA ) ;

/*==============================================================================
! . Public procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
LJParameterContainer *LJParameterContainer_Allocate ( const Integer ntypes )
{
    LJParameterContainer *self = NULL ;
    if ( ntypes > 0 )
    {
        auto Integer i, ij, j, n, size ;
        self = ( LJParameterContainer * ) malloc ( sizeof ( LJParameterContainer ) ) ;
        self->ntypes     = ntypes ;
        size = ( ntypes * ( ntypes + 1 ) ) / 2 ;
        self->epsilon    = ( Real * ) calloc ( ntypes, sizeof ( Real ) ) ;
        self->sigma      = ( Real * ) calloc ( ntypes, sizeof ( Real ) ) ;
        self->tableA     = ( Real * ) calloc ( size,   sizeof ( Real ) ) ;
        self->tableB     = ( Real * ) calloc ( size,   sizeof ( Real ) ) ;
        self->tableindex = ( Integer    * ) calloc ( ( ntypes * ntypes ), sizeof ( Integer ) ) ;
        /* . Fill tableindex. */
        n = 0 ;
        for ( i = 0, n = 0 ; i < ntypes ; i++ )
        {
            for ( j = 0 ; j < ntypes ; j++, n++ )
            {
                if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
                else          ij = ( j * ( j + 1 ) ) / 2 + i ;
                self->tableindex[n] = ij ;
            }
        }
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
LJParameterContainer *LJParameterContainer_Clone ( const LJParameterContainer *self )
{
    LJParameterContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = LJParameterContainer_Allocate ( self->ntypes ) ;
        for ( i = 0 ; i < self->ntypes ; i++ ) new->epsilon[i] = self->epsilon[i] ;
        for ( i = 0 ; i < self->ntypes ; i++ ) new->sigma[i]   = self->sigma[i]   ;
        for ( i = 0 ; i < ( self->ntypes * ( self->ntypes + 1 ) ) / 2 ; i++ ) { new->tableA[i] = self->tableA[i] ; new->tableB[i] = self->tableB[i] ; }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void LJParameterContainer_Deallocate ( LJParameterContainer **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self)->epsilon    ) ;
        Memory_Deallocate ( (*self)->sigma      ) ;
        Memory_Deallocate ( (*self)->tableA     ) ;
        Memory_Deallocate ( (*self)->tableB     ) ;
        Memory_Deallocate ( (*self)->tableindex ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Different representation procedures.
!-----------------------------------------------------------------------------*/
void LJParameterContainer_MakeEpsilonSigmaAMBER ( LJParameterContainer *self ) { LJParameterContainer_MakeEpsilonSigma ( self, False ) ; }
void LJParameterContainer_MakeEpsilonSigmaOPLS  ( LJParameterContainer *self ) { LJParameterContainer_MakeEpsilonSigma ( self, True  ) ; }
void LJParameterContainer_MakeTableAMBER        ( LJParameterContainer *self ) { LJParameterContainer_MakeTable        ( self, False, False ) ; }
void LJParameterContainer_MakeTableOPLS         ( LJParameterContainer *self ) { LJParameterContainer_MakeTable        ( self, True,  True  ) ; }

/*------------------------------------------------------------------------------
! . Merging - epsilon and sigma only.
!-----------------------------------------------------------------------------*/
LJParameterContainer *LJParameterContainer_MergeEpsilonSigma ( const LJParameterContainer *self, const LJParameterContainer *other )
{
    LJParameterContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = LJParameterContainer_Allocate ( self->ntypes + other->ntypes ) ;
        for ( i = 0 ; i < self->ntypes ; i++ )
        {
            new->epsilon[i] = self->epsilon[i] ;
            new->sigma  [i] = self->sigma  [i] ;
        }
        for ( i = 0 ; i < other->ntypes ; i++ )
        {
            new->epsilon[i+self->ntypes] = other->epsilon[i] ;
            new->sigma  [i+self->ntypes] = other->sigma  [i] ;
        }
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Scaling.
!-----------------------------------------------------------------------------*/
void LJParameterContainer_Scale ( LJParameterContainer *self, const Real scale )
{
   if ( self != NULL )
   {
      auto Integer i ;
      for ( i = 0 ; i < self->ntypes ; i++ ) self->epsilon[i] *= scale ;
      for ( i = 0 ; i < ( self->ntypes * ( self->ntypes + 1 ) ) / 2 ; i++ )
      {
         self->tableA[i] *= scale ;
         self->tableB[i] *= scale ;
      }
   }
}

/*==============================================================================
! . Private procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Convert table values to epsilon/sigma.
! . Diagonal values only used.
!-----------------------------------------------------------------------------*/
static void LJParameterContainer_MakeEpsilonSigma ( LJParameterContainer *self, const Boolean QSIGMA )
{
    if ( self != NULL )
    {
        auto Real es12, es6 ;
        auto Integer i, n ;
        n = 0 ;
        for ( i = 0 ; i < self->ntypes ; i++ )
        {
            n    = self->tableindex[i+i*self->ntypes] ;
            es12 = fabs ( self->tableA[n] ) ;
            es6  = fabs ( self->tableB[n] ) ;
            if ( QSIGMA ) { es12 /= 4.0e+00 ; es6 /= 4.0e+00 ; }
            else          { es6  /= 2.0e+00 ; }
            if ( ( es6 < NULL_INTERACTION ) || ( es12 < NULL_INTERACTION ) )
            {
                self->epsilon[i] = 0.0e+00 ;
                self->sigma[i]   = 0.0e+00 ;
            }
            else
            {
                self->epsilon[i] = ( es6 * es6 ) / es12 ;
                self->sigma[i]   = exp ( log ( es12 / es6 ) / 6.0e+00 ) ;
            }
        }
    }
}

/*------------------------------------------------------------------------------
! . Convert epsilon/sigma values to table values.
!-----------------------------------------------------------------------------*/
static void LJParameterContainer_MakeTable ( LJParameterContainer *self, const Boolean QGEOMETRIC, const Boolean QSIGMA )
{
    if ( self != NULL )
    {
        auto Real eij, sij, sij12, sij6 ;
        auto Integer i, j, n ;
        n = 0 ;
        for ( i = 0 ; i < self->ntypes ; i++ )
        {
            for ( j = 0 ; j <= i ; j++ )
            {
                eij = sqrt ( self->epsilon[i] * self->epsilon[j] ) ;
                if ( QGEOMETRIC ) sij =      sqrt ( self->sigma[i] * self->sigma[j] ) ;
	        else              sij = 0.5e+00 * ( self->sigma[i] + self->sigma[j] ) ;
                sij6  = pow ( sij, 6 ) ;
	        sij12 = sij6 * sij6 ;
                if ( QSIGMA ) { eij  *= 4.0e+00 ; }
	        else          { sij6 *= 2.0e+00 ; }
	        self->tableA[n] = eij * sij12 ;
	        self->tableB[n] = eij * sij6  ;
	        self->tableindex[j+i*self->ntypes] = n ;
	        self->tableindex[i+j*self->ntypes] = n ;
                n ++ ;
            }
        }
    }
}

# ifndef _DFTGRIDWEIGHTS
# define _DFTGRIDWEIGHTS

# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The grid weights type. */
typedef struct {
          Real         *aij            ;
          Real         *rij            ;
    const Coordinates3 *qcCoordinates3 ; /* . Pointer only - object doesn't belong to type. */
} DFTGridWeights ;

/* . The grid weights derivatives work type. */
typedef struct {
    Real *A    ;
    Real *R    ;
    Real *dAdm ;
} DFTGridWeightsDerivativesWork ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DFTGridWeights *DFTGridWeights_Allocate    ( const Coordinates3                  *qcCoordinates3   ,
                                                    const Real                          *radii            ,
                                                          Status                        *status           ) ;
extern void            DFTGridWeights_Deallocate  (       DFTGridWeights               **self             ) ;
extern void            DFTGridWeights_Derivatives ( const DFTGridWeights                *self             ,
                                                    const Integer                        gridAtom         ,
                                                    const Integer                        numberOfPoints   ,
                                                    const Coordinates3                  *gridCoordinates3 ,
                                                    const RealArray1D                   *gridWeights      ,
                                                    const RealArray1D                   *eXC              ,
                                                          Coordinates3                  *gradients3       ,
                                                          DFTGridWeightsDerivativesWork *work             ) ;
extern Real            DFTGridWeights_Weight      (       DFTGridWeights                *self             ,
                                                    const Integer                        iqm              ,
                                                    const Real                          *rg               ,
                                                          Real                          *psmu             ,
                                                          Real                          *rtemp            ) ;

extern DFTGridWeightsDerivativesWork *DFTGridWeightsDerivativesWork_Allocate   ( const DFTGridWeights                 *gridWeights ,
                                                                                       Status                         *status      ) ;
extern void                           DFTGridWeightsDerivativesWork_Deallocate (       DFTGridWeightsDerivativesWork **self        ) ;

# endif

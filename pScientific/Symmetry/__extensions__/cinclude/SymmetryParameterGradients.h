# ifndef _SYMMETRYPARAMETERGRADIENTS
# define _SYMMETRYPARAMETERGRADIENTS

# include "Coordinates3.h"
# include "Matrix33.h"
# include "Real.h"
# include "Status.h"
# include "SymmetryParameters.h"
# include "Transformation3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The symmetry parameter gradients type. */
typedef struct {
    Real      dEda      ;
    Real      dEdb      ;
    Real      dEdc      ;
    Real      dEdalpha  ;
    Real      dEdbeta   ;
    Real      dEdgamma  ;
    Matrix33 *dEdH      ;
} SymmetryParameterGradients ;


/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SymmetryParameterGradients *SymmetryParameterGradients_Allocate              (       void                                            ) ;
extern SymmetryParameterGradients *SymmetryParameterGradients_AllocateWithMatrix    ( const Matrix33                    *dEdH               ,
                                                                                            Status                      *status             ) ;
extern void                        SymmetryParameterGradients_Deallocate            (       SymmetryParameterGradients **self               ) ;
extern void                        SymmetryParameterGradients_CrystalDerivatives    (       SymmetryParameterGradients  *self               ,
                                                                                      const SymmetryParameters          *symmetryParameters ) ;
extern void                        SymmetryParameterGradients_FractionalDerivatives (       SymmetryParameterGradients  *self               ,
                                                                                      const SymmetryParameters          *symmetryParameters ,
                                                                                      const Coordinates3                *coordinates3       ,
                                                                                            Coordinates3                *gradients3         ) ;
extern void                        SymmetryParameterGradients_ImageDerivatives      (       SymmetryParameterGradients  *self               ,
                                                                                      const SymmetryParameters          *symmetryParameters ,
                                                                                      const Transformation3             *transformation3    ,
                                                                                      const Coordinates3                *coordinates3       ,
                                                                                      const Coordinates3                *gradients3         ) ;
extern void                        SymmetryParameterGradients_Initialize            (       SymmetryParameterGradients  *self               ) ;

# endif

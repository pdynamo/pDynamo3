# ifndef _TRANSFORMATION3
# define _TRANSFORMATION3

# include "Boolean.h"
# include "Matrix33.h"
# include "Status.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
   Matrix33 *rotation    ;
   Vector3  *translation ;
} Transformation3 ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Transformation3 *Transformation3_Allocate                (       Status           *status      ) ;
extern Transformation3 *Transformation3_AllocateFull            (       Status           *status      ) ;
extern void             Transformation3_ApplyToVector3          ( const Transformation3  *self        ,
                                                                        Vector3          *a           ) ;
extern Transformation3 *Transformation3_Clone                   ( const Transformation3  *self        ,
                                                                        Status           *status      ) ;
extern void             Transformation3_CopyTo                  ( const Transformation3  *self        ,
                                                                        Transformation3  *other       ) ;
extern void             Transformation3_Deallocate              (       Transformation3 **self        ) ;
extern Transformation3 *Transformation3_FromRotationTranslation (       Matrix33         *rotation    ,
                                                                        Vector3          *translation ,
                                                                        Status           *status      ) ;
extern Boolean          Transformation3_IsEqual                 ( const Transformation3  *self        ,
                                                                  const Transformation3  *other       ) ;
extern Boolean          Transformation3_IsIdentity              ( const Transformation3  *self        ) ;
extern void             Transformation3_Orthogonalize           (       Transformation3  *self        ,
                                                                  const Matrix33         *A           ,
                                                                  const Matrix33         *B           ) ;

# endif

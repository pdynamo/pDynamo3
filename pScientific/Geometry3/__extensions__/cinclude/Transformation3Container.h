# ifndef _TRANSFORMATION3CONTAINER
# define _TRANSFORMATION3CONTAINER

# include "Boolean.h"
# include "Integer.h"
# include "Integer.h"
# include "Status.h"
# include "Transformation3.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean           isOwner  ;
    Integer           capacity ;
    Integer           identity ;
    Integer          *inverses ;
    Transformation3 **items    ;
} Transformation3Container ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Transformation3Container *Transformation3Container_Allocate                      ( const Integer                    capacity    ,
                                                                                                Status                    *status      ) ;
extern void                      Transformation3Container_Deallocate                    (       Transformation3Container **self        ) ;
extern void                      Transformation3Container_FindIdentity                  (       Transformation3Container  *self        ) ;
extern void                      Transformation3Container_FindInverseIntegerTranslation ( const Transformation3Container  *self        ,
                                                                                          const Integer                    t           ,
                                                                                          const Integer                    a           ,
                                                                                          const Integer                    b           ,
                                                                                          const Integer                    c           ,
                                                                                                Vector3                   *translation ,
                                                                                                Integer                   *aInverse    ,
                                                                                                Integer                   *bInverse    ,
                                                                                                Integer                   *cInverse    ) ;
extern void                      Transformation3Container_FindInverses                  (       Transformation3Container  *self        ,
                                                                                                Status                    *status      ) ;
extern Integer                   Transformation3Container_NumberOfNonSelfInverses       ( const Transformation3Container  *self        ) ;
extern Integer                   Transformation3Container_NumberOfSelfInverses          ( const Transformation3Container  *self        ) ;

# endif

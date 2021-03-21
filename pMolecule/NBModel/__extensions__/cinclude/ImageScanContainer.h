# ifndef _IMAGESCANCONTAINER
# define _IMAGESCANCONTAINER

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "Status.h"
# include "SymmetryParameters.h"
# include "Transformation3Container.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The image scan type. */
typedef struct {
    Boolean  doSkip ;
    Integer  t      ;
    Integer  a      ;
    Integer  b      ;
    Integer  c      ;
    Real     scale  ;
} ImageScan ;

/* . The image scan container type. */
typedef struct {
    Integer    capacity ;
    Integer    count    ;
    Real       cutOff   ;
    ImageScan *records  ;
} ImageScanContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern ImageScanContainer *ImageScanContainer_Allocate    (       Status                   *status             ) ;
extern void                ImageScanContainer_Append      (       ImageScanContainer       *self               ,
                                                                  Boolean                   doSkip             ,
                                                                  Integer                   t                  ,
                                                                  Integer                   a                  ,
                                                                  Integer                   b                  , 
                                                                  Integer                   c                  ,
                                                                  Real                      scale              ,
                                                                  Status                   *status             ) ;
extern ImageScanContainer *ImageScanContainer_Constructor ( const Real                      cutOff             ,
                                                            const Coordinates3             *coordinates3       ,
                                                            const SymmetryParameters       *symmetryParameters ,
                                                            const Transformation3Container *transformations    ,
                                                            const Boolean                   checkForInverses   ,
                                                            const Integer                   expandFactor       ,
                                                                  Status                   *status             ) ;
extern void                ImageScanContainer_Deallocate  (       ImageScanContainer      **self               ) ;
extern Boolean             ImageScanContainer_Reallocate  (       ImageScanContainer       *self               ,
                                                            const Integer                   capacity           ,
                                                                  Status                   *status             ) ;

# endif

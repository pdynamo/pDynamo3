# ifndef _IMAGEPAIRLISTCONTAINER
# define _IMAGEPAIRLISTCONTAINER

# include "Boolean.h"
# include "Coordinates3.h"
# include "ImageScanContainer.h"
# include "Integer.h"
# include "PairList.h"
# include "PairListGenerator.h"
# include "Real.h"
# include "RegularGrid.h"
# include "RegularGridOccupancy.h"
# include "Selection.h"
# include "Status.h"
# include "SymmetryParameterGradients.h"
# include "SymmetryParameters.h"
# include "Transformation3.h"
# include "Transformation3Container.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The image pairlist type. */
typedef struct {
    Integer          a               ;
    Integer          b               ;
    Integer          c               ;
    Real             scale           ;
    PairList        *pairList        ; /* . Owned. */
    Transformation3 *transformation3 ; /* . Alias. */
} ImagePairList ;

/* . The image pairlist container type. */
typedef struct {
    Integer         capacity      ;
    Integer         count         ;
    Integer         numberOfPairs ;
    ImagePairList **records       ;
} ImagePairListContainer ;

/* . The image pairlist iterator type. */
typedef struct {
    Boolean                     doGradients                ;
    Integer                     current                    ;
    Real                        scale                      ;
    Coordinates3               *coordinates3               ; /* . Alias. */
    Coordinates3               *gradients3                 ;
    Coordinates3               *iCoordinates3              ; /* . Owned. */
    Coordinates3               *iGradients3                ;
    ImagePairListContainer     *target                     ; /* . Alias. */
    PairList                   *pairList                   ;
    SymmetryParameterGradients *symmetryParameterGradients ;
    SymmetryParameters         *symmetryParameters         ;
    Transformation3            *iTransformation3           ; /* . Owned. */
    Transformation3            *xTransformation3           ;
} ImagePairListIterator ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Image pairlist. */
extern ImagePairList          *ImagePairList_Allocate                (       Status                      *status                     ) ;
extern void                    ImagePairList_Deallocate              (       ImagePairList              **self                       ) ;
extern ImagePairList          *ImagePairList_FromItems               ( const Integer                      a                          ,
                                                                       const Integer                      b                          ,
                                                                       const Integer                      c                          ,
                                                                       const Real                         scale                      ,
                                                                             PairList                    *pairList                   ,
                                                                             Transformation3             *transformation3            ,
                                                                             Status                      *status                     ) ;

/* . Image pairlist container. */
extern ImagePairListContainer *ImagePairListContainer_Allocate       (       Status                      *status                     ) ;
extern void                    ImagePairListContainer_Append         (       ImagePairListContainer      *self                       ,
                                                                             ImagePairList               *record                     ,
                                                                             Status                      *status                     ) ;
extern ImagePairListContainer *ImagePairListContainer_Constructor    ( const PairListGenerator           *generator                  ,
                                                                             Selection                   *atomsA                     ,
                                                                             Selection                   *atomsB                     ,
                                                                             Selection                   *freeAtoms                  ,
                                                                       const Coordinates3                *coordinates3A              ,
                                                                       const Coordinates3                *coordinates3B              ,
                                                                       const SymmetryParameters          *symmetryParameters         ,
                                                                       const Transformation3Container    *transformations            ,
                                                                       const ImageScanContainer          *scanData                   ,
                                                                             RegularGrid                 *gridA                      ,
                                                                             RegularGridOccupancy        *occupancyA                 ,
                                                                       const Boolean                      checkForInverses           ,
                                                                             Status                      *status                     ) ;
extern void                    ImagePairListContainer_Deallocate     (       ImagePairListContainer     **self                       ) ;
extern Integer                 ImagePairListContainer_NumberOfImages ( const ImagePairListContainer      *self                       ) ;
extern Integer                 ImagePairListContainer_NumberOfPairs  ( const ImagePairListContainer      *self                       ) ;
extern Boolean                 ImagePairListContainer_Reallocate     (       ImagePairListContainer      *self                       ,
                                                                       const Integer                      capacity                   ,
                                                                             Status                      *status                     ) ;

/* . Image pairlist container iterator. */
extern void                    ImagePairListIterator_Finalize        (       ImagePairListIterator      *self                        ) ;
extern void                    ImagePairListIterator_Gradients       (       ImagePairListIterator      *self                        ) ;
extern void                    ImagePairListIterator_Initialize      (       ImagePairListIterator      *self                        ,
                                                                             ImagePairListContainer     *target                      ,
                                                                             Coordinates3               *coordinates3                ,
                                                                             SymmetryParameters         *symmetryParameters          ,
                                                                             Coordinates3               *gradients3                  ,
                                                                             SymmetryParameterGradients *symmetryParameterGradients  ,
                                                                             Status                     *status                      ) ;
extern Boolean                 ImagePairListIterator_Next            (       ImagePairListIterator      *self                        ) ;

# endif

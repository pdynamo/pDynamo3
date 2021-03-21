# ifndef _PAIRLIST
# define _PAIRLIST

# include "Boolean.h"
# include "Integer.h"
# include "Selection.h"
# include "SelectionContainer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The pair connections type. */
typedef struct {
    Integer  capacityI ;
    Integer  capacityJ ;
    Integer *itemsI ;
    Integer *itemsJ ;
} PairConnections ;

/* . The pair record type. */
typedef struct {
   Integer  index    ;
   Integer  capacity ;
   Integer *indices  ;
} PairRecord ;

/* . The pair excluded type. */
typedef struct {
    Integer     capacity ;
    Integer    *indices  ;
    Integer    *work     ;
    PairRecord  record   ;
} PairExcluded ;

/* . The pair list type. */
typedef struct {
    Boolean           isSelf        ;
    Boolean           isSorted      ;
    Integer           capacity      ;
    Integer           count         ;
    Integer           numberOfPairs ;
    PairConnections  *connections   ;
    PairExcluded     *excluded      ;
    PairRecord      **records       ;
} PairList ;

/* . The pair list iterator type. */
typedef struct {
    Integer   current ;
    PairList *target  ;
} PairListIterator ;


/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Pair connection functions. */
extern PairConnections    *PairConnections_Allocate             ( const Integer           capacityI     ,
                                                                  const Integer           capacityJ     ,
                                                                        Status           *status        ) ;
extern void                PairConnections_Deallocate           (       PairConnections **self          ) ;

extern PairExcluded       *PairExcluded_Allocate                ( const Integer           capacity      ,
                                                                        Status           *status        ) ;
extern void                PairExcluded_Deallocate              (       PairExcluded    **self          ) ;
extern PairExcluded       *PairExcluded_FromIndices             ( const Integer           capacity      ,
                                                                  const Integer          *indices       ,
                                                                        Status           *status        ) ;

/* . General pair-list functions. */
extern PairList           *PairList_Allocate                    ( const Integer           capacity      ,
                                                                        Status           *status        ) ;
extern void                PairList_Append                      (       PairList         *self          ,
                                                                        PairRecord       *record        ,
                                                                        Status           *status        ) ;
extern void                PairList_ClearRepresentations        (       PairList         *self          ) ;
extern void                PairList_Deallocate                  (       PairList        **self          ) ;
extern PairRecord         *PairList_GetRecord                   (       PairList         *self          ,
                                                                  const Integer           index         ) ;
extern void                PairList_Initialize                  (       PairList         *self          ) ;
extern Integer             PairList_MaximumRecordSize           ( const PairList         *self          ) ;
extern Integer             PairList_NumberOfPairs               ( const PairList         *self          ) ;
extern Integer             PairList_NumberOfRecords             ( const PairList         *self          ) ;
extern Boolean             PairList_Reallocate                  (       PairList         *self          ,
                                                                  const Integer           capacity      ,
                                                                        Status           *status        ) ;
extern void                PairList_Sort                        (       PairList         *self          ) ;
extern Integer             PairList_UpperBound                  ( const PairList         *self          ,
                                                                  const Boolean           isSelf        ) ;

/* . Cross pair-list functions. */
extern PairConnections    *CrossPairList_MakeConnections        (       PairList         *self          ,
                                                                  const Integer           upperBound    ,
                                                                        Status           *status        ) ;
extern PairList           *CrossPairList_MakeFull               ( const Integer           capacity1     ,
                                                                        Selection        *andSelection1 ,
                                                                  const Integer           capacity2     ,
                                                                        Selection        *andSelection2 ,
                                                                        Status           *status        ) ;
extern PairList           *CrossPairList_MakeFullExcluded       ( const Integer           capacity1     ,
                                                                        Selection        *andSelection1 ,
                                                                  const Integer           capacity2     ,
                                                                        Selection        *andSelection2 ,
                                                                        Status           *status        ) ;

/* . Self pair-list functions. */
extern SelectionContainer *SelfPairList_GetConnectedComponents  (       PairList         *self          ,
                                                                        Integer           upperBound    ,
                                                                        Status           *status        ) ;
extern PairConnections    *SelfPairList_MakeConnections         (       PairList         *self          ,
                                                                  const Integer           upperBound    ,
                                                                        Status           *status        ) ;
extern void                SelfPairList_Renumber                (       PairList         *self          ,
                                                                        Selection        *mapping       ,
                                                                        Status           *status        ) ;
extern PairList           *SelfPairList_ToCrossPairList         (       PairList         *self          ,
                                                                        Selection        *andSelection1 ,
                                                                        Selection        *andSelection2 ,
                                                                        Selection        *orSelection   ,
                                                                        Status           *status        ) ;
extern PairList           *SelfPairList_ToCrossPairListExcluded (       PairList         *self          ,
                                                                        Selection        *andSelection1 ,
                                                                        Selection        *andSelection2 ,
                                                                        Selection        *orSelection   ,
                                                                        Status           *status        ) ;
extern PairList           *SelfPairList_ToSelfPairList          (       PairList         *self          ,
                                                                        Selection        *andSelection  ,
                                                                        Selection        *orSelection   ,
                                                                        Status           *status        ) ;
extern PairList           *SelfPairList_ToSelfPairListExcluded  (       PairList         *self          ,
                                                                        Integer           capacity      ,
                                                                        Selection        *andSelection  ,
                                                                        Selection        *orSelection   ,
                                                                        Status           *status        ) ;

/* . Pair-list iterator functions. */
extern void                PairListIterator_Initialize          (       PairListIterator *self          ,
                                                                        PairList         *target        ) ;
extern PairRecord         *PairListIterator_Next                (       PairListIterator *self          ) ;

/* . Pair record functions. */
extern PairRecord         *PairRecord_Allocate                  ( const Integer           capacity      ,
                                                                        Status           *status        ) ;
extern void                PairRecord_Deallocate                (       PairRecord      **self          ) ;
extern PairRecord         *PairRecord_FromIndices               ( const Integer           index         ,
                                                                  const Integer           capacity      ,
                                                                  const Integer          *indices       ,
                                                                        Status           *status        ) ;
extern void                PairRecord_Initialize                (       PairRecord       *self          ) ;
extern void                PairRecord_Sort                      (       PairRecord       *self          ) ;

# endif

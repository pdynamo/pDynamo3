# ifndef _SLICE
# define _SLICE

# include "Boolean.h"
# include "Integer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The slice type. */
typedef struct {
    Boolean isScalar ;
    Integer extent   ;
    Integer start    ;
    Integer stop     ;
    Integer stride   ;
} Slice ;

/* . The multislice type. */
typedef struct {
    Integer  capacity ;
    Integer  rank     ;
    Slice   *items    ;
} MultiSlice ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Multislice. */
extern MultiSlice *MultiSlice_Allocate        ( const Integer      capacity  ,
                                                      Status      *status    ) ;
extern void        MultiSlice_Deallocate      (       MultiSlice **self      ) ;
extern Integer     MultiSlice_GetCapacity     ( const MultiSlice  *self      ) ;
extern Integer     MultiSlice_GetExtent       ( const MultiSlice  *self      ,
                                                const Integer      dimension ,
                                                      Status      *status    ) ;
extern Integer     MultiSlice_GetRank         ( const MultiSlice  *self      ) ;
extern Integer     MultiSlice_GetSize         ( const MultiSlice  *self      ) ;
extern void        MultiSlice_SetRank         (       MultiSlice  *self      ) ;
extern void        MultiSlice_SetSliceIndices (       MultiSlice  *self      ,
                                                const Integer      dimension ,
                                                const Boolean      isScalar  ,
                                                const Integer      start     ,
                                                const Integer      stop      ,
                                                const Integer      stride    ,
                                                const Integer      n         ) ;

/* . Utility functions. */
extern Integer     SliceIndices_CheckScalar   ( const Integer      index     ,
                                                const Integer      extent    ,
                                                      Status      *status    ) ;
extern Integer     SliceIndices_CheckSlice    ( const Integer     *pStart    ,
                                                const Integer     *pStop     ,
                                                const Integer     *pStride   ,
                                                const Integer      extent    ,
                                                      Integer     *qStart    ,
                                                      Integer     *qStop     ,
                                                      Integer     *qStride   ,
                                                      Status      *status    ) ;

# endif

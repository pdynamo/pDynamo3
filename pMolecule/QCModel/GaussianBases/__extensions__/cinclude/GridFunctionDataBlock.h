# ifndef _GRIDFUNCTIONDATABLOCK
# define _GRIDFUNCTIONDATABLOCK

# include "Boolean.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The grid function data block type. */
typedef struct {
    Integer numberOfFunctions ;
    Integer numberOfPoints    ;
    Integer order             ;
    IntegerArray1D *indices   ;
    RealArray2D    *f         ;
    RealArray2D    *fX        ;
    RealArray2D    *fY        ;
    RealArray2D    *fZ        ;
    RealArray2D    *fXX       ;
    RealArray2D    *fXY       ;
    RealArray2D    *fXZ       ;
    RealArray2D    *fYY       ;
    RealArray2D    *fYZ       ;
    RealArray2D    *fZZ       ;
    RealArray2D    *fXXX      ;
    RealArray2D    *fXXY      ;
    RealArray2D    *fXXZ      ;
    RealArray2D    *fXYY      ;
    RealArray2D    *fXYZ      ;
    RealArray2D    *fXZZ      ;
    RealArray2D    *fYYY      ;
    RealArray2D    *fYYZ      ;
    RealArray2D    *fYZZ      ;
    RealArray2D    *fZZZ      ;
} GridFunctionDataBlock ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern GridFunctionDataBlock *GridFunctionDataBlock_Allocate     ( const Integer                 numberOfFunctions ,
                                                                   const Integer                 numberOfPoints    ,
                                                                   const Integer                 order             ,
                                                                         Status                 *status            ) ;
extern Real                   GridFunctionDataBlock_ByteSize     ( const GridFunctionDataBlock  *self              ) ;
extern void                   GridFunctionDataBlock_Deallocate   (       GridFunctionDataBlock **self              ) ;
extern void                   GridFunctionDataBlock_FilterValues (       GridFunctionDataBlock  *self              ,
                                                                   const Integer                 fStart            ,
                                                                   const Real                  *tolerance          ) ;
extern void                   GridFunctionDataBlock_Initialize   (       GridFunctionDataBlock  *self              ) ;
extern void                   GridFunctionDataBlock_Resize       (       GridFunctionDataBlock  *self              ,
                                                                   const Integer                 numberOfFunctions ,
                                                                         Status                 *status            ) ;
# endif

/* . Need _ArrayDataType. */

/* . Need to include relevant data type and array 1D and 2D types. */

# include "Boolean.h"
# include "Integer.h"
# include "Slice.h"
# include "Status.h"
# include "TemplateMacros.h"
# include "ViewND.h"

# define _Array1DType Array1D
# define _Array2DType Array2D
# define _ArrayType   ArrayND
# define _Array1DName _MakeToken ( _ArrayDataType, _Array1DType )
# define _Array2DName _MakeToken ( _ArrayDataType, _Array2DType )
# define _ArrayName   _MakeToken ( _ArrayDataType, _ArrayType   )
# define _BlockType   _MakeToken ( _ArrayDataType,  Block       )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Include iterator here (cached)? */

/* . The array type. */
typedef struct {
    ViewND         *view  ;
    _BlockType     *block ;
    _ArrayDataType *data  ;
} _ArrayName ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
! . These are both general (any array dimension) and specific.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "ArrayND_Macros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern _ArrayName     * _MakeToken ( _ArrayName, _Allocate              ) (        Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _AllocateWithRank      ) ( const  Integer        rank          ,
                                                                                   Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _AllocateWithShape     ) ( const  Integer        rank          ,
                                                                            const  Integer       *extents       ,
                                                                                   Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _CloneDeep             ) ( const _ArrayName     *self          ,
                                                                                   Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _CloneShallow          ) ( const _ArrayName     *self          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _CopyTo                ) ( const _ArrayName     *self          ,
                                                                                  _ArrayName     *other         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Deallocate            ) (       _ArrayName    **self          ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _FromViewBlock         ) (        ViewND        *view          ,
                                                                                  _BlockType     *block         ,
                                                                            const  Boolean        withReference ,
                                                                                   Status        *status        ) ;
extern _ArrayDataType   _MakeToken ( _ArrayName, _GetItem               ) ( const _ArrayName     *self          ,
                                                                            const  Integer       *indices       ,
                                                                                   Status        *status        ) ;
extern _ArrayDataType   _MakeToken ( _ArrayName, _GetItemMultiSlice     ) ( const _ArrayName     *self          ,
                                                                            const  MultiSlice    *multiSlice    ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Resize                ) (       _ArrayName     *self          ,
                                                                            const  Integer        extent0       ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ResizeWithInitializer ) (       _ArrayName     *self          ,
                                                                            const  Integer        extent0       ,
                                                                            const _ArrayDataType  initializer   ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Set                   ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _SetItem               ) (       _ArrayName     *self          ,
                                                                            const  Integer       *indices       ,
                                                                            const _ArrayDataType  value         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _SetItemMultiSlice     ) (       _ArrayName     *self          ,
                                                                            const  MultiSlice    *multiSlice    ,
                                                                            const _ArrayDataType  value         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ViewTail2D            ) ( const _ArrayName     *self          ,
                                                                                   Integer       *indices       ,
                                                                            const  Boolean        withReference ,
                                                                                  _Array2DName   *tail          ,
                                                                                   Status        *status        ) ;

# undef _Array1DName
# undef _Array1DType
# undef _Array2DName
# undef _Array2DType
# undef _ArrayName
# undef _ArrayType
# undef _BlockType

/* . Need _ArrayDataType. */

/* . Need to include relevant data type and block headers. */

# include "Integer.h"
# include "Selection.h"
# include "Slice.h"
# include "Status.h"
# include "TemplateMacros.h"
# include "View1D.h"
# include "View2D.h"

# define _Array1DType Array1D
# define _ArrayType   Array2D
# define _Array1DName _MakeToken ( _ArrayDataType, _Array1DType )
# define _ArrayName   _MakeToken ( _ArrayDataType, _ArrayType   )
# define _BlockType   _MakeToken ( _ArrayDataType,  Block       )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The block type. */
typedef struct {
    View2D_Fields ;
    _BlockType     *block ;
    _ArrayDataType *data  ;
} _ArrayName ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "Array2D_Macros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
extern _ArrayDataType   _MakeToken ( _ArrayName, _AbsoluteMaximum       ) ( const _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _Add                   ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ,
                                                                            const _ArrayName     *other         ,
                                                                                   Status        *status        ) ;
# endif
extern _ArrayName     * _MakeToken ( _ArrayName, _Allocate              ) (        Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _AllocateWithExtents   ) ( const  Integer        rows          ,
                                                                            const  Integer        columns       ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _AssignBlock           ) (       _ArrayName     *self          ,
                                                                                  _BlockType     *block         ,
                                                                            const  Boolean        withReference ,
                                                                                   Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _CloneDeep             ) ( const _ArrayName     *self          ,
                                                                                 Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _CloneShallow          ) ( const _ArrayName     *self          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ColumnView            ) ( const _ArrayName     *self          ,
                                                                            const  Integer        i             ,
                                                                            const  Boolean        withReference ,
                                                                                  _Array1DName   *view          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _CopyTo                ) ( const _ArrayName     *self          ,
                                                                                  _ArrayName     *other         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Deallocate            ) (       _ArrayName    **self          ) ;
extern _ArrayDataType   _MakeToken ( _ArrayName, _GetItem               ) ( const _ArrayName     *self          ,
                                                                            const  Integer        i             ,
                                                                            const  Integer        j             ,
                                                                                   Status        *status        ) ;
extern _ArrayDataType   _MakeToken ( _ArrayName, _GetItemMultiSlice     ) ( const _ArrayName     *self          ,
                                                                            const  MultiSlice    *multiSlice    ,
                                                                                   Status        *status        ) ;
extern _ArrayDataType * _MakeToken ( _ArrayName, _PointerToData         ) ( const _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _Print                 ) ( const _ArrayName     *self          ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _PruneByRow            ) ( const _ArrayName     *self          ,
                                                                            const  Selection     *toKeep        ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Resize                ) (       _ArrayName     *self          ,
                                                                            const  Integer        extent0       ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ResizeWithInitializer ) (       _ArrayName     *self          ,
                                                                            const  Integer        extent0       ,
                                                                            const _ArrayDataType  initializer   ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _RowView               ) ( const _ArrayName     *self          ,
                                                                            const  Integer        i             ,
                                                                            const  Boolean        withReference ,
                                                                                  _Array1DName   *view          ,
                                                                                   Status        *status        ) ;
# ifndef _NoNumeric
extern void             _MakeToken ( _ArrayName, _Scale                 ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ) ;
# endif
extern void             _MakeToken ( _ArrayName, _Set                   ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ) ;
extern void             _MakeToken ( _ArrayName, _SetItem               ) (       _ArrayName     *self          ,
                                                                            const  Integer        i             ,
                                                                            const  Integer        j             ,
                                                                            const _ArrayDataType  value         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _SetItemMultiSlice     ) (       _ArrayName     *self          ,
                                                                            const  MultiSlice    *multiSlice    ,
                                                                            const _ArrayDataType  value         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _View                  ) ( const _ArrayName     *self          ,
                                                                            const  Integer        start0        ,
                                                                            const  Integer        start1        ,
                                                                            const  Integer        extent0       ,
                                                                            const  Integer        extent1       ,
                                                                            const  Integer        stride0       ,
                                                                            const  Integer        stride1       ,
                                                                            const  Boolean        withReference ,
                                                                                  _ArrayName     *view          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _View1D                ) ( const _ArrayName     *self          ,
                                                                            const  Integer        dimension     ,
                                                                            const  Integer        start0        ,
                                                                            const  Integer        start1        ,
                                                                            const  Integer        extent        ,
                                                                            const  Integer        stride        ,
                                                                            const  Boolean        withReference ,
                                                                                  _Array1DName   *view          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ViewOfRaw             ) (       _ArrayName     *self          ,
                                                                            const  Integer        offset        ,
                                                                            const  Integer        extent0       ,
                                                                            const  Integer        extent1       ,
                                                                            const  Integer        stride0       ,
                                                                            const  Integer        stride1       ,
                                                                                  _ArrayDataType *data          ) ;
# undef _Array1DName
# undef _Array1DType
# undef _ArrayName
# undef _ArrayType
# undef _BlockType

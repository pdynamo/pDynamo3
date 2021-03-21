/* . Need _ArrayDataType and to include relevant data type and block headers. */

# include "Integer.h"
# include "IntegerBlock.h"
# include "Status.h"
# include "TemplateMacros.h"
# include "View1D.h"

# define _ArrayType Array1D
# define _BlockType _MakeToken ( _ArrayDataType,  Block     )
# define _ArrayName _MakeToken ( _ArrayDataType, _ArrayType )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The block type. */
typedef struct {
    View1D_Fields ;
    _BlockType     *block  ;
    _ArrayDataType *data   ;
} _ArrayName ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "Array1D_Macros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
extern _ArrayDataType   _MakeToken ( _ArrayName, _AbsoluteMaximum       ) ( const _ArrayName     *self          ) ;
extern Integer          _MakeToken ( _ArrayName, _AbsoluteMaximumIndex  ) ( const _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _Add                   ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ,
                                                                            const _ArrayName     *other         ,
                                                                                   Status        *status        ) ;
# endif
extern _ArrayName     * _MakeToken ( _ArrayName, _Allocate              ) (        Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _AllocateWithExtent    ) ( const  Integer        extent        ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _AssignBlock           ) (       _ArrayName     *self          ,
                                                                                  _BlockType     *block         ,
                                                                            const  Boolean        withReference ,
                                                                                   Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _CloneDeep             ) ( const _ArrayName     *self          ,
                                                                                   Status        *status        ) ;
extern _ArrayName     * _MakeToken ( _ArrayName, _CloneShallow          ) ( const _ArrayName     *self          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _CopyTo                ) ( const _ArrayName     *self          ,
                                                                                  _ArrayName     *other         ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Deallocate            ) (       _ArrayName    **self          ) ;
# ifndef _NoNumeric
extern _ArrayDataType   _MakeToken ( _ArrayName, _Dot                   ) ( const _ArrayName     *self          ,   
                                                                            const _ArrayName     *other         ,   
                                                                                   Status        *status        ) ; 
# endif
extern _ArrayDataType   _MakeToken ( _ArrayName, _GetItem               ) ( const _ArrayName     *self          ,   
                                                                            const  Integer        i             ,   
                                                                                   Status        *status        ) ; 
# ifndef _NoNumeric
extern void             _MakeToken ( _ArrayName, _Increment             ) (       _ArrayName     *self          ,   
                                                                            const _ArrayDataType  value         ) ; 
# endif
extern void             _MakeToken ( _ArrayName, _LeftCircularShift     ) (       _ArrayName     *self          ) ; 
# ifndef _NoNumeric
extern _ArrayDataType   _MakeToken ( _ArrayName, _Maximum               ) ( const _ArrayName     *self          ) ;
extern _ArrayDataType   _MakeToken ( _ArrayName, _Minimum               ) ( const _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _Multiply              ) (       _ArrayName     *self          ,
                                                                            const _ArrayName     *other         ,
                                                                                   Status        *status        ) ;
# endif
extern _ArrayDataType * _MakeToken ( _ArrayName, _PointerToData         ) ( const _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _Print                 ) ( const _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _Resize                ) (       _ArrayName     *self          ,
                                                                            const  Integer        extent        ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ResizeWithInitializer ) (       _ArrayName     *self          ,
                                                                            const  Integer        extent        ,
                                                                            const _ArrayDataType  initializer   ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _Reverse               ) (       _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _RightCircularShift    ) (       _ArrayName     *self          ) ;
# ifndef _NoNumeric
extern void             _MakeToken ( _ArrayName, _Scale                 ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ) ;
# endif
extern void             _MakeToken ( _ArrayName, _Set                   ) (       _ArrayName     *self          ,
                                                                            const _ArrayDataType  value         ) ;
extern void             _MakeToken ( _ArrayName, _SetItem               ) (       _ArrayName     *self          ,
                                                                            const  Integer        i             ,
                                                                            const _ArrayDataType  value         ,
                                                                                   Status        *status        ) ;
# ifndef _NoNumeric
extern void             _MakeToken ( _ArrayName, _Sort                  ) (       _ArrayName     *self          ) ;
extern void             _MakeToken ( _ArrayName, _SortIndex             ) ( const _ArrayName     *self          ,
                                                                                   IntegerBlock  *indices       ,
                                                                                   Status        *status        ) ;
extern Integer          _MakeToken ( _ArrayName, _SortUnique            ) (       _ArrayName     *self          ) ;
extern _ArrayDataType   _MakeToken ( _ArrayName, _Sum                   ) ( const _ArrayName     *self          ) ;
# endif
extern void             _MakeToken ( _ArrayName, _View                  ) ( const _ArrayName     *self          ,
                                                                            const  Integer        start         ,
                                                                            const  Integer        extent        ,
                                                                            const  Integer        stride        ,
                                                                            const  Boolean        withReference ,
                                                                                  _ArrayName     *view          ,
                                                                                   Status        *status        ) ;
extern void             _MakeToken ( _ArrayName, _ViewOfRaw             ) (       _ArrayName     *self          ,
                                                                            const  Integer        offset        ,
                                                                            const  Integer        extent        ,
                                                                            const  Integer        stride        ,
                                                                                  _ArrayDataType *data          ) ;

# undef _ArrayName
# undef _ArrayType
# undef _BlockType

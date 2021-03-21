/* . Need _CoreDataType. */

# include "Status.h"
# include "TemplateMacros.h"

# define _BlockType Block
# define _BlockName _MakeToken ( _CoreDataType, _BlockType )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The block type. */
typedef struct {
    Integer        capacity   ;
    Integer        references ;
    _CoreDataType *items      ; 
} _BlockName ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "TypedMemoryBlock_Macros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern _BlockName *  _MakeToken ( _BlockName, _Allocate    ) ( const  Integer       capacity ,
                                                                      Status       *status   ) ;
extern _BlockName *  _MakeToken ( _BlockName, _Clone       ) ( const _BlockName    *self     ,
                                                                      Status       *status   ) ;
extern void          _MakeToken ( _BlockName, _CopyTo      ) ( const _BlockName    *self     ,
                                                                     _BlockName    *other    ,
                                                                      Status       *status   ) ;
extern void          _MakeToken ( _BlockName, _Deallocate  ) (       _BlockName   **self     ) ;
extern void          _MakeToken ( _BlockName, _Dereference ) (       _BlockName    *self     ) ;
extern _CoreDataType _MakeToken ( _BlockName, _GetItem     ) ( const _BlockName    *self     ,
                                                               const  Integer       index    ,
                                                                      Status       *status   ) ;
# ifndef _ExcludeNegate
extern void          _MakeToken ( _BlockName, _Negate      ) (       _BlockName    *self     ) ;
# endif
extern void          _MakeToken ( _BlockName, _Reference   ) (       _BlockName    *self     ) ;
extern void          _MakeToken ( _BlockName, _Resize      ) (       _BlockName    *self     ,
                                                               const  Integer       capacity ,
                                                                      Status       *status   ) ;
extern void          _MakeToken ( _BlockName, _Set         ) (       _BlockName    *self     ,
                                                               const _CoreDataType  value    ) ;
extern void          _MakeToken ( _BlockName, _SetItem     ) (       _BlockName    *self     ,
                                                               const  Integer       index    ,
                                                               const _CoreDataType  value    ,
                                                                      Status       *status   ) ;

# undef _BlockName
# undef _BlockType

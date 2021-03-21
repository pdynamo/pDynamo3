/* . Need _CoreDataType. */

# include "Integer.h"
# include "Status.h"
# include "TemplateMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern _CoreDataType * _MakeToken ( _CoreDataType, _Allocate   ) ( const  Integer        capacity ,
                                                                          Status        *status   ) ;
extern _CoreDataType * _MakeToken ( _CoreDataType, _Clone      ) ( const _CoreDataType  *self     ,
                                                                   const  Integer        capacity ,
                                                                          Status        *status   ) ;
extern void            _MakeToken ( _CoreDataType, _CopyTo     ) ( const _CoreDataType  *self     ,
                                                                   const  Integer        capacity ,
                                                                         _CoreDataType  *other    ,
                                                                          Status        *status   ) ;
extern Integer         _MakeToken ( _CoreDataType, _Count      ) ( const _CoreDataType  *self     ,
                                                                   const  Integer        capacity ,
                                                                         _CoreDataType   value    ) ;
extern void            _MakeToken ( _CoreDataType, _Deallocate ) (       _CoreDataType **self     ) ;
extern _CoreDataType   _MakeToken ( _CoreDataType, _GetItem    ) ( const _CoreDataType  *self     ,
                                                                   const  Integer        capacity ,
                                                                   const  Integer        index    ,
                                                                          Status        *status   ) ;
extern void            _MakeToken ( _CoreDataType, _Set        ) (       _CoreDataType  *self     ,
                                                                   const  Integer        capacity ,
                                                                   const _CoreDataType   value    ) ;
extern void            _MakeToken ( _CoreDataType, _SetItem    ) (       _CoreDataType  *self     ,
                                                                   const  Integer        capacity ,
                                                                   const  Integer        index    ,
                                                                   const _CoreDataType   value    ,
                                                                          Status        *status   ) ;
# ifndef _ExcludeSort
extern void            _MakeToken ( _CoreDataType, _Sort       ) (       _CoreDataType  *self     ,
                                                                   const  Integer        capacity ) ;
extern Integer         _MakeToken ( _CoreDataType, _SortUnique ) (       _CoreDataType  *self     ,
                                                                   const  Integer        capacity ) ;
# endif

/*==================================================================================================================================
! . This module provides procedures for handling a list of objects.
!
! . Notes:
!
!   The rules for calculating indices are as follows:
!
!   index >= 0         the index is taken as is.
!   index <  0         the index is index + nelements.
!
!   For inserting, insertion with index 0 replaces the head of the list whereas
!   insertion with index -1 inserts as the penultimate element in the list. To
!   insert as the last element use Append or use index = nelements.
!
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "List.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
List *List_Allocate ( void )
{
    List *self ;
    self = ( List * ) malloc ( sizeof ( List ) ) ;
    List_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocate a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Deallocate ( List **list )
{
    if ( (*list) != NULL )
    {
        List_Empty        ( (*list) ) ;
        Memory_Deallocate ( (*list) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Append a non-NULL element to the end of the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Element_Append ( List *self, void *node )
{
    if ( ( self != NULL ) && ( node != NULL ) )
    {
        auto ListElement *new ;

        /* . Create the node. */
        new = ( ListElement * ) malloc ( sizeof ( ListElement ) ) ;
        new->next = NULL ;
        new->node = node ;

        /* . This is the only object on the list. */
        if ( self->nelements == 0 )
        {
            self->first = new ;
            self->last  = new ;
        }
        /* . There are objects on the list. */
        else
        {
            self->last->next = new ;
            self->last       = new ;
        }

        /* . Increment the number of objects. */
        self->nelements ++ ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Append an element by index.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Element_Append_By_Index ( List *list, void *node, const Integer index )
{
   Integer      i, target ;
   ListElement *current = NULL, *new = NULL, *next = NULL, *previous = NULL ;

   /* . Get the target index. */
   target = index ;
   if ( index < 0 ) target += list->nelements ;

   /* . Create the node. */
   new = ( ListElement * ) malloc ( sizeof ( ListElement ) ) ;
   new->next = NULL ;
   new->node = node ;

   /* . This is the only object on the list. */
   if ( list->nelements == 0 )
   {
      list->first = new ;
      list->last  = new ;
   }
   /* . There are objects on the list. */
   else
   {
      /* . The first element on the list. */
      if ( target <= 0 )
      {
         new->next   = list->first ;
         list->first = new ;
      }
      /* . The last element on the list. */
      else if ( target >= list->nelements )
      {
         list->last->next = new ;
         list->last       = new ;
      }
      /* . Somewhere in the middle of the list - insert before. */
      else
      {
         for ( i = 0 ; i < list->nelements ; i++ )
         {
            if ( i == 0 )
            {
               current = list->first       ;
               next    = list->first->next ;
            }
            else
            {
               previous = current    ;
               current  = next       ;
               next     = next->next ;
            }
            if ( i == target )
            {
               new->next      = current ;
               previous->next = new     ;
               break ;
            }
         }
      }
   }

   /* . Increment the number of objects. */
   list->nelements ++ ;

}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find an element on the list and return its node address if it exists.
!---------------------------------------------------------------------------------------------------------------------------------*/
void *List_Element_Find_By_Index ( List *list, const Integer index )
{
   Integer      i, target ;
   ListElement *current = NULL, *next = NULL ;
   void        *node    = NULL ;
   target = index ;
   if ( index < 0 ) target += list->nelements ;
   for ( i = 0 ; i < list->nelements ; i++ )
   {
      if ( i == 0 )
      {
         current = list->first       ;
         next    = list->first->next ;
      }
      else
      {
         current = next       ;
         next    = next->next ;
      }
      if ( i == target )
      {
	 node = current->node ;
         break ;
      }
   }
   return node ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find an element on the list and return its node address if it exists.
!---------------------------------------------------------------------------------------------------------------------------------*/
void *List_Element_Find_By_Match ( List *list, void *reference )
{
   Integer      i ;
   ListElement *current = NULL, *next = NULL ;
   void        *node    = NULL ;
   for ( i = 0 ; i < list->nelements ; i++ )
   {
      if ( i == 0 )
      {
         current = list->first       ;
         next    = list->first->next ;
      }
      else
      {
         current  = next       ;
         next     = next->next ;
      }
      if ( (*list->Element_Match) ( reference, current->node ) )
      {
	 node = current->node ;
         break ;
      }
   }
   return node ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find an element on the list by using its list index, remove the element and
! . return its node address. The node is not deallocated.
!---------------------------------------------------------------------------------------------------------------------------------*/
void *List_Element_Pop_By_Index ( List *list, const Integer index )
{
   void *node = NULL ;
   if ( list != NULL )
   {
      if ( ( list->nelements > 0 ) && ( index <= list->nelements ) )
      {
         auto Integer      i ;
         auto ListElement *current = NULL, *next = NULL, *previous = NULL ;
         for ( i = 0 ; i < list->nelements ; i++ )
         {
            if ( i == 0 )
            {
               current = list->first       ;
               next    = list->first->next ;
            }
            else
            {
               previous = current    ;
               current  = next       ;
               next     = next->next ;
            }
            if ( i == index )
            {
               /* . The lists are now empty. */
	       if ( list->nelements == 1 )
               {
                  list->first = list->last = NULL ;
               }
               /* . The object is the first on the list. */
               else if ( i == 0 )
               {
                  list->first = current->next ;
               }
               /* . The object is the last on the list. */
               else if ( i == ( list->nelements - 1 ) )
               {
                  list->last       = previous ;
                  list->last->next = NULL     ;
               }
               /* . Redirect the next pointer of previous. */
               else
               {
                  previous->next = current->next ;
               }

               /* . Save the data but remove the list element. */
	       node = current->node ;
               current->node = current->next = NULL ;
               Memory_Deallocate ( current ) ;

               /* . Reset nelements. */
               list->nelements -- ;

               break ;
            }
         }
      }
   }
   return node ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find an element on the list by using the list's match function, remove the
! . element and return its node address. The node is not deallocated.
!---------------------------------------------------------------------------------------------------------------------------------*/
void *List_Element_Pop_By_Match ( List *list, void *reference )
{
   Integer      i ;
   ListElement *current = NULL, *next = NULL, *previous = NULL ;
   void        *node    = NULL ;
   for ( i = 0 ; i < list->nelements ; i++ )
   {
      if ( i == 0 )
      {
         current = list->first       ;
         next    = list->first->next ;
      }
      else
      {
         previous = current    ;
         current  = next       ;
         next     = next->next ;
      }
      if ( (*list->Element_Match) ( reference, current->node ) )
      {
         /* . The lists are now empty. */
	 if ( list->nelements == 1 )
         {
            list->first = list->last = NULL ;
         }
         /* . The object is the first on the list. */
         else if ( i == 0 )
         {
            list->first = current->next ;
         }
         /* . The object is the last on the list. */
         else if ( i == ( list->nelements - 1 ) )
         {
            list->last       = previous ;
            list->last->next = NULL     ;
         }
         /* . Redirect the next pointer of previous. */
         else
         {
            previous->next = current->next ;
         }

         /* . Save the data but remove the list element. */
	 node = current->node ;
         current->node = current->next = NULL ;
         Memory_Deallocate ( current ) ;

         /* . Reset nelements. */
         list->nelements -- ;

         break ;
      }
   }
   return node ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Empty a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Empty ( List *self )
{
    if ( self != NULL )
    {
        auto Integer      i ;
        auto ListElement *current = NULL, *next = NULL ;
        for ( i = 0 ; i < self->nelements ; i++ )
        {
            if ( i == 0 )
            {
                current = self->first       ;
                next    = self->first->next ;
            }
            else
            {
                current = next       ;
                next    = next->next ;
             }
             (*self->Element_Deallocate) ( current->node ) ;
             Memory_Deallocate ( current ) ;
        }
        List_Initialize ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Initialize ( List *self )
{
   if ( self != NULL )
   {
       self->nelements =    0 ;
       self->first     = NULL ;
       self->iterator  = NULL ;
       self->last      = NULL ;
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iterate over the elements in the list.
!---------------------------------------------------------------------------------------------------------------------------------*/
void *List_Iterate ( List *list )
{
   if ( list->iterator == NULL ) list->iterator = list->first ;
   else list->iterator = list->iterator->next ;
   if ( list->iterator == NULL ) return NULL ;
   else                          return list->iterator->node ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the current position of the list iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
ListElement *List_Iterate_Current ( const List *list ) { return list->iterator ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize the list iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Iterate_Initialize ( List *list ) { list->iterator = NULL ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the current position of the list iterator.
! . This procedure requires careful use.
!---------------------------------------------------------------------------------------------------------------------------------*/
void List_Iterate_Set ( List *list, ListElement *iterator ) { list->iterator = iterator ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the size of a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
int List_Size ( const List *list ) { return list->nelements ; }

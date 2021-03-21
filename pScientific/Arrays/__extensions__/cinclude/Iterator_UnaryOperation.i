/* . Need:

_InnerLoopOperation
_Locals
_Operation
_Reduction
_ReductionResultInitializer
_ReductionResultType
_SelfIsModified

*/

# ifdef _Reduction
    auto _ReductionResultType result = _ReductionResultInitializer ;
# endif
    if ( self != NULL )
    {
        auto const IteratorType *type  = self->type  ;
        auto       void         *state = self->vSelf ;
# ifdef _Locals
        _Locals
# endif
        type->Reset ( state ) ;
/*
# ifdef _SelfIsModified
        if ( self->readOnly ) Status_Set ( status, Status_InvalidArrayOperation ) ; 
        else
# endif
*/
        if ( self->extent > 1 )
        {
            auto Boolean ( *next ) ( void*, Integer*, Integer*, Integer* ) = type->NextInnerLoop ;
            auto Integer extent, first, stride ;
# ifndef _InnerLoopOperation
            auto Integer e, i ;
# endif
            while ( next ( state, &first, &extent, &stride ) )
            {
# ifdef _InnerLoopOperation
                _InnerLoopOperation
# else
                for ( e = 0, i = first ; e < extent ; e++, i += stride ) _Operation
# endif
            }
        }
        else
        {
            auto Integer ( *next ) ( void* ) = type->NextIndex ;
            auto Integer e, i ;
            for ( e = 0, i = next ( state ) ; e < self->size ; e++, i = next ( state ) ) _Operation
        }
    }
# ifdef _Reduction
    return result ;
# endif

/* . Need:

_InnerLoopOperation
_Locals
_Operation
_OtherIsModified
_Reduction
_ReductionResultInitializer
_ReductionResultType
_SelfIsModified

*/

/*

  Loop types:

  1 - Inner loop 1 / inner loop 2 with e1 >= e2, e1 / e2 = m, and e1 % e2 = 0

     ( R1 and R2 ) and ( e1 >= e2 ) and ( e1 % e2 = 0 )

  2 - Inner loop 1 / inner loop 2 with e1 <  e2, e2 / e1 = m, and e2 % e1 = 0

     ( R1 and R2 ) and ( e1 < e2 ) and ( e2 % e1 = 0 )

  3 - Inner loop 1 / index 2

     ( I1 and I2 ) and ( e1 >= e2 )
     ( R1 and I2 )
     ( R1 and R2 ) and ( e1 > e2 ) and ( e1 % e2 != 0 )

  4 - Index 1 / inner loop 2

     ( I1 and I2 ) and ( e1 < e2 )
     ( I1 and R2 )
     ( R1 and R2 ) and ( e1 < e2 ) and ( e2 % e1 != 0 )

  Default (0) - Index 1 / index 2

     ( I1 and I2 ) and ( e1 == e2 == 1 )

*/

# ifdef _Reduction
    auto _ReductionResultType result = _ReductionResultInitializer ;
# endif
    if ( ( self != NULL ) && ( other != NULL ) && ( self->size > 0 ) )
    {
        auto const IteratorType *type1  = self->type , *type2  = other->type  ;
        auto       void         *state1 = self->vSelf, *state2 = other->vSelf ;
        type1->Reset ( state1 ) ; type2->Reset ( state2 ) ;
        if      ( self->size != other->size ) Status_Set ( status, Status_NonConformableArrays  ) ;
        else if ( self       == other       ) Status_Set ( status, Status_InvalidArrayOperation ) ;
/*
# ifdef _SelfIsModified
        else if ( self->readOnly  ) Status_Set ( status, Status_InvalidArrayOperation ) ; 
# endif
*/
/*
# ifdef _OtherIsModified
        else if ( other->readOnly ) Status_Set ( status, Status_InvalidArrayOperation ) ; 
# endif
*/
        else
        {
            auto div_t   e12 ;
            auto Integer extent1 = self->extent,  extent2 = other->extent, lType = 0, nInner = 0 ;
# ifdef _Locals
            _Locals
# endif
            if ( self->isRegular && other->isRegular )
            {
                if ( extent1 >= extent2 )
                {
                    e12 = div ( extent1, extent2 ) ;
                    if ( e12.rem == 0 ) { lType = 1 ; nInner = e12.quot ; }
                    else                { lType = 3 ; }
                }
                else
                {
                    e12 = div ( extent2, extent1 ) ;
                    if ( e12.rem == 0 ) { lType = 2 ; nInner = e12.quot ; }
                    else                { lType = 4 ; }
                }
            }
            else if (  self->isRegular ) { lType = 3 ; }
            else if ( other->isRegular ) { lType = 4 ; }
            else if ( ( extent1 > 1 ) && ( extent1 >= extent2 ) ) { lType = 3 ; }
            else if ( ( extent2 > 1 ) && ( extent2 >  extent1 ) ) { lType = 4 ; }
            /* . Inner loop 1 / inner loop 2 with e1 > e2, e1 / e2 = m, and e1 % e2 = 0. */
            if ( lType == 1 )
            {
                auto Integer extent, extent1, first1, first2, l1, l2, stride1, stride2 ;
# ifndef _InnerLoopOperation
                auto Integer e, i1, i2 ;
# endif
                for ( l1 = 0 ; l1 < self->numberOfLoops ; l1++ )
                {
                    type1->NextInnerLoop ( state1, &first1, &extent1, &stride1 ) ;
                    for ( l2 = 0 ; l2 < nInner ; l2++ )
                    {
                        type2->NextInnerLoop ( state2, &first2, &extent, &stride2 ) ;
# ifdef _InnerLoopOperation
                        _InnerLoopOperation
# else
                        for ( e = 0, i1 = first1, i2 = first2 ; e < extent ; e++, i1 += stride1, i2 += stride2 ) _Operation
# endif
                        first1 += ( extent * stride1 ) ;
                    }
                }
            }
            /* . Inner loop 1 / inner loop 2 with e1 < e2, e2 / e1 = m, and e2 % e1 = 0. */
            else if ( lType == 2 )
            {
                auto Integer extent, extent2, first1, first2, l1, l2, stride1, stride2 ;
# ifndef _InnerLoopOperation
                auto Integer e, i1, i2 ;
# endif
                for ( l2 = 0 ; l2 < other->numberOfLoops ; l2++ )
                {
                    type2->NextInnerLoop ( state2, &first2, &extent2, &stride2 ) ;
                    for ( l1 = 0 ; l1 < nInner ; l1++ )
                    {
                        type1->NextInnerLoop ( state1, &first1, &extent, &stride1 ) ;
# ifdef _InnerLoopOperation
                        _InnerLoopOperation
# else
                        for ( e = 0, i1 = first1, i2 = first2 ; e < extent ; e++, i1 += stride1, i2 += stride2 ) _Operation
# endif
                        first2 += ( extent * stride2 ) ;
                    }
                }
            }
            /* . Inner loop 1 / index 2. */
            else if ( lType == 3 )
            {
                auto Boolean ( *next1 ) ( void*, Integer*, Integer*, Integer* ) = type1->NextInnerLoop ;
                auto Integer ( *next2 ) ( void* ) = type2->NextIndex ;
                auto Integer e, extent1, first1, i1, i2 = next2 ( state2 ), stride1 ;
                while ( next1 ( state1, &first1, &extent1, &stride1 ) )
                {
                    for ( e = 0, i1 = first1 ; e < extent1 ; e++, i1 += stride1, i2 = next2 ( state2 ) ) _Operation
                }
            }
            /* . Index 1 / inner loop 2. */
            else if ( lType == 4 )
            {
                auto Boolean ( *next2 ) ( void*, Integer*, Integer*, Integer* ) = type2->NextInnerLoop ;
                auto Integer ( *next1 ) ( void* ) = type1->NextIndex ;
                auto Integer e, extent2, first2, i1 = next1 ( state1 ), i2, stride2 ;
                while ( next2 ( state2, &first2, &extent2, &stride2 ) )
                {
                    for ( e = 0, i2 = first2 ; e < extent2 ; e++, i1 = next1 ( state1 ), i2 += stride2 ) _Operation
                }
            }
            /* . Index 1 / index 2. */
            else
            {
                auto Integer ( *next1 ) ( void* ) = type1->NextIndex, ( *next2 ) ( void* ) = type2->NextIndex ;
                auto Integer e, i1 = next1 ( state1 ), i2 = next2 ( state2 ) ;
                for ( e = 0 ; e < self->size ; e++, i1 = next1 ( state1 ), i2 = next2 ( state2 ) ) _Operation
            }
        }
    }
# ifdef _Reduction
    return result ;
# endif

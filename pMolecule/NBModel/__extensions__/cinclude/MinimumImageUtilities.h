# ifndef _MINIMUMIMAGEUTILITIES
# define _MINIMUMIMAGEUTILITIES

/* . Utilities for minimum image pairwise interactions. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MinimumImage_Allocate( coordinates3I, coordinates3J ) \
    auto Integer       c ; \
    auto Real          d, fI[3], t ; \
    auto Coordinates3 *fractionalI, *fractionalJ ; \
    auto RealArray2D  *fDisplacements, *rDisplacements, *translations ; \
    auto RealArray2D   fView, rView ; \
    n              = PairList_MaximumRecordSize ( pairList ) ; \
    fractionalI    = SymmetryParameters_MakeFractionalCoordinates ( symmetryParameters, coordinates3I, status ) ; \
    fractionalJ    = SymmetryParameters_MakeFractionalCoordinates ( symmetryParameters, coordinates3J, status ) ; \
    fDisplacements = RealArray2D_AllocateWithExtents ( n, 3, status ) ; \
    rDisplacements = RealArray2D_AllocateWithExtents ( n, 3, status ) ; \
    translations   = RealArray2D_AllocateWithExtents ( n, 3, status ) ; \
    if ( ( fractionalI    != NULL ) && \
         ( fractionalJ    != NULL ) && \
         ( fDisplacements != NULL ) && \
         ( rDisplacements != NULL ) && \
         ( translations   != NULL ) ) \
    {

# define MinimumImage_Deallocate \
    } \
    Coordinates3_Deallocate ( &fractionalI    ) ; \
    Coordinates3_Deallocate ( &fractionalJ    ) ; \
    RealArray2D_Deallocate  ( &fDisplacements ) ; \
    RealArray2D_Deallocate  ( &rDisplacements ) ; \
    RealArray2D_Deallocate  ( &translations   ) ; \

# define MinimumImage_Displacements( i, j ) \
    for ( c = 0 ; c < 3 ; c++ ) fI[c] = Coordinates3_Item ( fractionalI, i, c ) ; \
    for ( n = 0 ; n < record->capacity ; n++ ) \
    { \
        j = record->indices[n] ; \
        for ( c = 0 ; c < 3 ; c++ ) \
        { \
            d = fI[c] - Coordinates3_Item ( fractionalJ, j, c ) ; \
            t = 0.0e+00 ; \
            if      ( d >=  0.5 ) { d -= 1.0 ; t =  1.0e+00 ; } \
            else if ( d <  -0.5 ) { d += 1.0 ; t = -1.0e+00 ; } \
            Coordinates3_Item ( fDisplacements, n, c ) = d ; \
            Coordinates3_Item ( translations  , n, c ) = t ; \
        } \
    } \
    RealArray2D_View ( fDisplacements, 0, 0, record->capacity, 3, 1, 1, False, &fView, status ) ; \
    RealArray2D_View ( rDisplacements, 0, 0, record->capacity, 3, 1, 1, False, &rView, status ) ; \
    RealArray2D_MatrixMultiply ( False, True, 1.0e+00, fDisplacements, symmetryParameters->H, 0.0e+00, rDisplacements, status ) ; \
    RealArray2D_Set ( fDisplacements, 0.0e+00 ) ;

# define MinimumImage_Gradients \
    RealArray2D_View ( fDisplacements, 0, 0, record->capacity, 3, 1, 1, False, &fView, status ) ; \
    RealArray2D_View ( translations  , 0, 0, record->capacity, 3, 1, 1, False, &rView, status ) ; \
    RealArray2D_MatrixMultiply ( True, False, 1.0e+00, &fView, &rView, 1.0e+00, symmetryParameterGradients->dEdH, status ) ;

# endif

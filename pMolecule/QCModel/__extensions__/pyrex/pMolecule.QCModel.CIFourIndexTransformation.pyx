"""TEI four-index transformation functions."""

# . Should probably check dimensions here.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_MOTEIs   = "moTwoElectronIntegrals" 
_MOTEI34  = "moTEI34"
_MOTEI234 = "moTEI234"

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def CIFourIndexTransformation_Make ( RealArray2D  activeMOs            not None ,
                                     BlockStorage twoElectronIntegrals not None ,
                                                  store                not None ):
    """Perform a four-index transformation."""
    cdef DoubleSymmetricMatrix moTEIs   = None
    cdef RealArray2D           moTEI34  = None
    cdef RealArrayND           moTEI234 = None
    ( nBasis, nActive ) = activeMOs.shape
    nATr   = ( nActive * ( nActive + 1 ) ) / 2 ;
    nBTr   = ( nBasis  * ( nBasis  + 1 ) ) / 2 ;
    moTEIs = store.Get ( _MOTEIs, None )
    if moTEIs is None:
        moTEIs = DoubleSymmetricMatrix.WithExtent ( nActive )
        store.Set ( _MOTEIs, moTEIs )
    moTEI34 = store.Get ( _MOTEI34, None )
    if moTEI34 is None:
        moTEI34 = RealArray2D.WithExtents ( nBTr, nATr )
        store.Set ( _MOTEI34, moTEI34 )
    moTEI234 = store.Get ( _MOTEI234, None )
    if moTEI234 is None:
        moTEI234 = RealArrayND.WithExtents ( nBasis, nActive, nATr )
        store.Set ( _MOTEI234, moTEI234 )
    CIFourIndexTransformation ( activeMOs.cObject            ,
                                twoElectronIntegrals.cObject ,
                                moTEI34.cObject              ,
                                moTEI234.cObject             ,
                                moTEIs.cObject               )
    return moTEIs

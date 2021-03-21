"""MNDO multipole evaluator."""

# . No weighted density matrix contributions as S = I.

from .QCDefinitions import BasisRepresentation

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOMultipoleEvaluator:

    def BondOrders ( self                                ,
                     target                              ,
                     SymmetricMatrix density    not None ,
                     SymmetricMatrix bondOrders not None ):
        """Bond orders."""
        cdef IntegerArray1D basisIndices
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        MNDO_BondOrders ( basisIndices.cObject ,
                          density.cObject      ,
                          bondOrders.cObject   )

    def FockMultipoleDerivatives ( self                          ,
                                   target                        , 
                                   multipoleOrder                ,
                                   RealArray1D     dXdQ not None ,
                                   SymmetricMatrix dXdP not None ):
        """Multipole density derivatives."""
        cdef IntegerArray1D          basisIndices
        cdef MNDOParametersContainer parameters
        if ( multipoleOrder < 0 ) or ( multipoleOrder > 2 ): raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        parameters   = target.qcState.mndoParameters
        MNDO_AtomicMultipolesFock ( parameters.cObject   ,
                                    basisIndices.cObject ,
                                    dXdQ.cObject         ,
                                    multipoleOrder       ,
                                    dXdP.cObject         )

    def FockMultipoles ( self                                ,
                         target                              ,
                         multipoleOrder                      ,
                         RealArray1D     multipoles not None ,
                         SymmetricMatrix density    = None   ,
                                         withNuclei = True   ):
        """Multipoles for interactions with external potentials."""
        cdef IntegerArray1D          basisIndices
        cdef MNDOParametersContainer parameters
        cdef SymmetricMatrix         onePDM
        if ( multipoleOrder < 0 ) or ( multipoleOrder > 2 ): raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        parameters   = target.qcState.mndoParameters
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        MNDO_AtomicMultipoles ( parameters.cObject                ,
                                basisIndices.cObject              ,
                                onePDM.cObject                    ,
                                MultipoleRepresentation_Cartesian ,
                                multipoleOrder                    ,
                                multipoles.cObject                )
        if withNuclei:
            Z = target.qcState.nuclearCharges
            multipoles[0:len ( Z )].Add ( Z )

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Actual

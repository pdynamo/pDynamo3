"""A container for CMAP dihedral energy terms."""

from pCore import Clone, logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CMAPDihedralContainer ( MMTerm ):

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef CMAPDihedralContainer new
        new               = self.__class__.Raw ( )
        new.label         = self.label
        new.cObject       = CMAPDihedralContainer_Clone ( self.cObject )
        new.isOwner       = True
        new.parameterKeys = Clone ( self.parameterKeys )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            CMAPDihedralContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i, iX, iY, lengthX, lengthY
        # . Parameters.
        parameterFields = [ "x", "y", "map" ]
        parameters      = []
        for i from 0 <= i < self.cObject.nParameters:
            x = []
            y = []
            f = []
            lengthX = self.cObject.parameters[i].lengthX
            lengthY = self.cObject.parameters[i].lengthY
            for iX from 0 <= iX < lengthX:
                x.append ( RealArray1D_GetItem ( self.cObject.parameters[i].x, iX, NULL ) )
            for iY from 0 <= iY < lengthY:
                y.append ( RealArray1D_GetItem ( self.cObject.parameters[i].y, iY, NULL ) )
            for iX from 0 <= iX < lengthX:
                for iY from 0 <= iY < lengthY:
                   f.append ( RealArray2D_GetItem ( self.cObject.parameters[i].f, iX, iY, NULL ) )
            parameters.append ( [ x, y, f ] )
        # . Terms.
        termFields = [ "atom1", "atom2", "atom3", "atom4", "atom5", "atom6", "atom7", "atom8", "parameter", "isActive" ]
        terms      = []
        for i from 0 <= i < self.cObject.nTerms:
            terms.append ( [ self.cObject.terms[i].atom1 ,
                             self.cObject.terms[i].atom2 ,
                             self.cObject.terms[i].atom3 ,
                             self.cObject.terms[i].atom4 ,
                             self.cObject.terms[i].atom5 ,
                             self.cObject.terms[i].atom6 ,
                             self.cObject.terms[i].atom7 ,
                             self.cObject.terms[i].atom8 ,
                             self.cObject.terms[i].type  ,
                             self.cObject.terms[i].isActive == CTrue ] )
        # . State.
        state = { "label"           : self.label      ,
                  "parameterFields" : parameterFields ,
                  "parameters"      : parameters      ,
                  "termFields"      : termFields      ,
                  "terms"           : terms           }
        if self.parameterKeys is not None: state["parameterKeys"] = self.parameterKeys
        return state

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "CMAP Dihedral" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfParameters, numberOfTerms )
        self.label = label

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger      i, iX, iY, lengthX, lengthY
        cdef CRealArray1D *x
        cdef CRealArray1D *y
        cdef CRealArray2D *f
        # . Allocate the object.
        parameters = state["parameters"]
        terms      = state["terms"     ]
        self._Allocate ( len ( parameters ), len ( terms ) )
        if "label" in state: self.label = state["label"]
        self.parameterKeys = state.get ( "parameterKeys", None )
        # . Fill the object.
        # . Parameters.
        for ( i, ( xP, yP, fP ) ) in enumerate ( parameters ):
            lengthX = len ( xP )
            lengthY = len ( yP )
            x = RealArray1D_AllocateWithExtent  ( lengthX, NULL )
            y = RealArray1D_AllocateWithExtent  ( lengthY, NULL )
            f = RealArray2D_AllocateWithExtents ( lengthX, lengthY, NULL )
            for iX from 0 <= iX < lengthX:
                RealArray1D_SetItem ( x, iX, xP[iX], NULL )
            for iY from 0 <= iY < lengthY:
                RealArray1D_SetItem ( y, iY, yP[iY], NULL )
            n = 0
            for iX from 0 <= iX < lengthX:
                for iY from 0 <= iY < lengthY:
                    RealArray2D_SetItem ( f, iX, iY, fP[n], NULL )
                    n += 1
            self.cObject.parameters[i] = BicubicSpline_MakeFromRealArray2D ( &x, &y, &f, BicubicSplineType_Periodic, NULL )
        # . Terms.
        for ( i, ( a1, a2, a3, a4, a5, a6, a7, a8, p, q ) ) in enumerate ( terms ):
            self.cObject.terms[i].atom1 = a1
            self.cObject.terms[i].atom2 = a2
            self.cObject.terms[i].atom3 = a3
            self.cObject.terms[i].atom4 = a4
            self.cObject.terms[i].atom5 = a5
            self.cObject.terms[i].atom6 = a6
            self.cObject.terms[i].atom7 = a7
            self.cObject.terms[i].atom8 = a8
            self.cObject.terms[i].type  = p
            if q: self.cObject.terms[i].isActive = CTrue
            else: self.cObject.terms[i].isActive = CFalse

    def _Allocate ( self, numberOfParameters, numberOfTerms ):
        """Allocation."""
        self.cObject = CMAPDihedralContainer_Allocate ( numberOfTerms, numberOfParameters )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject       = NULL
        self.isOwner       = False
        self.label         = "CMAP Dihedral"
        self.parameterKeys = None

    def ActivateTerms ( self ):
        """Activate terms."""
        CMAPDihedralContainer_ActivateTerms ( self.cObject )

    def DeactivateFixedAtomTerms ( self, Selection freeAtoms ):
        """Deactivate terms involving fixed atoms."""
        if freeAtoms is not None: CMAPDihedralContainer_DeactivateTerms ( self.cObject, freeAtoms.cObject )

    def DeactivateQCAtomTerms ( self, Selection mmAtoms, Selection pureMMAtoms ):
        """Deactivate terms involving QC atoms."""
        if pureMMAtoms is not None: CMAPDihedralContainer_DeactivateTerms ( self.cObject, pureMMAtoms.cObject )

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        if gradients3 is None: energy = CMAPDihedralContainer_Energy ( self.cObject, coordinates3.cObject, NULL               )
        else:                  energy = CMAPDihedralContainer_Energy ( self.cObject, coordinates3.cObject, gradients3.cObject )
        return energy

    def GetUniqueTermIndices ( self ):
        """Return the term indices."""
        cdef CInteger i
        indices = []
        for i from 0 <= i < self.cObject.nTerms:
            indices.append ( self.cObject.terms[i].atom1 )
            indices.append ( self.cObject.terms[i].atom2 )
            indices.append ( self.cObject.terms[i].atom3 )
            indices.append ( self.cObject.terms[i].atom4 )
            indices.append ( self.cObject.terms[i].atom5 )
            indices.append ( self.cObject.terms[i].atom6 )
            indices.append ( self.cObject.terms[i].atom7 )
            indices.append ( self.cObject.terms[i].atom8 )
        return ( self.cObject.nTerms, indices )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging - slow version."""
        # . Initialization.
        new = None
        # . Get the selection increments.
        increments = information.get ( "Index Increments", None )
        # . Do the merge.
        if increments is not None:
            # . Gather terms.
            np         = 0
            parameters = []
            terms      = []
            for ( increment, item ) in zip ( increments, items ):
                if item is not None:
                    state = item.__getstate__ ( )
                    if np == 0:
                        label           = state["label"          ]
                        parameterFields = state["parameterFields"]
                        termFields      = state["termFields"     ]
                    parameters0 = state["parameters"]
                    terms0      = state["terms"     ]
                    np0         = len ( parameters0 )
                    parameters.extend ( parameters0 )
                    for ( a1, a2, a3, a4, a5, a6, a7, a8, p, q ) in terms0:
                        terms.append ( ( a1 + increment ,
                                         a2 + increment ,
                                         a3 + increment ,
                                         a4 + increment ,
                                         a5 + increment ,
                                         a6 + increment ,
                                         a7 + increment ,
                                         a8 + increment ,
                                         p  + np , q ) )
                    np += np0
            # . Construct the object.
            state = { "label"           : label           ,
                      "parameterFields" : parameterFields ,
                      "parameters"      : parameters      ,
                      "termFields"      : termFields      ,
                      "terms"           : terms           }
            new = selfClass.Raw ( )
            new.__setstate__ ( state )
        return new

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef CMAPDihedralContainer new
        new         = CMAPDihedralContainer.Raw ( )
        new.cObject = CMAPDihedralContainer_Prune ( self.cObject, selection.cObject )
        new.isOwner = True
        new.label   = self.label
        if new.cObject == NULL: return None
        else:                   return new

    def Sort ( self ):
        """Sorting."""
        CMAPDihedralContainer_Sort ( self.cObject )

    def SummaryItems ( self ):
        """Summary entry."""
        cdef CInteger n
        n = CMAPDihedralContainer_NumberOfInactiveTerms ( self.cObject )
        items = [ ( self.label + " Terms",      "{:d}".format ( self.cObject.nTerms      ) ) ,
                  ( self.label + " Parameters", "{:d}".format ( self.cObject.nParameters ) ) ]
        if n > 0: items.append ( ( self.label + " Inactive", "{:d}".format ( n ) ) )
        return items

    def UpperBound ( self ):
        return CMAPDihedralContainer_UpperBound ( self.cObject )

    @property
    def numberOfParameters ( self ):
        """Return the number of parameters."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nParameters
    @property
    def numberOfTerms ( self ):
        """Return the number of terms."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nTerms

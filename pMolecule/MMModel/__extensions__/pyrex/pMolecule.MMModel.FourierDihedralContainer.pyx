"""A container for Fourier dihedrals."""

from pCore import Clone, logFile, LogFileActive

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class FourierDihedralContainer ( MMTerm ):

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef FourierDihedralContainer new
        new               = self.__class__.Raw ( )
        new.label         = self.label
        new.cObject       = FourierDihedralContainer_Clone ( self.cObject )
        new.isOwner       = True
        new.parameterKeys = Clone ( self.parameterKeys )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            FourierDihedralContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i
        # . Parameters.
        parameterFields = [ "forceConstant", "period", "phase" ]
        parameters      = []
        for i from 0 <= i < self.cObject.nParameters:
            parameters.append ( [ self.cObject.parameters[i].fc     ,
                                  self.cObject.parameters[i].period ,
                                  self.cObject.parameters[i].phase  ] )
        # . Terms.
        termFields = [ "atom1", "atom2", "atom3", "atom4", "parameter", "isActive" ]
        terms      = []
        for i from 0 <= i < self.cObject.nTerms:
            terms.append ( [ self.cObject.terms[i].atom1 ,
                             self.cObject.terms[i].atom2 ,
                             self.cObject.terms[i].atom3 ,
                             self.cObject.terms[i].atom4 ,
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

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "Fourier Dihedral" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfParameters, numberOfTerms )
        self.label = label

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger i
        # . Allocate the object.
        parameters = state["parameters"]
        terms      = state["terms"     ]
        self._Allocate ( len ( parameters ), len ( terms ) )
        if "label" in state: self.label = state["label"]
        self.parameterKeys = state.get ( "parameterKeys", None )
        # . Fill the object.
        # . Parameters.
        for ( i, ( fc, period, phase ) ) in enumerate ( parameters ):
            self.cObject.parameters[i].fc     = fc
            self.cObject.parameters[i].period = period
            self.cObject.parameters[i].phase  = phase
        # . Terms.
        for ( i, ( a1, a2, a3, a4, p, q ) ) in enumerate ( terms ):
            self.cObject.terms[i].atom1 = a1
            self.cObject.terms[i].atom2 = a2
            self.cObject.terms[i].atom3 = a3
            self.cObject.terms[i].atom4 = a4
            self.cObject.terms[i].type  = p
            if q: self.cObject.terms[i].isActive = CTrue
            else: self.cObject.terms[i].isActive = CFalse
        # . Finish processing.
        self.FillCosSinPhases ( )

    def _Allocate ( self, numberOfParameters, numberOfTerms ):
        """Allocation."""
        self.cObject = FourierDihedralContainer_Allocate ( numberOfTerms, numberOfParameters )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject       = NULL
        self.isOwner       = False
        self.label         = "Fourier Dihedral"
        self.parameterKeys = None

    def ActivateTerms ( self ):
        """Activate terms."""
        FourierDihedralContainer_ActivateTerms ( self.cObject )

    def DeactivateFixedAtomTerms ( self, Selection freeAtoms ):
        """Deactivate terms involving fixed atoms."""
        if freeAtoms is not None: FourierDihedralContainer_DeactivateTerms ( self.cObject, freeAtoms.cObject )

    def DeactivateQCAtomTerms ( self, Selection mmAtoms, Selection pureMMAtoms ):
        """Deactivate terms involving QC atoms."""
        if pureMMAtoms is not None: FourierDihedralContainer_DeactivateTerms ( self.cObject, pureMMAtoms.cObject )

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        if gradients3 is None: energy = FourierDihedralContainer_Energy ( self.cObject, coordinates3.cObject, NULL               )
        else:                  energy = FourierDihedralContainer_Energy ( self.cObject, coordinates3.cObject, gradients3.cObject )
        return energy

    def FillCosSinPhases ( self ):
        """Fill the cos and sin phases of the parameter array."""
        FourierDihedralContainer_FillCosSinPhases ( self.cObject )

    def GetUniqueTermIndices ( self ):
        """Return the term indices without redundancy."""
        cdef CInteger i, j, k, l, m
        indices = set ( )
        for m from 0 <= m < self.cObject.nTerms:
            i = self.cObject.terms[m].atom1
            j = self.cObject.terms[m].atom2
            k = self.cObject.terms[m].atom3
            l = self.cObject.terms[m].atom4
            indices.add ( ( i, j, k, l ) )
        indices = list ( indices )
        indices.sort ( )
        newIndices = []
        for ijkl in indices: newIndices.extend ( ijkl )
        return ( len ( indices ), newIndices )

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
            doKeys        = True
            np            = 0
            parameterKeys = []
            parameters    = []
            terms         = []
            for ( increment, item ) in zip ( increments, items ):
                if item is not None:
                    state = item.__getstate__ ( )
                    if np == 0:
                        label           = state["label"          ]
                        parameterFields = state["parameterFields"]
                        termFields      = state["termFields"     ]
                    parameterKeys0 = state.get ( "parameterKeys", None )
                    parameters0    = state["parameters"]
                    terms0         = state["terms"     ]
                    np0            = len ( parameters0 )
                    parameters.extend ( parameters0 )
                    for ( a1, a2, a3, a4, p, q ) in terms0:
                        terms.append ( ( a1 + increment ,
                                         a2 + increment ,
                                         a3 + increment ,
                                         a4 + increment ,
                                         p  + np , q ) )
                    if parameterKeys0 is None: doKeys = False
                    else: parameterKeys.extend ( parameterKeys0 )
                    np += np0
            # . Construct the new keys and parameters.
            if doKeys:
                ( doKeys, oldToNew, parameterKeys, parameters ) = selfClass.MergeKeys ( parameterKeys, parameters )
            else:
                parameterKeys = None
            # . Construct the object.
            state = { "label"           : label           ,
                      "parameterFields" : parameterFields ,
                      "parameterKeys"   : parameterKeys   ,
                      "parameters"      : parameters      ,
                      "termFields"      : termFields      ,
                      "terms"           : terms           }
            new = selfClass.Raw ( )
            new.__setstate__ ( state )
            # . Reindex.
            if doKeys: new.ResetTermTypes ( oldToNew )
        return new

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef FourierDihedralContainer new
        new               = FourierDihedralContainer.Raw ( )
        new.cObject       = FourierDihedralContainer_Prune ( self.cObject, selection.cObject )
        new.isOwner       = True
        new.label         = self.label
        new.parameterKeys = self.parameterKeys
        if new.cObject == NULL: return None
        else:                   return new

    def ResetTermTypes ( self, mapping ):
        """Reset term types."""
        cdef CInteger i
        for i from 0 <= i < self.cObject.nTerms:
            iOld = self.cObject.terms[i].type
            self.cObject.terms[i].type = mapping[iOld]

    def Sort ( self ):
        """Sorting."""
        FourierDihedralContainer_Sort ( self.cObject )

    def SummaryItems ( self ):
        """Summary entry."""
        cdef CInteger n
        n = FourierDihedralContainer_NumberOfInactiveTerms ( self.cObject )
        items = [ ( self.label + " Terms",      "{:d}".format ( self.cObject.nTerms      ) ) ,
                  ( self.label + " Parameters", "{:d}".format ( self.cObject.nParameters ) ) ]
        if n > 0: items.append ( ( self.label + " Inactive", "{:d}".format ( n ) ) )
        return items

    def UpperBound ( self ):
        return FourierDihedralContainer_UpperBound ( self.cObject )

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

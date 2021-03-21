"""Containers for cosine MM terms."""

from pCore import Clone, logFile, LogFileActive

#===================================================================================================================================
# . Base class.
#===================================================================================================================================
cdef class CosineTermContainer ( MMTerm ):

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef CosineTermContainer new
        new               = self.__class__.Raw ( )
        new.label         = self.label
        new.cObject       = CosineTermContainer_Clone ( self.cObject )
        new.isOwner       = True
        new.parameterKeys = Clone ( self.parameterKeys )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            CosineTermContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i, t
        # . Parameters.
        parameterFields = [ "coefficient", "period" ]
        parameters      = []
        for i from 0 <= i < self.cObject.nParameters:
            data = []
            for t from 0 <= t < self.cObject.parameters[i].nTerms:
                data.append ( [ self.cObject.parameters[i].termCoefficients[t], self.cObject.parameters[i].periods[t] ] )
            parameters.append ( data )
        # . Terms.
        termFields = [ "indices", "parameter", "isActive" ]
        terms      = []
        for t from 0 <= t < self.cObject.nTerms:
            indices = []
            for i from 0 <= i < self.cObject.nIndices:
                indices.append ( self.cObject.terms[t].indices[i] )
            terms.append ( [ indices                                 ,
                             self.cObject.terms[i].type              ,
                             self.cObject.terms[i].isActive == CTrue ] )
        # . State.
        state = { "indices"         : self.cObject.nIndices ,
                  "label"           : self.label            ,
                  "parameterFields" : parameterFields       ,
                  "parameters"      : parameters            ,
                  "termFields"      : termFields            ,
                  "terms"           : terms                 }
        if self.parameterKeys is not None: state["parameterKeys"] = self.parameterKeys
        return state

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "Cosine Term" ):
        """Constructor."""
        pass

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger a, i, n, p, t
        cdef CReal c
        # . Allocate the object.
        n          = state["indices"   ]
        parameters = state["parameters"]
        terms      = state["terms"     ]
        self._Allocate ( n, len ( parameters ), len ( terms ) )
        if "label" in state: self.label = state["label"]
        self.parameterKeys = state.get ( "parameterKeys", None )
        # . Fill the object.
        # . Parameters.
        for ( i, data ) in enumerate ( parameters ):
            t = len ( data )
            CosineParameter_Allocate ( &(self.cObject.parameters[i]), t )
            for ( t, ( c, p ) ) in enumerate ( data ):
                self.cObject.parameters[i].termCoefficients[t] = c
                self.cObject.parameters[i].periods         [t] = p
        # . Terms.
        for ( i, ( indices, p, q ) ) in enumerate ( terms ):
            for ( t, a ) in enumerate ( indices ):
                self.cObject.terms[i].indices[t] = a
            self.cObject.terms[i].type  = p
            if q: self.cObject.terms[i].isActive = CTrue
            else: self.cObject.terms[i].isActive = CFalse
        # . Finish processing.
        CosineTermContainer_MakePowers ( self.cObject )

    def _Allocate ( self, numberOfIndices, numberOfParameters, numberOfTerms ):
        """Allocation."""
        self.cObject = CosineTermContainer_Allocate ( numberOfIndices, numberOfTerms, numberOfParameters )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject       = NULL
        self.isOwner       = False
        self.label         = "Cosine Term"
        self.parameterKeys = None

    def ActivateTerms ( self ):
        """Activate terms."""
        CosineTermContainer_ActivateTerms ( self.cObject )

    def DeactivateFixedAtomTerms ( self, Selection freeAtoms ):
        """Deactivate terms involving fixed atoms."""
        if freeAtoms is not None: CosineTermContainer_DeactivateTerms ( self.cObject, freeAtoms.cObject )

    def DeactivateQCAtomTerms ( self, Selection mmAtoms, Selection pureMMAtoms ):
        """Deactivate terms involving QC atoms."""
        if pureMMAtoms is not None: CosineTermContainer_DeactivateTerms ( self.cObject, pureMMAtoms.cObject )

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
                        numberOfIndices = state["indices"        ]
                        parameterFields = state["parameterFields"]
                        termFields      = state["termFields"     ]
                    parameterKeys0 = state.get ( "parameterKeys", None )
                    parameters0    = state["parameters"]
                    terms0         = state["terms"     ]
                    np0            = len ( parameters0 )
                    parameters.extend ( parameters0 )
                    for ( indices, p, q ) in terms0:
                        newIndices = [ a + increment for a in indices ]
                        terms.append ( [ newIndices, p + np, q ] )
                    if parameterKeys0 is None: doKeys = False
                    else: parameterKeys.extend ( parameterKeys0 )
                    np += np0
            # . Construct the new keys and parameters.
            if doKeys:
                ( doKeys, oldToNew, parameterKeys, parameters ) = selfClass.MergeKeys ( parameterKeys, parameters )
            else:
                parameterKeys = None
            # . Construct the object.
            state = { "indices"         : numberOfIndices ,
                      "label"           : label           ,
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
        cdef CosineTermContainer new
        new               = CosineTermContainer.Raw ( )
        new.cObject       = CosineTermContainer_Prune ( self.cObject, selection.cObject )
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

    def SummaryItems ( self ):
        """Summary entry."""
        cdef CInteger n
        n = CosineTermContainer_NumberOfInactiveTerms ( self.cObject )
        items = [ ( self.label + " Terms",      "{:d}".format ( self.cObject.nTerms      ) ) ,
                  ( self.label + " Parameters", "{:d}".format ( self.cObject.nParameters ) ) ]
        if n > 0: items.append ( ( self.label + " Inactive", "{:d}".format ( n ) ) )
        return items

    def UpperBound ( self ):
        return CosineTermContainer_UpperBound ( self.cObject )

    # . Properties.
    property numberOfIndices:
        def __get__ ( self ):
            if self.cObject == NULL: return 0
            else:                    return self.cObject.nIndices

    property numberOfParameters:
        def __get__ ( self ):
            if self.cObject == NULL: return 0
            else:                    return self.cObject.nParameters

    property numberOfTerms:
        def __get__ ( self ):
            if self.cObject == NULL: return 0
            else:                    return self.cObject.nTerms


#===================================================================================================================================
# . Specific classes.
#===================================================================================================================================
cdef class CosineAngleContainer ( CosineTermContainer ):

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "Cosine Angle" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( 3, numberOfParameters, numberOfTerms )
        self.label = label

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        cdef CRealArray2D *cGradients3 = NULL
        if gradients3 is not None: cGradients3 = gradients3.cObject
        energy = CosineTermEnergy_Angle ( self.cObject, coordinates3.cObject, cGradients3 )
        return energy

cdef class CosineDihedralContainer ( CosineTermContainer ):

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "Cosine Dihedral" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( 4, numberOfParameters, numberOfTerms )
        self.label = label

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        cdef CRealArray2D *cGradients3 = NULL
        if gradients3 is not None: cGradients3 = gradients3.cObject
        energy = CosineTermEnergy_Dihedral ( self.cObject, coordinates3.cObject, cGradients3 )
        return energy

cdef class CosineOutOfPlaneContainer ( CosineTermContainer ):

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "Cosine Out-Of-Plane" ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( 4, numberOfParameters, numberOfTerms )
        self.label = label

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        cdef CRealArray2D *cGradients3 = NULL
        if gradients3 is not None: cGradients3 = gradients3.cObject
        energy = CosineTermEnergy_OutOfPlane ( self.cObject, coordinates3.cObject, cGradients3 )
        return energy

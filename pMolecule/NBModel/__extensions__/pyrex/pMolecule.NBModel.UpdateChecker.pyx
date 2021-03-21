"""A class for determining whether distance-dependent NB lists need to be updated."""

from pCore                        import Clone                , \
                                         logFile              , \
                                         LogFileActive        , \
                                         RawObjectConstructor
from pMolecule.NBModel            import NBModelError
from pMolecule.NBModel.NBDefaults import _NumberOfCalls       , \
                                         _NumberOfUpdates     , \
                                         _PairListStatistics  , \
                                         _UpdatablePairLists  , \
                                         _UpdateChecker

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class UpdateChecker:
    """A NB list update checker."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            UpdateChecker_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        return { "buffer" : self.buffer }

    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions ( **options )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = UpdateChecker_Allocate ( &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise NBModelError ( "Error allocating update checker." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def Check ( self, target ):
        """Check for an update."""
        cdef Coordinates3            coordinates3
        cdef Coordinates3            rCoordinates3
        cdef ImagePairListContainer  images
        cdef Selection               freeAtoms
        cdef SymmetryParameters      rSymmetryParameters
        cdef SymmetryParameters      symmetryParameters
        cdef CReal                   maximumDisplacement = 0.0
        cdef CSelection             *cFreeAtoms          = NULL
        cdef CSymmetryParameters    *cSymmetryParameters = NULL
        coordinates3       = target.coordinates3
        scratch            = target.scratch
        symmetryParameters = target.symmetryParameters
        # . Update check data.
        uNodeExists        = hasattr ( scratch, _UpdateChecker )
        uNode              = scratch.GetSetNode ( _UpdateChecker )
        if uNodeExists:
            rCoordinates3       = uNode.coordinates3
            rSymmetryParameters = uNode.Get ( "symmetryParameters", None )
        else:
            rCoordinates3      = Coordinates3.WithExtent ( coordinates3.rows )
            uNode.coordinates3 = rCoordinates3
            if symmetryParameters is not None:
                rSymmetryParameters      = SymmetryParameters ( )
                uNode.symmetryParameters = rSymmetryParameters
        # . Pairlist data.
        doUpdate = ( not uNodeExists ) or ( not hasattr ( scratch, _UpdatablePairLists ) )
        pNode    = scratch.GetSetNode ( _UpdatablePairLists )
        if not doUpdate:
            freeAtoms = target.freeAtoms
            if freeAtoms is not None: cFreeAtoms = freeAtoms.cObject
            doUpdate  = ( UpdateChecker_CheckForUpdate ( coordinates3.cObject  ,
                                                         rCoordinates3.cObject ,
                                                         cFreeAtoms            ,
                                                         self.cObject.buffer   ,
                                                         &maximumDisplacement  ) == CTrue )
            if ( not doUpdate ) and ( symmetryParameters is not None ):
                images = pNode.Get ( "mmmmImage", None )
                if images is None:
                    doUpdate = True
                else:
                    doUpdate = ( UpdateChecker_CheckForImageUpdate ( symmetryParameters.cObject  ,
                                                                     rSymmetryParameters.cObject ,
                                                                     images.cObject              ,
                                                                     self.cObject.buffer         ,
                                                                     maximumDisplacement         ) == CTrue )
        # . Statistics and updating.
        sNode = scratch.Get ( _PairListStatistics )
        sNode[_NumberOfCalls] += 1
        if doUpdate:
            scratch.Delete ( _UpdatablePairLists )
            sNode[_NumberOfUpdates] += 1
            coordinates3.CopyTo ( rCoordinates3 )
            if symmetryParameters is not None: symmetryParameters.CopyTo ( rSymmetryParameters )

    @classmethod
    def WithOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "buffer" in options: self.cObject.buffer = options.pop ( "buffer" )
        return options

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Update Checker Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Buffer", "{:.3f}".format ( self.buffer ) ) ]

    @property
    def buffer ( self ):
        if self.cObject == NULL: return 0.0
        else:                    return self.cObject.buffer

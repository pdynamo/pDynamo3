"""A container for transformations of dimension 3."""

from  pCore          import Clone                , \
                            logFile              , \
                            LogFileActive        , \
                            RawObjectConstructor
from .Geometry3Error import Geometry3Error

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Transformation3Container:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef Transformation3Container new
        items = Clone ( self.items )
        new   = self.__class__.WithTransformations ( items )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Transformation3Container_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, CInteger  index ):
        """Get an item."""
        if index >= len ( self ): raise IndexError
        else:                     return self.items[index]

    def __getstate__ ( self ):
        """Return the state."""
        return { "items" : self.items }

    def __init__ ( self, transformations ):
        """Constructor given a set of transformations."""
        self._Initialize ( )
        self.CreateTransformationList ( transformations )
        self._Allocate ( len ( self.items ) )
        self.FillCObject ( )

    def __len__ ( self ):
        """Return the number of transformations."""
        return len ( self.items )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self.items = state["items"]
        self._Allocate ( len ( self.items ) )
        self.FillCObject ( )

    def _Allocate ( self, capacity ):
        """Allocation."""
        cdef CStatus status = CStatus_OK
        self.cObject = Transformation3Container_Allocate ( capacity, &status )
        if status != CStatus_OK: raise Geometry3Error ( "Unable to allocate Transformation3Container." )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.items   = None

    def CreateTransformationList ( self, transformations ):
        """Create a list of unique transformations."""
        cdef Transformation3 new, old
        self.items = []
        for new in transformations:
            isUnique = True
            for old in self.items:
                if ( Transformation3_IsEqual ( new.cObject, old.cObject ) == CTrue ):
                    isUnique = False
                    break
            if isUnique: self.items.append ( new )

    def FillCObject ( self ):
        """Fill the C object from items."""
        cdef CInteger        i
        cdef CStatus         status = CStatus_OK
        cdef Transformation3 t
        for ( i, t ) in enumerate ( self.items ):
            self.cObject.items[i] = t.cObject
        Transformation3Container_FindIdentity ( self.cObject )
        Transformation3Container_FindInverses ( self.cObject, &status )
        if status != CStatus_OK: raise Geometry3Error ( "Unable to fill Transformation3Container." )

    @classmethod
    def Identity ( selfClass ):
        """Construct a container containing only the identity transformation."""
        cdef Transformation3Container self
        identity = Transformation3.Identity ( )
        self     = Transformation3Container.WithTransformations ( [ identity ] )
        return self

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Transformation3 Container Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        # . Proper/improper rotation information also interesting but not necessarily valid unless symmetry parameters exist.
        nNonSelf = Transformation3Container_NumberOfNonSelfInverses ( self.cObject )
        nSelf    = Transformation3Container_NumberOfSelfInverses    ( self.cObject )
        nT3      = len ( self.items )
        return [ ( "Transformations"    , "{:d}".format ( nT3                    ) ) ,
                 ( "Identity"           , "{:d}".format ( self.cObject.identity  ) ) ,
                 ( "Inverse Pairs"      , "{:d}".format ( nNonSelf // 2          ) ) ,
                 ( "Self Inverses"      , "{:d}".format ( nSelf                  ) ) ,
                 ( "Absent Inverses"    , "{:d}".format ( nT3 - nNonSelf - nSelf ) ) ]

    @classmethod
    def WithTransformations ( selfClass, transformations ):
        """Constructor given a set of transformations."""
        return selfClass ( transformations )

    @property
    def isIdentity ( self ):
        """Does this container correspond to the identity?"""
        cdef Transformation3 item
        isIdentity = False
        if len ( self.items ) == 1:
            item       = self.items[0]
            isIdentity = ( Transformation3_IsIdentity ( item.cObject ) == CTrue )
        return isIdentity

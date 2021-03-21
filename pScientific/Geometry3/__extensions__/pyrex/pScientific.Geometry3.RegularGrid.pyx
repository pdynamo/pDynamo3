"""Handle regular grids."""

from pCore import logFile, LogFileActive, RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RegularGrid:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef RegularGrid new
        new = self.__class__.Raw ( )
        new.cObject = RegularGrid_Clone ( self.cObject, NULL )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            RegularGrid_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i
        dimensionFields = [ "isPeriodic", "bins", "stride", "binSize", "lower", "midPointLower", "period", "upper" ]
        dimensions      = []
        for i from 0 <= i < len ( self ):
            dimensions.append ( ( self.cObject.dimensions[i].isPeriodic == CTrue ,
                                  self.cObject.dimensions[i].bins                ,
                                  self.cObject.dimensions[i].stride              ,
                                  self.cObject.dimensions[i].binSize             ,
                                  self.cObject.dimensions[i].lower               ,
                                  self.cObject.dimensions[i].midPointLower       ,
                                  self.cObject.dimensions[i].period              ,
                                  self.cObject.dimensions[i].upper             ) )
        return { "dimensionFields" : dimensionFields, "dimensions" : dimensions }

    def __init__ ( self, rank ):
        """Constructor with rank."""
        self._Initialize ( )
        self._Allocate ( rank )

# . Is this good?
    def __len__ ( self ):
        """Length."""
        return self.rank

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger i
        dimensions = state["dimensions"]
        self._Allocate ( len ( dimensions ) )
        for ( i, datum ) in enumerate ( dimensions ):
            if datum[0]: self.cObject.dimensions[i].isPeriodic = CTrue
            else:        self.cObject.dimensions[i].isPeriodic = CFalse
            self.cObject.dimensions[i].bins          = datum[1]
            self.cObject.dimensions[i].stride        = datum[2]
            self.cObject.dimensions[i].binSize       = datum[3]
            self.cObject.dimensions[i].lower         = datum[4]
            self.cObject.dimensions[i].midPointLower = datum[5]
            self.cObject.dimensions[i].period        = datum[6]
            self.cObject.dimensions[i].upper         = datum[7]

    def _Allocate ( self, rank ):
        """Allocation."""
        self.cObject = RegularGrid_Allocate ( rank, NULL )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

#    def Lengths ( self ):
#        """Return the number of bins in each of the grid dimensions."""
#        return self.shape

    @classmethod
    def FromDimensionData ( selfClass, data ):
        """Constructor from dimension data."""
        cdef CInteger i
        cdef RegularGrid self
        # . Allocation.
        rank = len ( data )
        self = selfClass.Raw ( )
        self._Allocate ( rank )
# . Basic version - unfinished.
# . Possible combinations (first only done):
#
#   lower/bins/binSize
#   lower/upper/bins
#   lower/upper/binSize - need option to adjust lower and upper (making range either smaller or larger).
#
# . Also isPeriodic/period options.
#
        # . Initialization.
        for i from 0 <= i < len ( self ):
            # . Required data.
            bins    = data[i]["bins"   ]
            binSize = data[i]["binSize"]
            lower   = data[i]["lower"  ]
            self.cObject.dimensions[i].bins          = bins
            self.cObject.dimensions[i].binSize       = binSize
            self.cObject.dimensions[i].lower         = lower
            # . Derived data.
            self.cObject.dimensions[i].midPointLower = lower + 0.5 * binSize
            self.cObject.dimensions[i].upper         = lower + float ( bins ) * binSize
        # . Do the strides - reverse order.
        stride = 1
        for i from len ( self ) > i >= 0:
            self.cObject.dimensions[i].stride = stride
            stride *= data[i]["bins"]
        # . Finish up.
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
            log.SummaryOfItems ( [ ( "Rank", "{:d}".format ( self.rank ) ) ,
                                     ( "Size", "{:d}".format ( self.size ) ) ] ,
                                     title = "Regular Grid Summary" )

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    # . Properties.
    @property
    def rank ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.ndimensions
    @property
    def shape ( self ):
        cdef CInteger d
        shape = []
        if self.cObject != NULL:
            for d from 0 <= d < self.cObject.ndimensions:
                shape.append ( self.cObject.dimensions[d].bins )
        return shape
    @property
    def size ( self ):
        return RegularGrid_NumberOfGridPoints ( self.cObject )

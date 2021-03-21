"""Array iterators."""

from .AntisymmetricMatrix   import AntisymmetricMatrix
from .DoubleSymmetricMatrix import DoubleSymmetricMatrix
from .SymmetricMatrix       import SymmetricMatrix

# . Simple python implementation of array iterators - for debugging and printing.
# . 1-D, 2-D, N-D, antisymmetric, double symmetric and symmetric.

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def ArrayIterator ( array, conditionFunction = None, includeIndices = False ):
    """An iterator that returns array indices and values in default canonical order."""
    if   array.rank == 1:                             iteratorClass = _ArrayIterator1D
    elif isinstance ( array, AntisymmetricMatrix   ): iteratorClass = _ArrayIteratorAntisymmetric
    elif isinstance ( array, DoubleSymmetricMatrix ): iteratorClass = _ArrayIteratorDoubleSymmetric
    elif isinstance ( array, SymmetricMatrix       ): iteratorClass = _ArrayIteratorSymmetric
    else:                                             iteratorClass = _ArrayIteratorND
    return iteratorClass ( array, conditionFunction = conditionFunction, includeIndices = includeIndices )

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def _Next ( self ):
    """Next item."""
    indices = self._NextIndices ( )
    if indices is None: raise StopIteration
    else: return self.array.__getitem__ ( indices )

def _NextWithConditionFunction ( self ):
    """Next item."""
    while True:
        indices = self._NextIndices ( )
        if indices is None: raise StopIteration
        else:
            value = self.array.__getitem__ ( indices )
            if self.conditionFunction ( value, *indices ): return value

def _NextWithConditionFunctionAndIndices ( self ):
    """Next item."""
    while True:
        indices = self._NextIndices ( )
        if indices is None: raise StopIteration
        else:
            value = self.array.__getitem__ ( indices )
            if self.conditionFunction ( value, *indices ):
                if isinstance ( indices, int ): indices = ( indices, )
                return ( value, indices )

def _NextWithIndices ( self ):
    """Next item."""
    indices = self._NextIndices ( )
    if indices is None: raise StopIteration
    else:
        value = self.array.__getitem__ ( indices )
        if isinstance ( indices, int ): indices = ( indices, )
        return ( value, indices )

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class _BaseArrayIterator:

    def __init__ ( self, array, conditionFunction = None, includeIndices = False ):
        """Constructor."""
        self.array             = array
        self.conditionFunction = conditionFunction
        self.includeIndices    = includeIndices
        self._SetNext ( )
        self.Reset    ( )

    def __iter__ ( self ): return self

    def __len__ ( self ): return self.array.size

    def __next__ ( self ): return self._Next ( self )

    def _NextIndices ( self ):
        """Next set of indices."""
        return None

    def _SetNext ( self ):
        """Set the next method."""
        if self.conditionFunction is None:
            if self.includeIndices: self._Next = _NextWithIndices
            else:                   self._Next = _Next
        else:
            if self.includeIndices: self._Next = _NextWithConditionFunctionAndIndices
            else:                   self._Next = _NextWithConditionFunction

    def Reset ( self ):
        """Reset the iterator."""
        pass

    @property
    def size ( self ): return self.array.size

#-----------------------------------------------------------------------------------------------------------------------------------
class _ArrayIterator1D ( _BaseArrayIterator ):
    """1-D iterator."""

    def _NextIndices ( self ):
        """Next set of indices."""
        next = self.count
        if next < self.extent:
            self.count += 1
            return next
        else:
            return None

    def Reset ( self ):
        self.count  = 0
        self.extent = self.array.size

#-----------------------------------------------------------------------------------------------------------------------------------
class _ArrayIteratorAntisymmetric ( _BaseArrayIterator ):
    """2-D antisymmetric iterator."""

    def _NextIndices ( self ):
        """Next set of indices."""
        if self.counters is None:
            return None
        else:
            indices = tuple ( self.counters )
            ( i, j ) = self.counters
            j += 1
            if j < i:
                self.counters[1] = j
            else:
                self.counters[1] = 0
                i = self.counters[0] + 1
                if i < self.extent: self.counters[0] = i
                else:               self.counters    = None
            return indices

    def Reset ( self ):
        self.extent = self.array.shape[0]
        if self.extent > 1: self.counters = [ 1, 0 ]
        else:               self.counters = None

#-----------------------------------------------------------------------------------------------------------------------------------
class _ArrayIteratorDoubleSymmetric ( _BaseArrayIterator ):
    """4-D double symmetric iterator."""

    def _NextIndices ( self ):
        """Next set of indices."""
        if self.counters is None:
            return None
        else:
            indices = tuple ( self.counters )
            ( i, j, k, l ) = self.counters
            if k == i: lUpper = j+1
            else:      lUpper = k+1
            limits = [ self.extent, i+1, i+1, lUpper ]
            for d in ( 3, 2, 1, 0 ):
                n = self.counters[d] + 1
                if n < limits[d]:
                    self.counters[d] = n
                    break
                else:
                    if d == 0: self.counters    = None
                    else:      self.counters[d] = 0
            return indices

    def Reset ( self ):
        self.extent = self.array.shape[0]
        if self.extent > 0: self.counters = [ 0, 0, 0, 0 ]
        else:               self.counters = None

#-----------------------------------------------------------------------------------------------------------------------------------
class _ArrayIteratorSymmetric ( _BaseArrayIterator ):
    """2-D symmetric iterator."""

    def _NextIndices ( self ):
        """Next set of indices."""
        if self.counters is None:
            return None
        else:
            indices  = tuple ( self.counters )
            ( i, j ) = self.counters
            j += 1
            if j <= i:
                self.counters[1] = j
            else:
                self.counters[1] = 0
                i = self.counters[0] + 1
                if i < self.extent: self.counters[0] = i
                else:               self.counters    = None
            return indices

    def Reset ( self ):
        self.extent = self.array.shape[0]
        if self.extent > 0: self.counters = [ 0, 0 ]
        else:               self.counters = None

#-----------------------------------------------------------------------------------------------------------------------------------
class _ArrayIteratorND ( _BaseArrayIterator ):
    """N-D iterator."""

    def _NextIndices ( self ):
        """Next set of indices."""
        if self.counters is None:
            return None
        else:
            indices = tuple ( self.counters )
            for d in range ( self.rank-1, -1, -1 ):
                n = self.counters[d] + 1
                if n < self.extents[d]:
                    self.counters[d] = n
                    break
                else:
                    if d == 0: self.counters    = None
                    else:      self.counters[d] = 0
            return indices

    def Reset ( self ):
        self.extents  = self.array.shape
        self.rank     = self.array.rank
        if self.size > 0: self.counters = [ 0 for d in range ( self.rank ) ]
        else:             self.counters = None

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

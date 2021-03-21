"""Handle transformations of dimension 3."""

from   pCore          import logFile              , \
                             LogFileActive        , \
                             RawObjectConstructor
from ..Arrays         import ArrayPrint2D
from  .Geometry3Error import Geometry3Error

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Transformation3:

    def __copy__ ( self ):
        """Copying."""
        cdef Transformation3 new
        new = self.__class__.Uninitialized ( )
        Transformation3_CopyTo ( self.cObject, new.cObject )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Transformation3_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "rotation" : self.rotation, "translation" : self.translation }

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._FromRotationTranslation ( state["rotation"], state["translation"] )

    def _Allocate ( self ):
        """Allocation."""
        cdef Matrix33 r
        cdef Vector3  t
        r = Matrix33.Null ( )
        t = Vector3.Null  ( )
        self._FromRotationTranslation ( r, t )

    def _FromRotationTranslation ( self, Matrix33 rotation not None, Vector3 translation not None ):
        """Make the object from a rotation and translation."""
        cdef CStatus  cStatus = CStatus_OK
        self.rotation    = rotation
        self.translation = translation
        self.cObject     = Transformation3_FromRotationTranslation ( rotation.cObject, translation.cObject, &cStatus )
        self.isOwner     = True
        if cStatus != CStatus_OK: raise Geometry3Error ( "Unable to allocate Transformation3." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject     = NULL
        self.isOwner     = False
        self.rotation    = None
        self.translation = None

    @classmethod
    def FromSymmetryOperationString ( selfClass, oString ):
        """Construct a transformation from a string encoding a crystallographic symmetry operation."""
        cdef Transformation3 self
        self = None
        if isinstance ( oString, str ):
            self = selfClass.Null ( )
            # . Remove spaces and enclosing parentheses.
            t = oString.replace ( " ", "" )
            l = list ( t )
            if l[ 0] == "(" : l.pop ( 0 )
            if l[-1] == ")" : l.pop (   )
            t = "".join ( l )
            # . Loop over the three different specifications.
            for ( i, s ) in enumerate ( t.split ( "," ) ):
                # . Ensure all numbers end in a decimal point.
                ns = [ s[0:1] ]
                for ( j, c ) in enumerate ( s[1:] ):
                    if s[j:j+1].isdigit ( ) and not ( c.isdigit ( ) or c == "." ): ns.append ( "." )
                    ns.append ( c.lower ( ) )
                if ns[-1].isdigit ( ): ns.append ( "." )
                s = "".join ( ns )
                # . Determine the transformation.
                x  = 0.0
                y  = 0.0
                z  = 0.0
                t  = eval ( s, {}, { "x" : x, "y" : y, "z" : z } )
                x  = 1.0
                rX = eval ( s, {}, { "x" : x, "y" : y, "z" : z } ) - t
                y  = 1.0
                rY = eval ( s, {}, { "x" : x, "y" : y, "z" : z } ) - t - rX
                z  = 1.0
                rZ = eval ( s, {}, { "x" : x, "y" : y, "z" : z } ) - t - rX - rY
                self.translation[i] = t
                self.rotation [i,0] = rX
                self.rotation [i,1] = rY
                self.rotation [i,2] = rZ
        return self

    @classmethod
    def Identity ( selfClass ):
        """A constructor for the identity transformation."""
        self = selfClass     (     )
        self.rotation.Set    ( 0.0 )
        self.translation.Set ( 0.0 )
        for i in range ( 3 ):
            self.rotation[i,i] = 1.0
        return self

    def IsIdentity ( self ):
        """Is this transformation the identity?"""
        return ( Transformation3_IsIdentity ( self.cObject ) == CTrue )

    @classmethod
    def Null ( selfClass ):
        """Constructor."""
        cdef Transformation3 self
        self = selfClass.Uninitialized ( )
        self.rotation.Set    ( 0.0 )
        self.translation.Set ( 0.0 )
        return self

    def Orthogonalize ( self, Matrix33 A, Matrix33 B ):
        """Orthogonalization."""
        Transformation3_Orthogonalize ( self.cObject, A.cObject, B.cObject )

    def Print ( self, log = logFile, title = None ):
        """Printing."""
        if LogFileActive ( log ):
            if title is not None: log.Heading ( title, includeBlankLine = True )
            ArrayPrint2D ( self.rotation   , itemFormat = "{:.5f}", log = log, title = "Rotation"    )
            ArrayPrint2D ( self.translation, itemFormat = "{:.5f}", log = log, title = "Translation" )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

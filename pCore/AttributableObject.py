"""Picklable objects that have attributes, some of which can be set directly as options."""

# . _attributable is a dictionary of attribute/initializer name/values.
# . _unpicklable is a set of attributes that cannot be pickled.
# . _unsetable is a set of attributes that cannot be set directly.

from .Serialization import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AttributableObject:
    """An object with predefined attributes and defaults."""

    _attributable = {}
    _unpicklable  = set ( ) # . Worth it?
    _unsetable    = set ( ) # . Worth it?

    def __getstate__ ( self ):
        """Get state."""
        state = {}
        for name in self.__class__._attributable:
            if name not in self.__class__._unpicklable:
                value = getattr ( self, name, None )
                if value is not None: state[name] = value
        return state

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

#
# . Use the following instead and remove WithDefaults and WithOptions?
#    def __init__ ( self, **options ):
#        """Constructor with options."""
#        self._Initialize ( )
#        self.SetOptions ( **options )
#

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set state."""
        for ( name, value ) in state.items ( ):
            if ( name in self.__class__._attributable ) and ( value is not None ): setattr ( self, name, value )

    def _CheckOptions ( self ):
        """Check options."""
        # . Raise ValueError if arguments have illegal values.
        # . Rederive non-optionable attributes.
        pass

    def _Initialize ( self ):
        """Initialization."""
        for ( name, value ) in self.__class__._attributable.items ( ):
            if callable ( value ): value = value ( )
            setattr ( self, name, value )

    def _SetOptions ( self, **options ):
        """Set options."""
        unknown   = set ( )
        unsetable = set ( )
        for ( name, value ) in options.items ( ):
            if name in self.__class__._attributable:
                if name in self.__class__._unsetable:
                    unsetable.add ( name )
                else:
                    if callable ( value ): value = value ( )
                    setattr ( self, name, value )
            else:
                unknown.add ( name )
        if len ( unknown   ) > 0: raise TypeError ( "Unknown options for \"{:s}\": {:s}.".format      ( self.__class__.__name__, ", ".join ( sorted ( unknown   ) ) ) )
        if len ( unsetable ) > 0: raise TypeError ( "Unsetable attributes for \"{:s}\": {:s}.".format ( self.__class__.__name__, ", ".join ( sorted ( unsetable ) ) ) )
       
    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **options ):
        """Set and check options."""
        self._SetOptions ( **options )
        self._CheckOptions ( )

    # . Redundant?
    @classmethod
    def WithDefaults ( selfClass ):
        """Constructor with defaults."""
        self = selfClass ( )
        self._CheckOptions ( )
        return self

    # . Redundant?
    @classmethod
    def WithOptions ( selfClass, **options ):
        """Constructor with options."""
        self = selfClass ( )
        self.SetOptions ( **options )
        return self

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

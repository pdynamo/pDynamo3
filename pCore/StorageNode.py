"""A basic storage object."""

#===================================================================================================================================
# . A simple storage object.
#===================================================================================================================================
class StorageNode:
    """An object for storing attributes."""

    def Clear ( self ):
        """Clear the node."""
        self.__dict__.clear ( )

    def Delete ( self, attribute ):
        """Delete an attribute."""
        self.__dict__.pop ( attribute, None )

    def Get ( self, attribute, *arguments ):
        """Get the value of an attribute."""
        value = self.__dict__.get ( attribute, None )
        if ( value is None ) and ( len ( arguments ) > 0 ): value = arguments[0]
        return value

    def GetSet ( self, attribute, default, *arguments ):
        """Get the value of an attribute and set it if it does not exist."""
        value = self.__dict__.get ( attribute, None )
        if value is None:
            if callable ( default ): value = default ( *arguments )
            else:                    value = default
            self.__dict__[attribute] = value
        return value

    def GetSetNode ( self, attribute ):
        """Get the value of an attribute and set it as a node if it does not exist."""
        return self.GetSet ( attribute, self.__class__ )

    def Pop ( self, attribute ):
        """Pop an attribute."""
        return self.__dict__.pop ( attribute, None )

    def Set ( self, attribute, value ):
        """Set the value of an attribute."""
        self.__dict__[attribute] = value

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

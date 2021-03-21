"""Base class for MM parameter containers."""

from pCore import AttributableObject

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The character for parameter keys.
_ParameterKeyFieldSeparator = ":"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMParameterContainer ( AttributableObject ):
    """Base class for MM parameter containers."""

    _attributable     = dict ( AttributableObject._attributable )
    _strictAssignment = True
    _termLabel        = "MM Term"
    _attributable.update ( { "items"               : None                        ,
                             "keySeparator"        : _ParameterKeyFieldSeparator ,
                             "label"               : None                        ,
                             "parameterFactory"    : None                        ,
                             "properties"          : None                        ,
                             "rawItems"            : None                        ,
                             "termLabel"           : None                        ,
                             "useStrictAssignment" : True                        } )

    def _Initialize ( self ):
        """Initialization."""
        super ( MMParameterContainer, self )._Initialize ( )
        self.termLabel           = self.__class__._termLabel
        self.useStrictAssignment = self.__class__._strictAssignment

    @staticmethod
    def MakeKey ( *arguments ):
        """Make a key."""
        return tuple ( arguments )

    def MakeMMTerms ( self, atomTypeLabels, termIndices, connectivity ):
        """Make the appropriate MM terms given a list of term indices."""
        missingParameters = set ( )
        mmTerms           = []
        return ( mmTerms, missingParameters )

    def MakeMMTermsFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make the appropriate MM terms given a connectivity."""
        termIndices = connectivity.GetTermIndices ( ) # . Doesn't exist.
        return self.MakeMMTerms ( atomTypeLabels, termIndices, connectivity )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        pass

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""Attributable objects that can be summarized."""

from  enum               import Enum
from .AttributableObject import AttributableObject
from .LogFileWriter      import logFile            , \
                                LogFileActive

# . _summarizable is a dictionary of attribute/summary label or attribute/(summary label, summary format) name/values.

# . Need a SummarizableMixin for cases where need to add Summary to AttributableObject (e.g. Sequence).

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SummarizableObject ( AttributableObject ):
    """An attributable object that can be summarized."""

    _attributable = dict ( AttributableObject._attributable )
    _classLabel   = "Summarizable Object"
    _summarizable = {}
    _unsetable    = set  ( AttributableObject._unsetable    )
    _withSections = False
    _attributable.update ( { "label" : None } )

    def __str__ ( self ): return self.__class__._classLabel # . __repr__ is a formal representation, __str__ is an informal one.

    @staticmethod
    def _GetString ( value ):
        """Get a string representation of a value."""
        if   isinstance ( value, float ): string = "{:g}".format ( value )
        elif isinstance ( value, str   ): string = value
        elif isinstance ( value, Enum  ): string = value.name
        else:                             string = str ( value )
        return string

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            title = "Summary of {:s}".format ( self.__class__._classLabel )
            if self.label is not None: title += " \"{:s}\"".format ( self.label )
            log.SummaryOfItems ( self.SummaryItems ( ), title = title, withSections = self.__class__._withSections )

    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( self.__class__._classLabel, True ) ]
        for ( name, label ) in self.__class__._summarizable.items ( ):
            value = getattr ( self, name, None )
            if value is not None: # . None is ignored.
                if label is None:
                    items.extend ( value.SummaryItems ( ) ) # . A None label implies another object.
                else:
                    if isinstance ( label, str ):           # . Automatic format.
                        string = self._GetString ( value )
                    else:                                   # . Format given.
                        ( label, format ) = label
                        string = format.format   ( value )
                    items.append ( ( label, string ) )
        return items

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

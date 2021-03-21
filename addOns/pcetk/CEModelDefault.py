from  pCore   import logFile, LogFileActive
from .CEModel import CEModel

class CEModelDefault (CEModel):
    """A class to represent a built-in continuum electrostatic model."""
    _attributable = {
        }
    _attributable.update (CEModel._attributable)

    defaultAttributeNames = {
        }
    defaultAttributeNames.update (CEModel.defaultAttributeNames)


    def __init__ (self, system, customFiles=None, log=logFile, **options):
        """Constructor."""
        super (CEModelDefault, self).__init__ (system, customFiles=customFiles, log=log, **options)

    @property
    def label (self):
        return "Default"


    #-------------------------------------------------------------------------------
    def _CreateSite (self, **options):
        """Create a site and its instances specific to the built-in CE model."""
        pass


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

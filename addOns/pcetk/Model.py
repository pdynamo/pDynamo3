import os, glob

from  pCore       import logFile, LogFileActive
from .CEModelMEAD import CEModelMEAD

class MEADModel (CEModelMEAD):
    """A class to maintain backward compatibility."""

    # . This is troublesome as yaml files are often used as configuration files for other purposes.
    # . It is best not to use this class and to pass list of custom files directly to CEModelMEAD.
    def __init__ (self, system, log=logFile, **options):
        """Constructor."""
        workDir     = os.getcwd ()
        customFiles = []
        for extension in ("*.yaml", "*.YAML", "*.est", "*.EST"):
            customFiles.extend (glob.glob (os.path.join (workDir, extension)))

        # . Call the constructor from the parent class
        super (MEADModel, self).__init__ (system, customFiles=customFiles, log=log, **options)

#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

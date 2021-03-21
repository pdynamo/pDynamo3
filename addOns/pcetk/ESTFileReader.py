"""Reading EST parameter files."""

import os

from pCore import logFile, LogFileActive, TextFileReader

class ESTFileReader (TextFileReader):
    """A class for reading EST files."""

    _classLabel = "EST File Reader"

    def _CheckOptions ( self ):
        """Constructor."""
        super ( ESTFileReader, self )._CheckOptions ( )
        self.siteLabel = os.path.splitext (os.path.basename (name))[0]

    def Parse (self, log=logFile):
        """Parse the data on the file."""
        if not self.isParsed:
            if LogFileActive (log): self.log = log

            self.Open ()
            atoms  = []
            line   = None

            # In some files center atom is absent?
            center = None

            try:
                while True:
                    line   = self.GetLine (signalWarning = False)
                    tokens = line.split ()
                    if len (tokens) > 0:
                        if tokens[0] == "Gmodel" : Gmodels = map (float, tokens[1:])
                        if tokens[0] == "proton" : protons = map (int, tokens[1:])
                        if tokens[0] == "label"  : labels  = tokens[1:]
                        if tokens[0] == "center" : center  = tokens[1]

                        if tokens[0] == self.siteLabel:
                            label   = tokens[1]
                            charges = map (float, tokens[2:])
                            atoms.append ((label, charges))
            except EOFError:
                pass
            self.WarningStop ()
            self.Close ()
            self.log      = None
            self.isParsed = True

            instances = []
            for instanceIndex, (instanceLabel, instanceGmodel, instanceProtons) in enumerate (zip (labels, Gmodels, protons)):
                instanceAtoms   = []
                instanceCharges = []
                for atomLabel, atomCharges in atoms:
                    instanceAtoms.append (atomLabel)
                    instanceCharges.append (atomCharges[instanceIndex])
                instances.append ({"label" : instanceLabel, "Gmodel" : instanceGmodel, "protons" : instanceProtons, "charges" : instanceCharges})

            self.siteAtoms     = instanceAtoms
            self.siteInstances = instances
            self.siteCenter    = center


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass

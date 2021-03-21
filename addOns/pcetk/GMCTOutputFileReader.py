"""Reading output files from GMCT."""

from pCore       import logFile        , \
                        LogFileActive  , \
                        TextFileReader
from pScientific import Constants      , \
                        Units

class GMCTOutputFileReader (TextFileReader):
    """A class for reading output files from GMCT."""

    _classLabel = "GMCT Output File Reader"

    def Parse (self, temperature=300., log=logFile):
        """Parse the data on the file."""
        if not self.isParsed:
            if LogFileActive (log):
                self.log = log
            self.Open ()

            convert       = -1.0 / ( Constants.Molar_Gas * Units.Energy_Kilojoules_Per_Mole_To_Kilocalories_Per_Mole * Constants.Ln10 * temperature )
            line          = None
            pHtable       = []
            probabilities = {}
            try:
                while True:
                    if not line:
                        line = self.GetLine (signalWarnings=False)
                    if line.startswith ("chemical potential"):
                        tokens = line.split ()
                        mu     = float (tokens[2])
                        pHtable.append (mu * convert)

                        for i in range (2):
                            self.GetLine ()
                        while True:
                            line = self.GetLine (signalWarnings=False)
                            if line.startswith ("chemical potential"):
                                break
                            tokens = line.split ()
                            label, probability, mu, protons, vlabel, vmemb = tokens

                            if label in probabilities:
                                entries = probabilities[label]
                            else:
                                entries = []
                            entries.append (float (probability))
                            probabilities[label] = entries
            except EOFError:
                pass
            self.log           = None
            self.isParsed      = True
            self.pHtable       = pHtable
            self.probabilities = probabilities

            self.WarningStop ()
            self.Close ()


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

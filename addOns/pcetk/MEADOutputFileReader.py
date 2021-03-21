"""Reading output files from MEAD."""

from pCore import logFile, LogFileActive, TextFileReader

class MEADOutputFileReader (TextFileReader):
    """A class for reading output files from my_2diel_solver and my_3diel_solver."""

    _classLabel = "MEAD Output File Reader"

    def Parse (self, log=logFile):
        """Parse the data on the file."""
        if not self.isParsed:
            if LogFileActive (log):
                self.log = log

            self.interactions = []
            self.Open ()

            try:
                while True:
                    line = self.GetLine (signalWarnings = False)

                    if line.startswith ("Self energy of"):
                        tokens    = line.split ()
                        self.born = float (tokens[-1])

                    elif line.startswith ("Interaction energy of"):
                        tokens    = line.split ()
                        self.back = float (tokens[-1])

                    elif line.startswith ("Interaction energies of"):
                        self.GetLine ()
                        self.GetLine ()

                        while True:
                            # tokens = self.GetTokens (converters = [int, int, None, float])
                            line   = self.GetLine ()
                            if line.startswith ("Total runtime"):
                                break
                            else:
                                tokens   = line.split ()
                                site     = int   (tokens[0])
                                instance = int   (tokens[1])
                                energy   = float (tokens[3])
                                self.interactions.append ([site, instance, energy])
            except EOFError:
                pass

            self.WarningStop ()
            self.Close ()
            self.log      = None
            self.isParsed = True


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass

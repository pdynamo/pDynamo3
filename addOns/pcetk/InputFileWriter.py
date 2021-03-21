from pCore  import TextFileWriter

class InputFileWriter (TextFileWriter):
    """A class for writing text files."""

    _classLabel = "Input File Writer"

    def Write (self, listOfLines):
        """Write a list of lines."""
        self.Open  ()
        self.file.writelines (listOfLines)
        self.Close ()

#===============================================================================
# . Helper functions
#===============================================================================
def WriteInputFile (filename, listOfLines, addLineBreaks=False):
    """Write a file."""
    output = InputFileWriter.FromPath (filename)
    if addLineBreaks:
        templ = []
        for line in listOfLines:        
            templ.append ( "{:s}\n".format ( line ) )
        listOfLines = templ
    output.Write (listOfLines)

#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

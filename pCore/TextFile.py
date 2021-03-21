#===================================================================================================================================
# . Classes and functions to handle text files.
#===================================================================================================================================

import os, os.path, string

from .LogFileWriter      import logFile       , \
                                LogFileActive
from .PrintObjects       import Align
from .SummarizableObject import SummarizableObject

#===================================================================================================================================
# . Error classes.
#===================================================================================================================================
class TextFileReaderError ( Exception ):
    """Text file reader errors."""
    pass

#===================================================================================================================================
# . Error classes.
#===================================================================================================================================
class TextFileWriterError ( Exception ):
    """Text file writer errors."""
    pass

#===================================================================================================================================
# . Text file class.
#===================================================================================================================================
class TextFile ( SummarizableObject ):
    """TextFile is the base class for text files."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Text File"
    _attributable.update ( { "file"     : None  ,
                             "isActive" : False ,
                             "mode"     : None  ,
                             "path"     : None  } )

    def _CheckOptions ( self ):
        """Check options."""
        if self.mode  is None: raise IOError ( "Text file with undefined mode." )
        if self.path  is None: raise IOError ( "Text file with undefined path." )
        if self.label is None: self.label = os.path.split ( self.path )[-1]

    def Close ( self ):
        """Close the file."""
        if self.isActive:
            try:    self.file.close ( )
            except: raise IOError ( "Cannot close file: " + self.path + "." )
            self.isActive = False

    @classmethod
    def FromPath ( selfClass, path, **options ):
        """Constructor from path."""
        options         = dict ( options )
        options["path"] = path
        return selfClass.WithOptions ( **options )

    def Open ( self ):
        """Open the file."""
        if not self.isActive:
            try:    self.file = open ( self.path, self.mode )
            except: raise IOError ( "Cannot open file: " + self.path + "." )
            self.isActive = True

#===================================================================================================================================
# . Text file reader class.
#===================================================================================================================================
class TextFileReader ( TextFile ):
    """TextFileReader is the base class for text files that are to be read."""

    _attributable = dict ( TextFile._attributable )
    _classLabel   = "Text File Reader"
    _attributable.update ( { "fatalErrors"     : 0     ,
                             "isParsed"        : False ,
                             "linesParsed"     : 0     ,
                             "log"             : None  ,
                             "maximumWarnings" : 100   ,
                             "mode"            : "r"   ,
                             "warnings"        : 0     ,
                             "warningTable"    : None  } )

    def _CheckOptions ( self ):
        """The file must exist and be readable."""
        super ( TextFileReader, self )._CheckOptions ( )
        if not os.access ( self.path, os.R_OK ): raise IOError ( "Unreadable file: \"{:s}\".".format ( str ( self.path ) ) )

    def GetFixedFormatArray ( self, numberOfItems, itemsPerLine, itemWidth, converter = None, default = None, signalWarnings = True ):
        """Parse a fixed format array that may span several lines."""
        items = None
        if ( numberOfItems > 0 ) and ( itemsPerLine > 0 ) and ( itemWidth > 0 ):
            items = []
            while numberOfItems > 0:
                itemsOnLine = min ( numberOfItems, itemsPerLine )
                line        = self.GetFixedFormatLine ( itemsOnLine * itemWidth, signalWarnings = signalWarnings )
                for i in range ( itemsOnLine ):
                    word = line[i*itemWidth:(i+1)*itemWidth].strip ( )
                    if converter is None:
                        items.append ( word )
                    else:
                        item = default
                        if len ( word ) > 0:
                            try:    item = converter ( word )
                            except: self.Warning ( "Unable to convert token " + repr ( i ) + ".", True )
                        items.append ( item )
                numberOfItems -= itemsOnLine
        return items

    def GetFixedFormatLine ( self, lineLength, signalWarnings = True ):
        """Get a fixed format line with (at least) the given length."""
        try:
            line = next ( self.file ).rstrip ( ).ljust ( lineLength )
            self.linesParsed += 1
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def GetFixedFormatTokens ( self, *arguments, **options ):
        """Get tokens in fixed format from a line."""
        signalWarnings = options.get ( "signalWarnings", True )
        tokens = []
        length = 0
        for arg in arguments: length = max ( length, arg[1] )
        line   = self.GetFixedFormatLine ( length, signalWarnings = signalWarnings )
        for ( i, ( start, stop, converter, default ) ) in enumerate ( arguments ):
            word  = line[start:stop].strip ( )
            token = default
            if len ( word ) > 0:
                if converter is None:
                    token = word
                else:
                    try:    token = converter ( word )
                    except: self.Warning ( "Unable to convert token " + repr ( i ) + ".", True )
            tokens.append ( token )
        return tokens

    def GetLine ( self, signalWarnings = True ):
        """Get a line."""
        try:
            line = next ( self.file ).strip ( )
            self.linesParsed += 1
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def GetTokens ( self, converters = None, separator = None, signalWarnings = True ):
        """Get and convert tokens on a line."""
        line = self.GetLine ( signalWarnings = signalWarnings )
        return self.TokenizeLine ( line, converters = converters, separator = separator )

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Parse all the lines.
            try:
                while True:
                    line = self.GetLine ( )
                    if line is None: break
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.isParsed = True
            self.log      = None

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( TextFileReader, self ).SummaryItems ( )
        if self.isParsed:
            items.extend ( [ ( "Lines Parsed", "{:d}".format ( self.linesParsed ) ) ,
                             ( "Fatal Errors", "{:d}".format ( self.fatalErrors ) ) ,
                             ( "Warnings",     "{:d}".format ( self.warnings - self.fatalErrors ) ) ] )
        return items

    def TokenizeLine ( self, line, converters = None, separator = None ):
        """Tokenize a line with optional converters and separator."""
        tokens = None
        if line is not None:
            if separator is None: tokens = line.split ( )
            else:                 tokens = line.split ( separator )
            if converters is not None:
                for ( i, ( token, converter ) ) in enumerate ( zip ( tokens, converters ) ):
                    if converter is None:
                        new = token
                    else:
                        try:
                            new = converter ( token )
                        except:
                            new = converter ( )
                            self.Warning ( "Unable to convert token " + repr ( i ) + ".", True )
                    tokens[i] = new
        return tokens

    def Warning ( self, message, isFatal ):
        """Print a warning."""
        if ( self.log is not None ) and ( self.maximumWarnings > 0 ):
            if self.warnings == 0: self.WarningStart ( )
            self.warnings += 1
            if self.warnings <= self.maximumWarnings:
                self.warningTable.Entry ( "{:d}".format ( self.linesParsed ) )
                self.warningTable.Entry ( "" )
                self.warningTable.Entry ( message, align = Align.Left )
                if isFatal:
                    self.warningTable.Entry ( "Fatal" )
                    self.fatalErrors += 1
                else:
                    self.warningTable.Entry ( "Warning" )

    def WarningStart ( self ):
        """Start warning printing."""
        if ( self.log is not None ):
            self.warningTable = self.log.GetTable ( columns = [ 10, 5, 80, 10 ] )
            self.warningTable.Start ( )
            self.warningTable.Title ( "Text File Reader Warnings" )
            self.warningTable.Heading ( "Line"     )
            self.warningTable.Heading ( ""         )
            self.warningTable.Heading ( "Message"  )
            self.warningTable.Heading ( "Severity" )

    def WarningStop ( self ):
        """Terminate warning printing."""
        if ( self.warningTable is not None ):
            self.warningTable.Stop ( )
            self.warningTable = None
            if self.fatalErrors > 0:
                self.log.Paragraph ( "There have been fatal errors!" )
                raise TextFileReaderError ( "Fatal errors reading text file: " + self.path + "." )
            else:
                self.log.Paragraph ( "There have been warnings. Proceed with caution!" )

#===================================================================================================================================
# . Text file writer class.
#===================================================================================================================================
class TextFileWriter ( TextFile ):
    """TextFileWriter is the base class for text files that are to be written."""

    _attributable = dict ( TextFile._attributable )
    _classLabel   = "Text File Writer"
    _attributable.update ( { "mode" : "w" } )

    def _CheckOptions ( self ):
        """The file must be writable."""
        super ( TextFileWriter, self )._CheckOptions ( )
        ( head, tail ) = os.path.split ( self.path )
        if len ( head ) == 0: head = "."
        if not os.access ( head, os.W_OK ): raise IOError ( "Unwritable file: \"{:s}\".".format ( str ( self.path ) ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

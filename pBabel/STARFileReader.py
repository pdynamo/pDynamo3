"""Read data from a STAR file."""

# . Current restrictions:
#
# - Only data_ and loop_ control words are handled
# - No nested loops.
# - No save frames.

import os.path

from pCore import logFile, LogFileActive, TextFileReader

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Comment character.
_CommentCharacter = "#"

# . Continuation character (first character on line only).
_ContinuationCharacter = ";"

# . Data token.
_DataToken = "data_"

# . Empty token.
_EmptyToken = ""

# . Key character.
_KeyCharacter = "_"

# . Loop token.
_LoopToken = "loop_"

# . New line.
_NewLine = "\n"

# . Quotes characters.
_DoubleQuotes = "\""
_SingleQuotes = "'"

# . Undefined character.
_UndefinedCharacter = "."

# . Unknown character.
_UnknownCharacter = "?"

#===================================================================================================================================
# . Token classes.
#===================================================================================================================================
class STARFileDataBlockToken:
    """A data block token."""
    def __init__ ( self, string ):
        self.label = string[ len ( _DataToken ):].lower ( )

class STARFileDataNameToken:
    """A data name token."""
    def __init__ ( self, string ):
        self.label = string[1:].lower ( )

class STARFileDataValueToken:
    """A data value token."""
    def __init__ ( self, string ):
        self.value = string

class STARFileLoopToken:
    """A loop token."""
    def __init__ ( self, string ): pass

#===================================================================================================================================
# . Table class for loop data.
#===================================================================================================================================
class STARFileTable:
    """A class for storing loop data from a STAR file as a table."""

    def __init__ ( self, columnNames ):
        """Constructor."""
        self.columnNames = columnNames
        self.rows        = []
        label = os.path.commonprefix ( columnNames )
        if ( len ( label ) > 0 ) and ( not label[-1].isalnum ( ) ): label = label[:-1] # . Remove trailing non-alphanumeric characters.
        self.label = label

    def AppendRowValue ( self, value ):
        """Append a value to a row."""
        if self.IsComplete ( ): self.rows.append ( [] )
        self.rows[-1].append ( value )

    def ColumnValues ( self, columnName ):
        """Return the values in a column."""
        values = None
        if self.HasColumn ( columnName ):
            index  = self.columnNames.index ( columnName )
            values = [ self.rows[i][index] for i in range ( self.numberOfRows ) ]
        return values

    def HasColumn ( self, columnName ):
        """Has the table a column with this name?"""
        return ( columnName in self.columnNames )

    def IsComplete ( self ):
        """Is the table complete."""
        return ( len ( self.rows ) == 0 ) or ( len ( self.rows[-1] ) == len ( self.columnNames ) )

    def RealColumnValues ( self, columnName ):
        """Return the real values in a column."""
        values = self.ColumnValues ( columnName )
        if values is not None:
            for ( i, v ) in enumerate ( values ):
                try:
                    if v.endswith ( ")" ): v = v.split ( "(", 1 )[0]
                    f = float ( v )
                except:
                    f = None
                values[i] = f
        return values

    @property
    def numberOfRows ( self ):
        """The number of rows in the table."""
        return len ( self.rows )

#===================================================================================================================================
# . File reader class.
#===================================================================================================================================
class STARFileReader ( TextFileReader ):
    """STARFileReader is the class for reading STAR files."""

    _attributable = dict ( TextFileReader._attributable )
    _classLabel   = "STAR File Reader"
    _attributable.update ( { "AppendToken"      : None ,
                             "currentDataBlock" : None ,
                             "currentDataKey"   : None ,
                             "dataBlocks"       : None ,
                             "loopKeys"         : None ,
                             "loopTable"        : None } )
    _withSections = True

    def AppendControl ( self, token ):
        """Append a control token."""
        if   isinstance ( token, STARFileDataNameToken  ):
            self.AppendToken      = self.AppendDataValue
            self.currentDataKey   = token.label
        elif isinstance ( token, STARFileDataBlockToken ):
            if token.label in self.dataBlocks: self.Warning ( "Duplicate data block labels: " + token.label + ".", True )
            self.AppendToken      = self.AppendControl
            self.currentDataBlock = {}
            self.dataBlocks[token.label] = self.currentDataBlock
        elif isinstance ( token, STARFileLoopToken ):
            self.AppendToken      = self.AppendLoopHeader
            self.loopKeys         = []
        else:
            self.Warning ( "Control token missing.", False )

    def AppendDataValue ( self, token ):
        """Append a data value token."""
        if isinstance ( token, STARFileDataValueToken ):
            self.currentDataBlock[self.currentDataKey] = token.value
            self.currentDataKey = None
            self.AppendToken    = self.AppendControl
        else: self.Warning ( "Data value token missing.", False )

    def AppendInitialize ( self, token ):
        """Append the first token."""
        if isinstance ( token, STARFileDataBlockToken ):
           self.currentDataBlock = {}
           self.dataBlocks       = { token.label : self.currentDataBlock }
           self.AppendToken      = self.AppendControl
        else: self.Warning ( "File does not start with a data block token.", False )

    def AppendLoopBody ( self, token ):
        """Append a loop body."""
        if isinstance ( token, STARFileDataValueToken ):
            self.loopTable.AppendRowValue ( token.value )
            self.AppendToken = self.AppendLoopBody
        else:
            if not self.loopTable.IsComplete ( ): self.Warning ( "Incomplete loop table: " + self.loopTable.label + ".", True )
            self.loopTable = None
            self.AppendControl ( token )

    def AppendLoopHeader ( self, token ):
        """Append a loop header."""
        if isinstance ( token, STARFileDataNameToken ):
            self.loopKeys.append ( token.label )
        elif ( self.loopKeys is None ) or ( len ( self.loopKeys ) <= 0 ):
            self.Warning ( "Loop with no data column headers.", True )
            self.AppendToken = self.AppendControl
        else:
            table = STARFileTable ( self.loopKeys )
            if table.label in self.dataBlocks: self.Warning ( "Duplicate loop labels: " + table.label + ".", True )
            self.currentDataBlock[table.label] = table
            self.loopTable                     = table
            self.loopKeys                      = None
            self.AppendLoopBody ( token )

    def ExtractReal ( self, dataBlock, name ):
        """Extract the value corresponding to the name."""
        v = dataBlock.get ( name, None )
        if v is not None:
            try:
                if v.endswith ( ")" ): v = v.split ( "(", 1 )[0]
                v = float ( v )
            except:
                v = None
        return v

    def ExtractString ( self, dataBlock, name ):
        """Extract the value corresponding to the name."""
        return dataBlock.get ( name, None )

    def Finalize ( self ):
        """Finalize after tokenizing."""
        # . Data key still exists.
        if self.currentDataKey is not None:
            self.Warning ( "Incomplete data name/value pair: " + self.currentDataKey + ".", False )
        # . Table still exists.
        elif self.loopTable is not None:
            if not self.loopTable.IsComplete ( ): self.Warning ( "Incomplete loop table: " + self.loopTable.label + ".", True )

    def GetDataBlock ( self, blockCode ):
        """Get the data block with a given code."""
        dataBlock = None
        if len ( self.dataBlocks ) > 0:
            if blockCode is None:
                keys = list ( self.dataBlocks.keys ( ) )
                keys.sort ( )
                blockCode = keys[0]
            dataBlock = self.dataBlocks.get ( blockCode, None )
        if dataBlock is None: raise ValueError ( "Unrecognized block code: " + blockCode + "." )
        return ( blockCode, dataBlock )

    def GetLine ( self, signalWarnings = False ):
        """Get a line of non-zero length."""
        try:
            while True:
                line = next ( self.file ).strip ( )
                self.linesParsed += 1
                if len ( line ) > 0: break
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parse data from the file."""
        if not self.isParsed:
            if LogFileActive ( log ): self.log = log
            # . Start parsing.
            self.Open ( )
            # . Tokenization.
            self.Tokenize ( )
            # . Close the file and warning table.
            self.WarningStop ( )
            self.Close       ( )
            # . Everything is now parsed.
            self.log     = None
            self.isParsed = True

    def SummaryItems ( self ):
        """Summary items."""
        # . Header.
        items = [ (  None        , False                                     ) ,
                  ( "Data Blocks", "{:d}".format ( len ( self.dataBlocks ) ) ) ,
                  ( "Warnings"   , "{:d}".format ( self.warnings           ) ) ]
        # . Data blocks and models.
        for dataBlockKey in sorted ( self.dataBlocks.keys ( ) ):
            dataBlock = self.dataBlocks[dataBlockKey]
            d = t = 0
            for value in dataBlock.values ( ):
                if isinstance ( value, STARFileTable ): t += 1
                else:                                   d += 1
            items.extend ( [ ( "Data Block {:s}".format ( dataBlockKey ), False               ) ,
                             ( "Data Values"                            , "{:d}".format ( d ) ) ,
                             ( "Tables"                                 , "{:d}".format ( t ) ) ] )
        return items

    def TokenFactory ( self, string ):
        """Create a token given a string."""
        if   string.startswith ( _DataToken    ): token = STARFileDataBlockToken ( string )
        elif string ==  _LoopToken              : token = STARFileLoopToken      ( string )
        elif string.startswith ( _KeyCharacter ): token = STARFileDataNameToken  ( string )
        else                                    : token = STARFileDataValueToken ( string )
        return token

    def Tokenize ( self ):
        """Tokenize the file by line."""
        continuations      = []
        insideContinuation = False
        self.AppendToken   = self.AppendInitialize
        try:
            while True:
                line = self.GetLine ( )
                # . Inside a continuation.
                if insideContinuation:
                    if line.startswith ( _ContinuationCharacter ):
                        self.AppendToken ( STARFileDataValueToken ( _NewLine.join ( continuations ) ) )
                        continuations      = []
                        insideContinuation = False
                        if len ( line ) > 1: self.Warning ( "Invalid terminal continuation line.", False )
                    else:
                        continuations.append ( line )
                # . Start a continuation.
                elif line.startswith ( _ContinuationCharacter ):
                    insideContinuation = True
                    partialToken       = line[1:].strip ( )
                    if len ( partialToken ) > 0: continuations = [ partialToken ]
                    else:                        continuations = [ ]
                # . Normal line.
                else:
                    self.TokenizeLine ( line )
        except EOFError:
            pass
        self.Finalize ( )

    def TokenizeLine ( self, line ):
        """Tokenize a line."""
        characters     = []
        insideString   = False
        quoteCharacter = None
        for c in line.strip ( ):
            # . Inside a string.
            if insideString:
                if c == quoteCharacter:
                    self.AppendToken ( STARFileDataValueToken ( "".join ( characters ) ) )
                    characters     = []
                    insideString   = False
                    quoteCharacter = None
                else:
                    characters.append ( c )
            # . Start a string.
            elif ( c == _DoubleQuotes ) or ( c == _SingleQuotes ):
                if len ( characters ) > 0:
                    self.AppendToken ( self.TokenFactory ( "".join ( characters ) ) )
                    characters = []
                insideString   = True
                quoteCharacter = c
            # . A comment ends the current line.
            elif c == _CommentCharacter:
                if len ( characters ) > 0:
                    self.AppendToken ( self.TokenFactory ( "".join ( characters ) ) )
                    characters = []
                break
            # . White space ends the current token.
            elif c.isspace ( ):
                if len ( characters ) > 0:
                    self.AppendToken ( self.TokenFactory ( "".join ( characters ) ) )
                    characters = []
            # . Other characters.
            else:
                characters.append ( c )
        # . Error if a string is not closed.
        if insideString: self.Warning ( "Unmatched quotes (" + quoteCharacter + ") on line.", False )
        # . Append final token.
        if len ( characters ) > 0: self.AppendToken ( self.TokenFactory ( "".join ( characters ) ) )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

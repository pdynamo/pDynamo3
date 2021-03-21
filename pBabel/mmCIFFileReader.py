"""Read data from a mmCIF file."""

from  pCore                 import logFile, LogFileActive, TextFileReader
from  pMolecule             import System
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer
from .mmCIFModel            import mmCIFModel

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Continuation character.
_COMMENTCHARACTER = "#"

# . Continuation character.
_CONTINUATIONCHARACTER = ";"

# . Data token.
_DATATOKEN = "data_"

# . Default model number.
_DEFAULTMODELNUMBER = "1"

# . Key character.
_KEYCHARACTER = "_"

# . Loop token.
_LOOPTOKEN = "loop_"

# . Quotes characters.
_DOUBLEQUOTES = "\""
_SINGLEQUOTES = "'"

# . Undefined character.
_UNDEFINEDCHARACTER = "."

# . Unknown character.
_UNKNOWNCHARACTER = "?"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFFileReader ( TextFileReader ):
    """mmCIFFileReader is the class for reading mmCIF files."""

    _attributable = dict ( TextFileReader._attributable )
    _classLabel   = "mmCIF File Reader"
    _withSections = True
    _attributable.update ( { "nextline" : None } )

    def GenerateModels ( self, log = logFile ):
        """Generate mmCIF models from the parsed data."""
        if self.isParsed:
            self.mmcifmodels = {}
            for ( dataname, dataBlock ) in self.dataBlocks.items ( ):
                # . Get the various tables to be processed.
                atomsite = dataBlock.get ( "_atom_site", None )
                entity   = dataBlock.get ( "_entity",    None )
                # . Do nothing if there are no atoms.
                if ( atomsite is None ) or ( atomsite.numberOfRows <= 0 ): continue
                # . Define the basic model.
                model = mmCIFModel.WithOptions ( label = dataname, log = log )
                # . Add elements to the model (order is important).
                for ( elementname, attributename ) in ( ( "_entity", "EntityDefinition" ), ( "_struct_asym", "AsymmetricUnitDefinition" ), ( "_atom_site", "AtomSite" ), ( "_struct_conn", "Connection" ) ):
                    element  = dataBlock.get ( elementname, None )
                    function = getattr ( model, "Add" + attributename, None )
                    if ( element is not None ) and ( function is not None ):
                        for row in element.rows: function ( **row )
                # . Finish up.
                model.Finalize ( )
                self.mmcifmodels[dataname] = model

    def GetLine ( self, signalWarnings = False ):
        """Get a line of non-zero length."""
        try:
            if self.nextline is None:
                while True:
                    line = next ( self.file ).strip ( )
                    self.linesParsed += 1
                    if len ( line ) > 0: break
            else:
                line = self.nextline
                self.nextline = None
                self.linesParsed  += 1
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def GetModel ( self, dataBlockname = None ):
        """Get a model."""
        model = None
        if self.isParsed:
            mmcifmodels = getattr ( self, "mmcifmodels", {} )
            if dataBlockname is None:
                if   len ( mmcifmodels ) == 1: model = list ( mmcifmodels.values ( ) )[0]
                elif len ( mmcifmodels ) == 0: raise ValueError ( "There are no models in the data block." )
                else: raise ValueError ( "A data block name must be specified when there are multiple models." )
            else:
                model = mmcifmodels.get ( dataBlockname, None )
                if model is None: raise ValueError ( "Unable to find model with the name " + dataBlockname + "." )
        return model

    def GetTokens ( self, converters = None, separator = None, signalWarnings = False ):
        """Get the tokens on a line (and continuation lines if there are any)."""
        # . Initialization.
        tokens = None
        while True:
            # . First line.
            line   = self.GetLine ( )
            tokens = self.TokenizeLine ( line )
            # . Continuation lines.
            if len ( tokens ) > 0:
                line = self.GetLine ( )
                if line.startswith ( _CONTINUATIONCHARACTER ):
                    partialtoken = line[1:].strip ( )
                    if len ( partialtoken ) > 0: continuations = [ partialtoken ]
                    else:                        continuations = [ ]
                    while True:
                        line = self.GetLine ( )
                        if line == _CONTINUATIONCHARACTER:
                            if len ( continuations ) > 0: tokens.append ( " ".join ( continuations ) )
                            break
                        else:
                            partialtoken = line.strip ( )
                            if len ( partialtoken ) > 0: continuations.append ( partialtoken )
                else:
                    self.nextline = line
                    self.linesParsed  -= 1
                    break
        return tokens

    def Parse ( self, log = logFile ):
        """Parse data from the file."""
        if not self.isParsed:
            # . Initialization.
            self.dataBlocks = {}
            self.dataBlock  = {}
            if LogFileActive ( log ): self.log = log
            # . Start parsing.
            self.Open ( )
            try:
                # . Initialization.
                QLOOPHEADER = False
                QLOOPBODY   = False
                loopkeys    = None
                looptable   = None
                # . Start looping.
                while True:
                    tokens = self.GetTokens ( )
                    token  = tokens[0]
                    # . New data set.
                    if token.startswith ( _DATATOKEN ):
                        QLOOPHEADER = False
                        QLOOPBODY   = False
                        loopkeys    = None
                        looptable   = None
                        self.ProcessDataToken ( token )
                    # . Start loop.
                    elif token == _LOOPTOKEN:
                        QLOOPHEADER = True
                        QLOOPBODY   = False
                        loopkeys    = []
                        looptable   = None
                    # . In a loop header.
                    elif QLOOPHEADER:
                        if token.startswith ( _KEYCHARACTER ):
                            looptable = self.ProcessLoopHeader ( looptable, loopkeys, token )
                        else:
                            QLOOPHEADER = False
                            QLOOPBODY   = True
                            self.ProcessLoopBody ( looptable, loopkeys, tokens )
                    # . In a loop body.
                    elif QLOOPBODY:
                        if token.startswith ( _KEYCHARACTER ):
                            QLOOPHEADER = False
                            QLOOPBODY   = False
                            loopkeys    = []
                            looptable   = None
                            self.ProcessKeyValue ( token, tokens[1:] )
                        else:
                            self.ProcessLoopBody ( looptable, loopkeys, tokens )
                    # . Non-loop data.
                    elif token.startswith ( _KEYCHARACTER ):
                        self.ProcessKeyValue ( token, tokens[1:] )
                    # . Unrecognized tokens.
                    else:
                        self.Warning ( "Unrecognized data line.", False )
            except EOFError:
                pass
            # . Close the file and warning table.
            self.WarningStop ( )
            self.Close ( )
            # . Everything is now parsed.
            self.log     = None
            self.isParsed = True

    def ParseKey ( self, key ):
        """Parse a key for table and column names."""
        items = key.split ( ".", 1 )
        if len ( items ) == 1: return ( items[0].lower ( ), ""                 )
        else:                  return ( items[0].lower ( ), items[1].lower ( ) )

    @classmethod
    def PathToSystem ( selfClass                                 ,
                       path                                      ,
                       altLoc              = _UNDEFINEDCHARACTER ,
                       dataBlockname       = None                ,
                       log                 = logFile             ,
                       modelNumber         = _DEFAULTMODELNUMBER ,
                       embeddedHydrogens   = False               ,
                       useComponentLibrary = False               ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse          ( log = log )
        inFile.GenerateModels ( log = log )
        inFile.Summary        ( log = log )
        # . Get the model.
        model = inFile.GetModel ( dataBlockname = dataBlockname )
        # . Make the system.
        system = model.ToSystem ( altLoc = altLoc, modelNumber = modelNumber, embeddedHydrogens = embeddedHydrogens )
        # . Print out whether there are undefined coordinates.
        if LogFileActive ( log ) and ( system.coordinates3.numberUndefined > 0 ):
            # . Determine the types of atoms with undefined coordinates.
            nheavy    = 0
            nhydrogen = 0
            for i in system.coordinates3.undefined:
                if system.atoms[i].atomicNumber == 1: nhydrogen += 1
                else:                                 nheavy    += 1
            # . Output a summary.
            log.SummaryOfItems ( [ ( "Heavy Atoms", "{:d}".format ( nheavy    ) ) ,
                                   ( "Hydrogens",   "{:d}".format ( nhydrogen ) ) ] ,
                                   title = "Undefined Coordinates" )
        # . Finish up.
        return system

    def ProcessDataToken ( self, token ):
        """Process a data token."""
        # . Get the data token key.
        key = token[len(_DATATOKEN):].upper ( )
        if key in self.dataBlocks: self.Warning ( "Duplicate data block keys: " + key + ".", True )
        self.dataBlocks[key] = self.dataBlock = {}

    def ProcessKeyValue ( self, key, values ):
        """Process a unique key value pair."""
        ( tablename, columnname ) = self.ParseKey ( key )
        table = self.dataBlock.get ( tablename, None )
        if table is None:
            table = mmCIFFileTable ( )
            self.dataBlock[tablename] = table
        if table.HasColumn ( columnname ): self.Warning ( "Table in single-entry mode already has column " + columnname + ".", True )
        table.SetValue ( key, values )

    def ProcessLoopBody ( self, looptable, loopkeys, values ):
        """Process loop body."""
        if looptable is None: self.Warning ( "Loop defined without data definitions.", True )
        else: looptable.AddRow ( loopkeys, values )

    def ProcessLoopHeader ( self, looptable, loopkeys, key ):
        """Process loop header."""
        # . Treat the key.
        ( tablename, columnname ) = self.ParseKey ( key )
        # . Treat the tables.
        table = self.dataBlock.get ( tablename, None )
        if ( table is None ) and ( looptable is None ):
            table = mmCIFFileTable ( )
            self.dataBlock[tablename] = table
        elif table is looptable: pass
        else: self.Warning ( "Invalid loop header definition.", True )
        # . Treat the columns.
        if table is not None:
            if table.HasColumn ( columnname ): self.Warning ( "Loop header with duplicate columns " + key + ".", True )
            table.AddColumn    ( columnname )
        loopkeys.append ( columnname )
        return table

    def SummaryItems ( self ):
        """Summary items."""
        # . Header.
        items = [ ( None         , False                                     ) ,
                  ( "Data Blocks", "{:d}".format ( len ( self.dataBlocks ) ) ) ,
                  ( "Warnings"   , "{:d}".format ( self.warnings           ) ) ]
        # . Data blocks and models.
        mmcifmodels = getattr ( self, "mmcifmodels", {} )
        for key in sorted ( self.dataBlocks.keys ( ) ):
            items.append ( ( "Data Block " + key, True ) )
            items.extend ( [ ( tableName, "{:d}".format ( self.dataBlocks[key][tableName].numberOfRows ) ) for tableName in self.dataBlocks[key].keys ( ) ] )
            if key in mmcifmodels:
                model = mmcifmodels[key]
                if model.label is None: title = "Model"
                else:                   title = "Model {:s}".format ( model.label )
                items.append ( ( title, True ) )
                items.extend ( model.SummaryItems ( ) )
        return items

    def TokenizeLine ( self, line ):
        """Tokenize a line."""
        basictokens    = line.split ( )
        quotecharacter = None
        stringtokens   = []
        tokens         = []
        for basictoken in basictokens:
            # . Inside string.
            if len ( stringtokens ) > 0:
                # . End of string.
                if basictoken.endswith ( quotecharacter ):
                    stringtokens.append ( basictoken[0:-1].strip ( ) )
                    tokens.append ( " ".join ( stringtokens ) )
                    quotecharacter = None
                    stringtokens   = []
                # . Middle of string.
                else:
                    stringtokens.append ( basictoken )
            # . Outside string.
            # . Start (and possibly the end) of string.
            # . Make sure to check for end as well.
            elif basictoken.startswith ( _DOUBLEQUOTES ):
                if basictoken.endswith ( _DOUBLEQUOTES ):
                    token = basictoken[1:-1].strip ( )
                    if len ( token ) > 0: tokens.append ( token         )
                    else:                 tokens.append ( _UNKNOWNTOKEN )
                else:
                    quotecharacter = _DOUBLEQUOTES
                    stringtokens.append ( basictoken[1:].strip ( ) )
            elif basictoken.startswith ( _SINGLEQUOTES ):
                if basictoken.endswith ( _SINGLEQUOTES ):
                    token = basictoken[1:-1].strip ( )
                    if len ( token ) > 0: tokens.append ( token         )
                    else:                 tokens.append ( _UNKNOWNTOKEN )
                else:
                    quotecharacter = _SINGLEQUOTES
                    stringtokens.append ( basictoken[1:].strip ( ) )
            # . Comment.
            elif basictoken.find ( _COMMENTCHARACTER ) > -1:
                break
            else:
                tokens.append ( basictoken )
        # . Error if a string is not closed.
        if len ( stringtokens ) > 0:
            self.Warning ( "Unmatched quotes (" + quotecharacter + ") on line.", False )
            tokens.append ( " ".join ( stringtokens ) )
        return tokens

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class mmCIFFileTable:
    """A class for storing data from a mmCIF file."""

    def __init__ ( self ):
        """Constructor."""
        self.columnnames = set ( )
        self.rows        = []

    def AddColumn ( self, columnname ):
        """Add a column to the table."""
        self.columnnames.add ( columnname )

    def AddRow ( self, keys, values ):
        """Add a row to the table."""
        row = {}
        for ( key, value ) in zip ( keys, values ):
            if key in self.columnnames: row[key] = value
        self.rows.append ( row )

    def GetUniqueRowValues ( self, columnname ):
        """Get the unique row values for a column."""
        values = set ( )
        for row in self.rows:
            if columnname in row: values.add ( row[columnname] )
        return values

    def HasColumn ( self, columnname ):
        """Has the table a column with this name?"""
        return ( columnname in self.columnnames )

    def SetValue ( self, columnname, value ):
        """Assign a value to the table."""
        self.columnnames.add ( columnname )
        if len ( self.rows ) <= 0: self.rows.append ( {} )
        self.rows[-1][columnname] = value

    @property
    def numberOfRows ( self ):
        """The number of rows in the table."""
        return len ( self.rows )

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { System : mmCIFFileReader.PathToSystem } , [ "mmcif", "MMCIF" ], "Macromolecular Crystallographic Information File" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

"""Classes and functions for test data."""

import math, os

from .AttributableObject import AttributableObject
from .LogFileWriter      import logFile           , \
                                LogFileActive     , \
                                TextLogFileWriter
from .PrintObjects       import Align
from .Serialization      import Pickle            , \
                                Unpickle

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Output options.
_Indent = 4

# . TestReal options.
_TestRealFieldMappings = { "Absolute Error Tolerance" : "absoluteErrorTolerance" ,
                           "Label"                    : "label"                  ,
                           "Percent Error Tolerance"  : "percentErrorTolerance"  ,
                           "Tolerance Format"         : "toleranceFormat"        ,
                           "Units"                    : "units"                  ,
                           "Value"                    : "value"                  ,
                           "Value Format"             : "valueFormat"            }

#===================================================================================================================================
# . Exception class.
#===================================================================================================================================
class TestDataException ( Exception ):
    """Test errors."""
    pass

#===================================================================================================================================
# . Test data.
#===================================================================================================================================
class TestDatum ( AttributableObject ):
    """Base class for an item of test data."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "label"  : None ,
                             "parent" : None ,
                             "value"  : None } )

    def __len__ ( self ): return 1

    def DeviationAsString ( self, value ): return ""

    def ToleranceAsString ( self ): return ""

    def ValueAsString     ( self ): return ""

    def VerifyAgainst ( self, value ):
        """Verify the datum against an input value."""
        return False

#-----------------------------------------------------------------------------------------------------------------------------------
class TestReal ( TestDatum ):
    """A real item of test data."""

    _attributable = dict ( TestDatum._attributable )
    _attributable.update ( { "absoluteErrorTolerance" : None ,
                             "percentErrorTolerance"  : None ,
                             "toleranceFormat"        : None ,
                             "units"                  : None ,
                             "valueFormat"            : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( TestReal, self )._CheckOptions ( )
        # . Initialization.
        self.usePercentError = False
        # . Get the absoluteErrorTolerance.
        if self.percentErrorTolerance is not None:
            if self.value != 0.0:
                absoluteErrorTolerance = ( self.percentErrorTolerance * self.value ) / 100.0
                if ( self.absoluteErrorTolerance is None ) or ( absoluteErrorTolerance < self.absoluteErrorTolerance ):
                    self.absoluteErrorTolerance = absoluteErrorTolerance
                    self.usePercentError        = True
        # . Overall checks.
        isOK = isinstance ( self.value, float ) and ( self.absoluteErrorTolerance is not None )
        if not isOK: raise TestDataException ( "Invalid real datum input values." )
        # . Formats.
        if self.toleranceFormat is None: self.toleranceFormat = "{:.4e}"
        if self.valueFormat     is None: self.valueFormat     = "{:.4e}"

    def Deviation ( self, value ):
        """Get the absolute deviation."""
        return ( value - self.value )

    def DeviationAsString ( self, value ):
        """Get the deviation as a string."""
        deviation = self.Deviation ( value )
        if self.usePercentError:
            deviation *= ( 100.0 / self.value )
            string = "{:.1%}".format ( deviation )
        else:
            string = self.toleranceFormat.format ( deviation )
        return string

    def ToleranceAsString ( self ):
        """Get the tolerance as a string."""
        if self.usePercentError: string = "{:.1%}".format ( self.percentErrorTolerance )
        else:                    string = self.toleranceFormat.format ( self.absoluteErrorTolerance )
        return string

    def ValueAsString ( self, value = None ):
        """Get the value as a string."""
        if value is None: value = self.value
        return self.valueFormat.format ( value )

    def VerifyAgainst ( self, value, results ):
        """Verify the datum against an input value."""
        if math.fabs ( self.Deviation ( value ) ) <= self.absoluteErrorTolerance:
            results.AddSuccess ( self.label, self.parent, value )
        else:
            results.AddFailure ( self.label, self.parent, value )

#===================================================================================================================================
# . Test data results.
#===================================================================================================================================
class TestDataResult ( AttributableObject ):
    """An object to hold the results of data set verifications."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "dataSet"    : None ,
                             "failures"   : dict ,
                             "missing"    : dict ,
                             "successes"  : dict ,
                             "unverified" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        if self.dataSet is None: raise ValueError ( "Unassigned data set." )
        self.unverified = len ( self.dataSet )

    def AddFailure ( self, label, parent, value ):
        """Add a failure."""
        items = self.failures.get ( parent, {} )
        items[label] = value
        self.failures[parent] = items
        self.unverified -= 1

    def AddMissing ( self, label, parent ):
        """Add a missing result."""
        items = self.missing.get ( parent, [] )
        items.append ( label )
        self.missing[parent] = items

    def AddSuccess ( self, label, parent, value ):
        """Add a success."""
        items = self.successes.get ( parent, {} )
        items[label] = value
        self.successes[parent] = items
        self.unverified -= 1

    @staticmethod
    def Counter ( dictionary ):
        """Count up the number of items in the dictionary."""
        length = 0
        for items in dictionary.values ( ): length += len ( items )
        return length

    def Summary ( self, fullSummary = False, log = logFile ):
        """Write a summary of the results."""
        if LogFileActive ( log ):
            if fullSummary:
                log.SummaryOfItems ( self.SummaryItems ( ), title = "Verification Summary for " + self.dataSet.label )
                self.dataSet.ResultsSummary ( log = log, failures = self.failures, missing = self.missing, successes = self.successes )
            else:
                if self.WasSuccessful ( ):
                    log.Paragraph ( "The data set " + self.dataSet.label + " was successfully verified." )
                else:
                    self.dataSet.ResultsSummary ( log = log, failures = self.failures, missing = self.missing )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Failures"  , "{:d}".format ( self.numberFailures  ) ) ,
                 ( "Missing"   , "{:d}".format ( self.numberMissing   ) ) ,
                 ( "Successes" , "{:d}".format ( self.numberSuccesses ) ) ,
                 ( "Unverified", "{:d}".format ( self.unverified      ) ) ]

    def WasSuccessful ( self ):
        """Was the verification successful?"""
        return ( ( len ( self.failures ) == 0 ) and ( len ( self.missing ) == 0 ) )

    @property
    def numberFailures  ( self ): return TestDataResult.Counter ( self.failures  )
    @property
    def numberMissing   ( self ): return TestDataResult.Counter ( self.missing   )
    @property
    def numberSuccesses ( self ): return TestDataResult.Counter ( self.successes )

#===================================================================================================================================
# . Test data sets.
#===================================================================================================================================
class TestDataSet ( AttributableObject ):
    """A collection of test data items."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "children" : dict ,
                             "data"     : dict ,
                             "label"    : None ,
                             "parent"   : None } )

    # . It is assumed that all data are real for the moment.
    def __getstate__ ( self ):
        """Get state."""
        state = { "Label" : self.label }
        if len ( self.children ) > 0:
            state["Children"] = [ value.__getstate__ ( ) for ( key, value ) in self.children.items ( ) ]
        if len ( self.data ) > 0:
            names  = list ( sorted ( _TestRealFieldMappings.keys ( ) ) )
            keys   = [ _TestRealFieldMappings[name] for name in names ]
            data   = []
            for dLabel in sorted ( self.data.keys ( ) ):
                datum = self.data[dLabel]
                data.append ( [ "-" if getattr ( datum, key ) is None else getattr ( datum, key ) for key in keys ] )
            state["Real Data Fields"] = names
            state["Real Data"       ] = data
        return state

    def __len__ ( self ):
        """Length."""
        length = 0
        for item in self.children.values ( ): length += len ( item )
        for item in self.data.values     ( ): length += len ( item )
        return length

    def __setstate__ ( self, state ):
        """Set state."""
        self.label = state["Label"]
        if "Children" in state:
            for cState in state["Children"]:
                child = self.__class__.Raw ( )
                child.__setstate__ ( cState )
                child.parent = self
                self.children[child.label] = child
        if "Real Data Fields" in state:
            attributes = [ _TestRealFieldMappings[name] for name in state["Real Data Fields"] ]
        if "Real Data" in state:
            for values in state["Real Data"]:
                options = {}
                for ( key, value ) in zip ( attributes, values ):
                    if value == "-": options[key] = None
                    else:            options[key] = value
                options["parent"]           = self
                self.data[options["label"]] = TestReal.WithOptions ( **options )

    def AddDatum ( self, datum ):
        """Add a datum."""
        uniqueLabel = ( datum.label not in self.children ) and ( datum.label not in self.data )
        if   isinstance ( datum, TestDataSet ) and uniqueLabel: self.children[datum.label] = datum
        elif isinstance ( datum, TestDatum   ) and uniqueLabel: self.data    [datum.label] = datum
        else: raise TestDataException ( "Invalid datum added to data set." )

    def Path ( self ):
        """Return the path for the data set."""
        items  = [ self.label ]
        parent = self.parent
        while ( parent is not None ):
            items.append ( parent.label )
            parent = parent.parent
        items.reverse ( )
        return ":".join ( items )

    def ResultsSummary ( self, failures = {}, log = logFile, missing = {}, successes = None ):
        """Write a results summary of the data set."""
        if LogFileActive ( log ):
            if len ( self ) > 0:
                # . Initialization.
                addFailureColumn = False
                # . Gather data.
                localFailures = failures.get ( self, {} )
                items         = dict ( localFailures )
                if successes is not None:
                    localSuccesses = successes.get ( self, {} )
                    items.update ( localSuccesses )
                    addFailureColumn = ( len ( localSuccesses ) < len ( items ) )
                # . Do the summary of successes and failures.
                if len ( items ) > 0:
                    lLength = max ( [ len ( datum.label ) for datum in self.data.values ( ) ] ) + 1
                    columns = [ max ( 20, lLength ), 20, 20, 20 ]
                    if addFailureColumn: columns.append ( 8 )
                    table = log.GetTable ( columns = columns )
                    table.Start   ( )
                    table.Title   ( "Results for " + self.Path ( ) )
                    table.Heading ( "Label"     )
                    table.Heading ( "Observed"  )
                    table.Heading ( "Reference" )
                    table.Heading ( "Deviation" )
                    if addFailureColumn: table.Heading ( "Fail" )
                    keys = list ( items.keys ( ) )
                    keys.sort ( )
                    for key in keys:
                        datum = self.data[key]
                        value = items[key]
                        table.Entry ( datum.label, align = Align.Left )
                        table.Entry ( datum.ValueAsString     ( value ) )
                        table.Entry ( datum.ValueAsString     ( ) )
                        table.Entry ( datum.DeviationAsString ( value ) )
                        if addFailureColumn:
                            if key in localFailures: table.Entry ( "*", align = Align.Center )
                            else:                    table.Entry ( ""  )
                    table.Stop ( )
                # . Missing.
                localMissing = missing.get ( self, [] )
                if len ( localMissing ) > 0:
                    localMissing.sort ( )
                    maximumLabelLength = 0
                    for label in localMissing: maximumLabelLength = max ( maximumLabelLength, len ( label ) )
                    columnWidth   = max ( 20, maximumLabelLength + 2 )
                    numberColumns = min ( len ( localMissing ), 4 )
                    table = log.GetTable ( columns = numberColumns * [ columnWidth ] )
                    table.Start   ( )
                    table.Title   ( "Missing Items in " + self.Path ( ) )
                    for label in localMissing: table.Entry ( label )
                    table.Stop ( )
                # . Header.
#                else:
#                    log.Heading ( self.Path ( ), includeBlankLine = True )
                # . Children.
                if len ( self.children ) > 0:
                    keys = list ( self.children.keys ( ) )
                    keys.sort ( )
                    for key in keys:
                        self.children[key].ResultsSummary ( failures = failures, log = log, missing = missing, successes = successes )

    def Summary ( self, log = logFile ):
        """Write a summary of the data set."""
        if LogFileActive ( log ):
            if len ( self ) == 0:
                log.Paragraph ( "Data set " + self.label + " is empty." )
            else:
                # . Data.
                if len ( self.data ) > 0:
                    data   = [ ( datum.label, datum.ValueAsString ( ), datum.ToleranceAsString ( ) ) for datum in self.data.values ( ) ]
                    widths = [ 20, 20, 20 ]
                    for datum in data:
                        for ( n, item ) in enumerate ( datum ):
                            widths[n] = max ( widths[n], len ( item ) )
                    for n in range ( len ( datum ) ): widths[n] += 1
                    table = log.GetTable ( columns = widths )
                    table.Start   ( )
                    table.Title   ( self.Path ( ) )
                    table.Heading ( "Label"     )
                    table.Heading ( "Reference" )
                    table.Heading ( "Tolerance" )
                    for datum in sorted ( data ):
                        table.Entry ( datum[0] , align = Align.Left )
                        table.Entry ( datum[1] )
                        table.Entry ( datum[2] )
                    table.Stop ( )
                # . Header.
#                else:
#                    log.Heading ( self.Path ( ), includeBlankLine = True )
                # . Children.
                if len ( self.children ) > 0:
                    keys = list ( self.children.keys ( ) )
                    keys.sort ( )
                    for key in keys:
                        self.children[key].Summary ( log = log )

    def VerifyAgainst ( self, dataSet, results = None ):
        """Verify the data set against a data set input as a dictionary."""
        if results is None: results = TestDataResult.WithOptions ( dataSet = self )
        for ( label, value ) in dataSet.items ( ):
            if   label in self.children: results = self.children[label].VerifyAgainst ( value, results = results )
            elif label in self.data    : self.data[label].VerifyAgainst ( value, results )
            else: results.AddMissing ( label, self )
        return results

#===================================================================================================================================
# . Utility functions.
#===================================================================================================================================
def PrettyYAMLChild ( yFile, self, indent = 0 ):
    """Output a test set node as pretty YAML."""
    spaces  = ( _Indent * indent     ) * " "
    spacesL = ( _Indent * indent - 2 ) * " " + "- "
    if indent == 0:
        yFile.write ( "Label : {:s}\n".format ( self.label ) )
    else:
        yFile.write ( "{:s}Label : {:s}\n".format ( spacesL, self.label ) )
    if len ( self.children ) > 0:
        yFile.write ( "{:s}Children :\n".format ( spaces ) )
        for key in sorted ( self.children.keys ( ) ):
            child   = self.children[key]
            PrettyYAMLChild ( yFile, child, indent = indent+1 )
    if len ( self.data ) > 0:
        # . Gather all non-null attributes.
        nonNull = set ( )
        for ( name, attribute ) in _TestRealFieldMappings.items ( ):
            if ( attribute == "label" ) or ( attribute == "value" ): continue
            for datum in self.data.values ( ):
                if getattr ( datum, attribute ) is not None:
                    nonNull.add ( name )
                    break
        names = [ "Label", "Value" ] + list ( sorted ( nonNull ) )
        # . Gather data.
        data   = []
        widths = [ 0 for i in range ( len ( names ) ) ]
        for datum in self.data.values ( ):
            items = []
            for ( n, name ) in enumerate ( names ):
                attribute = _TestRealFieldMappings[name]
                value     = getattr ( datum, attribute )
                if value is None:
                    item = "-"
                elif attribute == "label":
                    item = value
                elif attribute == "value":
                    item = datum.valueFormat.format ( value )
                elif attribute == "absoluteErrorTolerance":
                    item = datum.toleranceFormat.format ( value )
                elif attribute == "percentErrorTolerance":
                    item = datum.toleranceFormat.format ( value )
                elif attribute.endswith ( "Format" ):
                    item = "\"{:s}\"".format ( value )
                else:
                    item = value
                widths[n] = max ( widths[n], len ( item ) )
                items.append ( item )
            data.append ( items )
        # . Write out fields.
        yFile.write ( "{:s}Real Data Fields :\n".format ( spaces ) )
        for name in names: yFile.write ( "{:s}  - {:s}\n".format ( spaces, name ) )
        # . Write out data.
        yFile.write ( "{:s}Real Data :\n".format ( spaces ) )
        for datum in sorted ( data ):
            items = [ item.ljust ( width ) for ( item, width ) in zip ( datum, widths ) ]
            yFile.write ( "{:s}  - [ {:s} ]\n".format ( spaces, " , ".join ( items ) ) )

def PrettyYAML ( path, title, self ):
    """Output a complete test set as pretty YAML."""
    yFile  = open ( path, "w" )
    yFile.write ( "# . {:s}\n".format ( title ) )
    yFile.write ( "---\n" )
    PrettyYAMLChild ( yFile, self, indent = 0 )
    yFile.write ( "...\n" )
    yFile.close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

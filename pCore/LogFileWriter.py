"""Classes for LogFileWriters."""

import os, os.path, sys, time

from .AttributableObject  import AttributableObject
from .PrintObjects        import TextHeading    , \
                                 TextParagraph  , \
                                 TextSummary    , \
                                 TextTable      , \
                                 Verbatim       , \
                                 XHTMLHeading   , \
                                 XHTMLParagraph , \
                                 XHTMLSummary   , \
                                 XHTMLTable     , \
                                 XHTMLVerbatim
from .Time                import CPUTime

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Page width.
_TextPageWidth = 100

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LogFileWriter ( AttributableObject ):
    """Base class for all LogFileWriters."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "bufferSize"       : 1    , # . bufferSize sets the buffering level: 0, unbuffered; 1, line buffered; other, use a buffer of that size.
                             "cpuTime"          : None ,
                             "file"             : None ,
                             "isActive"         : True ,
                             "path"             : None ,
                             "printObjectStack" : list } )

    # . Public methods.
    def __del__ ( self ):
        """Destructor."""
        self.Close ( )

    def _CheckOptions ( self ):
        """Check options."""
        if self.path is None:
            self.file = sys.stdout
        else:
            ( head, tail ) = os.path.split ( self.path )
            if len ( head ) == 0: head = "."
            if not os.access ( head, os.W_OK ): raise IOError ( "Unwriteable file: " + self.path + "." )
            try:    self.file = open ( self.path, "w", self.bufferSize )
            except: raise IOError ( "Cannot open file: " + self.path + "." )

    def Close ( self ):
        """Close the file."""
        self.PrintObjectsClose ( )
        if self.isActive:
            self.isActive = False
            if self.path is not None:
                try:
                    self.file.close ( )
                except:
                    raise IOError ( "Cannot close file: " + self.path + "." )

    # . Methods that return print objects.
    def GetHeading ( self, **options ):
        pass

    def GetObjectTree ( self, **options ):
        pass

    def GetParagraph ( self, **options ):
        pass

    def GetSummary ( self, **options ):
        pass

    def GetTable ( self, **options ):
        pass

    def GetVerbatim ( self, **options ):
        pass

    # . Print object methods.
    def PrintObjectsClose ( self ):
        """Close all print objects."""
        if self.isActive:
            for item in self.printObjectStack:
                item.Stop ( )

    def PrintObjectPop ( self ):
        """Remove a print object from the stack."""
        if len ( self.printObjectStack ) > 0: self.printObjectStack.pop ( )

    def PrintObjectPush ( self, item ):
        """Add a print object to the stack."""
        self.printObjectStack.append ( item )

    # . Output methods.
    def Footer ( self ):
        """Write a footer and close the file."""
        if self.cpuTime is not None: self.Paragraph ( "CPU Time: " + self.cpuTime.CurrentAsString ( ) )
        self.Paragraph ( "Stop Time: " + time.asctime ( time.localtime ( None ) ) )
        self.Close ( )

    def Header ( self, title = None ):
        """Write a header."""
        self.cpuTime = CPUTime ( )
        if title is not None: self.Heading   ( title )
        self.Paragraph ( "Start Time: " + time.asctime ( time.localtime ( None ) ) )

    def Heading ( self, text, **options ):
        """Heading."""
        item = self.GetHeading ( **options )
        item.Start ( )
        item.Text  ( text )
        item.Stop  ( )

    def LineBreak ( self ):
        """Line break."""
        pass

    def ObjectTree ( self, object, **options ):
        """Object tree."""
        item = self.GetObjectTree ( **options )
        item.Start ( )
        item.Entry ( object )
        item.Stop  ( )

    def Paragraph ( self, text, **options ):
        """Paragraph."""
        item = self.GetParagraph ( **options )
        item.Start ( )
        item.Text  ( text )
        item.Stop  ( )

    def Separator ( self ):
        """Separator."""
        pass

    def SummaryOfItems ( self, items, **options ):
        """Summary of items.

        Items consist of list of section headers ( label, boolean ) and entries ( label, string ).
        """
        if self.isActive:
            title   = options.pop ( "title", None )
            summary = self.GetSummary ( **options ) # . Include order and withSections.
            ( title, sections ) = summary.SetUp ( title, items )
            summary.Print ( title, sections )

#### . MJF
    def TableOfRecords ( self, records, **options ):
        """Table of records."""
        # . Headers must be present.
        if self.isActive and ( len ( records ) > 0 ):
            alignments = options.pop ( "alignments"        , None )
            columns    = options.pop ( "columns"           , None )
            headers    = options.pop ( "headers"           , None )
            padding    = options.pop ( "columnPadding"     ,    2 )
            title      = options.pop ( "title"             , "Table Of Records" )
            width      = options.pop ( "minimumColumnWidth",   20 )
            if columns is None:
                columns = []
                for c in range ( len ( records[0] ) ):
                    columns.append ( max ( width, max ( [ len ( record[c] ) + padding for record in records ] ) ) )
                options["columns"] = columns
            table = self.GetTable ( **options )
            table.Start ( )
            table.Title ( title )
            for header in headers: table.Heading ( header )
            for record in records:
                for ( alignment, item ) in zip ( alignments, record ):
                   table.Entry ( item, align = alignment )
            table.Stop ( )
#### . MJF

    def Text ( self, text ):
        """Text."""
        if self.isActive and ( text is not None ): self.file.write ( text )

    def Verbatim ( self, text, **options ):
        """Verbatim."""
        item = self.GetVerbatim ( **options )
        item.Start ( )
        item.Text  ( text )
        item.Stop  ( )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class TextLogFileWriter ( LogFileWriter ):
    """Text LogFileWriters."""

    def GetHeading ( self, **options ):
        return TextHeading ( self, **options )

    def GetObjectTree ( self, **options ):
        return TextObjectTree ( self, **options )

    def GetParagraph ( self, **options ):
        return TextParagraph ( self, **options )

    def GetSummary ( self, **options ):
        return TextSummary ( self, **options )

    def GetTable ( self, **options ):
        return TextTable ( self, **options )

    def GetVerbatim ( self, **options ):
        return Verbatim ( self, **options )

    def LineBreak ( self ):
        """Line break."""
        self.Text ( "\n" )

    def Separator ( self ):
        """Separator."""
        self.Text ( _TextPageWidth * "=" )
        self.LineBreak ( )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class XHTMLLogFileWriter ( LogFileWriter ):
    """XHTML LogFileWriters."""

    _attributable = dict ( LogFileWriter._attributable )
    _attributable.update ( { "styleFile" : None ,
                             "title"     : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( XHTMLLogFileWriter, self )._CheckOptions ( )
        # . Basic header.
        self.Text ( "<?xml version = \"1.0\" encoding = \"utf-8\" ?>\n" )
        self.Text ( "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">\n" )
        self.Text ( "<html xmlns = \"http://www.w3.org/1999/xhtml\" xml:lang = \"en\" lang = \"en\">\n" )
        self.Text ( "<head>\n" )
        if self.title is not None: self.Text ( "    <title>" + self.title + "</title>\n" )
        # . Get the style file.
        if self.styleFile is None: self.styleFile = os.getenv ( "PDYNAMO3_STYLE" )
        if self.styleFile is None:
            self.Text ( "    <style type = \"text/css\">\n" )
            self.Text ( "    body               { background : #ffffff ; color : #000000 }\n"  )
            self.Text ( "    h1, h2, h3, h4, h5 { text-align : center }\n"                   )
            self.Text ( "    table { margin-left : auto ; margin-right : auto ; margin-top : 15px ; margin-bottom : 15px }\n" )
            self.Text ( "    </style>\n"      )
        else:
            self.Text ( "    <link href = \"" + self.styleFile + "\" rel = \"stylesheet\" type = \"text/css\" />\n" )
        self.Text ( "</head>\n<body>" )

    def Close ( self ):
        """Close the file."""
        self.PrintObjectsClose ( )
        if self.isActive:
            self.Text ( "</body>\n</html>\n" )
            self.isActive = False
            if self.path is not None:
                try:
                    self.file.close ( )
                except:
                    raise IOError ( "Cannot close file: " + self.path + "." )

    def GetHeading ( self, **options ):
        return XHTMLHeading ( self, **options )

    def GetParagraph ( self, **options ):
        return XHTMLParagraph ( self, **options )

    def GetSummary ( self, **options ):
        return XHTMLSummary ( self, **options )

    def GetTable ( self, **options ):
        return XHTMLTable ( self, **options )

    def GetVerbatim ( self, **options ):
        return XHTMLVerbatim ( self, **options )

    def LineBreak ( self ):
        """Line break."""
        self.Text ( "<br/>\n" )

    def Separator ( self ):
        """Separator."""
        self.Text ( "<hr/>\n" )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def LogFileActive ( log ):
    """Check to see if printing is to be done on a log file."""
    return ( ( log is not None ) and isinstance ( log, LogFileWriter ) )

#===================================================================================================================================
# . Default log file.
#===================================================================================================================================
logFile = TextLogFileWriter.WithDefaults ( )
#logFile = XHTMLLogFileWriter ( title = "pDynamo Output" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    pass

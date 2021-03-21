"""Classes for print objects."""

from enum import Enum

#===================================================================================================================================
# . Various definitions and default parameters.
#===================================================================================================================================
# . Alignment choices.
class Align ( Enum ):
    """Text align."""
    Center = 1
    Left   = 2
    Right  = 3

# . Heading.
_HeadingPageWidth    = 100

# . Summary.
_SummaryPageWidth    =  80 # . Must be even.
_SummaryValueWidth   =  14
_SummaryKeyWidth     = ( _SummaryPageWidth - 8 ) // 2 - _SummaryValueWidth

# . Table.
_TableColumns        =   2
_TableColumnWidth    =  20

# . XHTML.
_XHTMLBlankCharacter = "&nbsp;"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PrintObject:
    """Base class for all print objects."""

    _attributable = { "isActive" : None ,
                      "owner"    : None }

    def __init__ ( self, owner, **options ):
        """Constructor."""
        self.__dict__.update ( self.__class__._attributable )
        for ( key, value ) in options.items ( ):
            if key in self.__class__._attributable: setattr ( self, key, value )
        self.owner = owner

    def Start ( self ):
        """Start a print object."""
        if self.isActive is None:
            self.owner.PrintObjectPush ( self )
            self.isActive = True

    def Stop ( self ):
        """Stop a print object."""
        if self.isActive:
            self.owner.PrintObjectPop ( )
            self.isActive = False

    def LineBreak ( self ):
        """Print a line break."""
        if self.isActive: self.owner.LineBreak ( )

    def Text ( self, text ):
        """Print some text."""
        if self.isActive: self.owner.Text ( text )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class XHTMLObject:
    """Base class for XHTML objects."""

    def ClassAttributeString ( self ):
        """Return a class attribute sting."""
        if hasattr ( self, "classAttribute" ) and ( self.classAttribute is not None ): return " class = \"" + self.classAttribute + "\""
        else:                                                                          return None

#===================================================================================================================================
# . Heading classes.
#===================================================================================================================================
class Heading ( PrintObject ):
    """Base class for heading writing."""

    _attributable = dict ( PrintObject._attributable )

class TextHeading ( Heading ):
    """Text heading writing."""

    _attributable = { "includeBlankLine" : False             ,
                      "pageWidth"        : _HeadingPageWidth }
    _attributable.update ( Heading._attributable )

    def __init__ ( self, owner, **options ):
        """Constructor."""
        super ( TextHeading, self ).__init__ ( owner )
        self.includeBlankLine = options.get ( "includeBlankLine", False             )
        self.pageWidth        = options.get ( "pageWidth"       , _HeadingPageWidth )

    def Start ( self ):
        """Activation."""
        super ( TextHeading, self ).Start ( )
        if self.isActive:
            if self.includeBlankLine: self.owner.LineBreak ( )
            self.owner.Text ( self.pageWidth * "-" + "\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive: self.owner.Text ( self.pageWidth * "-" + "\n" )
        super ( TextHeading, self ).Stop ( )

    def Text ( self, text ):
        """Text."""
        length = min ( self.pageWidth, len ( text ) )
        nspaces = ( self.pageWidth - length + 1 ) // 2
        self.owner.Text ( nspaces * " " + text[0:length] + "\n" )

class XHTMLHeading ( Heading, XHTMLObject ):
    """XHTML heading writing."""

    _attributable = { "headingTag" : "h1" }
    _attributable.update ( Heading._attributable )

    def __init__ ( self, owner, **options ):
        """Constructor."""
        super ( XHTMLHeading, self ).__init__ ( owner )
        self.headingTag = options.get ( "headingTag", "h1" )

    def Start ( self ):
        """Activation."""
        super ( XHTMLHeading, self ).Start ( )
        if self.isActive:
            self.owner.Text ( "<" + self.headingTag )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive: self.owner.Text ( "\n</" + self.headingTag + ">\n" )
        super ( XHTMLHeading, self ).Stop ( )

#===================================================================================================================================
# . Object tree printing classes - __dict__ attributes only.
#===================================================================================================================================
class ObjectTree ( PrintObject ):
    """Base class for object tree writing."""

    def __init__ ( self, owner, **options ):
        """Constructor."""
        super ( ObjectTree, self ).__init__ ( owner )
        self.nested     = options.get ( "nested"    , False )
        self.withValues = options.get ( "withValues", False )

class TextObjectTree ( ObjectTree ):
    """Text object printer tree writing."""

    def __init__ ( self, owner, **options ):
        """Constructor."""
        super ( TextObjectTree, self ).__init__ ( owner )
        self.indent = options.get ( "indent", 0 )
        self.step   = options.get ( "step"  , 4 )

    def Entry ( self, object, level = 0 ):
        """Entry."""
        if self.isActive:
            indent = self.indent + level * self.step
            if hasattr ( object, "label" ): string = "{:s}{:s} - {:s}:\n".format ( " " * indent, object.__class__.__name__, object.label )
            else:                           string = "{:s}{:s}:\n".format        ( " " * indent, object.__class__.__name__               )
            self.owner.Text ( string )
            indent += self.step
            spacing = " " * indent
            items   = [ ( "{:s}{:s}".format ( spacing, key ), value ) for ( key, value ) in object.__dict__.items ( ) if value is not None ]
            if self.withValues:
                width  = 0
                for ( string, value ) in items: width = max ( width, len ( string ) )
                width += 2
            for ( string, value ) in sorted ( items ):
                doMore = self.nested and hasattr ( value, "__dict__" )
                if doMore: string + ":"
                if withValues: self.owner.Text ( "{:s} - {:s}\n".format ( string.ljust ( width ), repr ( value ) ) )
                else:          self.owner.Text ( "{:s}\n".format ( string ) )
                if doMore: self.Entry ( value, level = level+2 )

#===================================================================================================================================
# . Paragraph classes.
#===================================================================================================================================
class Paragraph ( PrintObject ):
    """Base class for paragraph writing."""
    pass

class TextParagraph ( Paragraph ):
    """Text paragraph writing."""

    def Start ( self ):
        """Activation."""
        super ( TextParagraph, self ).Start ( )
        if self.isActive: self.owner.LineBreak ( )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive: self.owner.LineBreak ( )
        super ( TextParagraph, self ).Stop ( )

class XHTMLParagraph ( Paragraph, XHTMLObject ):
    """XHTML paragraph writing."""

    def __init__ ( self, owner ):
       super ( XHTMLParagraph, self ).__init__ ( owner )

    def Start ( self ):
        """Activation."""
        super ( XHTMLParagraph, self ).Start ( )
        if self.isActive:
            self.owner.Text ( "<p" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive: self.owner.Text ( "\n</p>\n" )
        super ( XHTMLParagraph, self ).Stop ( )

#===================================================================================================================================
# . Summary classes.
#===================================================================================================================================
class Summary ( PrintObject ):
    """Base class for summary writing."""

    _attributable = { "numberOfItems" : 0     ,
                      "order"         : True  ,
                      "withSections"  : False }
    _attributable.update ( PrintObject._attributable )

    def Print ( self, title, sections ):
        """Print the processed summary."""
        if len ( sections ) > 0:
            self.Start ( )
            if title is not None: self.Title ( title )
            for ( section, entries ) in sections:
                if section is not None: self.Section ( section )
                for ( label, value ) in entries: self.Entry ( label, value )
            self.Stop  ( )

    def SetUp ( self, title, items ):
        """Set up the summary."""
        # . Process items list.
        entries  = []
        order    = self.order
        sections = []
        if self.withSections:
            section = None
        else:
            section = title # . Use title as section header for entries-only summary.
            title   = None
        for ( label, value ) in items:
            # . Entry.
            if isinstance ( label, str ) and isinstance ( value, str ):
                entries.append ( ( label, value ) )
            # . Section but only if sections are to be treated.
            elif self.withSections and isinstance ( value, bool ):
                # . Save old section.
                if len ( entries ) > 0:
                    if order: entries.sort ( )
                    sections.append ( ( section, entries ) )
                    entries = []
                # . Start new section.
                order   = value
                section = label
                # . A None section is only possible directly after a title if it is present to avoid duplicating separator.
                if ( section is None ) and ( ( title is not None ) or ( len ( sections ) > 0 ) ): section = ""
        # . Save last section.
        if len ( entries ) > 0:
           if order: entries.sort ( )
           sections.append ( ( section, entries ) )
        return ( title, sections )

    def Start ( self ):
        """Activation."""
        super ( Summary, self ).Start ( )
        self.owner.LineBreak ( )

class TextSummary ( Summary ):
    """Text summary writing."""

    _attributable = { "keyWidth"   : _SummaryKeyWidth   ,
                      "minimize"   : False              ,
                      "pageWidth"  : _SummaryPageWidth  ,
                      "valueWidth" : _SummaryValueWidth }
    _attributable.update ( Summary._attributable )

    # . P = 2 * ( K + V ) + 8.
    def DetermineFieldSizes ( self, title, sections ):
        """Determine field sizes."""
        # . Get maximum field lengths.
        kMax = 0
        tMax = 0
        vMax = 0
        for ( section, entries ) in sections:
            if section is not None:
                tMax = max ( tMax, len ( section ) )
            for ( key, value ) in entries:
                kMax = max ( kMax, len ( key   ) )
                vMax = max ( vMax, len ( value ) )
        if title is not None:
            tMax = max ( tMax, len ( title ) )
        # . Use defaults for zero key or value widths.
        if ( kMax > 0 ) and ( vMax > 0 ):
            # . Minimal pagewidth.
            pWidth = max ( tMax + 2, 2 * ( kMax + vMax ) + 8 )
            if pWidth % 2 != 0: pWidth += 1 # . Should be even.
            # . Minimized representation.
            if self.minimize:
                self.keyWidth   = kMax
                self.pageWidth  = pWidth
                self.valueWidth = vMax
            # . Have a page at least as wide as the default.
            else:
                self.pageWidth = max ( self.pageWidth, pWidth )
                kv             = ( self.pageWidth - 8 ) // 2
                if kMax > self.keyWidth:
                    self.keyWidth   = kMax
                    self.valueWidth = kv - kMax
                elif vMax > self.valueWidth:
                    self.valueWidth = vMax
                    self.keyWidth   = kv - vMax
                elif ( kMax + vMax ) < kv:
                    self.keyWidth   = kMax + ( ( kv - ( kMax + vMax ) ) // 2 )
                    self.valueWidth = kv - self.keyWidth

    def Entry ( self, key, value ):
        """Entry."""
        if self.isActive:
            isNullEntry = ( key   is None ) or ( len ( key   ) == 0 ) or \
                          ( value is None ) or ( len ( value ) == 0 )
            if not isNullEntry:
                if ( self.numberOfItems == 1 ): self.owner.Text ( 2 * " " )
                self.owner.Text ( key.ljust ( self.keyWidth ) + " = " + value.rjust ( self.valueWidth ) )
            if ( self.numberOfItems == 1 ): self.owner.Text ( "\n" )
            self.numberOfItems = ( self.numberOfItems + 1 ) % 2

    def Section ( self, title ):
        """Section divider."""
        if self.isActive:
            if self.numberOfItems > 0: self.Entry ( None, None )
            length = min ( self.pageWidth - 2, len ( title ) )
            n  = length + 2
            n1 = ( self.pageWidth - n + 1 ) // 2
            n2 = self.pageWidth - n - n1
            self.owner.Text ( n1 * "-" + " " + title[0:length] + " " + n2 * "-" + "\n" )

    def SetUp ( self, title, items ):
        """Set up the summary."""
        ( title, sections ) = super ( TextSummary, self ).SetUp ( title, items )
        self.DetermineFieldSizes ( title, sections )
        return ( title, sections )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive:
            if self.numberOfItems > 0: self.Entry ( None, None )
            self.owner.Text ( self.pageWidth * "-" + "\n" )
        super ( TextSummary, self ).Stop ( )

    def Title ( self, title ):
        """Title."""
        if self.isActive:
            if self.numberOfItems > 0: self.Entry ( None, None )
            length = min ( self.pageWidth, len ( title ) )
            before = ( self.pageWidth - length + 1 ) // 2
            self.owner.Text ( self.pageWidth * "-" + "\n" )
            self.owner.Text ( before * " " + title[0:length]+ "\n" )
            self.owner.Text ( self.pageWidth * "-" + "\n" )

class XHTMLSummary ( Summary, XHTMLObject ):
    """XHTML summary writing."""

    def Entry ( self, key, value ):
        """Entry."""
        if self.isActive:
            isNullEntry = ( key   is None ) or ( len ( key   ) == 0 ) or \
                          ( value is None ) or ( len ( value ) == 0 )
            if ( self.numberOfItems == 0 ): self.owner.Text ( "<tr>\n" )
            self.owner.Text ( "<td" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">" )
            if ( isNullEntry ): self.owner.Text ( _XHTMLBlankCharacter )
            else:               self.owner.Text ( key )
            self.owner.Text ( "</td>\n<td" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">" )
            if ( isNullEntry ): self.owner.Text ( _XHTMLBlankCharacter )
            else:               self.owner.Text ( value )
            self.owner.Text ( "</td>\n" )
            if ( self.numberOfItems == 1 ): self.owner.Text ( "</tr>\n" )
            self.numberOfItems = ( self.numberOfItems + 1 ) % 2

    def Section ( self, title ):
        """Section divider."""
        if self.isActive:
            if self.numberOfItems > 0: self.Entry ( None, None )
            self.owner.Text ( "<div" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n<table border = \"1\" cellpadding = \"10\" cellspacing = \"0\">\n<tr><th" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( " colspan = \"4\">" + title + "</th></tr>\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive:
            if self.numberOfItems > 0: self.Entry ( None, None )
            self.owner.Text ( "</table></div>\n" )
        super ( XHTMLSummary, self ).Stop ( )

    def Title ( self, title ):
        """Title."""
        self.Section ( title )

#===================================================================================================================================
# . Table classes.
#===================================================================================================================================
class Table ( PrintObject ):
    """Base class for table writing."""

    _attributable = { "numberOfColumns" : 0 ,
                      "numberOfItems"   : 0 }
    _attributable.update ( PrintObject._attributable )

    def Heading ( self, text, columnSpan = 1 ):
        """Heading."""
        self.Entry ( text, align = Align.Center, columnSpan = columnSpan, isHeader = True )

    def Title ( self, text ):
        """Title spanning the table."""
        self.EndRow ( )
        self.Entry  ( text, align = Align.Center, columnSpan = self.numberOfColumns, isHeader = True )

class TextTable ( Table ):
    """Text table writing."""

    _attributable = { "columnWidths" : 0 ,
                      "pageWidth"    : 0 }
    _attributable.update ( Table._attributable )

    def __init__ ( self, owner, **options ):
        """Constructor."""
        super ( TextTable, self ).__init__ ( owner )
        if "columns" in options: self.columnWidths = options["columns"]
        else:                             self.columnWidths = [ _TableColumnWidth for i in range ( _TableColumns) ]
        self.numberOfColumns = len ( self.columnWidths )
        self.numberOfItems   = 0
        self.pageWidth       = sum ( self.columnWidths )

    def EndRow ( self ):
        """Terminate a row."""
        if self.isActive:
            if self.numberOfItems > 0:
                self.numberOfItems = 0

    def Entry ( self, text, align = Align.Right, columnSpan = 1, fillCharacter = " ", isHeader = False ):
        """Entry."""
        if self.isActive and ( columnSpan > 0 ):
            start = self.numberOfItems + 1
            stop  = min ( self.numberOfColumns, self.numberOfItems + columnSpan )
            n     = sum ( self.columnWidths[start-1:stop] )
            if ( start == 1 ): self.owner.Text ( "\n" )
            if ( text is None ) or ( len ( text ) == 0 ): self.owner.Text ( n * fillCharacter )
            else:
                after  = 0
                before = 0
                blanks = 0
                if len ( text ) <= n: blanks = n - len ( text )
                if   align == Align.Left : after  = blanks
                elif align == Align.Right: before = blanks
                else:
                    before = ( blanks + 1 ) // 2
                    after  = blanks - before
                if ( fillCharacter == " " ) or ( before <= 1 ):
                    self.owner.Text (   before     * " "           +       text )
                else:
                    self.owner.Text ( ( before-1 ) * fillCharacter + " " + text )
                if ( stop < self.numberOfColumns ) or ( fillCharacter != " " ):
                    if after > 1: self.owner.Text ( " " + ( after-1 ) * fillCharacter )
                    else:         self.owner.Text (         after     * fillCharacter )
            if ( stop == self.numberOfColumns ) and isHeader: self.owner.Text ( "\n" + self.pageWidth * "-" )
            self.numberOfItems = stop
            if self.numberOfItems == self.numberOfColumns: self.numberOfItems = 0

    def Section ( self, text ):
        """Section spanning the table."""
        self.EndRow ( )
        self.Entry  ( text, align = Align.Center, columnSpan = self.numberOfColumns, fillCharacter = "-" )

    def Start ( self ):
        """Activation."""
        super ( TextTable, self ).Start ( )
        self.owner.Text ( "\n" + self.pageWidth * "-" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive:
            self.owner.Text ( "\n" + self.pageWidth * "-" + "\n" )
        super ( TextTable, self ).Stop ( )

class XHTMLTable ( Table, XHTMLObject ):
    """XHTML table writing."""

    def __init__ ( self, owner, **options ):
        """Constructor."""
        super ( XHTMLTable, self ).__init__ ( owner )
        if "columns" in options: self.numberOfColumns = len ( options["columns"] )
        else:                             self.numberOfColumns = _TableColumns
        self.numberOfItems = 0

    def EndRow ( self ):
        """Terminate a row."""
        if self.isActive:
            if self.numberOfItems > 0:
                for i in range ( self.numberOfItems, self.numberOfColumns ): self.Entry ( None )
                self.numberOfItems = 0

    def Entry ( self, text, align = Align.Right, columnSpan = 1, isHeader = False ):
        """Entry."""
        if self.isActive and ( columnSpan > 0 ):
            start = self.numberOfItems + 1
            stop  = min ( self.numberOfColumns, self.numberOfItems + columnSpan )
            if start == 1: self.owner.Text ( "\n<tr>" )
            if isHeader  : self.owner.Text ( "<th" )
            else:          self.owner.Text ( "<td" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( " colspan = \"{:d}\">".format ( columnSpan ) )
            if ( text is None ) or ( len ( text ) == 0 ): self.owner.Text ( _XHTMLBlankCharacter )
            else:                                         self.owner.Text ( text )
            if ( isHeader ): self.owner.Text ( "</th>" )
            else:            self.owner.Text ( "</td>" )
            self.owner.Text ( "\n" )
            if stop == self.numberOfColumns: self.owner.Text ( "</tr>\n" )
            self.numberOfItems = stop
            if self.numberOfItems == self.numberOfColumns: self.numberOfItems = 0

    def Section ( self, title ):
        """Section divider."""
        self.EndRow ( )
        self.Entry  ( text, align = Align.Center, columnSpan = self.numberOfColumns, isHeader = True )

    def Start ( self ):
        """Activation."""
        super ( XHTMLTable, self ).Start ( )
        if self.isActive:
            self.owner.Text ( "<div" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n<table border = \"1\" cellpadding = \"10\" cellspacing = \"0\">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive:
            self.EndRow ( )
            self.owner.Text ( "</table></div>\n" )
        super ( XHTMLTable, self ).Stop ( )

#===================================================================================================================================
# . Verbatim classes.
#===================================================================================================================================
class Verbatim ( PrintObject ):
    """Base class for verbatim writing.

    This can be used directly for text output.
    """

    pass

class XHTMLVerbatim ( Verbatim, XHTMLObject ):
    """XHTML verbatim writing."""

    def __init__ ( self, owner ):
       super ( XHTMLVerbatim, self ).__init__ ( owner )

    def Start ( self ):
        """Activation."""
        super ( XHTMLVerbatim, self ).Start ( )
        if self.isActive:
            self.owner.Text ( "<pre" )
            self.owner.Text ( self.ClassAttributeString ( ) )
            self.owner.Text ( ">\n" )

    def Stop ( self ):
        """Deactivation."""
        if self.isActive: self.owner.Text ( "</pre>\n" )
        super ( XHTMLVerbatim, self ).Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":

    pass

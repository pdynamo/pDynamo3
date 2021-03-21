"""Classes and functions for importing and exporting objects to and from external files."""

import os.path

from  pCore                 import Align           , \
                                   logFile         , \
                                   LogFileActive   , \
                                   Pickle          , \
                                   Unpickle        , \
                                   YAMLPickle      , \
                                   YAMLUnpickle
from  pMolecule             import System      
from  pScientific.Geometry3 import Coordinates3
from .TrajectoryMixin       import TrajectoryMixin

# . Need a more general way of handling subclasses of particular object types (especially system).

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class ExportImportError ( Exception ): pass

class ExportImportHandler:
    """Class for a handler that imports or exports objects."""

    def __init__ ( self, functionIndex, identifiers, label, defaultFunction = None ):
        """Constructor."""
        self.functionIndex = functionIndex
        self.identifiers   = identifiers
        self.label         = label
        if defaultFunction is None: self.defaultFunction = functionIndex.get ( None, None )
        else:                       self.defaultFunction = defaultFunction

    # . A fudge until trajectory composition/inheritance sorted out.
    def ActivateTrajectory ( self, path, owner, **options ):
        """Activate a trajectory for export or import."""
        function = self.GetFunction ( TrajectoryMixin )
        if options is None: item = function ( path, owner            )
        else:               item = function ( path, owner, **options )
        if not getattr ( item, "isTrajectory", False ): #not isinstance ( item, TrajectoryMixin ):
            raise ExportImportError ( "Invalid class for activated trajectory: " + item.__class__.__name__ + "." )
        return item

    def ExportObject ( self, path, object, **options ):
        """Export an object."""
        function = self.GetFunction ( object.__class__ )
        if options is None: function ( path, object            )
        else:               function ( path, object, **options )

    def GetFunction ( self, objectClass ):
        """Get the appropriate function to handle a given object class."""
        function = self.functionIndex.get ( objectClass, self.defaultFunction )
        if function is None:
            if objectClass is None: tag = "default class"
            else:                   tag = "class \"" + objectClass.__name__ + "\""
            raise ExportImportError ( "Handler \"{:s}\" cannot treat objects of {:s}.".format ( self.label, tag ) )
        return function

    def ImportObject ( self, path, objectClass, **options ):
        """Import an object."""
        function = self.GetFunction ( objectClass )
        if options is None: item = function ( path            )
        else:               item = function ( path, **options )
        if ( objectClass is not None ) and ( not isinstance ( item, objectClass ) ):
            raise ExportImportError ( "Invalid class for imported object: " + item.__class__.__name__ + "." )
        return item

    @property
    def identifierString ( self ):
        """Return an identifier string."""
        items = set ( self.identifiers )
        items.remove ( self.label )
        items = list ( items )
        items.sort ( )
        return ", ".join ( items )

    @property
    def objectString ( self ):
        """Return an object string."""
        items = []
        for item in self.functionIndex.keys ( ):
            if item is None: items.append ( "All"         )
            else:            items.append ( item.__name__ )
        items.sort ( )
        return ", ".join ( items )

class ExportImportHandlerContainer:
    """A container for import/export handlers."""

    def __init__ ( self ):
        """Constructor."""
        self.handlerIndex = {}
        self.handlers     = []
        self.identifiers  = set ( )

    def AddHandler ( self, functionIndex, identifiers, label, defaultFunction = None ):
        """Add a handler."""
        # . Check for duplicate identifiers.
        identifiers = set ( identifiers )
        identifiers.add ( label )
        common      = self.identifiers & identifiers
        if len ( common ) > 0:
            common = list ( common )
            common.sort ( )
            raise ExportImportError ( "Duplicate handler identifiers: {:s}.".format ( str ( common ) ) )
        # . Create the object.
        handler = ExportImportHandler ( functionIndex, identifiers, label, defaultFunction = defaultFunction )
        for identifier in identifiers:
            self.handlerIndex[identifier] = handler
        self.handlers.append    ( handler     )
        self.identifiers.update ( identifiers )

    def FileFormats ( self ):
        """Get a list of file formats and their identifiers."""
        data = []
        for handler in self.handlers:
            label       = handler.label
            identifiers = list ( handler.identifiers )
            identifiers.remove ( label )
            identifiers.sort ( )
            data.append ( ( label, identifiers ) )
        data.sort ( )
        return data

    def GetHandler ( self, path, format = None ):
        """Get a handler given a path and format.

        The format is determined from the path extension if no format is given.
        """
        if format is None:
            format = os.path.splitext ( path )[-1]
            if format.startswith ( "." ): format = format[1:]
        try   : return self.handlerIndex[format]
        except: raise ExportImportError ( "Unrecognized format: " + format + "." )

    def HandlerSummary ( self, log = logFile, title = "Export/Import Handler Container" ):
        """Handler summary."""
        if LogFileActive ( log ):
            # . Gather data.
            data = []
            l0 = l1 = l2 = 18
            for handler in self.handlers:
                label       = handler.label
                identifiers = handler.identifierString
                objects     = handler.objectString
                l0 = max ( l0, len ( label       ) )
                l1 = max ( l1, len ( identifiers ) )
                l2 = max ( l2, len ( objects     ) )
                data.append ( ( label, identifiers, objects ) )
            data.sort ( )
            # . Output.
            table = log.GetTable ( columns = [ l0+2, l1+2, l2 ] )
            table.Start   ( )
            table.Title   ( title )
            table.Heading ( "Format Label"      )
            table.Heading ( "Other Identifiers" )
            table.Heading ( "Objects"           )
            for ( label, identifiers, objects ) in data:
                table.Entry ( label      , align = Align.Left )
                table.Entry ( identifiers, align = Align.Left )
                table.Entry ( objects    , align = Align.Left )
            table.Stop ( )

#===================================================================================================================================
# . Instantiate an exporter and an importer.
#===================================================================================================================================
# . Instantiation.
_Exporter = ExportImportHandlerContainer ( )
_Importer = ExportImportHandlerContainer ( )

# . Add formats from outside pBabel.
# . None means that objects of all classes are handled.
_Exporter.AddHandler ( { None :       Pickle }, [ "pkl" , "PKL"  ], "Pickle"      )
_Exporter.AddHandler ( { None :   YAMLPickle }, [ "yaml", "YAML" ], "YAML Pickle" )
_Importer.AddHandler ( { None :     Unpickle }, [ "pkl" , "PKL"  ], "Pickle"      )
_Importer.AddHandler ( { None : YAMLUnpickle }, [ "yaml", "YAML" ], "YAML Pickle" )

# . Surface files.
from pScientific.Surfaces import OOGLOffFileReader , \
                                 OOGLOffFileWriter , \
                                 STLFileReader     , \
                                 STLFileWriter     , \
                                 PolygonalSurface
_Exporter.AddHandler ( { PolygonalSurface : OOGLOffFileWriter.PathFromPolygonalSurface } , [ "off", "oogloff" ], "OOGL Off" )
_Exporter.AddHandler ( { PolygonalSurface : STLFileWriter.PathFromPolygonalSurface     } , [ "stl", "STL"     ], "STL"      )
_Importer.AddHandler ( { PolygonalSurface : OOGLOffFileReader.PathToPolygonalSurface   } , [ "off", "oogloff" ], "OOGL Off" )
_Importer.AddHandler ( { PolygonalSurface : STLFileReader.PathToPolygonalSurface       } , [ "stl", "STL"     ], "STL"      )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
# . In these functions the keyword arguments include "format", "log" and the arguments
# . applicable to the export/import function found by the handler. When function-specific
# . arguments are required it is probably best to use the appropriate function directly
# . and avoid these generic export/import functions.
def ExportImportPathOutput ( path ):
    """Path shortener."""
    home = os.path.expanduser ( "~" )
    if path.startswith ( home ): return path.replace ( home, "~", 1 )
    else:                        return path

def ExportFileFormats ( ):
    """Return a list of the export file formats."""
    return _Exporter.FileFormats ( )

def ExportOptions ( log = logFile ):
    """List the export options."""
    _Exporter.HandlerSummary ( log = log, title = "Export Options" )

def ExportSystem ( path, system, **options ):
    """Export a system."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Exporter.GetHandler ( path, format = format )
    handler.ExportObject ( path, system, **options )
    if LogFileActive ( log ):
        log.Paragraph ( "System exported to \"{:s}\" in {:s} format.".format ( ExportImportPathOutput ( path ), handler.label ) )

def ExportTrajectory ( path, owner, **options ):
    """Activate a trajectory for export."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Exporter.GetHandler ( path, format = format )
    item    = handler.ActivateTrajectory ( path, owner, **options )
    if LogFileActive ( log ):
        log.Paragraph ( "Trajectory activated for export to \"{:s}\" in {:s} format.".format ( ExportImportPathOutput ( path ), handler.label ) )
    return item

def ImportCoordinates3 ( path, **options ):
    """Import a set of coordinates."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ImportObject ( path, Coordinates3, **options )
    if LogFileActive ( log ):
        log.Paragraph( "Coordinates3 imported from \"{:s}\" in {:s} format.".format ( ExportImportPathOutput ( path ), handler.label ) )
    return item

def ImportFileFormats ( ):
    """Return a list of the import file formats."""
    return _Importer.FileFormats ( )

def ImportObjects ( path, **options ):
    """Import objects from the file."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ImportObject ( path, None, **options )
    if LogFileActive ( log ):
        log.Paragraph( "Objects imported from \"{:s}\" in {:s} format.".format ( ExportImportPathOutput ( path ), handler.label ) )
    return item

def ImportOptions ( log = logFile ):
    """List the import options."""
    _Importer.HandlerSummary ( log = log, title = "Import Options" )

def ImportSystem ( path, **options ):
    """Import a system."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ImportObject ( path, System, **options )
    if LogFileActive ( log ):
        log.Paragraph( "System imported from \"{:s}\" in {:s} format.".format ( ExportImportPathOutput ( path ), handler.label ) )
    return item

def ImportTrajectory ( path, owner, **options ):
    """Activate a trajectory for import."""
    format  = options.pop ( "format", None    )
    log     = options.pop ( "log"   , logFile )
    handler = _Importer.GetHandler ( path, format = format )
    item    = handler.ActivateTrajectory ( path, owner, **options )
    if LogFileActive ( log ):
        log.Paragraph( "Trajectory activated for import from \"{:s}\" in {:s} format.".format ( ExportImportPathOutput ( path ), handler.label ) )
    return item

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

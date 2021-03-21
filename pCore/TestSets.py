"""Classes and functions for running sets of test scripts."""

import glob, os, os.path, subprocess, sys, time

from .AttributableObject import AttributableObject
from .Serialization      import YAMLUnpickle
from .Time               import CPUTime
from .Version            import __version__

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Definitions.
_ErrorExtension      = ".err"
_LabelWidth          = 50
_LineWidth           = 85
_LogExtension        = ".log"
_ScriptExtension     = ".py"

# . Path names.
_DataPath            = "data"
_ErrorPath           = "errors"
_ExamplesPath        = "examples"
_GeneratedFilesPath  = "generatedFiles"
_LogPath             = "logs"
_OutputSubPaths      = [ _ErrorPath, _GeneratedFilesPath, _LogPath ]
_ScriptsFileName     = "__scripts__.yaml"

# . Status flags.
_StatusCode_Pass           = 0 ; _StatusLabel_Pass           = "Pass"
_StatusCode_Error          = 1 ; _StatusLabel_Error          = "Error"
_StatusCode_Fail           = 2 ; _StatusLabel_Fail           = "Fail"
_StatusCode_NotImplemented = 3 ; _StatusLabel_NotImplemented = "Not Implemented"
_StatusCode_NotInstalled   = 4 ; _StatusLabel_NotInstalled   = "Not Installed"
_StatusCodes               = { _StatusCode_Pass           : _StatusLabel_Pass           ,
                               _StatusCode_Error          : _StatusLabel_Error          ,
                               _StatusCode_Fail           : _StatusLabel_Fail           ,
                               _StatusCode_NotImplemented : _StatusLabel_NotImplemented ,
                               _StatusCode_NotInstalled   : _StatusLabel_NotInstalled   }
_StatusLabel_Crash         = "Crash"
_StatusLabel_Missing       = "Missing"

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class TestSetException ( Exception ):
    """Exception class."""
    pass

class TestSet ( AttributableObject ):
    """The test set class."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "doLong"        : False ,
                             "inputRoot"     : None  ,
                             "inputPath"     : None  ,
                             "logExtensions" : None  ,
                             "name"          : None  ,
                             "outputPath"    : None  ,
                             "outputRoot"    : None  ,
                             "pythonCommand" : None  ,
                             "scripts"       : None  } )

    def __len__ ( self ):
        if self.scripts is None: return 0
        else:                    return len ( self.scripts )

    def FinalizeOutputPaths ( self ):
        """Finalize output paths."""
        for tail in _OutputSubPaths:
            path = os.path.join ( self.outputPath, tail )
            if os.path.exists ( path ) and ( len ( os.listdir ( path ) ) == 0 ): os.rmdir ( path )

    @staticmethod
    def GenerateTestScriptsPath ( items ):
        """Generate a path given items from a test scripts file."""
        if   isinstance ( items, str  ): return items
        elif isinstance ( items, list ): return os.path.join ( *items )
        else: raise TestSetException ( "Invalid path items: {!r}.".format ( items ) )

    def InitializeOutputPaths ( self ):
        """Initialize output paths."""
        if not os.path.exists ( self.outputPath ): os.makedirs ( self.outputPath )
        for tail in _OutputSubPaths:
            path = os.path.join ( self.outputPath, tail )
            if not os.path.exists ( path ): os.mkdir ( path )

    def ProcessScriptsFile ( self ):
        """Process the test scripts file."""
        data               = YAMLUnpickle ( os.path.join ( self.inputPath, _ScriptsFileName ) )
        self.logExtensions = data.get ( "Log Extensions", {} )
        self.scripts       = [ TestSet.GenerateTestScriptsPath ( datum ) for datum in data["Scripts"] ]
        if not self.doLong:
            longPaths = [ TestSet.GenerateTestScriptsPath ( datum ) for datum in data.get ( "Long Scripts", [] ) ]
            for longPath in longPaths: self.scripts.remove ( longPath )

    def Run ( self ):
        """Run the test set."""
        # . Setup paths.
        self.InitializeOutputPaths ( )
        # . Local header.
        print ( "\n" + _LineWidth * "-"   )
        print ( "Example script set: {:s}".format ( self.name ) )
        print ( _LineWidth * "-" + "\n"   )
        # . Loop over scripts.
        numberOfSuccesses = 0
        totalTime         = 0.0
        for script in self.scripts:
            # . Paths.
            logExtension       = self.logExtensions.get ( script, _LogExtension )
            ( head, tail )     = os.path.split ( script )
            errorPath          = os.path.join ( self.outputPath, _ErrorPath         , tail   + _ErrorExtension  )
            generatedFilesPath = os.path.join ( self.outputPath, _GeneratedFilesPath                            )
            logPath            = os.path.join ( self.outputPath, _LogPath           , tail   + logExtension     )
            scriptPath         = os.path.join ( self.inputPath ,                      script + _ScriptExtension )
            for path in ( errorPath, logPath ):
                if os.path.exists ( path ): os.remove ( path )
            # . Check that the script file exists.
            if os.path.exists ( scriptPath ):
                # . Run the script.
                time0 = time.time ( )
                eFD   = open ( errorPath, "w" )
                oFD   = open ( logPath  , "w" )
                try:
                    process = subprocess.Popen ( [ self.pythonCommand, scriptPath ], stderr = eFD, stdout = oFD )
                    process.wait ( )
                    label = _StatusCodes.get ( process.returncode, _StatusLabel_Crash )
                except Exception as e:
                    eFD.write ( e.args[0] ) # . Better here (also traceback).
                    label = _StatusLabel_Error
                scriptTime  = time.time ( ) - time0
                totalTime  += scriptTime
                # . Close files.
                eFD.close ( )
                oFD.close ( )
                # . Check for empty files.
                for path in ( errorPath, logPath ):
                    if os.path.exists ( path ):
                        fd    = open ( path, "r" )
                        lines = fd.readlines ( )
                        fd.close ( )
                        isEmpty  = ( len ( lines ) == 0 )
                        if isEmpty: os.remove ( path )
            # . Script file not found.
            else:
                label      = _StatusLabel_Missing
                scriptTime = 0.0
            # . Check the result.
            if label == _StatusLabel_Pass: numberOfSuccesses += 1
            # . Printing.
            print ( "{:s}{:15s}{:>20s}".format ( script.ljust ( _LabelWidth ), label, CPUTime.TimeToString ( scriptTime ) ) )
        # . Finalization.
        if len ( self.scripts ) > 1:
            label = "{:d}/{:d}".format ( numberOfSuccesses, len ( self.scripts ) )
            print ( "{:s}{:15s}{:>20s}".format ( "** Totals **".ljust ( _LabelWidth ), label, CPUTime.TimeToString ( totalTime ) ) )
        self.FinalizeOutputPaths ( )
        # . Finish up.
        return totalTime

    def SetUpPaths ( self ):
        """Set up basic paths."""
        directory       = self.name.replace ( ".", os.sep )
        self.inputPath  = os.path.join ( self.inputRoot  , directory )
        self.outputPath = os.path.join ( self.outputRoot , directory )

    @classmethod
    def WithOptions ( selfClass, **options ):
        """Constructor with options."""
        self = super ( TestSet, selfClass ).WithOptions ( **options )
        self.SetUpPaths         ( )
        self.ProcessScriptsFile ( )
        return self

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def _FindScriptsFiles ( rootPath ):
    """Find all example directories that possess a test scripts file."""
    examplesPath = os.path.join ( rootPath, _ExamplesPath )
    n            = len ( examplesPath ) + len ( os.sep )
    paths        = [ ]
    for path in glob.iglob ( os.path.join ( examplesPath, "**", _ScriptsFileName ), recursive = True ):
        ( head, tail ) = os.path.split ( path )
        paths.append ( head[n:].replace ( os.sep, "." ) )
    return sorted ( paths )

def _GetEnvironmentVariable ( variable, label ):
    """Get an environment variable."""
    result = os.getenv ( variable )
    if result is None: raise TestSetException ( label + " not found." )
    return result

def _InputRoot ( ):
    """Set up the input root path."""
    root      = _GetEnvironmentVariable ( "PDYNAMO3_HOME", "pDynamo root path" )
    inputRoot = os.path.join ( root, _ExamplesPath )
    return inputRoot

def _OutputRoot ( ):
    """Set up the output root path."""
    # . Recursively create all directories (including scratch) if they do not exist - useful for tests after installation.
    scratch    = _GetEnvironmentVariable ( "PDYNAMO3_SCRATCH", "pDynamo scratch path" )
    outputRoot = os.path.join ( scratch, _ExamplesPath + "-" + __version__ )
    if not os.path.exists ( outputRoot ): os.makedirs ( outputRoot ) 
    return outputRoot

def TestScript_InputDataPath ( name ):
    """Get the input data path for a set of examples."""
    return os.path.join ( TestScript_InputPath ( name ), _DataPath )

def TestScript_InputPath ( name ):
    """Get the input path for a set of examples."""
    inputRoot = _InputRoot  ( )
    return os.path.join ( inputRoot, name )

def TestScript_OutputDataPath ( name ):
    """Get and set up the output data path for a set of examples."""
    outputRoot = TestScript_OutputPath ( name )
    path       = os.path.join ( outputRoot, _GeneratedFilesPath )
    if not os.path.exists ( path ): os.mkdir ( path )
    return path

def TestScript_OutputPath ( name ):
    """Get and set up the output path for a set of examples."""
    outputRoot = _OutputRoot ( )
    path       = os.path.join ( outputRoot, name )
    if not os.path.exists ( path ): os.mkdir ( path )
    return path

def TestScriptExit_Fail ( ):
    """Signal a failure on exit from a test script."""
    sys.exit ( _StatusCode_Fail )

def TestScriptExit_NotImplemented ( ):
    """Signal a not implemented exit from a test script."""
    sys.exit ( _StatusCode_NotImplemented )

def TestScriptExit_NotInstalled ( ):
    """Signal a not installed exit from a test script."""
    sys.exit ( _StatusCode_NotInstalled )

def TestSets_Get ( names ):
    """Get the names of known and unknown test sets."""
    rootPath    = _GetEnvironmentVariable ( "PDYNAMO3_HOME", "pDynamo root path" )
    allTestSets = _FindScriptsFiles   ( rootPath )
    if ( names is None ) or ( len ( names ) == 0 ):
        return ( allTestSets, [] )
    else:
        possible  = set ( allTestSets )
        specified = set ( names             )
        found     = ( specified & possible )
        notFound  = ( specified - possible )
        return ( sorted ( found ), sorted ( notFound ) )

def TestSets_Run ( names, doLong = False ):
    """Run a list of test sets."""
    # . Options.
    inputRoot     = _InputRoot  ( )
    outputRoot    = _OutputRoot ( )
    pythonCommand = _GetEnvironmentVariable ( "PDYNAMO3_PYTHONCOMMAND", "pDynamo python command" )
    # . Global header.
    print ( "\nRunning test script sets ..." )
    print ( "Please be patient ... this may take some time ..." )
    print ( "\nLogs and other files are in " + outputRoot + "." )
    # . Loop over sets.
    numberRun = 0
    totalTime = 0.0
    for name in names:
        testSet = TestSet.WithOptions ( doLong        = doLong        ,
                                        inputRoot     = inputRoot     ,
                                        name          = name          ,
                                        outputRoot    = outputRoot    ,
                                        pythonCommand = pythonCommand )
        if len ( testSet ) == 0:
            print ( "\nThere were no scripts for {:s}".format ( testSet.name ) )
        else:
            totalTime += testSet.Run ( )
            numberRun += 1
    # . Finish up.
    print ( "\n" + _LineWidth * "-" )
    if numberRun > 1:
        print ( "\nTotal time used: " + CPUTime.TimeToString ( totalTime, compact = False ) + "." )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    rootPath = _RootPath ( )
    paths    = _FindScriptsFiles ( rootPath )
    if len ( paths ) > 0:
        print ( "Test script sets found:" )
        for path in paths: print ( path )
    else:
        print ( "No test script sets were found." )

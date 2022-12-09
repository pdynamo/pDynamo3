"""Utilities for compilation and installation."""

# . Note that distutils is deprecated so this will have to be replaced by setuptools!
# . In any case, paths need to be sorted out!

import functools, glob, imp, os, os.path, shutil, sys, sysconfig

from distutils.core import setup, Extension
#from distutils.util import get_platform

# . Attempts have been made to use unqualified pyx names.
#   Everything compiles fine but runtime errors arise
#   due to module names not being found.

# . Note that the order of dependencies in the dependency files matters (least dependent first).

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Directory names.
_CIncludeDirectory     = "cinclude"
_CLibraryDirectory     = "clibrary"
_ExtensionsDirectory   = "__extensions__"
_PSourceDirectory      = "psource"
_PyrexDirectory        = "pyrex"
_ThirdPartyDirectory   = "thirdParty"

# . Environment variables.
_PDynamoHome           = "PDYNAMO3_HOME"

# . Extension names.
_CExtension            = ".c"
_LibraryExtension      = ".a"
_PyrexExtension        = ".pyx"
_SharedObjectExtension = ".so"
_SharedObjectBuildPathExtension = ( "." + imp.get_tag ( ) + "*" + _SharedObjectExtension )
#_SharedObjectBuildPathExtension = ( "." + imp.get_tag ( ) + "m*" + _SharedObjectExtension )

# . File names.
_DependencyFile        = "__dependencies__"
_InitFile              = "__init__.py"
_LibraryFile           = "__libraries__"

# . Language level (2 or 3).
_LanguageLevel = 3

# . Libraries.
# . May need to add "irc" to libraries for Intel C compiler.
# . Default library data.
_SystemLibraries    = [ "m" ]
_SystemLibraryPaths = [ "/usr/local/lib" , "/usr/lib" ]

# . Atlas library data.
# . This will depend on compiler and operating system so it needs to be generalized.
_AtlasCLibrariesToOmit  = [ "dcblas", "df2c", "df2cblas", "df2cdlamch", "df2clapack" ]
_AtlasLibraryPaths      = [ "/usr/local/lib/gcc/x86_64-linux-gnu/4.6" , "/usr/local/lib" , "/usr/lib" ]
_AtlasSerialLibraries   = [ "lapack"   , "f77blas"   , "cblas"   , "atlas" , "m" , "gfortran" ]
_AtlasThreadedLibraries = [ "ptlapack" , "ptf77blas" , "ptcblas" , "atlas" , "m" , "gfortran" ]

# . Compiler flags and options.
_FullOptimizationFlag   = "-O3"
_NoOptimizationFlag     = "-O0"
_NumberOfThreadsFlag    = "-DMAXIMUMNUMBEROFTHREADS={:d}"
_OpenMPFlag             = "-fopenmp"
_OpenMPPreprocessorFlag = "-DUSEOPENMP"
_OptimizationFlagPrefix = "-O"

# . Shell script paths.
_ShellScriptPath        = "shellScripts"
_TemplatePath           = os.path.join ( "shellScripts", "templates" )

#===================================================================================================================================
# . Class to allow compilation of C-libraries with additional compiler options.
# . Originally designed for removing all optimization (necessary for dlamch in clapack!).
#===================================================================================================================================
# . Specific imports.
from distutils                    import log
from distutils.command.build_clib import build_clib
from distutils.errors             import DistutilsSetupError

# . Class definition.
class BuildClibWithCompileOptions ( build_clib ):
    """Subclass of distutils.command.build_clib which enables compiler options to be played with."""

    def build_libraries ( self, libraries ):
        """Build libraries."""
        # . Loop over libraries.
        for ( lib_name, build_info ) in libraries:
            # . Get sources.
            sources = build_info.get ( 'sources' )
            if ( sources is None ) or ( type ( sources ) not in ( list, tuple ) ):
                raise DistutilsSetupError ( "in 'libraries' option (library '{:s}') 'sources' must be present and must be a list of source filenames".format ( lib_name ) )
            sources = list ( sources )
            # . Remove optimization strings from compiler options if requested.
            if build_info.get ( "removeOptimizationFlags", False ):
                new = []
                for token in self.compiler.compiler_so:
                    if not token.startswith ( _OptimizationFlagPrefix ): new.append ( token )
                self.compiler.compiler_so = new
            # . Do the work.
            log.info ( "building '{:s}' library".format ( lib_name ) )
            # . Compile the source code to object files in the library directory.
            extra_preargs = build_info.get ( "extra_preargs" )
            include_dirs  = build_info.get ( "include_dirs"  )
            macros        = build_info.get ( "macros"        )
            objects       = self.compiler.compile ( sources, output_dir = self.build_temp, macros = macros, include_dirs = include_dirs, debug = self.debug, extra_preargs = extra_preargs )
            #  Create a static library.
            self.compiler.create_static_lib ( objects, lib_name, output_dir = self.build_clib, debug = self.debug )

#===================================================================================================================================
# . Error class.
#===================================================================================================================================
class InstallationError ( Exception ):
    pass

#===================================================================================================================================
# . Installation options class.
#===================================================================================================================================
class InstallationOptions:
    """Class for determining installation options."""

    def __init__ ( self ):
        """Constructor."""
        pass

    @classmethod
    def FromCommandLine ( selfClass ):
        """Constructor by getting arguments from a command line."""
        # . Initialization.
        self = selfClass ( )
        # . Parse the command line.
        from argparse import ArgumentParser
        parser = ArgumentParser ( epilog = "The package names are optional as all packages are installed by default." )
        parser.add_argument (        "--atlas"          , action = "store_true"  , dest = "useSerialAtlas"   , default = False , help = "link with the serial ATLAS libraries"           )
        parser.add_argument ( "-c" , "--cLibraries"     , action = "store_true"  , dest = "doCLibraries"     , default = False , help = "install C libraries                  [default]" )
        parser.add_argument ( "-e" , "--extensions"     , action = "store_true"  , dest = "doExtensions"     , default = False , help = "install Cython/Pyrex extensions      [default]" )
        parser.add_argument ( "-f" , "--full"           , action = "store_true"  , dest = "doFull"           , default = False , help = "do a full installation"                         )
        parser.add_argument ( "-l" , "--list"           , action = "store_true"  , dest = "listPackages"     , default = False , help = "list packages only"                             )
        parser.add_argument (        "--noClearUp"      , action = "store_false" , dest = "doClearUp"        , default = True  , help = "do not clear up after installation"             )
        parser.add_argument ( "-n" , "--noThirdParty"   , action = "store_false" , dest = "doThirdParty"     , default = True  , help = "do not install thirdparty libraries"            )
        parser.add_argument (        "--openMP"         , action = "store_true"  , dest = "useOpenMP"        , default = False , help = "compile with OpenMP"                            )
        parser.add_argument (        "packageName"      , nargs  = "*"           ,                                               help = "the name of the package to install"             )
        parser.add_argument (        "--ptAtlas"        , action = "store_true"  , dest = "useThreadedAtlas" , default = False , help = "link with the threaded ATLAS libraries"         )
        parser.add_argument ( "-p" , "--pyrex"          , action = "store_true"  , dest = "doPyrex"          , default = False , help = "convert Cython/Pyrex to C"                      )
        parser.add_argument ( "-q" , "--quiet"          , action = "store_false" , dest = "verbose"          ,                   help = "quiet mode"                                     )
        parser.add_argument ( "-s" , "--shell"          , action = "store_true"  , dest = "doShellFiles"     , default = False , help = "write shell files                    [default]" )
        parser.add_argument (        "--threads"        , type   = int           , dest = "threads"          , default = -1    , help = "the number of threads"                          )
        parser.add_argument (        "--thirdPartyOnly" , action = "store_true"  , dest = "doThirdPartyOnly" , default = False , help = "install thirdparty libraries only"              )
        parser.add_argument ( "-u" , "--update"         , action = "store_true"  , dest = "update"           , default = False , help = "update package and all those that depend on it" )
        parser.add_argument ( "-v" , "--verbose"        , action = "store_true"  , dest = "verbose"          , default = True  , help = "verbose mode                         [default]" )
        arguments = parser.parse_args ( )
        # . Check what is to be done.
        # . A full installation overrides everything else.
        if arguments.doFull:
            doCLibraries = doExtensions = doPyrex = doShellFiles = True
        # . Get flags.
        else:
            doCLibraries = arguments.doCLibraries
            doExtensions = arguments.doExtensions
            doPyrex      = arguments.doPyrex
            doShellFiles = arguments.doShellFiles
            # . If all are False do a standard installation.
            if ( not doCLibraries ) and ( not doExtensions ) and ( not doPyrex ) and ( not doShellFiles ): doCLibraries = doExtensions = doShellFiles = True
        # . Get the maximum number of threads for the current processor.
        maximumThreads   = selfClass.MaximumNumberOfThreads ( )
        # . Remaining arguments.
        # . Threaded Atlas overrides serial Atlas on threaded machines.
        doClearUp        = arguments.doClearUp
        doThirdParty     = arguments.doThirdParty
        doThirdPartyOnly = arguments.doThirdPartyOnly
        numberOfThreads  = 0
        packageNames     = set ( arguments.packageName )
        useOpenMP        = arguments.useOpenMP        and ( maximumThreads > 1 )
        useThreadedAtlas = arguments.useThreadedAtlas and ( maximumThreads > 1 )
        useSerialAtlas   = arguments.useSerialAtlas and ( not useThreadedAtlas )
        update           = arguments.update
        verbose          = arguments.verbose
        # . Process the number of threads when OpenMP is to be used.
        # . Is this needed?
        if useOpenMP:
            if maximumThreads > 2: numberOfThreads = max ( min ( arguments.threads, maximumThreads ), 2 )
            else:                  numberOfThreads = 2
        # . Set all attributes.
        self.doClearUp        = doClearUp
        self.doCLibraries     = doCLibraries
        self.doExtensions     = doExtensions
        self.listPackages     = arguments.listPackages
        self.doPyrex          = doPyrex
        self.doShellFiles     = doShellFiles
        self.doThirdParty     = doThirdParty
        self.doThirdPartyOnly = doThirdPartyOnly
        self.numberOfThreads  = numberOfThreads
        self.packageNames     = packageNames
        self.useOpenMP        = useOpenMP
        self.useSerialAtlas   = useSerialAtlas
        self.useThreadedAtlas = useThreadedAtlas
        self.update           = update
        self.verbose          = verbose
        if self.doThirdPartyOnly:
            self.doExtensions = False
            self.doPyrex      = False
            self.doShellFiles = False
        # . Finish up.
        return self

    def IsCLibraryToBeCompiled ( self, libraryName, isThirdParty ):
        """Check whether a C library is to be compiled."""
        doCompilation = ( isThirdParty and ( self.doThirdParty or self.doThirdPartyOnly ) ) or \
                        ( ( not isThirdParty ) and ( not self.doThirdPartyOnly ) )
        if doCompilation and ( self.useSerialAtlas or self.useThreadedAtlas ):
            doCompilation = ( libraryName not in _AtlasCLibrariesToOmit )
        return doCompilation

    def IsCLibraryToBeLinked ( self, libraryName ):
        """Check whether a C library is to be linked."""
        if ( self.useSerialAtlas or self.useThreadedAtlas ):
            doLink = ( libraryName not in _AtlasCLibrariesToOmit )
        else:
            doLink = True
        return doLink

    @staticmethod
    def MaximumNumberOfThreads ( ):
        """Find the maximum number of threads."""
        try:
            import multiprocessing
            return multiprocessing.cpu_count ( )
        except:
            return 1

    # . Properties.
    @property
    def compilerFlags ( self ):
        item = self.__dict__.get ( "_compilerFlags", None )
        if item is None:
            if self.useOpenMP: item = [ _NumberOfThreadsFlag.format ( self.numberOfThreads ), _OpenMPPreprocessorFlag, _OpenMPFlag ]
            else: item = []
            setattr ( self, "_compilerFlags", item )
        return item
    @property
    def linkerFlags ( self ):
        item = self.__dict__.get ( "_linkerFlags", None )
        if item is None:
            if self.useOpenMP: item = [ _OpenMPFlag ]
            else: item = []
            setattr ( self, "_linkerFlags", item )
        return item
    @property
    def systemLibraries ( self ):
        item = self.__dict__.get ( "_systemLibraries", None )
        if item is None:
            if   self.useSerialAtlas   : item = _AtlasSerialLibraries
            elif self.useThreadedAtlas : item = _AtlasThreadedLibraries
            else                       : item = _SystemLibraries 
            setattr ( self, "_systemLibraries", item )
        return item
    @property
    def systemLibraryPaths ( self ):
        item = self.__dict__.get ( "_systemLibraryPaths", None )
        if item is None:
            if   self.useSerialAtlas   : item = _AtlasLibraryPaths
            elif self.useThreadedAtlas : item = _AtlasLibraryPaths
            else                       : item = _SystemLibraryPaths 
            setattr ( self, "_systemLibraryPaths", item )
        return item

#===================================================================================================================================
# . Directory install class.
#===================================================================================================================================
class PackageToInstall:
    """Class for a package that is to be installed."""

    def __init__ ( self, **options ):
        """Constructor."""
        for ( key, value ) in options.items ( ): setattr ( self, key, value )

    def CheckForCLibraries ( self ):
        """Check for C libraries."""
        # . Read the library file.
        items  = []
        target = os.path.join ( self.extensionsPath, _LibraryFile )
        if os.path.exists ( target ):
            tFile = open ( target, "r" )
            for line in tFile:
                tokens = line.split ( )
                if ( len ( tokens ) > 0 ) and ( not tokens[0].startswith ( "#" ) ):
                    directoryName    = tokens[0]
                    libraryName      = directoryName
                    isThirdParty     = False
                    optimizationFlag = None
                    if len ( tokens ) > 1:
                        libraryName      = tokens[1]
                        if len ( tokens ) > 2:
                            isThirdParty = ( "ThirdParty" in tokens[2:] )
                            if   "FullOptimization" in tokens[2:]: optimizationFlag = _FullOptimizationFlag
                            elif "NoOptimization"   in tokens[2:]: optimizationFlag = _NoOptimizationFlag
                    items.append ( ( directoryName, libraryName, isThirdParty, optimizationFlag ) )
            tFile.close ( )
        setattr ( self, "cLibraries", items )

    def CheckForDependencies ( self ):
        """Check for package dependencies."""
        items  = []
        target = os.path.join ( self.path, _ExtensionsDirectory, _DependencyFile )
        if os.path.exists ( target ):
            tFile = open ( target, "r" )
            for line in tFile:
                line = line.strip ( )
                if len ( line ) > 0: items.append ( line )
            tFile.close ( )
        setattr ( self, "dependencyNames", items )

    def CheckForExtensions ( self ):
        """Check for package extensions."""
        # . Extension directory.
        target = os.path.join ( self.path, _ExtensionsDirectory )
        if os.path.exists ( target ):
            setattr ( self, "extensionsPath", target )
            # . Check for a C include directory.
            target = os.path.join ( self.extensionsPath, _CIncludeDirectory )
            if os.path.exists ( target ): setattr ( self, "cIncludePath", target )
            # . Check for libraries.
            self.CheckForCLibraries ( )
            # . Check for a C library directory and make it if necessary - can be zero length.
            if len ( self.cLibraries ) >= 0:
                target = os.path.join ( self.extensionsPath, _CLibraryDirectory )
                if not os.path.exists ( target ): os.mkdir ( target )
                setattr ( self, "cLibraryPath" , target )
                setattr ( self, "hasCLibraries", True   )
            # . Check for a Pyrex directory with files.
            target = os.path.join ( self.extensionsPath, _PyrexDirectory )
            if os.path.exists ( target ):
                # . Get all "pyx" files no matter whether they begin with the package name or not.
                #   Store both the actual name and the module-qualified name.
                pyrexFiles = glob.glob ( os.path.join ( target, "*" + _PyrexExtension ) )
                if len ( pyrexFiles ) > 0:
                    setattr ( self, "pyrexPath", target )
                    # . Process the file names.
                    pyrexRoots = []
                    for pyrexFile in pyrexFiles:
                        ( head, tailext ) = os.path.split    ( pyrexFile )
                        ( tail, ext     ) = os.path.splitext ( tailext   )
                        if tail.find ( "." ) >= 0: tailM = tail
                        else:                      tailM = self.name + "." + tail
                        pyrexRoots.append ( ( tail, tailM ) )
                    pyrexRoots.sort ( )
                    setattr ( self, "pyrexFileRoots", pyrexRoots )
                    setattr ( self, "hasPyrexFiles",  True       )
                    # . Check for a psource directory.
                    target = os.path.join ( self.extensionsPath, _PSourceDirectory )
                    if not os.path.exists ( target ): os.mkdir ( target )
                    setattr ( self, "pSourcePath", target )

    def CheckPyrexFiles ( self ):
        """Check that the Pyrex files have been converted."""
        if getattr ( self, "hasPyrexFiles", False ):
            errors = []
            for ( name, nameM ) in self.pyrexFileRoots:
                pPath = os.path.join ( self.pyrexPath,   name  + _PyrexExtension )
                cPath = os.path.join ( self.pSourcePath, nameM + _CExtension     )
                if os.path.exists ( pPath ):
# . Could check for modification time here but no good if just unpacked distribution as these are not kept.
#                    isOK = os.path.exists ( cPath ) and ( os.path.getmtime ( cPath ) > os.path.getmtime ( pPath ) )
                    if not os.path.exists ( cPath ): errors.append ( name + _PyrexExtension )
            if len ( errors ) > 0:
                if self.options.verbose:
                    print ( "\nInvalid Pyrex-derived C files:\n" )
                    for name in errors: print ( name )
                    print ( "" )
                raise InstallationError ( "There are missing or out-of-date Pyrex-derived C files." )

    def CompileCLibraries ( self ):
        """Compile the C-libraries."""
        if getattr ( self, "hasCLibraries", False ):
            # . Get the build directory name.
            buildPath = BuildPath ( "temp" ) #os.path.join ( "build", "temp" + ".{:s}-cpython-{:s}".format ( get_platform ( ), VersionString ( ) ) )
            # . Get the include directories.
            includeDirectories = [ self.cIncludePath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "cIncludePath", None )
                if path is not None: includeDirectories.append ( path )
            # . Loop over libraries.
            for ( directoryName, libraryName, isThirdParty, optimizationFlag ) in self.cLibraries:
                # . Check to see whether this library is to be compiled.
                if self.options.IsCLibraryToBeCompiled ( libraryName, isThirdParty ):
                    # . Get the source file list.
                    sourceFiles = glob.glob ( os.path.join ( self.extensionsPath, directoryName, "*" + _CExtension ) )
                    if len ( sourceFiles ) > 0:
                        # . Make build_info.
                        build_info = { }
                        build_info["sources"      ] = sourceFiles
                        build_info["extra_preargs"] = list ( self.options.compilerFlags )
                        build_info["include_dirs" ] = includeDirectories
                        build_info["macros"       ] = None
                        if optimizationFlag is not None:
                            build_info["removeOptimizationFlags"] = True
                            build_info["extra_preargs"          ].append ( optimizationFlag )
                        # . Compile the library.
                        setup ( name        = libraryName ,
                                libraries   = [ ( libraryName, build_info )   ] ,
                                script_args = [ "BuildClibWithCompileOptions" ] ,
                                cmdclass    = { "BuildClibWithCompileOptions" : BuildClibWithCompileOptions } )
                        self.report["C Files"    ] += len ( sourceFiles )
                        self.report["C Libraries"] += 1
                        # . Move the library to the appropriate place.
                        os.rename ( os.path.join ( buildPath, "lib" + libraryName + _LibraryExtension ), os.path.join ( self.cLibraryPath, "lib" + libraryName + _LibraryExtension ) )

    def ConvertPyrexToC ( self ):
        """Convert Pyrex to C files."""
         # . Check for Pyrex files.
        if getattr ( self, "hasPyrexFiles", False ):
            # . Import the correct modules.
            tag = "Cython"
            try:
                import Cython.Compiler.Errors as Errors, Cython.Compiler.Main as Main, Cython.Compiler.Version as Version
            except:
                raise InstallationError ( "Error locating the " + tag + " compiler." )
            # . Compile.
            # . Get the pxdDirectories (self and then dependencies in reverse order).
            pxdDirectories = [ self.pyrexPath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "pyrexPath", None )
                if path is not None: pxdDirectories.append ( path )
            # . Convert the files.
            if self.options.verbose: print ( "\nConverting files in " + self.pyrexPath + " with " + tag + " version " + Version.version + ":\n" )
            for ( name, nameM ) in self.pyrexFileRoots:
                if self.options.verbose: print ( " -> " + name + _PyrexExtension )
                PyrexCompile ( Errors, Main, os.path.join ( self.pyrexPath, name + _PyrexExtension ), pxdDirectories = pxdDirectories )
                os.rename ( os.path.join ( self.pyrexPath, name + _CExtension ), os.path.join ( self.pSourcePath, nameM + _CExtension ) )
                self.report["Pyrex Files"] += 1
            if self.options.verbose: print ( "" )

    @classmethod
    def FromPathName ( selfclass, name, path ):
        """Constructor from path name."""
        return selfclass ( name = name, path = path )

    def GetCLibraryNames ( self, options ):
        """Return a list of C library names in reverse order."""
        names = []
        for data in reversed ( getattr ( self, "cLibraries", [] ) ):
            name = data[1]
            if options.IsCLibraryToBeLinked ( name ): names.append ( name )
        return names

    def HasExtensions ( self ):
        """Does this package have extensions?"""
        return ( getattr ( self, "hasCLibraries", False ) or getattr ( self, "hasPyrexFiles", False ) )

    def MakeExtensions ( self ):
        """Make extensions."""
        if getattr ( self, "hasPyrexFiles", False ):
            # . Get the build and destination path names.
            nameAsPath      = self.name.replace ( ".", os.sep )
            buildPath       = BuildPath ( "lib" ) #os.path.join ( "build", "lib" + ".{:s}-cpython-{:s}".format ( get_platform ( ), VersionString ( ) ), nameAsPath )
            destinationPath = self.path
            # . Get the include directories.
            includeDirectories = [ self.cIncludePath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "cIncludePath", None )
                if path is not None: includeDirectories.append ( path )
            # . Get the library directories.
            libraryDirectories = [ self.cLibraryPath ]
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                path = getattr ( item, "cLibraryPath", None )
                if path is not None: libraryDirectories.append ( path )
            libraryDirectories.extend ( self.options.systemLibraryPaths )
            # . Get the library names.
            cLibraries = self.GetCLibraryNames ( self.options )
            for item in reversed ( getattr ( self, "dependencyObjects", [] ) ):
                cLibraries.extend ( item.GetCLibraryNames ( self.options ) )
            cLibraries.extend ( self.options.systemLibraries )
            # . Compile and link arguments.
            compileArguments = list ( self.options.compilerFlags )
            linkArguments    = list ( self.options.linkerFlags   )
            # . Make the list of extension modules.
            extensions = []
            for ( name, nameM ) in self.pyrexFileRoots:
                # . The full extension name is needed here with package qualifiers.
                # . We use the same name for the psource file to ensure that there
                # . are not two completely different C files with the same name.
                extensions.append ( Extension ( nameM, [ os.path.join ( self.pSourcePath, nameM + _CExtension ) ] ,
                                                                        extra_compile_args   = compileArguments   ,
                                                                        extra_link_args      = linkArguments      ,
                                                                        include_dirs         = includeDirectories ,
                                                                        libraries            = cLibraries         ,
                                                                        library_dirs         = libraryDirectories ,
                                                                        runtime_library_dirs = libraryDirectories ) )
            # . Compile the extension modules.
            setup ( name = nameAsPath, ext_modules = extensions, script_args = [ "build_ext" ] )
            self.report["Extension Files"] += len ( extensions )
            # . Move the files to the appropriate place.
            # . Here the unqualified name is required.
            for ( name, nameM ) in self.pyrexFileRoots:
                ( head, tail ) = name.rsplit ( ".", 1 ) # . How robust is this?
                head = head.replace ( ".", os.sep )
#                if name.find ( "." ) >= 0: 
#                else:                      ( head, tail ) = ( "", name )
#                print ( "XXX>", buildPath, head, tail, name, nameM )
                paths = glob.glob ( os.path.join ( buildPath, head, tail + _SharedObjectBuildPathExtension ) )
                if len ( paths ) != 1: raise InstallationError ( "Wrong number of shared object paths ({:d}) for {:s}.".format ( len ( paths ), tail ) )
                os.rename ( paths[0], os.path.join ( destinationPath, tail + _SharedObjectExtension ) )

    def ProcessWithOptions ( self, options, report ):
        """Process the package given a set of installation options."""
        # . Reporting.
        report["Packages"] += 1
        # . Assign options and report.
        self.options = options
        self.report  = report
        # . C libraries.
        if options.doCLibraries: self.CompileCLibraries ( )
        # . Handle pyrex.
        if options.doPyrex: self.ConvertPyrexToC ( )
        # . Create the extension modules.
        if options.doExtensions:
            if not options.doPyrex: self.CheckPyrexFiles ( )
            self.MakeExtensions ( )

    def ResolveDependencies ( self, itemDictionary ):
        """Resolve the package dependencies."""
        names = getattr ( self, "dependencyNames", [] )
        items = []
        for name in names:
            try:    items.append ( itemDictionary[name] )
            except: raise InstallationError ( "Unable to resolve dependency \"" + name + "\" for package " + self.name + "." )
        setattr ( self, "dependencyObjects", items )

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def BuildPath ( name ):
    """Get the build path."""
    return os.path.join ( "build", "{:s}.{:s}-{:s}".format ( name, sysconfig.get_platform ( ), sys.implementation.cache_tag ) )

def FindRootDirectory ( ):
    """Find or guess a root directory."""
    # . Initialization.
    rootDirectory = None
    # . Guess a value from the current directory.
    ( head, tail ) = os.path.split ( os.getcwd ( ) )
    # . If this is the installation directory assume everything is OK.
    if tail == "installation":
        rootDirectory = head
    # . Look for a directory that has "pdynamo3" in it.
    else:
        while len ( tail ) > 0:
            if "pdynamo3" in tail.lower ( ):
                rootDirectory = os.path.join ( head, tail )
                break
            ( head, tail ) = os.path.split ( head )
        # . Use the environment variable in case of failure.
        if rootDirectory is None:
            rootDirectory = os.getenv ( _PDynamoHome )
    # . Finish up.
    if rootDirectory is None:
        raise InstallationError ( "Unable to find pDynamo root directory." )
    else:
        os.environ[_PDynamoHome] = rootDirectory
    return rootDirectory

def GetRootDirectories ( names, update = False ):
    """Get the root directories recursively in increasing order of dependency.
    All directories are returned by default otherwise only those specified by |names|.

    A root directory is one that has both an "__extensions__" directory and an "__init__.py" file.
    """
    # . Get the root directory name.
    rootDirectory  = os.getenv ( _PDynamoHome )
    thirdPartyRoot = os.path.join ( rootDirectory, _ThirdPartyDirectory )
    # . Get all directories recursively in root.
    toProcess = {}
    for candidate in glob.iglob ( os.path.join ( rootDirectory, "**", _ExtensionsDirectory ), recursive = True ):
        if os.path.isdir ( candidate ):
            # . Candidate is of form "root|packagePath|__extensions__" so extract packagePath and replace "/" by ".".
            ( head, tail ) = os.path.split ( candidate )
            # . Check for an "__init__.py" file.
            if os.path.exists ( os.path.join ( head, _InitFile ) ):
                name = head[len(rootDirectory)+len(os.sep):].replace ( os.sep, "." )
                directory = PackageToInstall.FromPathName ( name = name, path = head )
                directory.CheckForDependencies ( )
                directory.CheckForExtensions   ( )
                toProcess[directory.name] = directory
    # . Satisfy the dependencies of the directories.
    for directory in toProcess.values ( ):
        directory.ResolveDependencies ( toProcess )
    # . Create the sorted directory list.
    directories = SortPackages ( toProcess, names, update = update )
    # . Finish up.
    return directories

def InstallationSummary ( report, totalTime, message ):
    """Summarize the results of the installation."""
    print ( "\n---------------------\nInstallation Summary:\n---------------------" )
    if len ( report ) > 0:
        for key in sorted ( report.keys ( ) ):
            print ( "{:<17s} {:5d}".format ( key, report[key] ) )
    print ( "Processing Time   {:s}".format ( TimeToString ( totalTime ) ) )
    print ( "\n{:s}".format ( message ) )

def PyrexCompile ( Errors, Main, source, pxdDirectories = None ):
    """Convert a Pyrex file to C."""
    options                   = Main.default_options
    options["language_level"] = _LanguageLevel
    # . Depending on the Cython or Pyrex version options can be a dictionary or a compilation object.
    if isinstance ( options, dict ): options = Main.CompilationOptions ( options )
    if pxdDirectories is not None: options.include_path.extend ( pxdDirectories )
    context = Main.Context ( options.include_path, {} )
    try:
        result = Main.compile ( source, options )
        failed = ( result.num_errors > 0 )
    except Errors.PyrexError as e:
        print ( e )
        failed = True
    if failed: raise InstallationError ( "There was a Pyrex compiler error." )

def RemoveBuildDirectory ( ):
    """Remove the build directory if it exists."""
    if os.path.exists ( "build" ): shutil.rmtree ( "build" )

def SortPackages ( toProcess, names, update = False ):
    """Sort the package list in order of dependency."""
    # . This could be improved by using the same graph structure for the DFS and topological sort.
    # . Possibilities (always ordered):
    #   Names Update Result
    #   Empty False  All packages
    #   Empty True   All packages
    #   True  False  Named packages
    #   True  True   Named packages + those that depend on them
    # . toProcess is a dictionary of package names : package instances.
    # . Generate the graph.
    doNames = ( names is not None ) and ( len ( names ) > 0 )
    if doNames: # . Remove packages from toProcess that are not reachable from the named packages.
        graph = { key : set ( ) for key in toProcess.keys ( ) }
        for ( key, value ) in toProcess.items ( ):
            for name in value.dependencyNames: graph[name].add ( key )
        visited = set ( )
        for name in names:
            if name in toProcess.keys ( ): # . To avoid packages which do not require treatment.
                stack = set ( [ name ] ) 
                while len ( stack ) > 0:
                    vertex = stack.pop ( )
                    if vertex not in visited:
                        visited.add ( vertex )
                        stack.update ( graph[vertex] )
        graph = { name : set ( [ dependency for dependency in toProcess[name].dependencyNames if dependency in visited ] ) for name in visited }
    else: # . Full graph.
        graph = { key : set ( value.dependencyNames ) for ( key, value ) in toProcess.items ( ) }
    # . Order the graph.
    ordered = list ( TopologicalSort ( graph ) )
    # . Prune the list with explicitly specified packages only.
    # . The update option keeps all descendants.
    if doNames and ( not update ):
        common  = [ group.intersection ( names ) for group in ordered ]
        ordered = [ group for group in common if len ( group ) > 0    ]
    #print ( ordered )
    #sys.exit ( )
    # . Finish up.
    return [ toProcess[name] for group in ordered for name in sorted ( group ) ]

def TimeToString ( time ):
    """Convert a floating point time to a string."""
    value  = time
    fields = []
    for ( size, tag ) in ( ( 60*60*24, "d" ), ( 60*60, "h" ), ( 60, "m" ) ):
        ( newf, value ) = divmod ( value, size )
        newi = int ( round ( newf ) )
        if not ( ( newi == 0 ) and ( len ( fields ) == 0 ) ):
            fields.append ( "{:2d}{:1s}".format ( newi, tag ) )
    tag = "s"
    fields.append ( "{:5.3f}{:1s}".format ( value, tag ) )
    return " ".join ( fields )

def TopologicalSort ( data ):
    """Topologically sort the dependency graph."""
    # . |data| is a dictionary whose keys are package names, and values are sets of dependencies.
    # . Adapted from toposort-1.5.
    # . Nothing to be done.
    if len ( data ) == 0: return
    # . Copy the input so as to leave it unmodified.
    data = data.copy()
    # . Ignore self dependencies.
    for ( key, value ) in data.items ( ): value.discard ( key )
    # . Find all items that don't depend on anything.
    extraItems = functools.reduce ( set.union, data.values ( ) ) - set ( data.keys ( ) )
    # . Add empty dependences where needed.
    data.update ( { item : set ( ) for item in extraItems } )
    # . Do the sorting.
    while True:
        ordered = set ( item for item, dep in data.items ( ) if len ( dep ) == 0 )
        if not ordered:
            break
        yield ordered
        data = { item: ( dep - ordered ) for item, dep in data.items ( ) if item not in ordered }
    # . There are circular dependencies.
    if len ( data ) != 0:
        width = max ( [ len ( key ) for key in data.keys ( ) ] )
        print ( "\nPackage circular dependencies:" )
        for key in sorted ( data.keys ( ) ): print ( "{:s} - {:s}", key.ljust ( width ), ", ".join ( sorted ( data[key] ) ) )
        raise InstallationError ( "The packages have circular dependencies." )

#def VersionString ( ):
#    """The current system version string for build paths."""
#    return sys.version.split ( " ", 1 )[0].rsplit ( ".", 1 )[0].replace ( ".", "" )

def WriteShellFile ( inPath, outPath, variables, report ):
    """Write a shell file."""
    inFile   = open ( inPath, "r" )
    inString = inFile.read ( ).format ( **variables )
    inFile.close  ( )
    outFile  = open ( outPath, "w" )
    outFile.write ( inString )
    outFile.close ( )
    report["Shell Files"] += 1

def WriteShellFiles ( report, parameters = None, scratch = None ):
    """Write shell files containing the environment variables needed by pDynamo."""
    # . Get all required paths.
    # . Root.
    pDynamoHome = os.getenv ( _PDynamoHome )
    # . Parameters and scratch.
    if parameters is None: parameters = "$PDYNAMO3_HOME/parameters"
    if scratch    is None: scratch    = "$PDYNAMO3_HOME/scratch"
    # . Get the new Python path.
    pythonPath = os.getenv ( "PYTHONPATH" )
    if pythonPath is None: pythonPath = "."
    tokens     = [ token for token in pythonPath.split ( ":" ) if token.upper ( ).find ( "PDYNAMO" ) < 0 ]
    tokens.append ( "$PDYNAMO3_HOME" )
    pythonPath = ":".join ( tokens )
    # . Write the shell files.
    variables       = { "orcaCommand"   : "$HOME/programs/orca_local/orca" , # . A guess.
                        "parameters"    : parameters  ,
                        "pDynamoHome"   : pDynamoHome ,
                        "pythonCommand" : "python3"   ,
                        "pythonPath"    : pythonPath  ,
                        "scratch"       : scratch     }
    shellScriptPath = os.path.join ( pDynamoHome, "installation", _ShellScriptPath )
    templatePath    = os.path.join ( pDynamoHome, "installation", _TemplatePath    )
    WriteShellFile ( os.path.join ( templatePath, "environment_bash.tmpl"   ), os.path.join ( shellScriptPath, "environment_bash.com"   ), variables, report )
    WriteShellFile ( os.path.join ( templatePath, "environment_cshell.tmpl" ), os.path.join ( shellScriptPath, "environment_cshell.com" ), variables, report )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

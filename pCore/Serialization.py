"""Serialization utilities."""

import glob, os, yaml

from pickle import dump, load

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
PickleFileExtension     = ".pkl"
YAMLPickleFileExtension = ".yaml"

#===================================================================================================================================
# . Standard pickling.
#===================================================================================================================================
def Pickle ( path, items ):
    """Pickle items to a file."""
    pFile = open ( path, "wb" )
    dump ( items, pFile )
    pFile.close ( )

def Unpickle ( path ):
    """Unpickle objects from a file."""
    pFile = open ( path, "rb" )
    items = load ( pFile )
    pFile.close ( )
    return items

#===================================================================================================================================
# . Raw object construction.
#===================================================================================================================================
def RawObjectConstructor ( selfClass ):
    """Return an empty, initialized object."""
    return selfClass.Raw ( )

#===================================================================================================================================
# . YAML pickling.
#===================================================================================================================================
def YAMLPickle ( path, items, defaultFlowStyle = True, title = None ):
    """Pickle items to a file."""
    pFile = open ( path, "w" )
    if title is not None: pFile.write ( "# . " + title + "\n" )
    yaml.dump ( items, pFile, default_flow_style = defaultFlowStyle )
    pFile.close ( )

def YAMLUnpickle ( path ):
    """Unpickle objects from a file."""
    pFile = open ( path, "r" )
    items = yaml.load ( pFile )
    pFile.close ( )
    return items

#===================================================================================================================================
# . YAML functions to dump and load object state.
#===================================================================================================================================
def YAMLMappingFile_FromObject ( path, title, object, defaultFlowStyle = True ):
    """Create a mapping file from an object."""
    pFile = open ( path, "w" )
    if title is not None: pFile.write ( "# . " + title + "\n" )
    pFile.write ( "---\n" )
    yaml.dump ( object.__getstate__ ( ), pFile, default_flow_style = defaultFlowStyle ,
                                                explicit_start     = False            ,
                                                explicit_end       = False            )
    pFile.write ( "...\n" )
    pFile.close ( )

def YAMLMappingFile_ToObject ( path, selfClass ):
    """Create an object from a mapping file."""
    pFile = open ( path, "r" )
    state = yaml.load ( pFile )
    pFile.close ( )
    self  = selfClass.Raw ( )
    self.__setstate__ ( state )
    return self

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

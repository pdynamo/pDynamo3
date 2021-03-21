"""Space group classes and functions."""

import glob, math, os, os.path

from  pCore         import AttributableObject        , \
                           YAMLUnpickle
from .CrystalSystem import CrystalSystemCubic        , \
                           CrystalSystemHexagonal    , \
                           CrystalSystemMonoclinic   , \
                           CrystalSystemOrthorhombic , \
                           CrystalSystemTetragonal   , \
                           CrystalSystemTriclinic    , \
                           CrystalSystemTrigonal
from .SymmetryError import SymmetryError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SpaceGroup ( AttributableObject ):
    """Space group class."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "crystalSystem"    : None ,
                             "hall"             : None ,
                             "hermannMauguin"   : None ,
                             "schoenflies"      : None ,
                             "spaceGroupNumber" : 0    ,
                             "spaceGroupSymbol" : None ,
                             "setting"          : 0    ,
                             "transformations"  : None } )

    @classmethod
    def FromYAML ( selfClass, path ):
        """Constructor from a YAML file."""
        data = YAMLUnpickle ( path )
        self = selfClass.WithOptions ( crystalSystem    = data["Crystal System"    ] ,
                                       hall             = data["Hall"              ] ,
                                       hermannMauguin   = data["Hermann-Mauguin"   ] ,
                                       schoenflies      = data["Schoenflies"       ] ,
                                       spaceGroupNumber = data["Space Group Number"] ,
                                       spaceGroupSymbol = data["Space Group Symbol"] ,
                                       setting          = data["Setting"           ] ,
                                       transformations  = data["Transformations"   ] )
        return self

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def SpaceGroup_CrystalSystemFromNumber ( number ):
    """Return a crystal class given a space group number."""
    if   ( number <    1 ) or  ( number >  230 ): raise SymmetryError ( "Invalid space group number: {:d}.".format ( number ) )
    elif ( number ==   1 ) or  ( number ==   2 ): return CrystalSystemTriclinic    ( )
    elif ( number >=   3 ) and ( number <=  15 ): return CrystalSystemMonoclinic   ( )
    elif ( number >=  16 ) and ( number <=  74 ): return CrystalSystemOrthorhombic ( )
    elif ( number >=  75 ) and ( number <= 142 ): return CrystalSystemTetragonal   ( )
    elif ( number >= 143 ) and ( number <= 167 ): return CrystalSystemTrigonal     ( )
    elif ( number >= 168 ) and ( number <= 194 ): return CrystalSystemHexagonal    ( )
    elif ( number >= 195 ) and ( number <= 230 ): return CrystalSystemCubic        ( )

def SpaceGroup_FromYAML ( name, inPath = None ):
    """Read a space group with the given Hermann-Mauguin name from a directory containing YAML files."""
    if inPath is None:
        try:    inPath = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "spaceGroups" )
        except: raise ( "Unable to find space group directory." )
    index         = YAMLUnpickle ( os.path.join ( inPath, "__index__.yaml" ) )
    settingNumber = index.get ( name, None )
    if settingNumber is None: raise SymmetryError ( "Space group \"{:s}\" not found.".format ( name ) )
    return SpaceGroup.FromYAML ( os.path.join ( inPath, "{:03d}.yaml".format ( settingNumber ) ) )

def SpaceGroups_FromYAML ( inPath = None ):
    """Read space groups from a directory containing YAML files.

    Return as a dictionary with Hermann-Mauguin names as keys.
    """
    if inPath is None:
        try:    inPath = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "spaceGroups" )
        except: raise SymmetryError ( "Unable to find space group directory." )
    spaceGroups = {}
    paths = glob.glob ( os.path.join ( inPath, "*.yaml" ) )
    for path in paths:
        if path.find ( "ReadMe" ) < 0:
            group = SpaceGroup.FromYAML ( path )
            key   = group.hermannMauguin
            spaceGroups[key] = group
    return spaceGroups

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

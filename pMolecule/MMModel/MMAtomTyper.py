"""MM atom typing."""

import os

from  pCore               import Align                    , \
                                 AttributableObject       , \
                                 logFile                  , \
                                 LogFileActive            , \
                                 Selection                , \
                                 YAMLMappingFile_ToObject , \
                                 YAMLPickleFileExtension
from .MMAtomTypeContainer import MMAtomTypeContainer
from .MMModelError        import MMModelError
from .MMPattern           import MMPatternContainer
from .MMSequenceLibrary   import MMSequenceLibrary

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Maximum width for table output.
_MaximumTableWidth = 100

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMAtomTyper ( AttributableObject ):
    """Type atoms."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "mmAtomTypes"       : None ,
                             "mmPatterns"        : None ,
                             "mmSequenceLibrary" : None ,
                             "path"              : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( MMAtomTyper, self )._CheckOptions ( )
        self.ReadData ( )

    def CheckUntypedAtoms ( self, connectivity, untypedAtoms, log ):
        """Check for untyped atoms."""
        # . Check for untyped atoms.
        if len ( untypedAtoms ) > 0:
            if LogFileActive ( log ):
                labels = []
                length = 0
                for atom in untypedAtoms:
                    label  = atom.path
                    labels.append ( label )
                    length = max ( length, len ( label ) )
                length  += 2
                number   = len ( labels )
                title    = "Untyped Atoms"
                width    = min ( _MaximumTableWidth, max ( number * length, len ( title ) + 3 ) )
                ncolumns = ( width + length - 1 ) // length
                table    = log.GetTable ( columns = ncolumns * [ length ] )
                table.Start ( )
                table.Title ( title )
                for label in labels: table.Entry ( label, align = Align.Left )
                table.Stop ( )
            raise MMModelError ( "There are {:d} untyped atoms.".format ( len ( untypedAtoms ) ), untypedAtoms )

    @classmethod
    def FromPath ( selfClass, path, **options ):
        """Constructor from path and other options."""
        options         = dict ( options )
        options["path"] = path
        return selfClass.WithOptions ( **options )

    def ReadData ( self ):
        """Read the atom typing data."""
        if self.path is not None:
            # . Atom types and patterns.
            for ( attribute, tag, pathClass ) in ( ( "mmAtomTypes", "atomTypes", MMAtomTypeContainer ), ( "mmPatterns", "patterns", MMPatternContainer ) ):
                path = os.path.join ( self.path, tag + YAMLPickleFileExtension )
                if os.path.exists ( path ): setattr ( self, attribute, YAMLMappingFile_ToObject ( path, pathClass ) )
            if self.mmPatterns is not None: self.mmPatterns.IndexAtomTypes ( self.mmAtomTypes )
            # . Sequence attributes.
            self.mmSequenceLibrary = MMSequenceLibrary.TestForLibrary ( self.path )

    def TypeAtoms ( self, connectivity, sequence, log ):
        """Assign atom types and charges to the atoms."""
        # . Initialization.
        atomCharges  = [ 0.0  for i in range ( len ( connectivity.atoms ) ) ]
        atomTypes    = [ None for i in range ( len ( connectivity.atoms ) ) ]
        untypedAtoms = set ( connectivity.atoms )
        # . Type the atoms.
        self.TypeBySequence ( sequence    , atomTypes, atomCharges, untypedAtoms )
        self.TypeByPattern  ( connectivity, atomTypes, atomCharges, untypedAtoms )
        # . Check for untyped atoms.
        self.CheckUntypedAtoms ( connectivity, untypedAtoms, log )
        # . Finish up.
        return ( atomTypes, atomCharges )

    def TypeByPattern ( self, connectivity, atomTypes, atomCharges, untypedAtoms ):
        """Type the atoms by pattern."""
        if ( self.mmAtomTypes is not None ) and ( self.mmPatterns is not None ):
            self.mmPatterns.TypeConnectivity ( connectivity, atomTypes, atomCharges, untypedAtoms )

    def TypeBySequence ( self, sequence, atomTypes, atomCharges, untypedAtoms ):
        """Type the atoms by sequence."""
        if self.mmSequenceLibrary is not None:
            self.mmSequenceLibrary.TypeSequence ( sequence, atomTypes, atomCharges, untypedAtoms )
 
#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

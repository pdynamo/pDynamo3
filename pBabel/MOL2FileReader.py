#===================================================================================================================================
# . Classes and functions to read MOL2 files.
#===================================================================================================================================

from  pCore                 import logFile                  , \
                                   LogFileActive            , \
                                   TextFileReader
from  pMolecule             import Atom                     , \
                                   Bond                     , \
                                   BondType                 , \
                                   Connectivity             , \
                                   ConvertInputConnectivity , \
                                   System
from  pScientific           import PeriodicTable
from  pScientific.Arrays    import Array
from  pScientific.Geometry3 import Coordinates3
from  pScientific.Symmetry  import PeriodicBoundaryConditions         , \
                                   SpaceGroup_CrystalSystemFromNumber , \
                                   SymmetryParameters
from .ExportImport          import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Bond type definitions.
_MOL2BondTypes = { "1"  : ( BondType.Single    , False ) ,
                   "2"  : ( BondType.Double    , False ) ,
                   "3"  : ( BondType.Triple    , False ) ,
                   "am" : ( BondType.Single    , False ) ,
                   "ar" : ( BondType.Undefined , True  ) ,
                   "du" : ( BondType.Undefined , False ) ,
                   "un" : ( BondType.Undefined , False ) ,
                   "nc" : ( None               , False ) }

# . Generic (non-element) atom types.
_GenericAtomTypes = ( "Any", "Du", "Hal", "Het", "Hev", "LP" )

# . Section header.
_RTIString  = "@<TRIPOS>"

#===================================================================================================================================
# . MOL2 file reader class.
#===================================================================================================================================
class MOL2FileReader ( TextFileReader ):
    """MOL2FileReader is the class for MOL2 files that are to be read."""

    _classLabel = "MOL2 File Reader"

    def AtomicNumberFromAtomType ( self, atomType ):
        """Get an atomic number from atom type."""
        atomicNumber = -1
        token        = atomType.split ( ".", 1 )[0]
        if token not in _GenericAtomTypes: atomicNumber = PeriodicTable.AtomicNumber ( token )
        return atomicNumber

    def GetLine ( self, signalWarnings = False ):
        """Get a line."""
        try:
            while True:
                line = next ( self.file ).strip ( )
                self.linesParsed += 1
                if ( len ( line ) > 0 ) and ( not line.startswith ( "#" ) ): break
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Parse all entries.
                while True:
                    line = self.GetLine ( )
                    if   line.startswith ( _RTIString + "ATOM"     ): self.ParseAtomSection     ( )
                    elif line.startswith ( _RTIString + "BOND"     ): self.ParseBondSection     ( )
                    elif line.startswith ( _RTIString + "CRYSIN"   ): self.ParseCrysinSection   ( )
                    elif line.startswith ( _RTIString + "MOLECULE" ): self.ParseMoleculeSection ( )
            except EOFError:
                pass
            # . Complete the connectivity.
            ConvertInputConnectivity ( self.connectivity, {} )
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseAtomSection ( self ):
        """Parse the ATOM section."""
        if hasattr ( self, "numberOfAtoms" ):
            self.connectivity = Connectivity ( )
            atomicCharges = Array.WithExtent        ( self.numberOfAtoms ) ; atomicCharges.Set ( 0.0 )
            xyz           = Coordinates3.WithExtent ( self.numberOfAtoms ) ; xyz.Set           ( 0.0 )
            for i in range ( self.numberOfAtoms ):
                items = self.GetTokens ( converters = [ int, None, float, float, float, None, None, None, float ] )
                if len ( items ) < 6:
                    self.Warning ( "Invalid ATOM line.", True )
                else:
                    atomicNumber = PeriodicTable.AtomicNumber ( items[1] )
                    if atomicNumber <= 0: atomicNumber = self.AtomicNumberFromAtomType ( items[5] )
                    self.connectivity.AddNode ( Atom.WithOptions ( atomicNumber = atomicNumber, label = items[1] ) )
                    xyz[i,0] = items[2]
                    xyz[i,1] = items[3]
                    xyz[i,2] = items[4]
                    if len ( items ) >= 9: atomicCharges[i] = items[8]
            self.atomicCharges = atomicCharges
            self.xyz           = xyz
        else:
            self.Warning ( "Unknown number of atoms in molecule.", True )

    def ParseBondSection ( self ):
        """Parse the BOND section."""
        if hasattr ( self, "numberOfBonds" ) and hasattr ( self, "connectivity" ):
            numberOfAtoms = getattr ( self, "numberOfAtoms", -1 )
            for i in range ( self.numberOfBonds ):
                items = self.GetTokens ( converters = [ int, int, int, None ] )
                if len ( items ) < 4:
                    self.Warning ( "Invalid BOND line.", True )
                else:
                    atom1 = items[1] - 1
                    atom2 = items[2] - 1
                    ( bondType, bondIsAromatic ) = _MOL2BondTypes.get ( items[3], ( BondType.Undefined, False ) )
                    if ( atom1 < 0 ) or ( atom1 >= self.numberOfAtoms ) or ( atom2 < 0 ) or ( atom2 >= self.numberOfAtoms ):
                        self.Warning ( "Bond atom indices out of range: {:d}, {:d}.".format ( atom1, atom2 ), True )
                    if bondType is not None:
                        self.connectivity.AddEdge ( Bond.WithNodes ( self.connectivity.nodes[atom1], self.connectivity.nodes[atom2], isAromatic = bondIsAromatic, type = bondType ) )
        else:
            self.Warning ( "Unknown number of bonds in molecule.", True )

    def ParseCrysinSection ( self ):
        """Parse the CRYSIN section."""
        # . Items are a, b, c, alpha, beta, gamma, space group number and crystal setting.
        items  = self.GetTokens ( converters = 6 * [ float ] + 2 * [ int ] )
        try:
            setting    = items.pop ( -1 )
            spaceGroup = items.pop ( -1 )
            self.crystalSystem      = SpaceGroup_CrystalSystemFromNumber ( spaceGroup )
            self.symmetryParameters = SymmetryParameters ( )
            self.symmetryParameters.SetCrystalParameters ( *items )
        except:
            self.Warning ( "Invalid CRYSIN line.", False )

    def ParseMoleculeSection ( self ):
        """Parse the MOLECULE section."""
        # . Just parse the first two lines for the moment.
        self.label = self.GetLine ( )                             # . Molecule name.
        items      = self.GetTokens ( converters = [ int, int ] ) # . Number of atoms and bonds.
        if len ( items ) > 0: self.numberOfAtoms = items[0]
        if len ( items ) > 1: self.numberOfBonds = items[1]

    @classmethod
    def PathToAtomNames ( selfClass, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToAtomNames ( )

    @classmethod
    def PathToCharges ( selfClass, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCharges ( )

    @classmethod
    def PathToCoordinates3 ( selfClass, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToSystem ( )

    def ToAtomNames ( self ):
        """Return atom names."""
        atomNames = None
        if self.isParsed and hasattr ( self, "connectivity" ):
            atomNames = [ atom.label for atom in self.connectivity.nodes ]
        return atomNames

    def ToCharges ( self ):
        """Return charges."""
        charges = None
        if self.isParsed: return getattr ( self, "atomicCharges", None )
        return charges

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        xyz = None
        if self.isParsed: return getattr ( self, "xyz", None )
        return xyz

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.isParsed:
            try:
                system              = System.FromConnectivity ( self.connectivity )
                system.label        = self.label
                system.coordinates3 = self.xyz
                if hasattr ( self, "crystalSystem" ) and hasattr ( self, "symmetryParameters" ):
                    parameters                = self.crystalSystem.GetUniqueSymmetryParameters ( self.symmetryParameters )
                    system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( self.crystalSystem )
                    system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )
            except:
                pass
        return system

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : MOL2FileReader.PathToCoordinates3 ,
                         System       : MOL2FileReader.PathToSystem       } ,
                       [ "mol2", "MOL2" ], "Tripos MOL2", defaultFunction = MOL2FileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

#===================================================================================================================================
# . Classes and functions to read MOL2 files.
#===================================================================================================================================

# . MOL2 files do not contain some information that is vital for setting up MM systems using pDynamo.
#   In particular the formal charges of the atoms are missing which means that many atomic derived
#   properties, such as geometry and oxidation state, will be meaningless.

from  pCore                 import DataType            , \
                                   logFile             , \
                                   LogFileActive       , \
                                   TextFileReader      , \
                                   TextFileReaderError
from  pMolecule             import Atom                , \
                                   Bond                , \
                                   BondType            , \
                                   Connectivity        , \
                                   Sequence            , \
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
# . Atom type definitions giving atomic numbers, geometries and formal charges.
_AtomTypes = { "H.spc" : (  1, "Ter" , 0 ) ,
               "H.t3p" : (  1, "Ter" , 0 ) ,
               "C.1"   : (  6, "Lin" , 0 ) ,
               "C.2"   : (  6, "Tri" , 0 ) ,
               "C.3"   : (  6, "Tet" , 0 ) ,
               "C.ar"  : (  6, "Res" , 0 ) ,
               "C.cat" : (  6, "Tri" , 1 ) ,
               "Du.C"  : (  6,  None , 0 ) ,
               "N.1"   : (  7, "Lin" , 0 ) ,
               "N.2"   : (  7, "Tri" , 0 ) ,
               "N.3"   : (  7, "Tet" , 0 ) ,
               "N.4"   : (  7, "Tet" , 0 ) ,
               "N.am"  : (  7, "Tri" , 0 ) ,
               "N.ar"  : (  7, "Res" , 0 ) ,
               "N.pl3" : (  7, "Tri" , 0 ) ,
               "O.2"   : (  8, "Tri" , 0 ) ,
               "O.3"   : (  8, "Tet" , 0 ) ,
               "O.co2" : (  8, "Tri" , 0 ) ,
               "O.spc" : (  8, "Tet" , 0 ) ,
               "O.t3p" : (  8, "Tet" , 0 ) ,
               "P.3"   : ( 15, "Tet" , 0 ) ,
               "P.3+"  : ( 15, "Tet" , 1 ) ,
               "S.2"   : ( 16, "Tri" , 0 ) ,
               "S.3"   : ( 16, "Tet" , 0 ) ,
               "S.3+"  : ( 16, "Tet" , 1 ) ,
               "S.o"   : ( 16, "Tet" , 0 ) ,
               "S.o2"  : ( 16, "Tet" , 0 ) ,
               "Cr.oh" : ( 24, "Oct" , 0 ) ,
               "Cr.th" : ( 24, "Tet" , 0 ) ,
               "Co.oh" : ( 27, "Oct" , 0 ) }

# . Bond type definitions.
_BondTypes = { "1"  : ( BondType.Single    , False ) ,
               "2"  : ( BondType.Double    , False ) ,
               "3"  : ( BondType.Triple    , False ) ,
               "am" : ( BondType.Single    , False ) ,
               "ar" : ( BondType.Undefined , True  ) ,
               "du" : ( BondType.Undefined , False ) ,
               "un" : ( BondType.Undefined , False ) ,
               "nc" : ( None               , False ) }

# . Generic (non-element) atom types.
_GenericAtomTypes = ( "Any", "Du", "Hal", "Het", "Hev", "LP" )

# . Null string.
_NullString = "****"

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
        ( atomicNumber, geometry, formalCharge ) = _AtomTypes.get ( atomType, ( -1, None, 0 ) )
        if atomicNumber == -1:
            token = atomType.split ( "." )[0]
            if token not in _GenericAtomTypes: atomicNumber = PeriodicTable.AtomicNumber ( token )
        return ( atomicNumber, geometry, formalCharge )

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
                    if   line.startswith ( _RTIString + "ATOM"         ): self.ParseAtomSection         ( )
                    elif line.startswith ( _RTIString + "BOND"         ): self.ParseBondSection         ( )
                    elif line.startswith ( _RTIString + "CRYSIN"       ): self.ParseCrysinSection       ( )
                    elif line.startswith ( _RTIString + "MOLECULE"     ): self.ParseMoleculeSection     ( )
                    elif line.startswith ( _RTIString + "SUBSTRUCTURE" ): self.ParseSubstructureSection ( )
            except EOFError:
                pass
            # . Some final parsing and checks.
            if not hasattr ( self, "atoms" ): self.Warning ( "No atoms specified in file.", True )
            self.PreprocessSubstructures ( )
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseAtomSection ( self ):
        """Parse the ATOM section."""
        if hasattr ( self, "numberOfAtoms" ):
            self.atomicCharges     = Array.WithExtent ( self.numberOfAtoms ) ; self.atomicCharges.Set ( 0.0 )
            self.atoms             = []
            self.atomSubstructures = []
            self.atomTypes         = []
            self.xyz               = Coordinates3.WithExtent ( self.numberOfAtoms ) ; self.xyz.Set ( 0.0 )
            for i in range ( self.numberOfAtoms ):
                items = self.GetTokens ( converters = [ int, None, float, float, float, None, int, None, float ] )
                if len ( items ) < 6:
                    self.Warning ( "Invalid ATOM line.", True )
                else:
                    ssID         = None
                    ssName       = None
                    ( atomicNumber, geometry, formalCharge ) = self.AtomicNumberFromAtomType ( items[5] )
                    if atomicNumber <= 0: atomicNumber = PeriodicTable.AtomicNumber ( items[1] )
                    self.atoms.append ( Atom.WithOptions ( atomicNumber = atomicNumber ,
                                                           geometry     = geometry     ,
                                                           formalCharge = formalCharge ,
                                                           label        = items[1]     ) )
                    self.xyz[i,0]     = items[2]
                    self.xyz[i,1]     = items[3]
                    self.xyz[i,2]     = items[4]
                    self.atomTypes.append ( items[5] )
                    if len ( items ) >= 7: ssID   = items[6]
                    if len ( items ) >= 8: ssName = items[7]
                    if len ( items ) >= 9: self.atomicCharges[i] = items[8]
                    self.atomSubstructures.append ( [ ssID, ssName ] )
        else:
            self.Warning ( "Unknown number of atoms in molecule.", True )

    def ParseBondSection ( self ):
        """Parse the BOND section."""
        if hasattr ( self, "numberOfBonds" ) and hasattr ( self, "atoms" ):
            numberOfAtoms = getattr ( self, "numberOfAtoms", -1 )
            self.bonds    = []
            for i in range ( self.numberOfBonds ):
                items = self.GetTokens ( converters = [ int, int, int, None ] )
                if len ( items ) < 4:
                    self.Warning ( "Invalid BOND line.", True )
                else:
                    atom1 = items[1] - 1
                    atom2 = items[2] - 1
                    ( bondType, bondIsAromatic ) = _BondTypes.get ( items[3], ( BondType.Undefined, False ) )
                    if ( atom1 < 0 ) or ( atom1 >= self.numberOfAtoms ) or ( atom2 < 0 ) or ( atom2 >= self.numberOfAtoms ):
                        self.Warning ( "Bond atom indices out of range: {:d}, {:d}.".format ( atom1, atom2 ), True )
                    if bondType is not None:
                        self.bonds.append ( Bond.WithNodes ( self.atoms[atom1], self.atoms[atom2], isAromatic = bondIsAromatic, type = bondType ) )
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
        self.label = self.GetLine ( )                                  # . Molecule name.
        items      = self.GetTokens ( converters = [ int, int, int ] ) # . Number of atoms, bonds and substructures.
        if len ( items ) > 0: self.numberOfAtoms         = items[0]
        if len ( items ) > 1: self.numberOfBonds         = items[1]
        if len ( items ) > 2: self.numberOfSubstructures = items[2]

    def ParseSubstructureSection ( self ):
        """Parse the SUBSTRUCTURE section."""
        if hasattr ( self, "numberOfSubstructures" ):
            self.substructures = []
            for i in range ( self.numberOfSubstructures ):
                # . subst_id subst_name root_atom subst_type dict_type chain sub_type inter_bonds
                items  = self.GetTokens ( converters = [ int, None, int, None, int, None, None, int ] )
                values = []
                for i in range ( 8 ):
                    if len ( items ) > 0:
                        value = items.pop ( 0 )
                        if ( ( i == 5 ) or ( i == 6 ) ) and ( value == _NullString ): value = None
                    else:
                        value = None
                    values.append ( value )
                self.substructures.append ( values )
        else:
            self.Warning ( "Unknown number of substructures in molecule.", True )

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

    def PreprocessSubstructures ( self ):
        """Preprocess substructures."""
        # . There are a number of assumptions made here to convert MOL2 substructures into a pDynamo sequence:
        #   - both ATOM and SUBSTRUCTURE sections are present.
        #   - the substructure IDs are assumed to be unique and the same in both sections.
        #   - component labels are constructed as either (i) ssSubtype + tail or (ii) alphabetic + integer.
        #   - atoms in each substructure are contiguous.
        #   - there are no overlapping or multiple substructures
        #     (this could perhaps be remedied by only processing residue substructures?)
        if self.isParsed:
            atomSubstructures = getattr ( self, "atomSubstructures", None )
            substructures     = getattr ( self,     "substructures", None )
            if ( atomStructures is None ) and ( substructures is None ):
                pass
            elif ( atomStructures is not None ) and ( substructures is not None ):
                # . Initialization.
                isOK           = True
                majorSeparator = Sequence._attributable["labelSeparator"]
                minorSeparator = Sequence._attributable["fieldSeparator"]
                # . Generate path heads.
                pathHeads = {}
                ssFirst   = {}
                ssIndices = []
                ssLast    = {}
                ssNames   = {}
                for i in range ( self.numberOfSubstructures ):
                    ( ssID, ssName, ssRoot, _, _, ssChain, ssSubtype, _ ) = self.substructures[i]
                    if ( ssSubtype is not None ) and ssName.startswith ( ssSubtype ):
                        head = ssSubtype
                        tail = ssName[len(head):]
                    else:
                        first = len ( ssName )
                        for c in range ( len ( ssName )-1, -1, -1 ):
                            if ssName[c].isdigit ( ): first = c
                            else:                     break
                        head = ssName[:first]
                        tail = ssName[first:]
                    componentLabel = head
                    if len ( tail ) > 0: componentLabel += ( minorSeparator + tail ) 
                    if ssChain is None: entityLabel = ""
                    else:               entityLabel = ssChain
                    pathHeads[ssID] = "".join ( [ entityLabel, majorSeparator, componentLabel, majorSeparator ] )
                    ssIndices.append ( ( ssRoot-1, ssID ) )
                    ssNames  [ssID] = ssName
                # . Get first and last indices of each ssID.
                ssIndices.sort ( )
                if ( ssIndices[0][0] == 0 ) and ( ssIndices[-1][0] < self.numberOfAtoms ):
                    ssOld = None
                    while len ( ssIndices ) > 0:
                        ( ssRoot, ssID ) = ssIndices.pop ( 0 )
                        if ssOld is not None: ssLast[ssOld] = ssRoot - 1
                        ssFirst[ssID] = ssRoot
                        ssOld         = ssID
                    ssLast[ssOld] = self.numberOfAtoms - 1
                else:
                    self.Warning ( "Substructures do not cover the complete system.", False )
                    isOK = False
                # . Generate atom paths.
                if isOK:
                    self.atomPaths = []
                    for i in range ( self.numberOfAtoms ):
                        atomLabel        = self.atoms[i].label
                        ( ssID, ssName ) = self.atomSubstructures[i]
                        if ( ssName != ssNames[ssID] ) or ( i < ssFirst[ssID] ) or ( i > ssLast[ssID] ):
                            self.Warning ( "Atom substructure name or index error.", False )
                            isOK = False
                            break
                        self.atomPaths.append ( pathHeads[ssID] + atomLabel )
            else:
                self.Warning ( "No sequence can be generated as atom and substructure data missing.", False )

    def ToAtomNames ( self ):
        """Return atom names."""
        atomNames = None
        if self.isParsed and hasattr ( self, "atoms" ):
            atomNames = [ atom.label for atom in self.atoms ]
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
                if hasattr ( self, "atomPaths" ): sequence = Sequence.FromAtomPaths ( atomPaths, atoms = self.atoms )
                else:                             sequence = Sequence.FromAtoms ( self.atoms )
                system = System.FromSequence ( sequence, bonds = ( self.bonds if hasattr ( self, "bonds" ) else None ) )
                system.label        = self.label
                system.coordinates3 = self.xyz
                if hasattr ( self, "crystalSystem" ) and hasattr ( self, "symmetryParameters" ):
                    parameters                = self.crystalSystem.GetUniqueSymmetryParameters ( self.symmetryParameters )
                    system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( self.crystalSystem )
                    system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( **parameters )
            except Error as e:
                raise TextFileReaderError ( "Error generating system from MOL2 file.", e.args )
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

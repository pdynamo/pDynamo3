"""Classes and functions for reading CHARMM topology files."""

import itertools, math

from pCore       import Clone, logFile, LogFileActive, TextFileReader
from pScientific import PeriodicTable

#
# . Notes:
#
#   In some topology files the ATOM lines contain extra terms. These can represent extra
#   exclusions (e.g. TOPH19) or fluctuating charge parameters (e.g. CHEQ).
#

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
def grouper ( n, iterable, fillValue = None ):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    arguments = [ iter ( iterable ) ] * n
# . Change to zip_longest in Python3.
    return itertools.zip_longest ( fillvalue = fillValue, *arguments )

#===================================================================================================================================
# . CHARMM topology classes.
#===================================================================================================================================
class CHARMMTopologyAtomType:
    """An atom type definition."""

    _attributable = { "atomicNumber" :   -1 ,
                      "index"        :   -1 ,
                      "label"        : None ,
                      "mass"         :  0.0 }

    def __init__ ( self, **options ):
        """Constructor."""
        for ( key, value ) in self.__class__._attributable.items ( ): setattr ( self, key, value )
        for ( key, value ) in                      options.items ( ): setattr ( self, key, value )

class CHARMMTopologyResidue:
    """A patch or residue definition."""

    _attributable = { "acceptors"           : None ,
                      "angles"              : None ,
                      "atoms"               : None ,
                      "bonds"               : None ,
                      "charge"              :  0.0 ,
                      "cmaps"               : None ,
                      "currentGroup"        :    0 ,
                      "defaultFirst"        : None ,
                      "defaultLast"         : None ,
                      "deleteableAcceptors" : None ,
                      "deleteableAngles"    : None ,
                      "deleteableAtoms"     : None ,
                      "deleteableBonds"     : None ,
                      "deleteableDihedrals" : None ,
                      "deleteableICs"       : None ,
                      "deleteableImpropers" : None ,
                      "dihedrals"           : None ,
                      "donors"              : None ,
                      "ics"                 : None ,
                      "impropers"           : None ,
                      "label"               : None }

    def __init__ ( self, **options ):
        """Constructor."""
        for ( key, value ) in self.__class__._attributable.items ( ): setattr ( self, key, value )
        for ( key, value ) in                      options.items ( ): setattr ( self, key, value )

    def AddAcceptor ( self, label ):
        """Add an acceptor."""
        if self.acceptors is None: self.acceptors = set ( )
        self.acceptors.add ( label )

    def AddAngle ( self, label1, label2, label3 ):
        """Add an angle."""
        if self.angles is None: self.angles = []
        self.angles.append ( ( label1, label2, label3 ) )

    def AddAtom ( self, label, atomType, charge ):
        """Add an atom."""
        if self.atoms is None: self.atoms = []
        self.atoms.append ( ( label, atomType, charge, self.currentGroup ) )

# Label unique?
# total charge = charge?
# Canonicalize all angle, bond, ... keys. Redundancy?

    def AddBond ( self, label1, label2 ):
        """Add a bond."""
        if self.bonds is None: self.bonds = []
        self.bonds.append ( ( label1, label2 ) )

    def AddCMap ( self, label1, label2, label3, label4 ):
        """Add a CMAP."""
        if self.cmaps is None: self.cmaps = []
        self.cmaps.append ( ( label1, label2, label3, label4 ) )

    def AddDeleteableAcceptor ( self, label ):
        """Mark an acceptor as deletable."""
        if self.deleteableAcceptors is None: self.deleteableAcceptors = set ( )
        self.deleteableAcceptors.add ( label )

    def AddDeleteableAngle ( self, label1, label2, label3 ):
        """Mark an angle as deletable."""
        if self.deleteableAngles is None: self.deleteableAngles = []
        self.deleteableAngles.append ( ( label1, label2, label3 ) )

    def AddDeleteableAtom ( self, label ):
        """Mark an atom as deletable."""
        if self.deleteableAtoms is None: self.deleteableAtoms = set ( )
        self.deleteableAtoms.add ( label )

    def AddDeleteableBond ( self, label1, label2 ):
        """Mark a bond as deletable."""
        if self.deleteableBonds is None: self.deleteableBonds = []
        self.deleteableBonds.append ( ( label1, label2 ) )

    def AddDeleteableDihedral ( self, label1, label2, label3, label4 ):
        """Mark a dihedral as deletable."""
        if self.deleteableDihedrals is None: self.deleteableDihedrals = []
        self.deleteableDihedrals.append ( ( label1, label2, label3, label4 ) )

    def AddDeleteableIC ( self, label1, label2, label3, label4 ):
        """Mark an IC as deletable."""
        if self.deleteableICs is None: self.deleteableICs = []
        self.deleteableICs.append ( ( label1, label2, label3, label4 ) )

    def AddDeleteableImproper ( self, label1, label2, label3, label4 ):
        """Mark an improper as deletable."""
        if self.deleteableImpropers is None: self.deleteableImpropers = []
        self.deleteableImpropers.append ( ( label1, label2, label3, label4 ) )

    def AddDihedral ( self, label1, label2, label3, label4 ):
        """Add a dihedral."""
        if self.dihedrals is None: self.dihedrals = []
        self.dihedrals.append ( ( label1, label2, label3, label4 ) )

    def AddDonor ( self, label ):
        """Add a donor."""
        if self.donors is None: self.donors = set ( )
        self.donors.add ( label )

    def AddIC ( self, labels, ics ):
        """Add an IC."""
        if self.ics is None: self.ics = []
        self.ics.append ( ( labels, ics ) )

    def AddImproper ( self, label1, label2, label3, label4 ):
        """Add an improper."""
        if self.impropers is None: self.impropers = []
        self.impropers.append ( ( label1, label2, label3, label4 ) )

    def AddPatch ( self, key, value ):
        """Add a patch."""
        if value == "NONE": value = None
        if   key[0:4] == "FIRS": self.defaultFirst = value
        elif key[0:4] == "LAST": self.defaultLast  = value

    def StartNewGroup ( self ):
        """Start a new group."""
        self.currentGroup += 1

#===================================================================================================================================
# . CHARMM topology file reader class.
#===================================================================================================================================
class CHARMMTopologyFileReader ( TextFileReader ):
    """CHARMMPSFFileReader is the class for CHARMM PSF files that are to be read."""

    _classLabel = "CHARMM Topology File Reader"

    def GetLine ( self, signalWarnings = False ):
        """Get a non-empty line removed of comments.

        Continuation lines are checked for and included in the line if they are present.
        """
        try:
            found      = ""
            isNotFound = True
            while isNotFound:
                line = next ( self.file ).strip ( ).upper ( )
                self.linesParsed += 1
                # . Check for comments.
                index = line.find ( "!" )
                if index >= 0: line = line[:index].strip ( )
                # . Check for continuations.
                if line.endswith ( "-" ): line = line[:-1].strip ( )
                else:                     isNotFound = ( len ( line ) <= 0 )
                found += line
            return found
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            if ( len ( found ) > 0 ): return found
            else:                     raise EOFError

    def Initialize ( self ):
        """Initialization for parsing."""
        self.atomTypes             = {}
        self.autogenerateAngles    = False
        self.autogenerateDihedrals = False
        self.backwardDeclarations  = set ( )
        self.forwardDeclarations   = set ( )
        self.defaultFirst          = None
        self.defaultLast           = None
        self.patches               = {}
        self.residues              = {}
        self.title                 = []
        self.version               = None

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                self.Initialize ( )
                # . Get the title and version number.
                while True:
                    line = self.GetLine ( )
                    if line.startswith ( "*" ):
                        line = line[1:].strip ( )
                        if len ( line ) > 0: self.title.append ( line )
                    else:
                        tokens       = line.split ( )
                        self.version = int ( tokens[0] )
                        break
                # . Get the remaining preamble lines (MASS, etc.).
                while True:
                    line    = self.GetLine ( )
                    tokens  = line.split   ( )
                    keyword = tokens.pop ( 0 )[0:4]
                    if   keyword == "AUTO": self.ParseAutogenerate ( tokens )
                    elif keyword == "DECL": self.ParseDeclaration  ( tokens )
                    elif keyword == "DEFA": self.ParseDefaults     ( tokens )
                    elif keyword == "MASS": self.ParseMass         ( tokens )
                    else: break
                # . Parse PRES and RESI definitions - keyword and tokens exist already.
                while True:
                    newLine = False
                    if keyword == "DEFA":
                        self.ParseDefaults ( tokens )
                        newLine = True
                    elif keyword == "END" : break
                    elif keyword == "PRES": ( keyword, tokens ) = self.ParseResidue ( tokens, self.patches , "PATCH"   )
                    elif keyword == "RESI": ( keyword, tokens ) = self.ParseResidue ( tokens, self.residues, "RESIDUE" )
                    else:
                        self.Warning ( "Unrecognized keyword: " + keyword + ".", True )
                        newLine = True
                    if newLine:
                        line    = self.GetLine ( )
                        tokens  = line.split ( )
                        keyword = tokens.pop ( 0 )[0:4]
            except EOFError:
                pass
            # . Processing here?

#check types, etc.
#check residues, etc.

            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseAutogenerate ( self, tokens ):
        """Autogenerate line."""
        for token in tokens:
            if   token[0:4] == "ANGL": self.autogenerateAngles    = True
            elif token[0:4] == "DIHE": self.autogenerateDihedrals = True

    def ParseDeclaration ( self, tokens ):
        """Declaration line."""
        isOK = False
        if ( len ( tokens ) == 1 ) and ( len ( tokens[0] ) > 1 ):
            isOK   = True
            sign   = tokens[0][0:1]
            label  = tokens[0][1:]
            if   sign == "+": self.forwardDeclarations.add  ( label )
            elif sign == "-": self.backwardDeclarations.add ( label )
            else: isOK = False
        if not isOK: self.Warning ( "Invalid DECLARATIONS line.", False )

    def ParseDefaults ( self, tokens ):
        """Defaults line."""
        if ( len ( tokens ) == 4 ) and ( tokens[0].startswith ( "FIRS" ) ) and ( tokens[2].startswith ( "LAST" ) ):
            if tokens[1] == "NONE": self.defaultFirst = None
            else:                   self.defaultFirst = tokens[1]
            if tokens[3] == "NONE": self.defaultLast  = None
            else:                   self.defaultLast  = tokens[3]
        else: self.Warning ( "Invalid DEFAULTS line.", False )

    def ParseMass ( self, tokens ):
        """Mass line."""
        if len ( tokens ) < 3: self.Warning ( "Invalid MASS line.", False )
        else:
            index = int   ( tokens[0] )
            label =         tokens[1]
            mass  = float ( tokens[2] )
            if len ( tokens ) > 3: atomicNumber = PeriodicTable.AtomicNumber         ( tokens[3] )
            else:                  atomicNumber = PeriodicTable.AtomicNumberFromMass ( mass      )
            self.atomTypes[label] = CHARMMTopologyAtomType ( atomicNumber = atomicNumber, index = index, label = label, mass = mass )

    def ParseResidue ( self, tokens, container, tag ):
        """Residue definition."""
        if len ( tokens ) > 0: label  = tokens[0]
        else:                  label  = None
        if len ( tokens ) > 1: charge = float ( tokens[1] )
        else:                  charge = 0.0
        if   label is None         : self.Warning ( "Invalid " + tag + " line.", False )
        elif label in self.residues: self.Warning ( "Duplicate " + tag + ": " + label + ".", False )
        else:
            item = CHARMMTopologyResidue ( charge = charge, defaultFirst = self.defaultFirst, defaultLast = self.defaultLast, label = label )
            container[label] = item 
        while True:
            line    = self.GetLine ( )
            tokens  = line.split   ( )
            keyword = tokens.pop ( 0 )[0:4]
            if keyword == "ACCE":
                for token in tokens: item.AddAcceptor ( token )
            elif ( keyword == "ANGL" ) or ( keyword == "THET" ):
                for ( label1, label2, label3 ) in grouper ( 3, tokens ): item.AddAngle ( label1, label2, label3 )
            elif keyword == "ATOM":
                label    = tokens[0]
                atomType = self.atomTypes.get ( tokens[1], None )
                charge   = tokens[2]
                if atomType is None: self.Warning ( "Unknown atom type - " + tokens[1] + " - for atom " + label + ".", False )
                item.AddAtom ( label, atomType, charge )
#                if len ( tokens ) > 3: print tokens
            elif ( keyword == "BILD" ) or ( keyword == "IC" ):
                ics = []
                for token in tokens[4:]: ics = float ( token )
                item.AddIC ( tokens[0:4], ics )
            elif ( keyword == "BOND" ) or ( keyword == "DOUB" ):
                for ( label1, label2 ) in grouper ( 2, tokens ): item.AddBond ( label1, label2 )
            elif keyword == "CMAP":
                for ( label1, label2, label3, label4 ) in grouper ( 4, tokens ): item.AddCMap ( label1, label2, label3, label4 )
            elif keyword == "DELE":
                option = tokens.pop ( 0 )[0:4]
                if option == "ACCE":
                    item.AddDeleteableAcceptor ( tokens[0] )
                elif ( option == "ANGL" ) or ( option == "THET" ):
                    item.AddDeleteableAngle ( tokens[0], tokens[1], tokens[2] )
                elif option == "ATOM":
                    item.AddDeleteableAtom ( tokens[0] )
                elif option == "BOND":
                    item.AddDeleteableBond ( tokens[0], tokens[1] )
                elif option == "DIHE":
                    item.AddDeleteableDihedral ( tokens[0], tokens[1], tokens[2], tokens[3] )
                elif option == "IC":
                    item.AddDeleteableIC ( tokens[0], tokens[1], tokens[2], tokens[3] )
                elif ( option == "IMPH" ) or ( option == "IMPR" ):
                    item.AddDeleteableImproper ( tokens[0], tokens[1], tokens[2], tokens[3] )
                else: self.Warning ( "Unable to interpret DELETE option: " + option + ".", False )
            elif keyword == "DIHE":
                for ( label1, label2, label3, label4 ) in grouper ( 4, tokens ): item.AddDihedral ( label1, label2, label3, label4 )
            elif keyword == "DONO":
                for token in tokens: item.AddDonor ( token )
            elif keyword == "GROU":
                item.StartNewGroup ( )
            elif ( keyword == "IMPH" ) or ( keyword == "IMPR" ):
                for ( label1, label2, label3, label4 ) in grouper ( 4, tokens ): item.AddImproper ( label1, label2, label3, label4 )
            elif keyword == "PATC":
                for ( key, value ) in grouper ( 2, tokens ): item.AddPatch ( key, value )
            else:
                break
        return ( keyword, tokens )

    @classmethod
    def PathsToTopologies ( selfClass, paths, log = logFile ):
        """Return the topologies from a list of topology files."""
        # . Initialization.
        topologies = None
        # . Get the individual topology sets.
        topologyData = []
        for path in paths:
            topologyFile = selfClass.FromPath ( path )
            topologyFile.Parse   ( log = log )
            topologyFile.Summary ( log = log )
            topologyData.append ( topologyFile.ToTopologies ( ) )
        # . Get a merged topology set.
        if len ( topologyData ) > 0:
            topologies = topologyData[0]
            for data in topologyData[1:]:
                topologies = topologies.Merge ( data )
        return topologies

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Autogenerate Angles"   , "{!r}".format ( self.autogenerateAngles           ) ) ,
                 ( "Autogenerate Dihedrals", "{!r}".format ( self.autogenerateDihedrals        ) ) ,
                 ( "Backward Declarations" , "{:d}".format ( len ( self.backwardDeclarations ) ) ) ,
                 ( "Forward Declarations"  , "{:d}".format ( len ( self.forwardDeclarations  ) ) ) ,
                 ( "Default First"         , "{:s}".format ( self.defaultFirst                 ) ) ,
                 ( "Default Last"          , "{:s}".format ( self.defaultLast                  ) ) ,
                 ( "Atom Types"            , "{:d}".format ( len ( self.atomTypes            ) ) ) ,
                 ( "Patches"               , "{:d}".format ( len ( self.patches              ) ) ) ,
                 ( "Residues"              , "{:d}".format ( len ( self.residues             ) ) ) ]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

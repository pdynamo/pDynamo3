"""Classes and functions for reading CHARMM parameter files."""

import math

from pCore       import Clone, logFile, LogFileActive, TextFileReader
from pScientific import Units

#
# . Notes:
#
#   This needs cleaning up. More methods need to be added to ParameterContainer for handling different parameter types. In
#   particular, this would remove the need for duplicate key-handling code in CHARMMPSFReader.
#

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Mapping from CHARMM section names to pMolecule section names.
_SectionMapping = { "ANGL" : "Angle"   , "ATOM" : "Atom" , "BOND" : "Bond"   , "CMAP" : "CMap"   , "DIHE" : "Dihedral", "END"  : "End"  , "EQUI" : "Equivalence", "HBON" : "HBond", "IMPH" : "Improper", \
                    "IMPR" : "Improper", "NBFI" : "NBFix", "NBON" : "NonBond", "NONB" : "NonBond", "PHI"  : "Dihedral", "PRIN" : "Print", "SPAS" : "SPAS" , "THET" : "Dihedral"  }

# . Section names excluding titles which start with "*".
_SectionNames = ( "ANGL", "ATOM", "BOND", "CMAP", "DIHE", "END", "EQUI", "HBON", "IMPH", "IMPR", "NBFI", "NBON", "NONB", "PHI", "PRIN", "SPAS", "THET" )

# . Wildcards.
# . For dihedrals it seems the only wildcard format used is "X A B X".
# . For impropers it seems the only wildcard format used is "A X X B".
_DihedralWildCard = "X"
_NonBondWildCards = ( "*", "%", "#", "+" )

#===================================================================================================================================
# . CHARMM parameter container.
#===================================================================================================================================
class CHARMMParameterContainer:
    """A container for CHARMM parameters."""

    def Merge ( self, other ):
        """Merge two sets of parameters."""
        new = CHARMMParameterContainer ( )
        new.atomtypes = self.atomtypes.union ( other.atomtypes )
        for attribute in ( "angles", "bonds", "cmaps", "dihedrals", "dihedralwilds", "impropers", "improperwilds", "nonbonds", "nonbondwilds", "nonbond14s", "nonbond14wilds", "ureybradleys" ):
            old1   = getattr ( self , attribute, {} )
            old2   = getattr ( other, attribute, {} )
            common = set ( old1.keys ( ) ).intersection ( set ( old2.keys ( ) ) )
            if len ( common ) > 0: print ( "There are {:d} common terms for attribute {:s}.".format ( len ( common ), attribute ) )
            newd = dict ( old1 )
            newd.update ( old2 )
            setattr ( new, attribute, newd )
        return new

#===================================================================================================================================
# . CHARMM parameter file reader class.
#===================================================================================================================================
class CHARMMParameterFileReader ( TextFileReader ):
    """CHARMMPSFFileReader is the class for CHARMM PSF files that are to be read."""

    _classLabel   = "CHARMM Parameter File Reader"
    _withSections = True

    def AddAtomType ( self, label ):
        """Add an atom type."""
        atomtype = label.upper ( )
        self.atomtypes.add ( atomtype )
        return atomtype

    def GetLine ( self, signalWarnings = False ):
        """Get a non-empty line removed of comments.

        Continuation lines are checked for and included in the line if they are present.
        """
        try:
            found      = ""
            isNotFound = True
            while isNotFound:
                line = next ( self.file ).strip ( )
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

    def IdentifySection ( self, line ):
        """Identify a section of the parameter file."""
        section = None
        if line.startswith ( "*" ):
            section = "Title"
            line    = line[1:].strip ( )
        else:
            tokens = line.split ( " ", 1 )
            head   = tokens[0].upper ( )
            for name in _SectionNames:
                if head.startswith ( name ):
                    section = _SectionMapping[name]
                    if len ( tokens ) == 1: line = ""
                    else:                   line = tokens[1]
                    break
        return ( section, line )

    def Initialize ( self ):
        """Initialization for parsing."""
        # . Read data.
        self.numberOfSectionLines = {}
        self.sectionBodies        = {}
        self.sectionHeaders       = {}
        # . Processed data.
        self.angles         = {}
        self.atomtypes      = set ( )
        self.bonds          = {}
        self.cmaps          = {}
        self.dihedrals      = {}
        self.dihedralwilds  = {}
        self.impropers      = {}
        self.improperwilds  = {}
        self.nonbonds       = {}
        self.nonbondwilds   = {}
        self.nonbond14s     = {}
        self.nonbond14wilds = {}
        self.ureybradleys   = {}

    def Parse ( self, convertUnits = True, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            self.SetUnitConversions ( convertUnits )
            # . Open the file.
            self.Open ( )
            try:
                self.Initialize ( )
                # . Loop over lines in the file.
                current = None
                while True:
                    line = self.GetLine ( )
                    ( section, line ) = self.IdentifySection ( line )
                    if section is None:
                        if current is None: self.Warning ( "Unrecognized section.", False )
                        else:
                            data = self.sectionBodies.get ( current, [] )
                            data.append ( line.split ( ) )
                            self.sectionBodies[current] = data
                    else:
                        if section == "End": break
                        current = section
                        if len ( line ) > 0:
                            data = self.sectionHeaders.get ( current, [] )
                            data.append ( line )
                            self.sectionHeaders[current] = data
                    if current is not None:
                        self.numberOfSectionLines[current] = self.numberOfSectionLines.get ( current, 0 ) + 1
            except EOFError:
                pass
            # . Process the data.
            for section in self.numberOfSectionLines.keys ( ):
                processmethod = getattr ( self, "Process" + section, None )
                if processmethod is not None: processmethod ( )
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    @classmethod
    def PathsToParameters ( selfClass, paths, convertUnits = True, log = logFile ):
        """Return the parameters from a list of parameter files."""
        # . Initialization.
        parameters = None
        # . Get the individual parameter sets.
        parameterData = []
        for path in paths:
            parameterFile = selfClass.FromPath ( path )
            parameterFile.Parse   ( convertUnits = convertUnits, log = log )
            parameterFile.Summary ( log = log )
            parameterData.append ( parameterFile.ToParameters ( ) )
        # . Get a merged parameter set.
        if len ( parameterData ) > 0:
            parameters = parameterData[0]
            for data in parameterData[1:]:
                parameters = parameters.Merge ( data )
        return parameters

    def ProcessAngle ( self ):
        """Process angles and Urey-Bradley terms."""
        data = self.sectionBodies.get ( "Angle", None )
        if data is not None:
            for datum in data:
                try:
                    # . Angle.
                    t1 = self.AddAtomType ( datum[0] )
                    t2 = self.AddAtomType ( datum[1] )
                    t3 = self.AddAtomType ( datum[2] )
                    key = ( max ( t1, t3 ), t2, min ( t1, t3 ) )
                    if key in self.angles: self.Warning ( "Duplicate angle for " + "-".join ( datum[0:3] ) + ".", False )
                    else: self.angles[key] = ( self.energyConverter * float ( datum[3] ), self.angleConverter * float ( datum[4] ) )
                    # . Urey-Bradley.
                    if len ( datum ) > 5:
                        if key in self.ureybradleys: self.Warning ( "Duplicate Urey-Bradley term for " + "-".join ( datum[0:3] ) + ".", False )
                        else: self.ureybradleys[key] = ( self.energyConverter * float ( datum[5] ), float ( datum[6] ) )
                except:
                    self.Warning ( "Unable to process angle/Urey-Bradley terms: " + " ".join ( datum ) + ".", False )

    def ProcessBond ( self ):
        """Process bonds."""
        data = self.sectionBodies.get ( "Bond", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = self.AddAtomType ( datum[0] )
                    t2 = self.AddAtomType ( datum[1] )
                    key = ( max ( t1, t2 ), min ( t1, t2 ) )
                    if key in self.bonds: self.Warning ( "Duplicate bond for " + "-".join ( datum[0:2] ) + ".", False )
                    else: self.bonds[key] = ( self.energyConverter * float ( datum[2] ), float ( datum[3] ) )
                except:
                    self.Warning ( "Unable to process bond: " + " ".join ( datum ) + ".", False )

    def ProcessCMap ( self ):
        """Process cmap terms."""
        data = self.sectionBodies.get ( "CMap", None )
        if data is not None:
            while len ( data ) > 0:
                datum = data.pop ( 0 )
                try:
                    t1 = self.AddAtomType ( datum[0] )
                    t2 = self.AddAtomType ( datum[1] )
                    t3 = self.AddAtomType ( datum[2] )
                    t4 = self.AddAtomType ( datum[3] )
                    t5 = self.AddAtomType ( datum[4] )
                    t6 = self.AddAtomType ( datum[5] )
                    t7 = self.AddAtomType ( datum[6] )
                    t8 = self.AddAtomType ( datum[7] )
                    n  = int ( datum[8] )
                    if   t2 >  t3: key1 = [ t1, t2, t3, t4 ]
                    elif t2 == t3: key1 = [ max ( t1, t4 ), t2, t3, min ( t1, t4 ) ]
                    else:          key1 = [ t4, t3, t2, t1 ]
                    if   t6 >  t7: key2 = [ t5, t6, t7, t8 ]
                    elif t6 == t7: key2 = [ max ( t5, t8 ), t6, t7, min ( t5, t8 ) ]
                    else:          key2 = [ t8, t7, t6, t5 ]
                    key = tuple ( key1 + key2 )
                    # . Calculate abscissa.
                    increment = 2.0 * math.pi / float ( n )
                    abscissa  = [ -math.pi + i * increment for i in range ( n + 1 ) ]
                    # . Gather data - one row at a time and expand so that upper boundary values are also included.
                    # . First dihedral changes most rapidly (i.e. columnwise storage).
                    cmap = []
                    for i in range ( n ):
                        tokens = []
                        while len ( tokens ) != n: tokens.extend ( data.pop ( 0 ) )
                        values = []
                        for token in tokens: values.append ( float ( token ) * self.energyConverter )
                        values.append ( values[0] )
                        cmap.extend ( values )
                    # . Duplicate first column.
                    cmap.extend ( cmap[0:n+1] )
                    # . Now transpose so have rowwise storage.
                    # . This is not necessary.
#                    n += 1
#                    for i in range ( 1, n ):
#                        for j in range ( 0, i ):
#                            cij = cmap[i*n+j]
#                            cji = cmap[j*n+i]
#                            cmap[i*n+j] = cji
#                            cmap[j*n+i] = cij
                    if key in self.cmaps: self.Warning ( "Duplicate cmap for " + "-".join ( key[0:4] ) + "/" + "-".join ( key[4:] ) + ".", False )
                    else: self.cmaps[key] = ( abscissa, cmap )
                except:
                    self.Warning ( "Unable to process cmap: " + " ".join ( datum ) + ".", False )

    def ProcessDihedral ( self ):
        """Process dihedrals."""
        data = self.sectionBodies.get ( "Dihedral", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = self.AddAtomType ( datum[0] )
                    t2 = self.AddAtomType ( datum[1] )
                    t3 = self.AddAtomType ( datum[2] )
                    t4 = self.AddAtomType ( datum[3] )
                    if   t2 >  t3: key = ( t1, t2, t3, t4 )
                    elif t2 == t3: key = ( max ( t1, t4 ), t2, t3, min ( t1, t4 ) )
                    else:          key = ( t4, t3, t2, t1 )
                    if _DihedralWildCard in key: parameters = self.dihedralwilds
                    else:                        parameters = self.dihedrals
                    # . Allow for multiple dihedrals.
                    kchi  = float ( datum[4] ) * self.energyConverter
                    n     = int   ( datum[5] )
                    delta = self.angleConverter * float ( datum[6] )
                    isOK  = True
                    terms = parameters.get ( key, [] )
                    for term in terms:
                        if term[0] == n:
                            isOK = False
                            self.Warning ( "Duplicate dihedral for " + "-".join ( datum ) + ".", False )
                            break
                    if isOK:
                        terms.append ( ( n, kchi, delta ) )
                        terms.sort ( )
                        parameters[key] = terms
                except:
                    self.Warning ( "Unable to process dihedral: " + " ".join ( datum ) + ".", False )

    def ProcessImproper ( self ):
        """Process impropers."""
        data = self.sectionBodies.get ( "Improper", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = self.AddAtomType ( datum[0] )
                    t2 = self.AddAtomType ( datum[1] )
                    t3 = self.AddAtomType ( datum[2] )
                    t4 = self.AddAtomType ( datum[3] )
                    if   t1 >  t4: key = ( t1, t2, t3, t4 )
                    elif t1 == t4: key = ( t1, max ( t2, t3 ), min ( t2, t3 ), t4 )
                    else:          key = ( t4, t3, t2, t1 )
                    if _DihedralWildCard in key: parameters = self.improperwilds
                    else:                        parameters = self.impropers
                    if key in parameters: self.Warning ( "Duplicate improper for " + "-".join ( datum[0:4] ) + ".", False )
                    else: parameters[key] = ( self.energyConverter * float ( datum[4] ), self.angleConverter * float ( datum[6] ) ) # . No check on central integral value (always "0").
                except:
                    self.Warning ( "Unable to process improper: " + " ".join ( datum ) + ".", False )

    def ProcessNonBond ( self ):
        """Process nonbonds."""
        data = self.sectionBodies.get ( "NonBond", None )
        if data is not None:
            for datum in data:
                try:
                    t = self.AddAtomType ( datum[0] )
                    # . Check for wildcards.
                    isWildCard = False
                    for wildcard in _NonBondWildCards:
                        if t.find ( wildcard ) >= 0:
                            isWildCard = True
                            break
                    # . Process both normal and 1-4 values.
                    increment = 0
                    for ( nonwilds, wilds, tag ) in ( ( self.nonbonds, self.nonbondwilds, "nonbond" ), ( self.nonbond14s, self.nonbond14wilds, "1-4 nonbond" ) ):
                        if isWildCard: parameters = wilds
                        else:          parameters = nonwilds
                        if t in parameters: self.Warning ( "Duplicate " + tag + " for " + t + ".", False )
                        else:
                            p = float ( datum[increment+1] )
                            e = float ( datum[increment+2] ) * self.energyConverter
                            r = float ( datum[increment+3] ) * 2.0
                            if ( p != 0.0 ) or ( e > 0.0 ): self.Warning ( "Unable to process non-epsilon/sigma " + tag + " for " + t + ".", False )
                            else: parameters[t] = ( math.fabs ( e ), r )
                        if len ( datum ) <= 4: break
                        increment += 3
                except:
                    self.Warning ( "Unable to process nonbond: " + " ".join ( datum ) + ".", False )

    def SetUnitConversions ( self, convertUnits ):
        """Determine unit conversion constants."""
        if convertUnits:
            self.angleConverter  = Units.Angle_Degrees_To_Radians
            self.energyConverter = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
        else:
            self.angleConverter  = 1.0
            self.energyConverter = 1.0

    def SummaryItems ( self ):
        """Summary items."""
        # . Parsed data.
        items = [ ( "Section Lines", True ) ]
        items.extend ( [ ( key, "{:d}".format ( self.numberOfSectionLines[key] ) ) for key in self.numberOfSectionLines.keys ( ) ] )
        # . Processed data.
        items.extend ( [ ( "Parameters"             , True                                          ) ,
                         ( "Atom Types"             , "{:d}".format ( len ( self.atomtypes      ) ) ) ,
                         ( "Angles"                 , "{:d}".format ( len ( self.angles         ) ) ) ,
                         ( "Bonds"                  , "{:d}".format ( len ( self.bonds          ) ) ) ,
                         ( "Cmaps"                  , "{:d}".format ( len ( self.cmaps          ) ) ) ,
                         ( "Dihedrals (Nonwild)"    , "{:d}".format ( len ( self.dihedrals      ) ) ) ,
                         ( "Dihedrals (Wild)"       , "{:d}".format ( len ( self.dihedralwilds  ) ) ) ,
                         ( "Impropers (Nonwild)"    , "{:d}".format ( len ( self.impropers      ) ) ) ,
                         ( "Impropers (Wild)"       , "{:d}".format ( len ( self.improperwilds  ) ) ) ,
                         ( "Nonbonds (Nonwild)"     , "{:d}".format ( len ( self.nonbonds       ) ) ) ,
                         ( "Nonbonds (Wild)"        , "{:d}".format ( len ( self.nonbondwilds   ) ) ) ,
                         ( "Nonbonds 1-4 (Nonwild)" , "{:d}".format ( len ( self.nonbond14s     ) ) ) ,
                         ( "Nonbonds 1-4 (Wild)"    , "{:d}".format ( len ( self.nonbond14wilds ) ) ) ,
                         ( "Urey-Bradleys"          , "{:d}".format ( len ( self.ureybradleys   ) ) ) ] )
        return items

    def ToParameters ( self ):
         """Return an object containing all the parameters."""
         parameters = None
         if self.isParsed:
             parameters = CHARMMParameterContainer ( )
             for attribute in ( "atomtypes", "angles", "bonds", "cmaps", "dihedrals", "dihedralwilds", "impropers", "improperwilds", "nonbonds", "nonbondwilds", "nonbond14s", "nonbond14wilds", "ureybradleys" ):
                 setattr ( parameters, attribute, getattr ( self, attribute ) )
         return parameters

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

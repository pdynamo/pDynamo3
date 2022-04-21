#===================================================================================================================================
#
# . File      : GromacsTopologyFileReader.py
#
# . Author    : Guilherme M. Arantes (University of Sao Paulo, Brazil, 2011)
#
# . Based on and to be used with the pDynamo library, copyright CEA, CNRS, Martin J. Field 
#
# . Web       : http://www.pdynamo.org
#
#===================================================================================================================================
# 
# . Notes:
#
# . CHARMM and AMBER force-field topologies may be read at the moment.
#
# . Pre-processed topologies are necessary, use "grompp -pp ...". 
#
# . pDynamo does not support SETTLE, use "define = -DFLEXIBLE" in your .mdp file. 
#
# . The CHARMM definition of TIP3P should be asked explicitly via -DCHARMM_TIP3P, if desired.
#
# . Online patches (i.e. parameter definition after moleculetype sections) will work only for harmonic bonds and angles.
#
# . This script is far from fully compatible with the Gromacs specification (see their manual), but it works for most topologies
#   generated with Gromacs 4.5.x.
#
#
#===================================================================================================================================
"""Classes and functions for reading Gromacs topologies."""

import itertools, math, sys 

from pCore              import Clone                     , \
                               DataType                  , \
                               logFile                   , \
                               LogFileActive             , \
                               SelfPairList              , \
                               TextFileReader
from pMolecule          import EnergyModelPriority       , \
                               Sequence                  , \
                               System
from pMolecule.MMModel  import CMAPDihedralContainer     , \
                               FourierDihedralContainer  , \
                               HarmonicAngleContainer    , \
                               HarmonicBondContainer     , \
                               HarmonicImproperContainer , \
                               LJForm                    , \
                               LJParameterContainer      , \
                               MMModelCHARMM             , \
                               MMModelAMBER              , \
                               MMModelOPLS               , \
                               MMModelState
from pScientific        import PeriodicTable             , \
                               Units
from pScientific.Arrays import Array

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The character for parameter keys.
_ParameterKeyFieldSeparator = ":"

# . Mapping from Gromacs to local names.
_SectionMapping = { "angles"                  : "Angles"       ,
                    "angletypes"              : "Angle"        ,
                    "atoms"                   : "Atoms"        ,
                    "atomtypes"               : "NonBond"      ,
                    "bonds"                   : "Bonds"        ,
                    "bondtypes"               : "Bond"         ,
                    "cmap"                    : "CMaps"        ,
                    "cmaptypes"               : "CMap"         ,
                    "constrainttypes"         : ""             ,        
                    "defaults"                : "ForceField"   ,
                    "exclusions"              : ""             ,
                    "implicit_genborn_params" : ""             ,
                    "impropers"               : "Impropers"    ,
                    "impropertypes"           : "Improper"     ,
                    "molecules"               : "Molecules"    ,
                    "moleculetype"            : "Molecule"     ,
                    "pairs"                   : "Exclusions14" ,
                    "pairtypes"               : "NonBond14"    ,
                    "position_restraints"     : ""             ,
                    "propertypes"             : "Dihedral"     ,
                    "propers"                 : "Dihedrals"    ,
                    "settles"                 : ""             ,
                    "system"                  : "Title"        }

_MoleculeSections = [ "angles", "atoms", "bonds", "cmap", "impropers", "pairs", "propers" ]

# . Wildcards.
_DihedralWildCard = "X"

#===================================================================================================================================
# . Gromacs parameter container.
#===================================================================================================================================
class GromacsParameterContainer:
    """A container for Gromacs parameters."""
    pass

#===================================================================================================================================
# . General Gromacs file reader class.
#===================================================================================================================================
class GromacsFileReader ( TextFileReader ):
    """Class for general Gromacs topology files that are to be read.""" # . Needs some work!!!

    _classLabel = "Gromacs File Reader"

    def GetDihedralSection ( self, section ):
        """Distinguish between proper and improper dihedrals."""
        line   = self.GetLine ( signalWarnings = False )
        tokens = line.split ( ) 
        if len( tokens )  < 5 : self.Warning ( "Invalid or unrecognized dihedral data line.", True )
        else                  : 
            angleType = int( tokens[4] )
            # . CHARMM {proper: 9, improper: 2}, AMBER {proper: 9, improper: 4}, OPLS {proper:3 ,improper:3 }
            if   angleType == 2 or angleType == 4 : section = "improper" + section[8:] 
            elif angleType == 3 or angleType == 9 : section =   "proper" + section[8:] 
        return ( section, line )

    def GetLine ( self, signalWarnings = False ):
        """Get a non-empty line removed of comments.
        Continuation lines are checked for and included in the line if they are present.
        """
        try:
            found      = ""
            isNotFound = True
            while isNotFound:
                line = next ( self.file ).strip ( )
                if line.startswith ( "*" ): continue
                self.linesParsed += 1
                # . Check for comments.
                index = line.find ( ";" )
                if index >= 0: line = line[:index].strip ( )
                # . Check for continuations.
                if line.endswith ( "\\" ): line = line[:-1].strip ( ) + " "
                else:                      isNotFound = ( len ( line ) <= 0 )
                found += line
            return found
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            if ( len ( found ) > 0 ): return found
            else:                     raise EOFError

    def IdentifySection ( self, line ):
        """Identify a section of the parameter file."""
        section = None
        tokens = line.split ( " ", 1 )
        if tokens[0].startswith ( "[" ):
            if len ( tokens[0] ) > 1 : head = tokens[0].split ( "[", 1 )[1].split ( "]", 1 )[0].strip()
            else                     : head = tokens[1].split ( "]", 1 )[0].strip()
            # . Proper and improper are both termed "dihedral", but have different function types.
            if head[0:8] == "dihedral" : ( head, line ) = self.GetDihedralSection ( head )
            for name in _SectionMapping.keys():
                if head == name :
                    section = _SectionMapping[name]
                    break
        return ( section, line )

    def Parse ( self, log = logFile, parameters = None ):
        """Parse the data on the file."""
        self.parameters = parameters
        if not self.isParsed:
            # . Initialization.
            self.Initialize ( )
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                current = None
                # . Loop over lines in the file.
                while True:
                    line = self.GetLine ( )
                    ( section, line ) = self.IdentifySection ( line )
                    # . Error or Data from previously found section.
                    if section is None:
                        if current is None: self.Warning ( "Unrecognized section.", False )
                        else:
                            data = self.topData.get ( current, [] )
                            data.append ( line.split ( ) )
                            self.topData[current] = data
                    # . New section found.
                    else:
                        if self.__class__.__name__ == "GromacsDefinitionsFileReader": self.ProcessMoleculeType ( current )
                        current = section
                        # . Include data used to identify dihedral type.
                        if "proper" in current or "ihedral" in current:
                            data = self.topData.get ( current, [] )
                            data.append ( line.split ( ) )
                            self.topData[current] = data
            except EOFError:
                pass
            # . Process the data.
            if self.__class__.__name__ == "GromacsParameterFileReader": 
                # . Determine the force field type to allow correct parsing.
                ffmethod = getattr ( self, "ProcessForceField", None )
                if ffmethod is not None: 
                    ffmethod ( )
                    del self.topData["ForceField"]
                else: self.Warning ( "Unable to determine force field type.", True )
                # . Then, parse the other sections.
                for section in self.topData.keys ( ):
                    if len ( self.topData[section] ) > 0:
                        processmethod = getattr ( self, "Process" + section, None )
                        if processmethod is not None: processmethod ( )
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and other options.
            self.log      = None
            self.isParsed = True

#===================================================================================================================================
# . Gromacs parameter file reader class.
#===================================================================================================================================
class GromacsParameterFileReader ( GromacsFileReader ):
    """Class for Gromacs parameter files that are to be read."""

    _classLabel   = "Gromacs Parameter File Reader"
    _withSections = True

    def Initialize ( self ):
        """Initialization for parsing."""
        # . Read data.
        self.topData        = {}
        # . Processed data.
        self.angle          = {}
        self.bond           = {}
        self.cmap           = {}
        self.ff             = None
        self.proper         = {}
        self.properwild     = {}
        self.improper       = {}
        self.improperwild   = {}
        self.molecules      = {}
        self.nonbond        = {}
        self.nonbond14      = {}
        self.title          = None
        self.ureybradley    = {}

    @classmethod
    def PathToParameters ( selfClass, path, log = logFile ):
        """Return the parameters from a file."""
        # . Get the parameter set.
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToParameters ( ) 

    def ProcessAngle ( self ):
        """Process angles and Urey-Bradley terms."""
        data = self.topData.get ( "Angle", None )
        if data is not None:
            for datum in data:
                try:
                    # . Angle.
                    t1 = datum[0]
                    t2 = datum[1]
                    t3 = datum[2]
                    key = ( max ( t1, t3 ), t2, min ( t1, t3 ) )
                    # . Harmonic term.
                    if key in self.angle: self.Warning ( "Duplicate angle for " + "-".join ( datum[0:3] ) + ".", False )
                    else: self.angle[key] = ( 0.5 * float ( datum[5] ), math.radians ( float ( datum[4] ) ) )
                    # . Urey-Bradley term.
                    if len ( datum ) > 6:
                        if float ( datum[7] ) != 0.0:
                            if key in self.ureybradley: self.Warning ( "Duplicate Urey-Bradley term for " + "-".join ( datum[0:3] ) + ".", False )
                            else: self.ureybradley[key] = ( 0.005 * float ( datum[7] ), 10 * float ( datum[6] ) )
                except:
                    self.Warning ( "Unable to process angle/Urey-Bradley terms: " + " ".join ( datum ) + ".", False )

    def ProcessBond ( self ):
        """Process bonds."""
        data = self.topData.get ( "Bond", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = datum[0]
                    t2 = datum[1]
                    key = ( max ( t1, t2 ), min ( t1, t2 ) )
                    if key in self.bond: self.Warning ( "Duplicate bond for " + "-".join ( datum[0:2] ) + ".", False )
                    # . Conversion => 0.01 / 2 : 0.01 because of nm=>\AA and /2 because of k_b/2.
                    else: self.bond[key] = ( 0.005 * float ( datum[4] ), 10 * float ( datum[3] ) )
                except:
                    self.Warning ( "Unable to process bond: " + " ".join ( datum ) + ".", False )

    def ProcessCMap ( self ):
        """Process cmap terms."""
        data = self.topData.get ( "CMap", None )[:]
        if data is not None:
            while len ( data ) > 0:
                datum = data.pop ( 0 )
                try:
                    t1 = datum[0]
                    t2 = datum[1]
                    t3 = datum[2]
                    t4 = datum[3] # could it have been t4 = t1 ?
                    t5 = t2
                    t6 = t3
                    t7 = t4
                    t8 = datum[4]
                    n  = int ( datum[7] )
                    if   t2 >  t3: key1 = [ t1, t2, t3, t4 ]
                    elif t2 == t3: key1 = [ max ( t1, t4 ), t2, t3, min ( t1, t4 ) ]
                    else:          key1 = [ t4, t3, t2, t1 ]
                    if   t6 >  t7: key2 = [ t5, t6, t7, t8 ]
                    elif t6 == t7: key2 = [ max ( t5, t8 ), t6, t7, min ( t5, t8 ) ]
                    else:          key2 = [ t8, t7, t6, t5 ]
                    key = tuple ( key1 + key2 )
                    # . Calculate abscissa. AKA the (dihedral) angle.
                    increment = 2.0 * math.pi / float ( n )
                    abscissa  = [ -math.pi + i * increment for i in range ( n + 1 ) ]
                    # . Gather data - one row at a time and expand so that upper boundary values are also included.
                    # . First dihedral changes most rapidly (i.e. columnwise storage).
                    cmap = []
                    numbers = datum[8:]
                    for i in range ( n ):
                        tokens = []
                        while len ( tokens ) != n: tokens.append ( numbers.pop ( 0 ) )
                        values = []
                        for token in tokens: values.append ( float ( token ) )
                        # . Duplicate first value for every abscissa angle. GMA: Why? 
                        values.append ( values[0] )
                        cmap.extend ( values )
                    # . Duplicate first column. GMA: Why?
                    cmap.extend ( cmap[0:n+1] )
                    if key in self.cmap: self.Warning ( "Duplicate cmap for " + "-".join ( key[0:4] ) + "/" + "-".join ( key[4:] ) + ".", False )
                    else: self.cmap[key] = ( abscissa, cmap )
                except:
                    self.Warning ( "Unable to process cmap: " + " ".join ( datum ) + ".", False )

    def ProcessDihedral ( self ):
        """Process dihedrals."""
        data = self.topData.get ( "Dihedral", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = datum[0]
                    t2 = datum[1]
                    t3 = datum[2]
                    t4 = datum[3]
                    if   t2 >  t3: key = ( t1, t2, t3, t4 )
                    elif t2 == t3: key = ( max ( t1, t4 ), t2, t3, min ( t1, t4 ) )
                    else:          key = ( t4, t3, t2, t1 )
                    if _DihedralWildCard in key: parameters = self.properwild
                    else:                        parameters = self.proper
                    # . Allow for multiple dihedrals.
                    kchi  = float ( datum[6] ) 
                    n     = int   ( datum[7] )
                    delta = math.radians ( float ( datum[5] ) )
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

    def ProcessForceField ( self ):
        """Process defaults and assign force-field."""
        data = self.topData.get ( "ForceField", None )
        if data is not None:
            if float ( data[0][3] ) < 1.0 :
                if float ( data[0][4] ) > 0.5 : self.ff = "AMBER"
                else                          : 
                    self.ff = "OPLS" 
                    self.Warning ( "Gromacs topologies with OPLS ForceField cannot be read.", True )
            else: self.ff = "CHARMM"

    def ProcessImproper ( self ):
        """Process impropers."""
        data = self.topData.get ( "Improper", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = datum[0]
                    t2 = datum[1]
                    t3 = datum[2]
                    t4 = datum[3]
                    if self.ff == "AMBER":
                        key = ( max ( t1, t4 ), max ( t2, t3 ), min ( t2, t3 ), min ( t1, t4 ) )
                    else:
                        if   t1 >  t4: key = ( t1, t2, t3, t4 )
                        elif t1 == t4: key = ( t1, max ( t2, t3 ), min ( t2, t3 ), t4 )
                        else:          key = ( t4, t3, t2, t1 )
                    if _DihedralWildCard in key: parameters = self.improperwild
                    else:                        parameters = self.improper
                    if key in parameters: self.Warning ( "Duplicate improper for " + "-".join ( datum[0:4] ) + ".", False )
                    else: 
                        if   self.ff == "CHARMM": parameters[key] = ( 0.5 * float ( datum[6] ), math.radians ( float ( datum[5] ) ) ) # . No check on central integer value (always "0").
                        elif self.ff == "AMBER" : parameters[key] = ( float ( datum[6] ), math.radians ( float ( datum[5] ) ), int( datum[7] ) )
                except:
                    self.Warning ( "Unable to process improper: " + " ".join ( datum ) + ".", False )

    def ProcessMolecules ( self ):
        """Process molecules."""
        for key, value in self.topData["Molecules"]: self.molecules[key] = value

    def ProcessNonBond ( self ):
        """Process nonbonds."""
        data = self.topData.get ( "NonBond", None )
        if data is not None:
            for datum in data:
                try:
                    t = datum[0]
                    if t in self.nonbond: self.Warning ( "Duplicate " + tag + " for " + t + ".", False )
                    else:
                        e    = float ( datum[6] )
                        r    = float ( datum[5] ) / 0.1781797 * 2.0
                        mass = float ( datum[2] )
                        self.nonbond[t] = ( math.fabs ( e ), r, mass )
                except:
                    self.Warning ( "Unable to process nonbond: " + " ".join ( datum ) + ".", False )

    def ProcessNonBond14 ( self ):
        """Process nonbonds for 1-4 exclusions."""
        data = self.topData.get ( "NonBond14", None )
        if data is not None:
            for datum in data:
                try:
                    t1 = datum[0]
                    t2 = datum[1]
                    if t1 == t2: 
                        if t1 in self.nonbond14: self.Warning ( "Duplicate nonbond14 for " + t1 + ".", False )
                        else                   : 
                            e = float ( datum[4] )
                            r = float ( datum[3] ) / 0.1781797 * 2.0 
                            self.nonbond14[t1] = ( math.fabs ( e ), r )
                except:
                    self.Warning ( "Unable to process nonbond14: " + " ".join ( datum ) + ".", False )

    def ProcessTitle ( self ):
        """Process the title."""
        self.title = self.topData["Title"][0][0]

    def SummaryItems ( self ):
        """Summary items."""
        # . Parsed data.
        inverted = {}
        items    = [ ( "Section Lines", True ) ]
        for ( key, value ) in _SectionMapping.items ( ): inverted[value] = key
        for section in self.topData.keys ( ):
            if len ( self.topData[section] ) > 0:
                items.append ( ( inverted[section], "{:d}".format ( len ( self.topData[section] ) ) ) )
        # . Processed data.
        items.extend ( [ ( "Parameters"          , True                                        ) ,
                         ( "Angles"              , "{:d}".format ( len ( self.angle        ) ) ) ,
                         ( "Bonds"               , "{:d}".format ( len ( self.bond         ) ) ) ,
                         ( "Cmaps"               , "{:d}".format ( len ( self.cmap         ) ) ) ,
                         ( "Dihedrals (Nonwild)" , "{:d}".format ( len ( self.proper       ) ) ) ,
                         ( "Dihedrals (Wild)"    , "{:d}".format ( len ( self.properwild   ) ) ) ,
                         ( "Impropers (Nonwild)" , "{:d}".format ( len ( self.improper     ) ) ) ,
                         ( "Impropers (Wild)"    , "{:d}".format ( len ( self.improperwild ) ) ) ,
                         ( "Nonbonds"            , "{:d}".format ( len ( self.nonbond      ) ) ) ,
                         ( "Nonbonds 1-4"        , "{:d}".format ( len ( self.nonbond14    ) ) ) ,
                         ( "Urey-Bradleys"       , "{:d}".format ( len ( self.ureybradley  ) ) ) ] )
        return items

    def ToParameters ( self ):
         """Return an object containing all the parameters."""
         parameters = None
         if self.isParsed:
             parameters = GromacsParameterContainer ( )
             for attribute in ( "angle", "bond", "cmap", "ff", "improper", "improperwild", "molecules", "nonbond", "nonbond14", "proper", "properwild", "title", "ureybradley" ):
                 setattr ( parameters, attribute, getattr ( self, attribute ) )
         return parameters

#===================================================================================================================================
# . Gromacs definitions that are to be read and processed.
#===================================================================================================================================
class GromacsDefinitionsFileReader ( GromacsFileReader ):
    """Class for Gromacs definitions that are to be read and processed."""

    _classLabel = "Gromacs Definition File Reader"

    def Initialize ( self ):
        """Initialization for parsing."""
        # . Processed data.
        self.firstMolecule = True
        self.topData       = {}
        # . Read data.
        self.angles        = []
        self.atoms         = []
        self.bonds         = []
        self.cmaps         = []
        self.dihedrals     = []
        self.exclusions    = []
        self.impropers     = []

    def ReInitialize ( self ):
        """Re-Initialization of items read."""
        for item in _MoleculeSections: 
            key = _SectionMapping[item]
            if key in self.topData: self.topData[key] = []

    def ProcessMoleculeType ( self, current ):
        """Process moleculetype section."""
        if current == "Molecule" or current == "Title":
            # . Do not process if this is the first moleculetype hit. Process only when another moleculetype or system is found.
            if self.firstMolecule: 
                self.firstMolecule = False
                self.natoms        = 0
            # . Process current (previously found) moleculetype. 
            else                 :
                self.molname = self.topData["Molecule"].pop ( 0 ) [ 0 ]
                if self.molname in self.parameters.molecules: 
                    nmol = int ( self.parameters.molecules [ self.molname ] )
                    for i in range ( nmol ):
                        self.counter = i
                        for section in _MoleculeSections:
                            key = _SectionMapping [ section ]
                            if key in self.topData:
                                if len ( self.topData [ key ] ) > 0:
                                    processmethod = getattr ( self, "Process" + key, None )
                                    if processmethod is not None: processmethod ( )
                        self.natoms += len ( self.topData["Atoms"] ) 
            
                    self.ReInitialize ( )

    def GetMass ( self, item ):
        """Include atomic mass in atom line."""
        atomType       = item[1]
        ( e, r, mass ) = self.parameters.nonbond[atomType]
        item += [ mass ]
        return item

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile, parameters = None ):
        """Return a system from a Gromacs file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log, parameters = parameters )
        inFile.Summary ( log = log )
        return inFile.ToSystem ( parameters = parameters )

    def ProcessAngles ( self ):
        """Process angles."""
        for item in self.topData["Angles"]: 
            for value in item[0:3]: self.angles.append ( int ( value ) + self.natoms )
            # . Extra patch for online terms.
            if len ( item ) > 4: 
                i1 = int ( item[0] ) - 1
                i2 = int ( item[1] ) - 1
                i3 = int ( item[2] ) - 1
                # . Find respective atomTypes.
                t1 = self.topData["Atoms"][i1][1]
                t2 = self.topData["Atoms"][i2][1]
                t3 = self.topData["Atoms"][i3][1]
                key = ( max ( t1, t3 ), t2, min ( t1, t3 ) )
                # . Check and include new data in parameters.
                if len ( item ) > 6 : 
                    nueq = math.radians ( float ( item[4] ) )
                    nufk = 0.5 *          float ( item[5] )
                    self.parameters.angle[ key ] = ( nufk, nueq )

    def ProcessAtoms ( self ):
        """Process atoms."""
        if hasattr ( self, "counter" )   : counter = self.counter
        else                             : counter = 0
        for item in self.topData["Atoms"]:
            # . Atomic masses for water and ions are not given on atoms section line in some topologies. 
            if len ( item ) == 7: item = self.GetMass ( item )
            newitem = [ int ( item[0] ) + self.natoms ] + item[1:2] + [ int ( item[2] ) + counter ] + item[3:] + [ self.molname ]
            self.atoms.append ( newitem )

    def ProcessBonds ( self ):
        """Process bonds."""
        for item in self.topData["Bonds"]: 
            for value in item[0:2]: self.bonds.append ( int ( value ) + self.natoms )
            # . Extra patch for online terms.
            if len ( item ) > 3: 
                i1 = int ( item[0] ) - 1
                i2 = int ( item[1] ) - 1
                # . Find respective atomTypes.
                t1 = self.topData["Atoms"][i1][1]
                t2 = self.topData["Atoms"][i2][1]
                key = ( max ( t1, t2 ), min ( t1, t2 ) )
                # . Include new data in parameters.
                if len ( item ) > 5 : 
                    nueq = 10    * float ( item[3] )
                    nufk = 0.005 * float ( item[4] )
                    self.parameters.bond[ key ] = ( nufk, nueq )

    def ProcessCMaps ( self ):
        """Process CMAP cross-terms."""
        for item in self.topData["CMaps"]: 
            for value in item[0:5]: self.cmaps.append ( int ( value ) + self.natoms )

    def ProcessDihedrals ( self ):
        """Process dihedrals."""
        for item in self.topData["Dihedrals"]: 
            for value in item[0:4]: self.dihedrals.append ( int ( value ) + self.natoms )

    def ProcessExclusions14 ( self ):
        """Process exclusions."""
        for item in self.topData["Exclusions14"]: 
            for value in item[0:2]: self.exclusions.append ( int ( value ) + self.natoms )

    def ProcessImpropers ( self ):
        """Process impropers."""
        for item in self.topData["Impropers"]: 
            for value in item[0:4]: self.impropers.append ( int ( value ) + self.natoms )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.isParsed:
            items.extend ( [ ( "Angles"         , "{:d}".format ( len ( self.angles     ) // 3 ) ) ,
                             ( "Atoms"          , "{:d}".format ( len ( self.atoms      )      ) ) ,
                             ( "Bonds"          , "{:d}".format ( len ( self.bonds      ) // 2 ) ) ,
                             ( "CMaps"          , "{:d}".format ( len ( self.cmaps      ) // 5 ) ) ,
                             ( "Dihedrals"      , "{:d}".format ( len ( self.dihedrals  ) // 4 ) ) ,
                             ( "Exclusions 1-4" , "{:d}".format ( len ( self.exclusions ) // 2 ) ) ,
                             ( "Impropers"      , "{:d}".format ( len ( self.impropers  ) // 4 ) ) ,
                             ( "MoleculeTypes"  , "{:d}".format ( len ( self.parameters.molecules ) ) ) ] )
        return items

    def ToCMAPDihedralContainer ( self, parameters = {} ):
        """Create a CMAP dihedral container."""
        mm     = None
        nterms = len ( self.cmaps ) // 5
        if nterms > 0 :
            # . Terms.
            atomIndices = []
            uniqueKeys  = set ( )
            termKeys    = []
            for n in range ( nterms ):
                i1  = self.cmaps[5*n  ] - 1
                j1  = self.cmaps[5*n+1] - 1
                k1  = self.cmaps[5*n+2] - 1
                l1  = self.cmaps[5*n+3] - 1
                i2  = j1
                j2  = k1
                k2  = l1 
                l2  = self.cmaps[5*n+4] - 1
                ti1 = self.typesOfAtoms[i1]
                tj1 = self.typesOfAtoms[j1]
                tk1 = self.typesOfAtoms[k1]
                tl1 = self.typesOfAtoms[l1]
                ti2 = self.typesOfAtoms[i2]
                tj2 = self.typesOfAtoms[j2]
                tk2 = self.typesOfAtoms[k2]
                tl2 = self.typesOfAtoms[l2]
                # . Specific key and then wildcard key.
                if   tj1 >  tk1: key1 = [ ti1, tj1, tk1, tl1 ]
                elif tj1 == tk1: key1 = [ max ( ti1, tl1 ), tj1, tk1, min ( ti1, tl1 ) ]
                else:            key1 = [ tl1, tk1, tj1, ti1 ]
                if   tj2 >  tk2: key2 = [ ti2, tj2, tk2, tl2 ]
                elif tj2 == tk2: key2 = [ max ( ti2, tl2 ), tj2, tk2, min ( ti2, tl2 ) ]
                else:            key2 = [ tl2, tk2, tj2, ti2 ]
                key = tuple ( key1 + key2 )
                uniqueKeys.add ( key )
                atomIndices.append ( [ i1, j1, k1, l1, i2, j2, k2, l2 ] )
                termKeys.append ( key )
            # . Parameters.
            uniqueKeys = list ( uniqueKeys )
            uniqueKeys.sort ( )
            keyIndex      = {}
            parameterKeys = []
            splines       = []
            n             =  0
            for key in uniqueKeys:
                ( abscissa, values ) = parameters[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                splines.append ( ( abscissa, abscissa, values ) )
                keyIndex[key] = n
                n += 1
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
            # . Construct the container.
            state = { "label"         : "CMAP Dihedral" ,
                      "parameterKeys" : parameterKeys   ,
                      "parameters"    : splines         ,
                      "terms"         : terms           }
            mm = CMAPDihedralContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToExclusionPairLists ( self ):
        """Create pairlists for the exclusions and 1-4 interactions.
        This is done from the bond list as well as taken from the Gromacs topology.
        """
        # . Initialization.
        exclusions     = None
        interactions14 = None
        nbonds         = len ( self.bonds ) // 2
        if nbonds > 0 :
            # . 1-2 lists.
            i12 = set ( )
            for n in range ( nbonds ):
                i   = self.bonds[2*n  ] - 1
                j   = self.bonds[2*n+1] - 1
                t12 = ( max ( i, j ), min ( i, j ) )
                i12.add ( t12 )
            # . Generate sorted connections for each atom.
            connections = [ [] for n in range ( self.natoms ) ]
            for ( i, j ) in i12:
                connections[i].append ( j )
                connections[j].append ( i )
            for n in range ( self.natoms ): connections[n].sort ( )
            # . 1-3 lists (i > k automatically).
            i13 = set ( )
            for j in range ( self.natoms ):
                jbonds = connections[j]
                nj     = len ( jbonds )
                for i in range ( 1, nj ):
                    for k in range ( 0, i ): i13.add ( ( jbonds[i], jbonds[k] ) )
            # . Generate 123 exclusions.
            i123 = i12.union  ( i13 )
            # . Add in explicitly defined exclusions.
            # . For strict reproduction of Gromacs topologies, only the 1-4 "Exclusions14" or self.exclusions list should be used.
            nexclu14 = len ( self.exclusions ) // 2
            ex14 = set ( )
            if nexclu14 > 0:
                for n in range ( nexclu14 ):
                    i = self.exclusions[2*n  ] - 1
                    j = self.exclusions[2*n+1] - 1
                    ex14.add ( ( max ( i, j ), min ( i, j ) ) )
            # . Possible 1-4 lists.
            i14p = set ( )
            for ( j, k ) in i12:
                jbonds = connections[j]
                kbonds = connections[k]
                nj     = len ( jbonds )
                nk     = len ( kbonds )
                for i in jbonds:
                    if ( i != k ):
                        for l in kbonds:
                            if ( l != i ) and ( l != j ): i14p.add ( ( max ( i, l ), min ( i, l ) ) )
            # . Get the final lists.
            i1234 = i123.union  ( ex14 )
            i1234 = i1234.union ( i14p )
            temp  = ex14.difference ( i123 )
            i14   = i14p.difference ( i123 )
            i14   = i14.union ( temp )
            # . Generate the lists.
            if len ( i1234 ) > 0:
                i1234 = list ( itertools.chain ( * list ( i1234 ) ) )
                exclusions = SelfPairList.FromIndexPairs ( i1234 )
                exclusions.label = "Exclusions"
            if len ( i14 ) > 0:
                i14 = list ( itertools.chain ( * list ( i14 ) ) )
                interactions14 = SelfPairList.FromIndexPairs ( i14 )
                interactions14.label = "1-4 Interactions"
            return ( exclusions, interactions14 )

    def ToFourierDihedralContainer ( self, parameters = {}, parameterswild = {} ):
        """Create a harmonic dihedral container."""
        mm     = None
        nterms = len ( self.dihedrals ) // 4
        if nterms >0 :
            # . Terms.
            atomIndices = []
            uniqueKeys  = set ( )
            termKeys    = []
            for n in range ( nterms ):
                i  = self.dihedrals[4*n  ] - 1
                j  = self.dihedrals[4*n+1] - 1
                k  = self.dihedrals[4*n+2] - 1
                l  = self.dihedrals[4*n+3] - 1
                ti = self.typesOfAtoms[i]
                tj = self.typesOfAtoms[j]
                tk = self.typesOfAtoms[k]
                tl = self.typesOfAtoms[l]
                # . Specific key and then wildcard key.
                if   tj >  tk: key = ( ti, tj, tk, tl )
                elif tj == tk: key = ( max ( ti, tl ), tj, tk, min ( ti, tl ) )
                else:          key = ( tl, tk, tj, ti )
                if key in parameters:
                    nperiods = len ( parameters[key] )
                else:
                    keyw0  = ( _DihedralWildCard, key[1], key[2], key[3] )
                    keyw3  = ( key[0], key[1], key[2], _DihedralWildCard )
                    keyw03 = ( _DihedralWildCard, max ( tj, tk ), min ( tj, tk ), _DihedralWildCard )
                    keys = [ keyw03, keyw3, keyw0 ]
                    isOK = False
                    for key in keys: 
                        if key in parameterswild: 
                            nperiods = len ( parameterswild[key] )
                            isOK = True
                            break
                    if not isOK: self.Warning ( "Unable to assign parameters for dihedral: " + key + " for centers " + [ i, j, k, l ] + ".", True )
                uniqueKeys.add ( key )
                for p in range ( nperiods ):
                    atomIndices.append ( [ i, j, k, l ] )
                    termKeys.append ( ( key, p ) )
            # . Parameters.
            # . Care is taken to unravel multiple dihedrals.
            uniqueKeys = list ( uniqueKeys )
            uniqueKeys.sort ( )
            parameterKeys   = []
            localParameters = []
            keyIndex        = {}
            n               =  0
            for key in uniqueKeys:
                data = parameters.get ( key, None )
                if data is None: data = parameterswild[key]
                baseKey = _ParameterKeyFieldSeparator.join ( key )
                for ( nperiod, datum ) in enumerate ( data ):
                    keyIndex[ ( key, nperiod )] = n
                    parameterKeys.append ( baseKey + _ParameterKeyFieldSeparator + "{:d}".format ( datum[0] ) )
                    localParameters.append ( ( datum[1], datum[0], datum[2] ) )
                    n += 1
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
            # . Construct the container.
            state = { "label"         : "Fourier Dihedral" ,
                      "parameterKeys" : parameterKeys      ,
                      "parameters"    : localParameters    ,
                      "terms"         : terms              }
            mm = FourierDihedralContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToFourierImproperContainer ( self, parameters = {}, parameterswild = {} ):
        """Create a periodic improper container for AMBER topologies."""
        mm     = None
        nterms = len ( self.impropers ) // 4
        if nterms > 0 :
            # . Terms.
            atomIndices = []
            termKeys    = []
            for n in range ( nterms ):
                i  = self.impropers[4*n  ] - 1
                j  = self.impropers[4*n+1] - 1
                k  = self.impropers[4*n+2] - 1 # . Center
                l  = self.impropers[4*n+3] - 1
                ti = self.typesOfAtoms[i]
                tj = self.typesOfAtoms[j]
                tk = self.typesOfAtoms[k]
                tl = self.typesOfAtoms[l]
                # . Specific key and then wildcard key.
                ttent = [ ti, tj, tl ]
                ttent.sort()
                tlist = ttent[0:2] + [ tk ] + ttent[2:]
                key = ( max ( ti, tl ), max ( tj, tk ), min ( tj, tk ), min ( ti, tl ) )
                if key not in parameters:
                    isOK = False
                    for n in [ ti, tl ]:
                        key = ( _DihedralWildCard, max ( tj, tk ), min ( tj, tk ), n )
                        if key in parameterswild: 
                            isOK = True
                            break
                        for m in [ tj, tk ]:
                            key = ( _DihedralWildCard, _DihedralWildCard, m, n )
                            if key in parameterswild: 
                                isOK = True
                                break
                        if isOK: break
                atomIndices.append ( [ i, j, k, l ] )
                termKeys.append ( key )
            # . Parameters.
            uniqueKeys = list ( set ( termKeys ) )
            uniqueKeys.sort ( )
            parameterKeys   = []
            localParameters = []
            keyIndex        = {}
            for ( i, key ) in enumerate ( uniqueKeys ):
                keyIndex[key] = i
                datum = parameters.get ( key, None )
                if datum is None: datum = parameterswild[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                localParameters.append ( ( datum[0], datum[2], datum[1] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
            # . Construct the container.
            state = { "label"         : "Fourier Improper" ,
                      "parameterKeys" : parameterKeys      ,
                      "parameters"    : localParameters    ,
                      "terms"         : terms              }
            mm = FourierDihedralContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicAngleContainer ( self, parameters = {} ):
        """Create a harmonic angle container."""
        mm     = None
        nterms = len ( self.angles ) // 3 
        if nterms > 0 :
            # . Terms.
            atomIndices = []
            termKeys    = []
            for n in range ( nterms ):
                i  = self.angles[3*n  ] - 1
                j  = self.angles[3*n+1] - 1
                k  = self.angles[3*n+2] - 1
                ti = self.typesOfAtoms[i]
                tj = self.typesOfAtoms[j]
                tk = self.typesOfAtoms[k]
                key = ( max ( ti, tk ), tj, min ( ti, tk ) )
                atomIndices.append ( [ i, j, k ] )
                termKeys.append ( key )
            # . Parameters.
            uniqueKeys = list ( set ( termKeys ) )
            uniqueKeys.sort ( )
            parameterKeys   = []
            localParameters = []
            keyIndex        = {}
            for ( i, key ) in enumerate ( uniqueKeys ):
                keyIndex[key] = i
                datum = parameters[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                localParameters.append ( ( datum[1], datum[0] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
            # . Construct the container.
            state = { "parameterKeys" : parameterKeys   ,
                      "parameters"    : localParameters ,
                      "terms"         : terms           }
            mm = HarmonicAngleContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicBondContainer ( self, parameters = {} ):
        """Create a harmonic bond container."""
        mm     = None
        nterms = len ( self.bonds ) // 2
        if nterms > 0 :
            # . Terms.
            atomIndices = []
            termKeys    = []
            for n in range ( nterms ):
                i  = self.bonds[2*n  ] - 1
                j  = self.bonds[2*n+1] - 1
                ti = self.typesOfAtoms[i]
                tj = self.typesOfAtoms[j]
                key = ( max ( ti, tj ), min ( ti, tj ) )
                atomIndices.append ( [ i, j ] )
                termKeys.append ( key )
            # . Parameters.
            uniqueKeys = list ( set ( termKeys ) )
            uniqueKeys.sort ( )
            localParameters = []
            parameterKeys   = []
            keyIndex        = {}
            for ( i, key ) in enumerate ( uniqueKeys ):
                keyIndex[key] = i
                datum = parameters[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                localParameters.append ( ( datum[1], datum[0] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
            # . Construct the container.
            state = { "parameterKeys" : parameterKeys   ,
                      "parameters"    : localParameters ,
                      "terms"         : terms           }
            mm = HarmonicBondContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicImproperContainer ( self, parameters = {}, parameterswild = {} ):
        """Create a harmonic improper container for CHARMM topologies."""
        mm     = None
        nterms = len ( self.impropers ) // 4
        if nterms > 0 :
            # . Terms.
            atomIndices = []
            termKeys    = []
            for n in range ( nterms ):
                i  = self.impropers[4*n  ] - 1
                j  = self.impropers[4*n+1] - 1
                k  = self.impropers[4*n+2] - 1
                l  = self.impropers[4*n+3] - 1
                ti = self.typesOfAtoms[i]
                tj = self.typesOfAtoms[j]
                tk = self.typesOfAtoms[k]
                tl = self.typesOfAtoms[l]
                # . Specific key and then wildcard key.
                if   ti >  tl: key = ( ti, tj, tk, tl )
                elif ti == tl: key = ( ti, max ( tj, tk ), min ( tj, tk ), tl )
                else:          key = ( tl, tk, tj, ti )
                if key not in parameters: key = ( max ( ti, tl ), _DihedralWildCard, _DihedralWildCard, min ( ti, tl ) )
                atomIndices.append ( [ i, j, k, l ] )
                termKeys.append ( key )
            # . Parameters.
            uniqueKeys = list ( set ( termKeys ) )
            uniqueKeys.sort ( )
            parameterKeys   = []
            localParameters = []
            keyIndex        = {}
            for ( i, key ) in enumerate ( uniqueKeys ):
                keyIndex[key] = i
                datum = parameters.get ( key, None )
                if datum is None: datum = parameterswild[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                localParameters.append ( ( datum[1], datum[0] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
            # . Construct the container.
            state = { "parameterKeys" : parameterKeys   ,
                      "parameters"    : localParameters ,
                      "terms"         : terms           }
            mm = HarmonicImproperContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicUreyBradleyContainer ( self, parameters = {} ):
        """Create a harmonic Urey-Bradley container.

        These are only assigned if the parameters are specified.
        """
        mm     = None
        nterms = len ( self.angles ) // 3
        if ( nterms > 0 ) and ( len ( parameters ) > 0 ):
            # . Terms.
            atomIndices = []
            termKeys    = []
            for n in range ( nterms ):
                i  = self.angles[3*n  ] - 1
                j  = self.angles[3*n+1] - 1
                k  = self.angles[3*n+2] - 1
                ti = self.typesOfAtoms[i]
                tj = self.typesOfAtoms[j]
                tk = self.typesOfAtoms[k]
                key = ( max ( ti, tk ), tj, min ( ti, tk ) )
                if key in parameters:
                    atomIndices.append ( [ i, k ] )
                    termKeys.append ( key )
            numterms = len ( termKeys )
            if numterms > 0:
                # . Parameters.
                uniqueKeys = list ( set ( termKeys ) )
                uniqueKeys.sort ( )
                parameterKeys   = []
                localParameters = []
                keyIndex        = {}
                for ( i, key ) in enumerate ( uniqueKeys ):
                    keyIndex[key] = i
                    datum = parameters[key]
                    parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                    localParameters.append ( ( datum[1], datum[0] ) )
                # . Terms.
                terms = []
                for ( i, indices ) in enumerate ( atomIndices ):
                    terms.append ( tuple ( indices + [ keyIndex[termKeys[i]], True ] ) )
                # . Construct the container.
                state = { "is12Interaction" : False           ,
                          "label"           : "Urey-Bradley"  ,
                          "parameterKeys"   : parameterKeys   ,
                          "parameters"      : localParameters ,
                          "terms"           : terms           }
                mm = HarmonicBondContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
        return mm

    def ToLJParameterContainers ( self, parameters = {}, parameters14 = {} ):
        """Create Lennard-Jones parameter containers."""
        ljparametes    = None
        ljParameters14 = None
        ntypes = len ( self.uniqueAtomTypes )
        if ntypes > 0:
            epsilons   = []
            epsilons14 = []
            sigmas     = []
            sigmas14   = []
            # . Loop over atom types.
            for atomType in self.uniqueAtomTypes:
                ( e, r, mass ) = parameters[atomType]
                data14         = parameters14.get ( atomType, None )
                epsilons.append ( e )
                sigmas.append   ( r )
                if data14 is None:
                    epsilons14.append ( e )
                    sigmas14.append   ( r )
                else:
                    epsilons14.append ( data14[0] )
                    sigmas14.append   ( data14[1] )
            # . This will have to be changed to OPLS, since its analytic form is different
            ljParameters   = LJParameterContainer.FromEpsilonSigma ( ntypes, epsilons  , sigmas  , analyticForm = LJForm.Amber )
            ljParameters14 = LJParameterContainer.FromEpsilonSigma ( ntypes, epsilons14, sigmas14, analyticForm = LJForm.Amber )
            ljParameters14.label = "1-4 LJ Parameters"
        return ( ljParameters, ljParameters14 )

    def ToMMAtomData ( self, mmState ):
        """Add atom data to a MM state."""
        if self.natoms > 0 :
            # . Gather some data.
            atomCharges = []
            atomTypes   = []
            for datum in self.atoms:
                atomCharges.append ( float ( datum[6] ) )
                atomTypes.append   ( datum[1] )
            # . Get unique types.
            uniqueTypes     = set ( atomTypes )
            uniqueAtomTypes = list ( uniqueTypes )
            uniqueAtomTypes.sort ( )
            # . Generate an index.
            uniqueAtomTypeIndex = {}
            for ( i, atomType ) in enumerate ( uniqueAtomTypes ): uniqueAtomTypeIndex[atomType] = i
            # . Save the atom types.
            self.typesOfAtoms    = atomTypes
            self.uniqueAtomTypes = uniqueAtomTypes
            # . Create the data.
            atomTypeIndices = Array.WithExtent ( self.natoms, dataType = DataType.Integer )
            charges         = Array.WithExtent ( self.natoms )
            for ( i, ( charge, label ) ) in enumerate ( zip ( atomCharges, atomTypes ) ):
                atomTypeIndices[i] = uniqueAtomTypeIndex[label]
                charges        [i] = charge
            # . Assign the data.
            mmState.atomTypeIndices = atomTypeIndices
            mmState.atomTypes       = uniqueAtomTypes
            mmState.charges         = charges
            mmState.ljTypeIndices   = Clone ( atomTypeIndices )

    def ToSequence ( self ):
        """Generate a sequence."""
        # . Create atom paths - chainID (datum[8]), resname (datum[3]), resID (datum[2]), atom name (datum[4]).
        majorSeparator = Sequence._attributable["labelSeparator"]
        minorSeparator = Sequence._attributable["fieldSeparator"]
        atomPaths      = []
        for datum in self.atoms:
            atomPaths.append ( datum[8] + majorSeparator + datum[3] + minorSeparator + str ( datum[2] ) + majorSeparator + datum[4] )
        # . Get the sequence.
        return Sequence.FromAtomPaths ( atomPaths )

    def ToSystem ( self, parameters = None ):
        """Create a system."""
        system = None
        if self.isParsed:
            # . Sequence data and then basic system.
            sequence     = self.ToSequence ( )
            system       = System.FromSequence ( sequence )
            system.label = " ".join ( self.parameters.title )
            # . Assign atomic numbers from masses.
            for ( datum, atom ) in zip ( self.atoms, system.atoms ):
                atom.atomicNumber = PeriodicTable.AtomicNumberFromMass ( float ( datum[7] ) )
            # . The MM model.
            if parameters is not None: self.ff = parameters.ff
            if   self.ff == "AMBER": mmModel = MMModelAMBER  ( )
            elif self.ff == "OPLS" : mmModel = MMModelOPLS   ( )
            else                   : mmModel = MMModelCHARMM ( )
            mmState = MMModelState ( )
            # . Build the model if parameters are available.
            if ( mmModel is not None ) and ( mmState is not None ) and ( parameters is not None ):
                # . The MM atom data.
                self.ToMMAtomData ( mmState )
                # . Various MM terms.
                mmTerms = mmState.mmTerms
                mm      = self.ToHarmonicBondContainer        ( parameters = parameters.bond  )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToHarmonicAngleContainer       ( parameters = parameters.angle )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToHarmonicUreyBradleyContainer ( parameters = parameters.ureybradley )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToFourierDihedralContainer     ( parameters = parameters.proper, parameterswild = parameters.properwild )
                if mm is not None: mmTerms.append ( mm )
                if self.ff == "CHARMM" : 
                   mm   = self.ToHarmonicImproperContainer    ( parameters = parameters.improper, parameterswild = parameters.improperwild )
                if self.ff == "AMBER"  : 
                   mm   = self.ToFourierImproperContainer     ( parameters = parameters.improper, parameterswild = parameters.improperwild )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToCMAPDihedralContainer        ( parameters = parameters.cmap )
                if mm is not None: mmTerms.append ( mm )
                # . Exclusions.
                ( exclusions, interactions14 ) = self.ToExclusionPairLists ( )
                if exclusions     is not None: mmState.exclusions     = exclusions
                if interactions14 is not None: mmState.interactions14 = interactions14
                # . LJ terms.
                ( ljParameters, ljParameters14 ) = self.ToLJParameterContainers ( parameters = parameters.nonbond, parameters14 = parameters.nonbond14 )
                if ljParameters is not None: mmState.ljParameters = ljParameters
                # . Apply proper FF _Scale factors for 1-4 pairs.
                if ljParameters14 is not None:
                    scale = getattr ( mmModel, "lennardJonesScale14", 1.0 )
                    if scale != 0.0:
                        ljParameters14.Scale ( scale )
                        mmState.ljParameters14 = ljParameters14
                # . Assign the model.
                mmState.target                   = system
                system._energyModels["_mmModel"] = EnergyModelPriority.MMModel
                system.__dict__     ["_mmModel"] = mmModel
                system.__dict__     [ "mmState"] = mmState
                system._UpdateEnergyClosures ( )
                # . Complete the system connectivity from the MM model.
                system.mmModel.CompleteConnectivity  ( system )
        return system

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass


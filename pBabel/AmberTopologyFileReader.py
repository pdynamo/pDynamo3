"""Classes and functions for reading Amber topology files."""

import math

from  pCore                import Clone                      , \
                                  DataType                   , \
                                  logFile                    , \
                                  LogFileActive              , \
                                  SelfPairList               , \
                                  TextFileReader
from  pMolecule            import EnergyModelPriority        , \
                                  Sequence                   , \
                                  System
from  pMolecule.MMModel    import FourierDihedralContainer   , \
                                  HarmonicAngleContainer     , \
                                  HarmonicBondContainer      , \
                                  LJParameterContainer       , \
                                  MMModelAMBER               , \
                                  MMModelState
from  pScientific          import PeriodicTable              , \
                                  Units
from  pScientific.Arrays   import Array
from  pScientific.Symmetry import CrystalSystemCubic         , \
                                  CrystalSystemMonoclinic    , \
                                  CrystalSystemOrthorhombic  , \
                                  CrystalSystemTetragonal    , \
                                  PeriodicBoundaryConditions
from .ExportImport         import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The constant used to scale charges in AMBER topology files.
_ELECTRONTOKCAL = 18.2223

# . The cutOff for checking system neutrality. This is the same as in AMBER (sander/ew_setup) for consistency.
_NEUTRALITYCUTOFF = 1.0e-2

# . The number of pointers in the topology file.
_NUMBEROFPOINTERS = 32 # . 31 before Amber12.

# . The index number of the various pointers.
_NATOMS =  0
_NTYPES =  1
_NBONH  =  2
_MBONA  =  3
_NTHETH =  4
_MTHETA =  5
_NPHIH  =  6
_MPHIA  =  7
_NHPARM =  8
_NPARM  =  9
_NNB    = 10
_NRES   = 11
_NBONA  = 12
_NTHETA = 13
_NPHIA  = 14
_NUMBND = 15
_NUMANG = 16
_NPTRA  = 17
_NATYP  = 18
_NPHB   = 19
_IFPERT = 20
_NBPER  = 21
_NGPER  = 22
_NDPER  = 23
_MBPER  = 24
_MGPER  = 25
_MDPER  = 26
_IFBOX  = 27
_NMXRS  = 28
_IFCAP  = 29
_NUMEXT = 30
_NCOPY  = 31

# . Names corresponding to the index numbers.
_NAMES = [ "natoms" ,
           "ntypes" ,
           "nbonh"  ,
           "mbona"  ,
           "ntheth" ,
           "mtheta" ,
           "nphih"  ,
           "mphia"  ,
           "nhparm" ,
           "nparm"  ,
           "nnb"    ,
           "nres"   ,
           "nbona"  ,
           "ntheta" ,
           "nphia"  ,
           "numbnd" ,
           "numang" ,
           "nptra"  ,
           "natyp"  ,
           "nphb"   ,
           "ifpert" ,
           "nbper"  ,
           "ngper"  ,
           "ndper"  ,
           "mbper"  ,
           "mgper"  ,
           "mdper"  ,
           "ifbox"  ,
           "nmxrs"  ,
           "ifcap"  ,
           "numext" ,
           "ncopy"  ]

# . Dimension functions.
def _YX         ( n ): return n
def _YXSquared  ( n ): return n * n
def _YXTriangle ( n ): return ( n * ( n + 1 ) ) // 2
def _YX3        ( n ): return 3 * n
def _YX4        ( n ): return 4 * n
def _YX5        ( n ): return 5 * n

# . In AMBER 12 - two additional sections after DIHEDRAL_PHASE are SCEE_SCALE_FACTOR
# . and SCNB_SCALE_FACTOR. One should also maybe take care of three sections (RADIUS_SET,
# . RADII and SCREEN) that occur after IROTAT but before the IFBOX option.

# . It might be better to allow sections in any order so that unknown sections could be ignored.

# . The section data (other than TITLE and POINTERS).
_SECTIONDATA = { "ATOM_NAME"                  : ( _NATOMS, _YX         ) ,
                 "CHARGE"                     : ( _NATOMS, _YX         ) ,
                 "ATOMIC_NUMBER"              : ( _NATOMS, _YX         ) ,
                 "MASS"                       : ( _NATOMS, _YX         ) ,
                 "ATOM_TYPE_INDEX"            : ( _NATOMS, _YX         ) ,
                 "NUMBER_EXCLUDED_ATOMS"      : ( _NATOMS, _YX         ) ,
                 "NONBONDED_PARM_INDEX"       : ( _NTYPES, _YXSquared  ) ,
                 "RESIDUE_LABEL"              : ( _NRES,   _YX         ) ,
                 "RESIDUE_POINTER"            : ( _NRES,   _YX         ) ,
                 "BOND_FORCE_CONSTANT"        : ( _NUMBND, _YX         ) ,
                 "BOND_EQUIL_VALUE"           : ( _NUMBND, _YX         ) ,
                 "ANGLE_FORCE_CONSTANT"       : ( _NUMANG, _YX         ) ,
                 "ANGLE_EQUIL_VALUE"          : ( _NUMANG, _YX         ) ,
                 "DIHEDRAL_FORCE_CONSTANT"    : ( _NPTRA,  _YX         ) ,
                 "DIHEDRAL_PERIODICITY"       : ( _NPTRA,  _YX         ) ,
                 "DIHEDRAL_PHASE"             : ( _NPTRA,  _YX         ) ,
                 "SCEE_SCALE_FACTOR"          : ( _NPTRA,  _YX         ) ,
                 "SCNB_SCALE_FACTOR"          : ( _NPTRA,  _YX         ) ,
                 "SOLTY"                      : ( _NATYP,  _YX         ) ,
                 "LENNARD_JONES_ACOEF"        : ( _NTYPES, _YXTriangle ) ,
                 "LENNARD_JONES_BCOEF"        : ( _NTYPES, _YXTriangle ) ,
                 "BONDS_INC_HYDROGEN"         : ( _NBONH,  _YX3        ) ,
                 "BONDS_WITHOUT_HYDROGEN"     : ( _NBONA,  _YX3        ) ,
                 "ANGLES_INC_HYDROGEN"        : ( _NTHETH, _YX4        ) ,
                 "ANGLES_WITHOUT_HYDROGEN"    : ( _NTHETA, _YX4        ) ,
                 "DIHEDRALS_INC_HYDROGEN"     : ( _NPHIH,  _YX5        ) ,
                 "DIHEDRALS_WITHOUT_HYDROGEN" : ( _NPHIA,  _YX5        ) ,
                 "EXCLUDED_ATOMS_LIST"        : ( _NNB,    _YX         ) ,
                 "HBOND_ACOEF"                : ( _NPHB,   _YX         ) ,
                 "HBOND_BCOEF"                : ( _NPHB,   _YX         ) ,
                 "HBCUT"                      : ( _NPHB,   _YX         ) ,
                 "AMBER_ATOM_TYPE"            : ( _NATOMS, _YX         ) ,
                 "TREE_CHAIN_CLASSIFICATION"  : ( _NATOMS, _YX         ) ,
                 "JOIN_ARRAY"                 : ( _NATOMS, _YX         ) ,
                 "IROTAT"                     : ( _NATOMS, _YX         ) ,
                 "RADIUS_SET"                 : ( -1     , _YX         ) ,
                 "RADII"                      : ( _NATOMS, _YX         ) ,
                 "SCREEN"                     : ( _NATOMS, _YX         ) ,
                 "IPOL"                       : ( -1,      _YX         ) ,
                 "SOLVENT_POINTERS"           : ( -1,      _YX         ) ,
                 "ATOMS_PER_MOLECULE"         : ( -1,      _YX         ) ,
                 "BOX_DIMENSIONS"             : ( -1,      _YX         ) }

#===================================================================================================================================
# . AmberTopology file reader class.
#===================================================================================================================================
class AmberTopologyFileReader ( TextFileReader ):
    """AmberTopologyFileReader is the class for AmberTopology files that are to be read."""

    _attributable = dict ( TextFileReader._attributable )
    _classLabel   = "Amber Topology File Reader"
    _attributable.update ( { "QUNTANGLE" : False ,
                             "QVERIFIED" : False } )

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Version line.
                self.ParseVersionLine ( )
                # . Title.
                self.title = self.ParseSection ( "TITLE", -1 )
                # . Pointers.
                self.pointers = self.ParseSection ( "POINTERS", _NUMBEROFPOINTERS )
                # . Loop over remaining sections.
                while True:
                    # . The FLAG line.
                    try:                  self.line = next ( self.file )
                    except StopIteration: break
                    words       = self.line.split ( )
                    sectionName = words[1]
                    if words[0] != "%FLAG": self.Warning ( "Invalid section header: " + self.line.rstrip ( ), True )
                    data = _SECTIONDATA.get ( sectionName, None )
                    if data is None: self.Warning ( "Unknown section name: " + sectionName, True )
                    ( pointerIndex, dimensioner ) = data
                    if sectionName == "ATOMS_PER_MOLECULE": numberOfItems = self.solvent_pointers[1]
                    elif pointerIndex < 0:                  numberOfItems = -1
                    else:                                   numberOfItems = self.pointers[pointerIndex]
                    setattr ( self, sectionName.lower ( ), self.ParseSectionData ( sectionName, dimensioner ( numberOfItems ) ) )
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseFormat ( self, format ):
        """Parse a format string of the form (20A4) or (10I8) or (5E16.8) or (5F16.8)."""
        try:
            words = format.strip ( ).upper ( ).split ( "." )
            for ( separator, converter, default ) in [ ( "A", None, None ), ( "E", float, 0.0 ), ( "F", float, 0.0 ), ( "I", int, 0 ) ]:
                integers = words[0].split ( separator )
                if len ( integers ) == 2:
                    number = int ( integers[0] )
                    width  = int ( integers[1] )
                    return ( number, width, converter, default )
        except:
            pass

    def ParseSection ( self, sectionName, numberOfItems ):
        """Parse a full section."""
        # . The FLAG line.
        self.line   = next ( self.file )
        words       = self.line.split ( )
        sectionName = words[1]
        if ( words[0] != "%FLAG" ) or ( words[1] != sectionName ): self.Warning ( "Invalid " + sectionName + " section header: " + self.line.rstrip ( ), True )
        # . Data.
        return self.ParseSectionData ( sectionName, numberOfItems )

    def ParseSectionData ( self, sectionName, numberOfItems ):
        """Parse a section."""
        items = None
        # . The FORMAT line.
        self.line = next ( self.file ).strip ( )
        if ( self.line[0:8] != "%FORMAT(" ) or ( self.line[-1] != ")" ): self.Warning ( "Invalid " + sectionName + " format line: " + self.line.rstrip ( ), True )
        try:
            ( number, width, converter, default ) = self.ParseFormat ( self.line[8:-1] )
        except:
            self.Warning ( "Invalid " + sectionName + " format: " + self.line.rstrip ( ), True )
        # . The data lines.
        if numberOfItems == 0:
            next ( self.file ) # . The line is empty for zero entries.
        else:
            if numberOfItems < 0: numberOfItems = number
            items = self.GetFixedFormatArray ( numberOfItems, number, width, converter = converter, default = default )
        return items

    def ParseVersionLine ( self ):
        """Parse the version line."""
        self.line = next ( self.file )
        words = self.line.split ( )
        if ( words[0] != "%VERSION" ): self.Warning ( "Invalid version line: " + self.line.rstrip ( ), True )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile, mmModel = None ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Verify  ( )
        inFile.Summary ( log = log )
        return inFile.ToSystem ( mmModel = mmModel )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.isParsed:
            items.extend ( [ ( name, "{:d}".format ( i ) ) for ( name, i ) in zip ( _NAMES, self.pointers ) ] )
        return items

    def ToExclusionsPairList ( self ):
        """Create an exclusions pairlist."""
        mm = None
        if self.pointers[_NNB] > 0:
            ijex = []
            nnb  = 0
            for i in range ( self.pointers[_NATOMS] ):
                n = self.number_excluded_atoms[i]
                for j in self.excluded_atoms_list[nnb:nnb+n]:
                    if j > 0:
                        ijex.append ( i )
                        ijex.append ( j - 1 )
                nnb += n
            mm = SelfPairList.FromIndexPairs ( ijex )
            mm.label = "Exclusions"
        return mm

    def ToHarmonicAngleContainer ( self ):
        """Create a harmonic angle container."""
        mm          = None
        nterms      = self.pointers[_NTHETH] + self.pointers[_NTHETA]
        nparameters = self.pointers[_NUMANG]
        if ( nterms > 0 ) and ( nparameters > 0 ):
            # . Parameters.
            parameters = []
            for ( eq, fC ) in zip ( self.angle_equil_value, self.angle_force_constant ):
                parameters.append ( ( eq, fC ) )
            # . Terms.
            data  = self.angles_inc_hydrogen + self.angles_without_hydrogen
            terms = []
            for i in range ( nterms ):
                terms.append ( ( data[4*i  ] // 3 ,
                                 data[4*i+1] // 3 ,
                                 data[4*i+2] // 3 ,
                                 data[4*i+3] - 1 ,
                                 True ) )
            state = { "parameters" : parameters, "terms" : terms }
            mm = HarmonicAngleContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicBondContainer ( self ):
        """Create a harmonic bond container."""
        mm          = None
        nterms      = self.pointers[_NBONH] + self.pointers[_NBONA]
        nparameters = self.pointers[_NUMBND]
        if ( nterms > 0 ) and ( nparameters > 0 ):
            # . Parameters.
            parameters = []
            for ( eq, fC ) in zip ( self.bond_equil_value, self.bond_force_constant ):
                parameters.append ( ( eq, fC ) )
            # . Terms.
            data  = self.bonds_inc_hydrogen + self.bonds_without_hydrogen
            terms = []
            for i in range ( nterms ):
                terms.append ( ( data[3*i  ] // 3 ,
                                 data[3*i+1] // 3 ,
                                 data[3*i+2] - 1 ,
                                 True ) )
            state = { "parameters" : parameters, "terms" : terms }
            mm = HarmonicBondContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToLJParameterContainer ( self, analyticForm ):
        """Create a Lennard-Jones parameter container."""
        mm     = None
        ntypes = self.pointers[_NTYPES]
        if ntypes > 0: mm = LJParameterContainer.FromTableCoefficients ( ntypes, self.lennard_jones_acoef, self.lennard_jones_bcoef, analyticForm = analyticForm )
        return mm

    def ToMMAtomData ( self, mmState ):
        """Add atom data to a MM state."""
        n = self.pointers[_NATOMS]
        if n > 0:
            # . Find the unique types.
            atomTypes = list ( set ( self.amber_atom_type ) )
            atomTypes.sort ( )
            # . Find the type index mapping.
            indices   = []
            for atomType in self.amber_atom_type:
                indices.append ( atomTypes.index ( atomType ) )
            # . Create the data.
            atomTypeIndices = Array.WithExtent ( n , dataType = DataType.Integer )
            ljTypeIndices   = Array.WithExtent ( n , dataType = DataType.Integer )
            charges         = Array.WithExtent ( n )
            for ( i, ( aType, charge, ljType ) ) in enumerate ( zip ( indices, self.charge, self.atom_type_index ) ):
                atomTypeIndices[i] = aType
                charges        [i] = charge
                ljTypeIndices  [i] = ljType - 1
            # . Assign the data.
            mmState.atomTypeIndices = atomTypeIndices
            mmState.atomTypes       = atomTypes
            mmState.charges         = charges
            mmState.ljTypeIndices   = ljTypeIndices

    def ToSequence ( self ):
        """Assign a sequence to a system."""
        # . Check for residues.
        nresidues = self.pointers[_NRES]
        if nresidues > 0:
            # . Intialization.
            majorSeparator = Sequence._attributable["labelSeparator"]
            minorSeparator = Sequence._attributable["fieldSeparator"]
            # . Get the atom names and residue labels.
            atomnames = self.atom_name
            labels    = self.residue_label
            # . Get pointers - suitably adjusted.
            pointers = self.residue_pointer
            pointers.append ( self.pointers[_NATOMS] + 1 )
            # . Initialization.
            atomPaths = []
            # . Loop over residues.
            for ires in range ( nresidues ):
                # . Generate a conforming residue/entity label.
                label  = labels[ires]
                tokens = label.split ( majorSeparator, 1 )
                componentLabel = tokens.pop ( )
                if len ( tokens ) > 0: entityLabel = tokens.pop ( ).strip ( )
                else:                  entityLabel = ""
                if componentLabel.find ( minorSeparator ) < 0: componentLabel += ( minorSeparator + "{:d}".format ( ires+1 ) )
                # . Set labels.
                pathHead = "".join ( [ entityLabel, majorSeparator, componentLabel, majorSeparator ] )
                # . Loop over atoms in residue.
                start = pointers[ires  ] - 1
                stop  = pointers[ires+1] - 1
                for iatom in range ( start, stop ): atomPaths.append ( pathHead + atomnames[iatom] )
            return Sequence.FromAtomPaths ( atomPaths )

    def ToSymmetry ( self, system ):
        """Assign symmetry to a system."""
        # . Check for symmetry.
        if hasattr ( self, "box_dimensions" ):
            # . Unpack data.
            beta = self.box_dimensions[0]
            a    = self.box_dimensions[1]
            b    = self.box_dimensions[2]
            c    = self.box_dimensions[3]
            # . Get the crystal class.
            if ( beta != 90.0 ):            crystalSystem = CrystalSystemMonoclinic   ( )
            elif ( a == b ) and ( a == c ): crystalSystem = CrystalSystemCubic        ( )
            elif ( a == b ):                crystalSystem = CrystalSystemTetragonal   ( )
            else:                           crystalSystem = CrystalSystemOrthorhombic ( )
            # . Define symmetry.
            system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystalSystem )
            system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( a = a, b = b, c = c, beta = beta )

    def ToSystem ( self, mmModel = None ):
        """Create a system."""
        system = None
        if self.isParsed:
            # . Sequence data and then basic system.
            sequence     = self.ToSequence ( )
            system       = System.FromSequence ( sequence )
            system.label = "".join ( self.title )
            # . Assign atomic numbers using atomic numbers if they exist, otherwise masses.
            atomicNumbers = getattr ( self, "atomic_number", None )
            if atomicNumbers is None:
                for ( mass, atom ) in zip ( self.mass, system.atoms ):
                    atom.atomicNumber = PeriodicTable.AtomicNumberFromMass ( mass )
            else:
                for ( n, atom ) in zip ( atomicNumbers, system.atoms ):
                    atom.atomicNumber = n
            # . The MM model (with some options).
            if mmModel is None: mmModel = MMModelAMBER ( )
            mmState = MMModelState ( )
            # . 1-4 scaling factors.
            if hasattr ( self, "scee_scale_factor" ): mmModel.electrostaticScale14 = 1.0 / self.scee_scale_factor[0]
            if hasattr ( self, "scnb_scale_factor" ): mmModel.lennardJonesScale14  = 1.0 / self.scnb_scale_factor[0]
            # . The MM atom data.
            self.ToMMAtomData ( mmState )
            # . Various MM terms.
            mmTerms = mmState.mmTerms
            mm      = self.ToHarmonicBondContainer ( )
            if mm is not None: mmTerms.append ( mm )
            mm      = self.ToHarmonicAngleContainer ( )
            if mm is not None: mmTerms.append ( mm )
            self.UntangleDihedralsAnd14s ( )
            if self.dihedrals is not None: mmTerms.append ( self.dihedrals )
            if self.impropers is not None: mmTerms.append ( self.impropers )
            # . Remaining MM data.
            if self.interactions14 is not None: mmState.interactions14 = self.interactions14
            mm = self.ToExclusionsPairList ( )
            if mm is not None: mmState.exclusions = mm
            mm = self.ToLJParameterContainer ( mmModel.lennardJonesStyle )
            if mm is not None:
                mmState.ljParameters = mm
                # . 1-4 parameters.
                scale = mmModel.lennardJonesScale14
                if scale != 0.0:
                    lj14       = Clone ( mm )
                    lj14.label = "1-4 Lennard-Jones"
                    lj14.Scale ( scale )
                    mmState.ljParameters14 = lj14
            # . Assign the model.
            mmState.target                   = system
            system._energyModels["_mmModel"] = EnergyModelPriority.MMModel
            system.__dict__     ["_mmModel"] = mmModel
            system.__dict__     [ "mmState"] = mmState
            system._UpdateEnergyClosures ( )
            # . Complete the system connectivity from the MM model.
            system.mmModel.CompleteConnectivity  ( system )
            # . Symmetry data.
            self.ToSymmetry ( system )
        return system

    def UntangleDihedralsAnd14s ( self ):
        """Create the dihedral, improper and 1-4 data structures."""
        if not self.QUNTANGLE:
            # . Initialization - old data.
# . Corrected by Fernando Bachega  January 20 - 2017.
# . It is likely that similar changes need to be made for other self section attributes.
            data        = getattr ( self, "dihedrals_inc_hydrogen", [] ) + getattr ( self, "dihedrals_without_hydrogen", [] )
            nparameters = len ( self.dihedral_force_constant )
            # . Initialization - new data.
            dIndices = [] ; dParameters = [] ; dTerms = []
            iIndices = [] ; iParameters = [] ; iTerms = []
            dpIndex  = [ 0 for i in range ( nparameters ) ]
            ipIndex  = [ 0 for i in range ( nparameters ) ]
            ij14     = []
            # . Determine the dihedral, improper and 1-4 terms.
            for i in range ( len ( data[0::5] ) ):
                atom1  = data[5*i  ] // 3
                atom2  = data[5*i+1] // 3
                atom3  = data[5*i+2] // 3
                atom4  = data[5*i+3] // 3
                pIndex = data[5*i+4] - 1
                # . 1-4 interactions.
                if ( atom3 >= 0 ) and ( atom4 >= 0 ):
                    ij14.append (       atom1   )
                    ij14.append ( abs ( atom4 ) )
                # . Dihedrals and impropers.
                if ( atom4 >= 0 ):
                    dIndices.append ( ( atom1, atom2, abs ( atom3 ), abs ( atom4 ), pIndex ) )
                    dpIndex[pIndex] += 1
                else:
                    iIndices.append ( ( atom1, atom2, abs ( atom3 ), abs ( atom4 ), pIndex ) )
                    ipIndex[pIndex] += 1
            # . Determine the dihedral and improper parameters allowing for duplications.
            ndparameters = niparameters = 0
            for i in range ( nparameters ):
                fc     =                    self.dihedral_force_constant[i]
                period = int ( math.floor ( self.dihedral_periodicity[i] + 0.5 ) ) # . Transform to the nearest integer.
                phase  =                    self.dihedral_phase[i]
                if dpIndex[i] > 0:
                    ndparameters += 1
                    dpIndex[i]    = ndparameters
                    dParameters.append ( ( fc, period, phase ) )
                if ipIndex[i] > 0:
                    niparameters += 1
                    ipIndex[i]    = niparameters
                    iParameters.append ( ( fc, period, phase ) )
            # . Create the terms.
            for ( a1, a2, a3, a4, p ) in dIndices:
                dTerms.append ( ( a1, a2, a3, a4, dpIndex[p] - 1, True ) )
            for ( a1, a2, a3, a4, p ) in iIndices:
                iTerms.append ( ( a1, a2, a3, a4, ipIndex[p] - 1, True ) )
            # . Create the data structures.
            # . Dihedrals.
            if ( len ( dParameters ) > 0 ) and ( len ( dTerms ) > 0 ):
                state = { "label" : "Fourier Dihedral", "parameters" : dParameters, "terms" : dTerms }
                mm = FourierDihedralContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
                self.dihedrals = mm
            else:
                self.dihedrals = None
            # . Impropers.
            if ( len ( iParameters ) > 0 ) and ( len ( iTerms ) > 0 ):
                state = { "label" : "Fourier Improper", "parameters" : iParameters, "terms" : iTerms }
                mm = FourierDihedralContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
                self.impropers = mm
            else:
                self.impropers = None
            # . Create the 1-4 pairlist.
            if ij14 is not None:
                self.interactions14 = SelfPairList.FromIndexPairs ( ij14 )
                self.interactions14.label = "1-4 Interactions"
            else:
                self.interactions14 = None
            # . Finish up.
            self.QUNTANGLE = True

    def Verify ( self ):
        """Verify the charges and the energy units of the input data."""
        if self.isParsed and not self.QVERIFIED:
            # . Charges.
            for i in range ( len ( self.charge ) ): self.charge[i] /= _ELECTRONTOKCAL
            q = sum ( self.charge )
            if math.fabs ( q ) <= _NEUTRALITYCUTOFF:
                q /= float ( len ( self.charge ) )
                for i in range ( len ( self.charge ) ): self.charge[i] -= q
            # . Units.
            for i in range ( len ( self.bond_force_constant     ) ): self.bond_force_constant[i]     *= Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
            for i in range ( len ( self.angle_force_constant    ) ): self.angle_force_constant[i]    *= Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
            for i in range ( len ( self.dihedral_force_constant ) ): self.dihedral_force_constant[i] *= Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
            for i in range ( len ( self.lennard_jones_acoef     ) ): self.lennard_jones_acoef[i]     *= Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
            for i in range ( len ( self.lennard_jones_bcoef     ) ): self.lennard_jones_bcoef[i]     *= Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
            self.QVERIFIED = True

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { System : AmberTopologyFileReader.PathToSystem } , [ "top", "TOP" ], "Amber Topology" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

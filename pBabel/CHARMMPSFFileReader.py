"""Classes and functions for reading CHARMM PSF files."""

import itertools

from  pCore              import Clone                     , \
                                DataType                  , \
                                logFile                   , \
                                LogFileActive             , \
                                SelfPairList              , \
                                Selection                 , \
                                TextFileReader
from  pMolecule          import EnergyModelPriority       , \
                                Sequence                  , \
                                System
from  pMolecule.MMModel  import CMAPDihedralContainer     , \
                                FourierDihedralContainer  , \
                                HarmonicAngleContainer    , \
                                HarmonicBondContainer     , \
                                HarmonicImproperContainer , \
                                LJForm                    , \
                                LJParameterContainer      , \
                                MMModelCHARMM             , \
                                MMModelState
from  pScientific        import PeriodicTable
from  pScientific.Arrays import Array
from .ExportImport       import _Importer

#
# . Notes:
#
#   A problem that is known but not yet corrected is with NAMD-generated PSF files that employ LES.
#   First, all LES atoms have the same sequence. Second, the NNB section employs the number of atoms
#   appropriate for the system without multiple copies of the given species. Therefore, an additional
#   number of entries need to be added.
#
#   The extended format appears to have changed between CHARMM versions 34 and 35 with the difference
#   being the replacement of A4 for A6 in the atom lines of XPLOR files. A6 is now used.
#
#   The DRUDE option is not fully implemented.
#

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The character for parameter keys.
_ParameterKeyFieldSeparator = ":"

# . Counter names for each section (if different from section name).
_SectionCounters = { "CheqLabels" : [ "Cheq Labels" ],    \
                     "CrossTerms" : [ "Cross Terms" ],    \
                     "Groups"     : [ "Groups", "ST2s" ], \
                     "LonePairs"  : [ "Lone Pairs", "Lone Pair Hydrogens"] }

# . Mapping from CHARMM to pMolecule names.
# . There are two alternative "Groups" sections as NAMD seems to omit "NST2".
_SectionMapping = { "MOLNT"  : "CheqLabels" , "NACC" : "Acceptors"  , "NATOM" : "Atoms"     , "NBOND"  : "Bonds"  , "NCRTERM" : "CrossTerms" , "NDON"         : "Donors"    , "NGRP" : "Groups", "NGRP NST2" : "Groups", \
                    "NIMPHI" : "Impropers"  , "NNB"  : "Exclusions" , "NPHI"  : "Dihedrals" , "NTHETA" : "Angles" , "NTITLE"  : "Title"      , "NUMLP NUMLPH" : "LonePairs" }

# . PSF file sections.
_SectionNames = ( "MOLNT", "NACC", "NATOM", "NBOND", "NCRTERM", "NDON", "NGRP", "NGRP NST2", "NIMPHI", "NNB", "NPHI", "NTHETA", "NTITLE", "NUMLP NUMLPH" )

# . Wildcards.
_DihedralWildCard = "X"

#===================================================================================================================================
# . CHARMM PSF file reader class.
#===================================================================================================================================
class CHARMMPSFFileReader ( TextFileReader ):
    """CHARMMPSFFileReader is the class for CHARMM PSF files that are to be read."""

    _classLabel = "CHARMM PSF File Reader"

    def GetSectionHeader ( self ):
        """Get a section header."""
        section = None
        while True:
            line = self.GetLine ( signalWarnings = False )
            if len ( line ) > 0:
                isOK   = False
                tokens = line.split ( "!", 1 )
#                print tokens
                if len ( tokens ) > 1:
                    name = tokens[1].split ( ":", 1 )[0]
#                    print name
                    if name in _SectionNames:
                        section = _SectionMapping[name]
                        head    = tokens[0].strip ( )
#                        print head
                        try:
                            for ( i, label ) in enumerate ( _SectionCounters.get ( section, [ section ] ) ):
                                self.counters[label] = int ( head[self.integerWidth*i:self.integerWidth*(i+1)] )
#                                print label
#                                print self.counters
                            isOK = True
                        except:
                            pass
                if isOK: break
                else:    self.Warning ( "Invalid or unrecognized section header.", True )
        return section

    def Initialize ( self ):
        """Initialization for parsing."""
        # . Read data.
        self.acceptors         = []
        self.angles            = []
        self.atoms             = []
        self.bonds             = []
        self.cheqlabels        = []
        self.crossterms        = []
        self.counters          = {}
        self.dihedrals         = []
        self.donors            = []
        self.exclusionsI       = []
        self.exclusionsJ       = []
        self.groups            = []
        self.impropers         = []
        self.lonepairhydrogens = []
        self.title             = []

    def Parse ( self, isXPLOR = False, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                self.Initialize ( )
                # . Header "PSF EXT CMAP CHEQ".
                tokens = set ( self.GetTokens ( ) )
                extendedFormat = ( "EXT"   in tokens ) ; tokens.discard ( "EXT"   )
                self.hasCHEQ   = ( "CHEQ"  in tokens ) ; tokens.discard ( "CHEQ"  )
                self.hasCMAP   = ( "CMAP"  in tokens ) ; tokens.discard ( "CMAP"  )
                self.hasDRUDE  = ( "DRUDE" in tokens ) ; tokens.discard ( "DRUDE" )
                self.isXPLOR   = isXPLOR
                if not ( ( len ( tokens ) == 1 ) and ( "PSF" in tokens ) ): self.Warning ( "Invalid header line.", False )
                self.SetFormatOptions ( extendedFormat )
                # . Process each section.
                while True:
                    section = self.GetSectionHeader ( )
                    if section is not None:
                        parsemethod = getattr ( self, "Parse" + section, None )
                        if parsemethod is not None: parsemethod ( )
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseAcceptors ( self ):
        """Parse acceptors."""
        self.acceptors = self.GetFixedFormatArray ( 2 * self.counters["Acceptors"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseAngles ( self ):
        """Parse angles."""
        self.angles = self.GetFixedFormatArray ( 3 * self.counters["Angles"], 9, self.integerWidth, converter = int, default = 0 )

    def ParseAtoms ( self ):
        """Parse atoms."""
        for i in range ( self.counters["Atoms"] ):
            self.atoms.append ( self.GetFixedFormatTokens ( *(self.atomLineFormat) ) )

    def ParseBonds ( self ):
        """Parse bonds."""
        self.bonds = self.GetFixedFormatArray ( 2 * self.counters["Bonds"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseCheqLabels ( self ):
        """Parse labels for CHEQ."""
        self.cheqlabels = self.GetFixedFormatArray ( self.counters["Atoms"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseCrossTerms ( self ):
        """Parse CMAP cross-terms."""
        self.crossterms = self.GetFixedFormatArray ( 8 * self.counters["Cross Terms"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseDihedrals ( self ):
        """Parse dihedrals."""
        self.dihedrals = self.GetFixedFormatArray ( 4 * self.counters["Dihedrals"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseDonors ( self ):
        """Parse donors."""
        self.donors = self.GetFixedFormatArray ( 2 * self.counters["Donors"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseExclusions ( self ):
        """Parse exclusions."""
        if self.counters["Exclusions"] > 0:
            self.exclusionsJ = self.GetFixedFormatArray ( self.counters["Exclusions"], 8, self.integerWidth, converter = int, default = 0 )
        else:
            self.GetLine ( )
        self.exclusionsI = self.GetFixedFormatArray ( self.counters["Atoms"]     , 8, self.integerWidth, converter = int, default = 0 )

    def ParseGroups ( self ):
        """Parse groups."""
        self.groups = self.GetFixedFormatArray ( 3 * self.counters["Groups"], 9, self.integerWidth, converter = int, default = 0 )

    def ParseImpropers ( self ):
        """Parse impropers."""
        self.impropers = self.GetFixedFormatArray ( 4 * self.counters["Impropers"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseLonePairs ( self ):
        """Parse lone pairs."""
        if self.counters["Lone Pairs"] > 0:
            for i in range ( self.counters["Lone Pairs"] ): self.GetLine ( ) # . Nothing done here.
            self.lonepairhydrogens = self.GetFixedFormatArray ( self.counters["Lone Pair Hydrogens"], 8, self.integerWidth, converter = int, default = 0 )

    def ParseTitle ( self ):
        """Parse the title."""
        for i in range ( self.counters["Title"] ):
            line = self.GetLine ( )
            if line.startswith ( "*" ): line = line[1:].strip ( )
            if len ( line ) > 0: self.title.append ( line )

    @classmethod
    def PathToSystem ( selfClass, path, isXPLOR = False, log = logFile, mmModel = None, parameters = None ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( isXPLOR = isXPLOR, log = log )
        inFile.Summary ( log = log )
        return inFile.ToSystem ( mmModel = mmModel, parameters = parameters )

    def SetFormatOptions ( self, extendedFormat ):
        """Set format options.

        CHARMM: I,SEGID(ISEG),RESID(IRES),RES(IRES),TYPE(I),    IAC(I) ,CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
        XPLOR : I,SEGID(ISEG),RESID(IRES),RES(IRES),TYPE(I),ATC(IAC(I)),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)

        Extended format:

        CHARMM - CHEQ: '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)'
        CHARMM + CHEQ: '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)'
        XPLOR  - CHEQ: '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8)'
        XPLOR  + CHEQ: '(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A6,1X,2G14.6,I8,2G14.6)'

        Normal format:

        CHARMM - CHEQ: '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8)'
        CHARMM + CHEQ: '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)'
        XPLOR  - CHEQ: '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8)'
        XPLOR  + CHEQ: '(I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)'
        """
        # . Item widths.
        iacWidth = 4
        if extendedFormat:
            self.integerWidth = 10
            stringWidth       =  8
            if self.isXPLOR: iacWidth = 6
        else:
            self.integerWidth =  8
            stringWidth       =  4
        # . Atom line format.
        p = 0
        i = self.integerWidth
        s = stringWidth
        formats =      [ ( p, p+i, int , 0  ) ] ; p += (i+1)
        formats.append ( ( p, p+s, None, "" ) ) ; p += (s+1)
        formats.append ( ( p, p+s, None, "" ) ) ; p += (s+1)
        formats.append ( ( p, p+s, None, "" ) ) ; p += (s+1)
        formats.append ( ( p, p+s, None, "" ) ) ; p += (s+1)
        if self.isXPLOR: formats.append ( ( p, p+iacWidth, None, "" ) )
        else:            formats.append ( ( p, p+iacWidth, int ,  0 ) )
        p += ( iacWidth + 1 )
        formats.append ( ( p, p+14, float, 0.0 ) ) ; p += 14
        formats.append ( ( p, p+14, float, 0.0 ) ) ; p += 14
        formats.append ( ( p, p+ 8, int  , 0   ) ) ; p +=  8
        # . Additional tokens for CHEQ and DRUDE.
        if   self.hasDRUDE: n = 4
        elif self.hasCHEQ : n = 2
        else: n = 0
        if n > 0:
            for i in range ( n ):
                formats.append ( ( p, p+14, float, 0.0 ) )
                p += 14
        self.atomLineFormat = formats

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( key, "{:d}".format ( self.counters[key] ) ) for key in self.counters.keys ( ) ]

    def ToCMAPDihedralContainer ( self, rawParameters = {} ):
        """Create a CMAP dihedral container."""
        mm     = None
        nterms = self.counters.get ( "Cross Terms", 0 )
        if nterms > 0:
            # . Terms.
            atomIndices = []
            uniquekeys  = set ( )
            termkeys    = []
            for n in range ( nterms ):
                i1  = self.crossterms[8*n  ] - 1
                j1  = self.crossterms[8*n+1] - 1
                k1  = self.crossterms[8*n+2] - 1
                l1  = self.crossterms[8*n+3] - 1
                i2  = self.crossterms[8*n+4] - 1
                j2  = self.crossterms[8*n+5] - 1
                k2  = self.crossterms[8*n+6] - 1
                l2  = self.crossterms[8*n+7] - 1
                ti1 = self.atomtypes[i1]
                tj1 = self.atomtypes[j1]
                tk1 = self.atomtypes[k1]
                tl1 = self.atomtypes[l1]
                ti2 = self.atomtypes[i2]
                tj2 = self.atomtypes[j2]
                tk2 = self.atomtypes[k2]
                tl2 = self.atomtypes[l2]
                # . Specific key and then wildcard key.
                if   tj1 >  tk1: key1 = [ ti1, tj1, tk1, tl1 ]
                elif tj1 == tk1: key1 = [ max ( ti1, tl1 ), tj1, tk1, min ( ti1, tl1 ) ]
                else:            key1 = [ tl1, tk1, tj1, ti1 ]
                if   tj2 >  tk2: key2 = [ ti2, tj2, tk2, tl2 ]
                elif tj2 == tk2: key2 = [ max ( ti2, tl2 ), tj2, tk2, min ( ti2, tl2 ) ]
                else:            key2 = [ tl2, tk2, tj2, ti2 ]
                key = tuple ( key1 + key2 )
                uniquekeys.add ( key )
                atomIndices.append ( [ i1, j1, k1, l1, i2, j2, k2, l2 ] )
                termkeys.append ( key )
            # . Parameters.
            # . Care is taken to unravel multiple dihedrals.
            uniquekeys    = list ( uniquekeys )
            uniquekeys.sort ( )
            keyindex      = {}
            parameterKeys = []
            parameters    = []
            n             =  0
            for key in uniquekeys:
                ( abscissa, values ) = rawParameters[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                parameters.append    ( ( abscissa, abscissa, values ) )
                keyindex[key] = n
                n += 1
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyindex[termkeys[i]], True ] ) )
            # . Construct the container.
            state = { "label"         : "CMAP Dihedral" ,
                      "parameterKeys" : parameterKeys   ,
                      "parameters"    : parameters      ,
                      "terms"         : terms           }
            mm = CMAPDihedralContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToExclusionPairLists ( self ):
        """Create pairlists for the exclusions and 1-4 interactions.

        This is done from the bond list as in CHARMM.
        """
        # . Initialization.
        exclusions     = None
        interactions14 = None
        natoms         = self.counters.get ( "Atoms", 0 )
        nbonds         = self.counters.get ( "Bonds", 0 )
        # . 1-2 lists.
        i12 = set ( )
        for n in range ( nbonds ):
            i   = self.bonds[2*n  ] - 1
            j   = self.bonds[2*n+1] - 1
            t12 = ( max ( i, j ), min ( i, j ) )
            i12.add ( t12 )
        # . Generate sorted connections for each atom.
        connections = [ [] for n in range ( natoms ) ]
        for ( i, j ) in i12:
            connections[i].append ( j )
            connections[j].append ( i )
        for n in range ( natoms ): connections[n].sort ( )
        # . 1-3 lists (i > k automatically).
        i13 = set ( )
        for j in range ( natoms ):
            jbonds = connections[j]
            nj     = len ( jbonds )
            for i in range ( 1, nj ):
                for k in range ( 0, i ): i13.add ( ( jbonds[i], jbonds[k] ) )
        # . Generate 123 exclusions.
        i123 = i12.union  ( i13 )
        # . Add in explicitly defined exclusions.
        if self.counters.get ( "Exclusions", 0 ) > 0:
            nnb = 0
            for i in range ( natoms ):
                n = self.exclusionsI[i]
                for j in self.exclusionsJ[nnb:n]: i123.add ( ( max ( i, j ), min ( i, j ) ) )
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
        i1234 = i123.union ( i14p )
        i14   = i14p.difference ( i123 )
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

    def ToFreeAtoms ( self, system ):
        """Define the free atoms for the system."""
        natoms = self.counters.get ( "Atoms", 0 )
        if natoms > 0:
            # . Get the fixed atom indices.
            indices = []
            for ( index, datum ) in enumerate ( self.atoms ):
                flag = datum[8]
                if flag == 0: indices.append ( index )
            # . Define the fixed atoms.
            if len ( indices ) < natoms: system.freeAtoms = Selection.FromIterable ( indices )

    def ToFourierDihedralContainer ( self, rawParameters = {}, rawParametersWild = {} ):
        """Create a harmonic dihedral container."""
        mm     = None
        nterms = self.counters.get ( "Dihedrals", 0 )
        if nterms > 0:
            # . Terms.
            atomIndices = []
            uniquekeys  = set ( )
            termkeys    = []
            for n in range ( nterms ):
                i  = self.dihedrals[4*n  ] - 1
                j  = self.dihedrals[4*n+1] - 1
                k  = self.dihedrals[4*n+2] - 1
                l  = self.dihedrals[4*n+3] - 1
                ti = self.atomtypes[i]
                tj = self.atomtypes[j]
                tk = self.atomtypes[k]
                tl = self.atomtypes[l]
                # . Specific key and then wildcard key.
                if   tj >  tk: key = ( ti, tj, tk, tl )
                elif tj == tk: key = ( max ( ti, tl ), tj, tk, min ( ti, tl ) )
                else:          key = ( tl, tk, tj, ti )
                if key in rawParameters:
                    nperiods = len ( rawParameters[key] )
                else:
                    key = ( _DihedralWildCard, max ( tj, tk ), min ( tj, tk ), _DihedralWildCard )
                    nperiods = len ( rawParametersWild[key] )
                uniquekeys.add ( key )
                for p in range ( nperiods ):
                    atomIndices.append ( [ i, j, k, l ] )
                    termkeys.append ( ( key, p ) )
            # . Parameters.
            # . Care is taken to unravel multiple dihedrals.
            uniquekeys = list ( uniquekeys )
            uniquekeys.sort ( )
            parameterKeys = []
            parameters    = []
            keyindex      = {}
            n             =  0
            for key in uniquekeys:
                data = rawParameters.get ( key, None )
                if data is None: data = rawParametersWild[key]
                baseKey = _ParameterKeyFieldSeparator.join ( key )
                for ( nperiod, datum ) in enumerate ( data ):
                    keyindex[( key, nperiod )] = n
                    parameterKeys.append ( baseKey + _ParameterKeyFieldSeparator + "{:d}".format ( datum[0] ) )
                    parameters.append ( ( datum[1], datum[0], datum[2] ) )
                    n += 1
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyindex[termkeys[i]], True ] ) )
            # . Construct the container.
            state = { "label"         : "Fourier Dihedral" ,
                      "parameterKeys" : parameterKeys      ,
                      "parameters"    : parameters         ,
                      "terms"         : terms              }
            mm = FourierDihedralContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicAngleContainer ( self, rawParameters = {} ):
        """Create a harmonic angle container."""
        mm     = None
        nterms = self.counters.get ( "Angles", 0 )
        if nterms > 0:
            # . Terms.
            atomIndices = []
            termkeys    = []
            for n in range ( nterms ):
                i  = self.angles[3*n  ] - 1
                j  = self.angles[3*n+1] - 1
                k  = self.angles[3*n+2] - 1
                ti = self.atomtypes[i]
                tj = self.atomtypes[j]
                tk = self.atomtypes[k]
                key = ( max ( ti, tk ), tj, min ( ti, tk ) )
                atomIndices.append ( [ i, j, k ] )
                termkeys.append ( key )
            # . Parameters.
            uniquekeys = list ( set ( termkeys ) )
            uniquekeys.sort ( )
            parameterKeys = []
            parameters    = []
            keyindex      = {}
            for ( i, key ) in enumerate ( uniquekeys ):
                keyindex[key] = i
                datum = rawParameters[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                parameters.append ( ( datum[1], datum[0] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyindex[termkeys[i]], True ] ) )
            # . Construct the container.
            state = { "parameterKeys" : parameterKeys ,
                      "parameters"    : parameters    ,
                      "terms"         : terms         }
            mm = HarmonicAngleContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicBondContainer ( self, rawParameters = {} ):
        """Create a harmonic bond container."""
        mm     = None
        nterms = self.counters.get ( "Bonds", 0 )
        if nterms > 0:
            # . Terms.
            atomIndices = []
            termkeys    = []
            for n in range ( nterms ):
                i  = self.bonds[2*n  ] - 1
                j  = self.bonds[2*n+1] - 1
                ti = self.atomtypes[i]
                tj = self.atomtypes[j]
                key = ( max ( ti, tj ), min ( ti, tj ) )
                atomIndices.append ( [ i, j ] )
                termkeys.append ( key )
            # . Parameters.
            uniquekeys = list ( set ( termkeys ) )
            uniquekeys.sort ( )
            parameterKeys = []
            parameters    = []
            keyindex      = {}
            for ( i, key ) in enumerate ( uniquekeys ):
                keyindex[key] = i
                datum = rawParameters[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                parameters.append ( ( datum[1], datum[0] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyindex[termkeys[i]], True ] ) )
            # . Construct the container.
            state = { "parameterKeys" : parameterKeys ,
                      "parameters"    : parameters    ,
                      "terms"         : terms         }
            mm = HarmonicBondContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicImproperContainer ( self, rawParameters = {}, rawParametersWild = {} ):
        """Create a harmonic improper container."""
        mm     = None
        nterms = self.counters.get ( "Impropers", 0 )
        if nterms > 0:
            # . Terms.
            atomIndices = []
            termkeys    = []
            for n in range ( nterms ):
                i  = self.impropers[4*n  ] - 1
                j  = self.impropers[4*n+1] - 1
                k  = self.impropers[4*n+2] - 1
                l  = self.impropers[4*n+3] - 1
                ti = self.atomtypes[i]
                tj = self.atomtypes[j]
                tk = self.atomtypes[k]
                tl = self.atomtypes[l]
                # . Specific key and then wildcard key.
                if   ti >  tl: key = ( ti, tj, tk, tl )
                elif ti == tl: key = ( ti, max ( tj, tk ), min ( tj, tk ), tl )
                else:          key = ( tl, tk, tj, ti )
                if key not in rawParameters: key = ( max ( ti, tl ), _DihedralWildCard, _DihedralWildCard, min ( ti, tl ) )
                atomIndices.append ( [ i, j, k, l ] )
                termkeys.append ( key )
            # . Parameters.
            uniquekeys = list ( set ( termkeys ) )
            uniquekeys.sort ( )
            parameterKeys = []
            parameters    = []
            keyindex      = {}
            for ( i, key ) in enumerate ( uniquekeys ):
                keyindex[key] = i
                datum = rawParameters.get ( key, None )
                if datum is None: datum = rawParametersWild[key]
                parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                parameters.append ( ( datum[1], datum[0] ) )
            # . Terms.
            terms = []
            for ( i, indices ) in enumerate ( atomIndices ):
                terms.append ( tuple ( indices + [ keyindex[termkeys[i]], True ] ) )
            # . Construct the container.
            state = { "parameterKeys" : parameterKeys ,
                      "parameters"    : parameters    ,
                      "terms"         : terms         }
            mm = HarmonicImproperContainer.Raw ( )
            mm.__setstate__ ( state )
            mm.Sort ( )
        return mm

    def ToHarmonicUreyBradleyContainer ( self, rawParameters = {} ):
        """Create a harmonic Urey-Bradley container.

        These are only assigned if the parameters are specified.
        """
        mm     = None
        nterms = self.counters.get ( "Angles", 0 )
        if ( nterms > 0 ) and ( len ( rawParameters ) > 0 ):
            # . Terms.
            atomIndices = []
            termkeys    = []
            for n in range ( nterms ):
                i  = self.angles[3*n  ] - 1
                j  = self.angles[3*n+1] - 1
                k  = self.angles[3*n+2] - 1
                ti = self.atomtypes[i]
                tj = self.atomtypes[j]
                tk = self.atomtypes[k]
                key = ( max ( ti, tk ), tj, min ( ti, tk ) )
                if key in rawParameters:
                    atomIndices.append ( [ i, k ] )
                    termkeys.append ( key )
            nterms = len ( termkeys )
            if nterms > 0:
                # . Parameters.
                uniquekeys = list ( set ( termkeys ) )
                uniquekeys.sort ( )
                parameterKeys = []
                parameters    = []
                keyindex      = {}
                for ( i, key ) in enumerate ( uniquekeys ):
                    keyindex[key] = i
                    datum = rawParameters[key]
                    parameterKeys.append ( _ParameterKeyFieldSeparator.join ( key ) )
                    parameters.append ( ( datum[1], datum[0] ) )
                # . Terms.
                terms = []
                for ( i, indices ) in enumerate ( atomIndices ):
                    terms.append ( tuple ( indices + [ keyindex[termkeys[i]], True ] ) )
                # . Construct the container.
                state = { "is12Interaction" : False          ,
                          "label"           : "Urey-Bradley" ,
                          "parameterKeys"   : parameterKeys  ,
                          "parameters"      : parameters     ,
                          "terms"           : terms          }
                mm = HarmonicBondContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
        return mm

    def ToLJParameterContainers ( self, parameters = {}, parametersWild = {}, parameters14 = {}, parameters14Wild = {} ):
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
            for atomtype in self.uniqueAtomTypes:
                ( e, r ) = parameters[atomtype]
                data14   = parameters14.get ( atomtype, None )
                epsilons.append ( e )
                sigmas.append   ( r )
                if data14 is None:
                    epsilons14.append ( e )
                    sigmas14.append   ( r )
                else:
                    epsilons14.append ( data14[0] )
                    sigmas14.append   ( data14[1] )
            ljParameters   = LJParameterContainer.FromEpsilonSigma ( ntypes, epsilons  , sigmas  , analyticForm = LJForm.Amber, parameterKeys = self.uniqueAtomTypes )
            ljParameters14 = LJParameterContainer.FromEpsilonSigma ( ntypes, epsilons14, sigmas14, analyticForm = LJForm.Amber, parameterKeys = self.uniqueAtomTypes )
            ljParameters14.label = "1-4 LJ Parameters"
        return ( ljParameters, ljParameters14 )

    def ToMMAtomData ( self, mmState ):
        """Add atom data to a MM state."""
        n = self.counters.get ( "Atoms", 0 )
        if n > 0:
            # . Gather some data.
            atomTypes   = []
            atomCharges = []
            for datum in self.atoms:
                atomTypes.append   ( datum[5].upper ( ) )
                atomCharges.append ( datum[6] )
            # . Get unique types.
            uniqueTypes     = set ( atomTypes )
            uniqueAtomTypes = list ( uniqueTypes )
            uniqueAtomTypes.sort ( )
            # . Generate an index.
            uniqueAtomTypeIndex = {}
            for ( i, atomType ) in enumerate ( uniqueAtomTypes ): uniqueAtomTypeIndex[atomType] = i
            # . Save the atomtypes.
            self.atomtypes       = atomTypes
            self.uniqueAtomTypes = uniqueAtomTypes
            # . Create the data.
            atomTypeIndices = Array.WithExtent ( n , dataType = DataType.Integer )
            charges         = Array.WithExtent ( n )
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
        # . Initialization.
        majorSeparator = Sequence._attributable["labelSeparator"]
        minorSeparator = Sequence._attributable["fieldSeparator"]
        # . Create atom paths - segID (datum[1]), resname (datum[3]), resID (datum[2]), atom name (datum[4]).
        atomPaths = []
        for datum in self.atoms:
            items = []
            for i in ( 4, 3, 2, 1 ):
                item = datum[i]
                if len ( item ) <= 0: item = " "
                items.append ( item )
            atomPaths.append ( items[3] + majorSeparator + items[1] + minorSeparator + items[2] + majorSeparator + items[0] )
        # . Get the sequence.
        return Sequence.FromAtomPaths ( atomPaths )

    def ToSystem ( self, mmModel = None, parameters = None ):
        """Create a system."""
        system = None
        if self.isParsed:
            # . Get the sequence.
            sequence = self.ToSequence ( )
            # . Basic system.
            system = System.FromSequence ( sequence )
            system.label = " ".join ( self.title )
            # . Assign atomic numbers from masses.
            for ( datum, atom ) in zip ( self.atoms, system.atoms ):
                atom.atomicNumber = PeriodicTable.AtomicNumberFromMass ( datum[7] )
            # . The MM model.
            if mmModel is None: mmModel = MMModelCHARMM ( )
            mmState = MMModelState ( )
            # . Build the model but only when PSF files are in XPLOR format.
            if ( mmModel is not None ) and ( mmState is not None ) and ( parameters is not None ) and ( self.isXPLOR ):
                # . The MM atom data.
                self.ToMMAtomData ( mmState )
                # . Various MM terms.
                mmTerms = mmState.mmTerms
                mm      = self.ToHarmonicBondContainer        ( rawParameters = parameters.bonds  )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToHarmonicAngleContainer       ( rawParameters = parameters.angles )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToHarmonicUreyBradleyContainer ( rawParameters = parameters.ureybradleys )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToFourierDihedralContainer     ( rawParameters = parameters.dihedrals, rawParametersWild = parameters.dihedralwilds )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToHarmonicImproperContainer    ( rawParameters = parameters.impropers, rawParametersWild = parameters.improperwilds )
                if mm is not None: mmTerms.append ( mm )
                mm      = self.ToCMAPDihedralContainer        ( rawParameters = parameters.cmaps )
                if mm is not None: mmTerms.append ( mm )
                # . Exclusions.
                ( exclusions, interactions14 ) = self.ToExclusionPairLists ( )
                if exclusions     is not None: mmState.exclusions     = exclusions
                if interactions14 is not None: mmState.interactions14 = interactions14
                # . LJ terms.
                ( ljParameters, ljParameters14 ) = self.ToLJParameterContainers ( parameters       = parameters.nonbonds       ,
                                                                                  parametersWild   = parameters.nonbondwilds   ,
                                                                                  parameters14     = parameters.nonbond14s     ,
                                                                                  parameters14Wild = parameters.nonbond14wilds )
                if ljParameters   is not None: mmState.ljParameters   = ljParameters
                if ljParameters14 is not None: mmState.ljParameters14 = ljParameters14
                # . Assign the model.
                mmState.target                   = system
                system._energyModels["_mmModel"] = EnergyModelPriority.MMModel
                system.__dict__     ["_mmModel"] = mmModel
                system.__dict__     [ "mmState"] = mmState
                system._UpdateEnergyClosures ( )
            # . Free atoms - do after MM model so MM terms are appropriately affected.
            self.ToFreeAtoms ( system )
        return system

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { System : CHARMMPSFFileReader.PathToSystem } , [ "psf", "psfx", "PSF", "PSFX" ], "Charmm PSF" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

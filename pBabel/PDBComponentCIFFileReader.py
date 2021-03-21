"""Classes and functions for reading PDB chemical component files in mmCIF format."""

from  pCore        import logFile, LogFileActive, TextFileReader
from  pMolecule    import BondType
from  pScientific  import PeriodicTable
from .PDBComponent import PDBComponent, PDBComponentAtom, PDBComponentBond

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Bond strings.
_BondStrings = { "AROM" : ( BondType.Single , True  ) ,
                 "DOUB" : ( BondType.Double , False ) ,
                 "SING" : ( BondType.Single , False ) ,
                 "TRIP" : ( BondType.Triple , False ) }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentCIFFileReader ( TextFileReader ):
    """PDBComponentCIFFileReader is the class for PDB chemical component files in mmCIF format that are to be read."""

    _classLabel = "PDB Component CIF File Reader"

    # . Some simplicity is assumed in the parsing as all data entries are of the form:
    # . data_XXX followed by "#".
    # . general info. followed by "#".
    # . atom info. followed by "#".
    # . bond info. followed by "#".
    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            self.active     = None
            self.components = []
            if LogFileActive ( log ): self.log = log
#            self.tags = set ( )
            # . Open the file.
            self.Open ( )
            try:
                QLOOP    = False
                signalWarnings = False
                while True:
                    # . Get a line.
                    line = self.GetLine ( signalWarnings = signalWarnings )
                    # . Terminator.
                    if line == "#":
                        QLOOP    = False
                        signalWarnings = False
                    # . Loop.
                    elif line == "loop_":
                        QLOOP    = True
                        signalWarnings = True
                    # . Header.
                    elif line.startswith ( "data_" ):
                        if self.active is not None: self.components.append ( self.active )
                        self.active = PDBComponent.WithOptions ( label = line[5:] )
                        QLOOP    = False
                        signalWarnings = True
                    # . Atom block.
                    elif line.startswith ( "_chem_comp_atom." ):
                        self.ParseAtomBlock ( line, QLOOP )
                        QLOOP    = False
                        signalWarnings = False
                    # . Bond block.
                    elif line.startswith ( "_chem_comp_bond." ):
                        self.ParseBondBlock ( line, QLOOP )
                        QLOOP    = False
                        signalWarnings = False
                    # . General block.
                    elif line.startswith ( "_chem_comp." ):
                        self.ParseGeneralBlock ( line )
                        QLOOP    = False
                        signalWarnings = False
                    # . Descriptor and identifier blocks.
                    elif line.startswith ( "_pdbx_chem_comp_descriptor." ) or \
                         line.startswith ( "_pdbx_chem_comp_identifier." ):
                        self.ParseDescriptorBlock ( line )
                        QLOOP    = False
                        signalWarnings = False
                    # . Other blocks.
                    else:
                        QLOOP    = False
                        signalWarnings = True
#                        print line
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Finish up.
            if self.active is not None: self.components.append ( self.active )
            del self.active
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True
#            tags = list ( self.tags )
#            tags.sort ( )
#            for tag in tags:
#                print tag

    def ParseAtomBlock ( self, line, QLOOP ):
        """Parse the block with atom information."""
        # . Initialization.
        atoms  = []
        QFIRST = True
        if QLOOP:
            # . Initialization.
            index      = {}
            converters = []
            QATOMS     = False
            # . Atom lines.
            while True:
                if not QFIRST: line = self.GetLine ( )
                if line.startswith ( "#" ): break
                # . Column headers.
                elif not QATOMS:
                    if   line.startswith ( "_chem_comp_atom.atom_id"       ):
                        index["atom_id"]     = len ( converters )
                        converters.append ( None )
                    elif line.startswith ( "_chem_comp_atom.charge"        ):
                        index["charge"]      = len ( converters )
                        converters.append ( int )
                    elif line.startswith ( "_chem_comp_atom.comp_id"       ):
                        index["comp_id"]     = len ( converters )
                        converters.append ( None )
                    elif line.startswith ( "_chem_comp_atom.pdbx_align"    ):
                        index["pdbx_align"]  = len ( converters )
                        converters.append ( int )
                    elif line.startswith ( "_chem_comp_atom.type_symbol"   ):
                        index["type_symbol"] = len ( converters )
                        converters.append ( PeriodicTable.AtomicNumber )
                    elif line.startswith ( "_chem_comp_atom." ):
                        converters.append ( None )
                    else:
                        QATOMS = True
                # . Atoms.
                if QATOMS:
                    tokens = self.ParseTokens ( line, converters = converters )
                    atoms.append ( PDBComponentAtom.WithOptions ( atomicNumber = tokens[index["type_symbol"]],  \
                                                                  formalCharge = tokens[index["charge"     ]],  \
                                                                  label        = tokens[index["atom_id"    ]],  \
                                                                  pdbAlign     = tokens[index["pdbx_align" ]] ) )
                QFIRST = False
        # . Single atom.
        else:
            while True:
                if not QFIRST: line = self.GetLine ( )
                tokens = self.ParseTokens ( line )
                if tokens[0] == "#": break
                elif len ( tokens ) > 2: self.Warning ( "Too many tokens.", True )
                elif tokens[0] == "_chem_comp_atom.atom_id"    : label        =       tokens[1]
                elif tokens[0] == "_chem_comp_atom.charge"     : formalCharge = int ( tokens[1] )
                elif tokens[0] == "_chem_comp_atom.pdbx_align" : pdbAlign     = int ( tokens[1] )
                elif tokens[0] == "_chem_comp_atom.type_symbol": atomicNumber = PeriodicTable.AtomicNumber ( tokens[1] )
                QFIRST = False
            atoms.append ( PDBComponentAtom.WithOptions ( atomicNumber = atomicNumber,  \
                                                          formalCharge = formalCharge,  \
                                                          label        = label,         \
                                                          pdbAlign     = pdbAlign     ) )
        # . Finish up.
        self.active.atoms = atoms

    def ParseBondBlock ( self, line, QLOOP ):
        """Parse the block with bond information."""
        # . Initialization.
        bonds  = []
        QFIRST = True
        # . Multiple bonds.
        if QLOOP:
            # . Initialization.
            index  = {}
            n      = 0
            QBONDS = False
            # . Bond lines.
            while True:
                if not QFIRST: line = self.GetLine ( )
                if line.startswith ( "#" ): break
                # . Column headers.
                elif not QBONDS:
                    if   line.startswith ( "_chem_comp_bond.atom_id_1"   ):
                        index["atom_id_1"]   = n
                        n += 1
                    elif line.startswith ( "_chem_comp_bond.atom_id_2"   ):
                        index["atom_id_2"]   = n
                        n += 1
                    elif line.startswith ( "_chem_comp_bond.comp_id"     ):
                        index["comp_id"]     = n
                        n += 1
                    elif line.startswith ( "_chem_comp_bond.value_order" ):
                        index["value_order"] = n
                        n += 1
                    elif line.startswith ( "_chem_comp_bond." ):
                        n += 1
                    else:
                        QBONDS = True
                # . Bonds.
                if QBONDS:
                    tokens                       = self.ParseTokens ( line )
                    ( bondType, bondIsAromatic ) = _BondStrings.get ( tokens[index["value_order"]], ( BondType.Undefined, False ) )
                    bonds.append ( PDBComponentBond.WithOptions ( atomLabel1 = tokens[index["atom_id_1"]],  \
                                                                  atomLabel2 = tokens[index["atom_id_2"]],  \
                                                                  bondType   = bondType       ,
                                                                  isAromatic = bondIsAromatic ) )
                    # . Checking.
                    # . Atoms in bond exist and type recognized.
                QFIRST = False
        # . Single bond.
        else:
            while True:
                if not QFIRST: line = self.GetLine ( )
                tokens = self.ParseTokens ( line )
                if tokens[0] == "#": break
                elif len ( tokens ) > 2: self.Warning ( "Too many tokens.", True )
                elif tokens[0] == "_chem_comp_bond.atom_id_1"  : atomLabel1 = tokens[1]
                elif tokens[0] == "_chem_comp_bond.atom_id_2"  : atomLabel2 = tokens[1]
                elif tokens[0] == "_chem_comp_bond.value_order": bondType   = _BondStrings.get ( tokens[1], BondType.Undefined )
                QFIRST = False
            bonds.append ( PDBComponentBond.WithOptions ( atomLabel1 = atomLabel1, atomLabel2 = atomLabel2, bondType = bondType ) )
        # . Finish up.
        self.active.bonds = bonds

    def ParseDescriptorBlock ( self, line ):
        """Parse the block with descriptor information."""
        # . Initialization.
        QFIRST = True
        # . Processing.
        while True:
            if not QFIRST: line = self.GetLine ( )
            if line == "#": break
            elif line.startswith ( "_pdbx_chem_comp_descriptor." ): pass
#                tokens = line.split ( )
#                self.tags.add ( tokens[0][27:] )
            elif line.startswith ( "_pdbx_chem_comp_identifier." ): pass
#                tokens = line.split ( )
#                self.tags.add ( tokens[0][27:] )
#            else: self.Warning ( "Unrecognized descriptor line type: " + line + ".", False )
# . Other tags are:
#comp_id
#descriptor
#identifier
#program
#program_version
#type
            QFIRST = False
        # . Finish up.

    def ParseGeneralBlock ( self, line ):
        """Parse the block with general information."""
        # . Initialization.
        # . Component data.
        componentClass = "Unknown"
        formalCharge   = None
        isHeteroatom    = True
        pdbClass       = "UNKNOWN"
        # . Local.
        namelines    = []
        synonymlines = []
        QFIRST       = True
        QNAME        = False
        QSYNONYM     = False
        # . Processing.
        while True:
            if not QFIRST: line = self.GetLine ( )
            if line == "#": break
            elif line.startswith ( "_chem_comp.id" ):
                tokens = line.split ( )
                if tokens[1] != self.active.label: self.Warning ( "Entry ID mismatch: " + self.active.label + "/" + tokens[1] + ".", False )
            elif line.startswith ( "_chem_comp.name" ):
                QNAME    = True
                QSYNONYM = False
                line     = line[15:].strip ( )
                if len ( line ) > 0: namelines.append ( line )
            elif line.startswith ( "_chem_comp.mon_nstd_flag" ):
                tokens = line.split ( )
                if len ( tokens ) > 1:
                    if   tokens[1].upper ( ) == "N": isHeteroatom = True
                    elif tokens[1].upper ( ) == "Y": isHeteroatom = False
            elif line.startswith ( "_chem_comp.pdbx_formal_charge" ):
                tokens = line.split ( )
                charge = tokens[1]
                sign   = 1
                if   charge[-1] == "+":
                    charge = charge[:-1]
                    sign   = 1
                elif charge[-1] == "-":
                    charge = charge[:-1]
                    sign   = -1
                else:
                    sign   = 1
                try:    formalCharge = int ( charge ) * sign
                except: formalCharge = None
            elif line.startswith ( "_chem_comp.pdbx_synonyms" ):
                QNAME    = False
                QSYNONYM = True
                line     = line[24:].strip ( )
                if len ( line ) > 0: synonymlines.append ( line )
            elif line.startswith ( "_chem_comp.pdbx_type" ):
                pdbClass = line.split ( )[1].title ( )
            elif line.startswith ( "_chem_comp.type" ):
                componentClass = line.split ( )[1].title ( )
            elif line.startswith ( "_chem_comp." ):
                QNAME    = False
                QSYNONYM = False
            elif QNAME:
                if line[0] == ";": line = line[1:].strip ( )
                if len ( line ) > 0: namelines.append ( line )
            elif QSYNONYM:
                if line[0] == ";": line = line[1:].strip ( )
                if len ( line ) > 0: synonymlines.append ( line )
            else: self.Warning ( "Unrecognized general line type: " + line + ".", False )
            QFIRST = False
        # . Finish up.
        if len ( namelines ) > 0: name = " ".join ( namelines )
        else:                     name = None
        # . Save data.
        self.active.componentClass = componentClass
        self.active.formalCharge   = formalCharge
        self.active.name           = name
        self.active.isHeteroatom   = isHeteroatom
        self.active.pdbClass       = pdbClass

    def ParseTokens ( self, line, converters = None ):
        """Parse tokens on a line.

        This method assumes that all tokens containing spaces or single or double quotes
        will be enclosed in single or double quotes. The enclosing quotes are removed."""
        tokens = None
        if line is not None:
            # . Do the splitting on whitespace as normal and then treat single and double quotes.
            # . This assumes that no more than one space occurs. Isolated quotes are not treated.
            # . Initialization - constants.
            doubleQuote  = '"'
            singleQuote  = "'"
            # . Initialization - variables.
            openQuotes   = False
            openToken    = None
            partialToken = None
            quoteToken   = None
            tokens       = []
            # . Loop over tokens.
            for token in line.split ( ):
                # . Open quotes.
                if openQuotes:
                    # . Closing.
                    if token[-1] == quoteToken:
                        partialToken.append ( token[:-1] )
                        tokens.append ( " ".join ( partialToken ) )
                        openQuotes   = False
                        partialToken = None
                        quoteToken   = None
                    # . Inside.
                    else:
                        partialToken.append ( token )
                # . Opening quotes.
                elif token[0] in ( doubleQuote, singleQuote ):
                    openQuotes = True
                    quoteToken = token[0]
                    # . Closing as well.
                    if token[-1] == quoteToken:
                        tokens.append ( token[1:-1] )
                        openQuotes = False
                        quoteToken = None
                    # . Inside.
                    else:
                        partialToken = [ token[1:] ]
                # . Outside quotes.
                else:
                    tokens.append ( token )
            if partialToken is not None:
#                print partial
#                print line
#                print line.split ( )
                self.Warning ( "Unpaired single or double quotes on line.", True )
            # . Conversion.
            if converters is not None:
                for ( i, ( token, converter ) ) in enumerate ( zip ( tokens, converters ) ):
                    if converter is None:
                        new = token
                    else:
                        try:
                            new = converter ( token )
                        except:
                            new = converter ( )
                            self.Warning ( "Unable to convert token {:d}.".format ( i ), True )
                    tokens[i] = new
        return tokens

    @classmethod
    def PathToComponents ( selfClass, cifPath, asDictionary = True, fullOutput = False, log = logFile ):
        """Read components from a PDB chemical component file in mmCIF format."""
        # . Parse file.
        inFile = selfClass.FromPath ( cifPath )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        components = inFile.ToComponents ( )
        # . Further output if necessary.
        if LogFileActive ( log ):
            # . Basic output.
            log.Paragraph ( "PDB component dictionary processed with {:d} components.".format ( len ( components ) ) )
            # . Full output.
            if fullOutput:
                table = log.GetTable ( columns = [ 10, 10, 10, 10 ] )
                table.Start ( )
                table.Title ( "PDB Components" )
                table.Heading ( "Index" )
                table.Heading ( "Label" )
                table.Heading ( "Atoms" )
                table.Heading ( "Bonds" )
                for ( i, component ) in enumerate ( components ):
                    table.Entry ( "{:d}".format ( i ) )
                    table.Entry ( component.label )
                    table.Entry ( "{:d}".format ( len ( component.atoms ) ) )
                    table.Entry ( "{:d}".format ( len ( component.bonds ) ) )
                table.Stop ( )
        # . Return component dictionary.
        if asDictionary:
            index = {}
            for component in components: index[component.label] = component
            return index
        # . Return component list.
        else:
            return components

    def ToComponents ( self ):
        """Return the components."""
        if self.isParsed: return self.components
        else:             return []

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

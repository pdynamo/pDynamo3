"""Classes and functions for writing CHARMM PSF files in XPLOR format."""

from  pCore             import Selection                 , \
                               TextFileWriter
from  pMolecule         import Sequence                  , \
                               System
from  pMolecule.MMModel import CMAPDihedralContainer     , \
                               FourierDihedralContainer  , \
                               HarmonicAngleContainer    , \
                               HarmonicBondContainer     , \
                               HarmonicImproperContainer
from .ExportImport      import _Exporter

#===================================================================================================================================
# . CHARMM PSF file writer class.
#===================================================================================================================================
class CHARMMPSFFileWriter ( TextFileWriter ):
    """CHARMMPSFFileWriter is the class for CHARMM PSF files that are to be written."""

    _classLabel = "CHARMM PSF File Writer"

    # . Mapping from pDynamo to CHARMM section header names.
    _SectionMapping = { "Acceptors"  : "NACC: acceptors"     ,
                        "Angles"     : "NTHETA: angles"      ,
                        "Atoms"      : "NATOM"               ,
                        "Bonds"      : "NBOND: bonds"        ,
                        "CheqLabels" : "MOLNT"               ,
                        "CrossTerms" : "NCRTERM: cross-terms",
                        "Dihedrals"  : "NPHI: dihedrals"     ,
                        "Donors"     : "NDON: donors"        ,
                        "Exclusions" : "NNB"                 ,
                        "Groups"     : "NGRP NST2"           ,
                        "Impropers"  : "NIMPHI: impropers"   ,
                        "LonePairs"  : "NUMLP NUMLPH"        ,
                        "Title"      : "NTITLE"              }

    # . Sections to write (order significant).
    _SectionNames = ( "Header"     ,
                      "Title"      ,
                      "Atoms"      ,
                      "Bonds"      ,
                      "Angles"     ,
                      "Dihedrals"  ,
                      "Impropers"  ,
                      "Donors"     ,
                      "Acceptors"  ,
                      "Exclusions" ,
                      "Groups"     ,
                      "CheqLabels" ,
                      "LonePairs"  ,
                      "CrossTerms" )

    # . Temporary attributes.
    _TemporaryAttributes = ( "extendedFormat",
                             "hasCHEQ"       ,
                             "hasCMAP"       ,
                             "isXPLOR"       ,
                             "sectionNames"  ,
                             "system"        ,
                             "systemData"    )

    # . Format-specific options.
    _Use2Counters = set ( [ "Groups", "LonePairs" ] )
    _Use9Columns  = set ( [ "Angles", "Groups"    ] )

    def Clear ( self ):
        """Clear temporary attributes."""
        for attribute in self.__class__._TemporaryAttributes:
            setattr ( self, attribute, None )

    def GetAtomSequenceData ( self, atom ):
        """Get atom sequence data given an atom from a sequence."""
        atomName            = atom.label
        sequence            = self.systemData["Sequence"]
        ( resName, resSeq ) = sequence.ParseLabel ( atom.parent.label       , fields = 2 )
        entityFields        = sequence.ParseLabel ( atom.parent.parent.label, fields = 2 )
        segName             = entityFields[1]
        if len ( segName ) <= 0: segName = entityFields[0]
        if not self.extendedFormat:
            atomName = atomName[0:4]
            resName  = resName [0:4]
            resSeq   = resSeq  [0:4]
            segName  = segName [0:4]
        return ( atomName, resName, resSeq, segName )

    def MakeGroupData ( self, sequence, charges, fixedAtoms ):
        """Make group data for output."""
        # . Assumes contiguous atoms within each component.
        accumulatedAtoms = 0
        numberOfGroups   = sequence.numberOfComponents
        indices          = []
        for entity in sequence.children:
            for component in entity.children:
                isFixed       = True
                localCharge   = 0.0
                startingIndex = numberOfGroups
                for atom in component.children:
                    i = atom.index
                    isFixed       = ( isFixed and ( i in fixedAtoms ) )
                    localCharge  += charges[i]
                    startingIndex = min ( i, startingIndex )
                igpbs = accumulatedAtoms + startingIndex - 1
                if localCharge == 0.0: igptyp = 0
                else:                  igptyp = 1
                if isFixed: imoveg =  0
                else:       imoveg = -1
                indices.extend ( [ igpbs, igptyp, imoveg ] )
                accumulatedAtoms += len ( component.children )
        return [ numberOfGroups, 0, indices ]

    def OutputSystem ( self, system, extendedFormat = True, hasCHEQ = True, isXPLOR = True ):
        """Write a system to the file."""
        # . System data.
        self.PrepareSystemData ( system, hasCHEQ )
        # . Options.
        self.extendedFormat = extendedFormat
        self.hasCHEQ        = hasCHEQ
        self.hasCMAP        = ( "CrossTerms" in self.systemData )
        self.isXPLOR        = isXPLOR
        self.SetupOutputOptions ( )
        # . Writing (with empty line after each written section).
        for sectionName in self.sectionNames:
            writeMethod = getattr ( self, "Write" + sectionName, None )
            if writeMethod is None: self.WriteGeneric ( sectionName )
            else:                   writeMethod ( )
            self.file.write ( "\n" )
        # . Finish up.
        self.Clear ( )

    @classmethod
    def PathFromSystem ( selfClass, path, system, extendedFormat = True, hasCHEQ = True, isXPLOR = True ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.Open  ( )
        outFile.OutputSystem ( system, extendedFormat = extendedFormat, hasCHEQ = hasCHEQ, isXPLOR = isXPLOR )
        outFile.Close ( )

    def PrepareSystemData ( self, system, hasCHEQ ):
        """Prepare data for output given a system."""
        # . Initialization.
        numberOfAtoms = len ( system.atoms )
        self.system   = system
        # . Extract basic system data.
        # . Fixed atoms.
        if system.freeAtoms is None: fixedAtoms = set ( )
        else:                        fixedAtoms = system.freeAtoms.Complement ( upperBound = len ( system.atoms ) )
        # . Atom data.
        try:
            mmState         = system.mmState
            mmTerms         = mmState.mmTerms
            atomTypeIndices = mmState.atomTypeIndices
            charges         = mmState.charges
            uniqueAtomTypes = mmState.atomTypes
            atomTypes       = [ uniqueAtomTypes[index] for index in atomTypeIndices ]
            masses          = [ atom.mass for atom in system.atoms ]
        except Exception as e:
            print ( e[0] )
            raise ValueError ( "Unable to obtain MM atom data for the system." )
        # . Sequence.
        if system.sequence is None: sequence = Sequence.FromAtoms ( system.atoms, componentLabel = "UNK.1", entityLabel = "UNK" )
        else                      : sequence = system.sequence
        # . Setup output terms.
        # . Atom data.
        systemData = { "AtomTypes" : atomTypes, "Charges" : charges, "FixedAtoms" : fixedAtoms, "Masses" : masses, "Sequence" : sequence }
        # . CHEQ - use sequence components rather than isolates.
        if hasCHEQ:
            index   = 0
            indices = []
            for entity in sequence.children:
                for component in entity.children:
                    indices.extend ( len ( component.children ) * [ index ] )
                    index += 1
            systemData["CheqLabels"] = [ numberOfAtoms, indices ]
        # . Exclusions.
        systemData["Exclusions"] = [ 0, None, numberOfAtoms * [ -1 ] ]
        # . Groups.
        systemData["Groups"] = self.MakeGroupData ( sequence, charges, fixedAtoms )
        # . MM terms.
        if mmTerms is not None:
            for mmTerm in mmTerms:
                tag = None
                if   isinstance ( mmTerm, CMAPDihedralContainer     ): tag = "CrossTerms"
                elif isinstance ( mmTerm, FourierDihedralContainer  ):
                    if ( "Improper"      in mmTerm.label ) or \
                       ( "Out-Of-Plane"  in mmTerm.label ): tag = "Impropers"
                    else:                                   tag = "Dihedrals"
                elif isinstance ( mmTerm, HarmonicAngleContainer    ): tag = "Angles"
                elif isinstance ( mmTerm, HarmonicBondContainer     ):
                    if ( "Harmonic Bond" in mmTerm.label ): tag = "Bonds"
                elif isinstance ( mmTerm, HarmonicImproperContainer ): tag = "Impropers"
                if tag is not None:
                    ( n, indices ) = mmTerm.GetUniqueTermIndices ( )
                    if n > 0: systemData[tag] = [ n, indices ]
        # . Finish up.
        self.systemData = systemData

    def SetupOutputOptions ( self ):
        """Setup the output formats and sections."""
        # . XPLOR format is the only one that makes sense at the moment.
        if not self.isXPLOR: raise ValueError ( "Unable to write CHARMM PSF files that are not in XPLOR format." )
        # . Item widths.
        if self.extendedFormat:
            self.integerFormat        = "{:10d}"
            self.sectionHeaderFormat1 = "{:10d} !{:s}\n"
            self.sectionHeaderFormat2 = "{:10d}{:10d} !{:s}\n"
            if self.isXPLOR: format = "{:10d} {:<8s} {:<8s} {:<8s} {:<8s} {:<6s} {:14.6E}{:14.6E}{:8d}" # . Note change of 6s from 4s (newer versions of CHARMM).
            else           : format = "{:10d} {:<8s} {:<8s} {:<8s} {:<8s} {:4d} {:14.6E}{:14.6E}{:8d}"
        else:
            self.integerFormat        = "{:8d}"
            self.sectionHeaderFormat1 = "{:8d} !{:s}\n"
            self.sectionHeaderFormat2 = "{:8d}{:8d} !{:s}\n"
            if self.isXPLOR: format = "{:8d} {:<4s} {:<4s} {:<4s} {:<4s} {:<4s} {:14.6E}{:14.6E}{:8d}"
            else           : format = "{:8d} {:<4s} {:<4s} {:<4s} {:<4s} {:4d} {:14.6E}{:14.6E}{:8d}"
        if self.hasCHEQ: format += "{:14.6E}{:14.6E}"
        self.atomLineFormat = format + "\n"
        # . Output sections.
        self.sectionNames = list ( self.__class__._SectionNames )
        if not self.hasCHEQ: self.sectionNames.remove ( "CheqLabels" )
        if not self.hasCMAP: self.sectionNames.remove ( "CrossTerms" )

    def WriteAtoms ( self ):
        """Write atoms."""
        # . Header.
        numberOfAtoms = len ( self.system.atoms )
        self.file.write ( self.sectionHeaderFormat1.format ( numberOfAtoms, self.__class__._SectionMapping["Atoms"] ) )
        if numberOfAtoms > 0:
            # . Prepare atom data.
            atomTypes  = self.systemData["AtomTypes" ]
            charges    = self.systemData["Charges"   ]
            fixedAtoms = self.systemData["FixedAtoms"]
            masses     = self.systemData["Masses"    ]
            sequence   = self.systemData["Sequence"  ]
            # . Output.
            for i in range ( numberOfAtoms ):
                atomType = atomTypes[i]
                charge   = charges  [i]
                mass     = masses   [i]
                if i in fixedAtoms: isFixed = 1
                else:               isFixed = 0
                ( atomName, resName, resSeq, segName ) = self.GetAtomSequenceData ( self.system.atoms[i] )
                if self.hasCHEQ:
                    self.file.write ( self.atomLineFormat.format ( i+1, segName, resSeq, resName, atomName, atomType, charge, mass, isFixed, 0.0, 0.0 ) )
                else:
                    self.file.write ( self.atomLineFormat.format ( i+1, segName, resSeq, resName, atomName, atomType, charge, mass, isFixed ) )

    def WriteGeneric ( self, sectionName ):
        """Write generic data."""
        # . Get all data.
        data = self.systemData.get ( sectionName, None )
        tag  = self.__class__._SectionMapping[sectionName]
        if sectionName in self.__class__._Use2Counters: counters = 2
        else:                                           counters = 1
        if sectionName in self.__class__._Use9Columns : columns  = 9
        else:                                           columns  = 8
        # . Default data.
        if data is None: data = counters * [ 0 ]
        # . Header.
        if counters == 1: self.file.write ( self.sectionHeaderFormat1.format ( data.pop ( 0 ),                 tag ) )
        else:             self.file.write ( self.sectionHeaderFormat2.format ( data.pop ( 0 ), data.pop ( 0 ), tag ) )
        # . Indices.
        while len ( data ) > 0:
            indices = data.pop ( 0 )
            if indices is None: self.file.write ( "\n" )
            else: self.WriteIntegerIndices ( indices, columns )

    def WriteHeader ( self ):
        """Write the header."""
        # . Same order as CHARMM.
        tokens = [ "PSF" ]
        if self.extendedFormat: tokens.append ( "EXT"  )
        if self.hasCMAP       : tokens.append ( "CMAP" )
        if self.hasCHEQ       : tokens.append ( "CHEQ" )
        self.file.write ( " ".join ( tokens ) + "\n" )

    def WriteIntegerIndices ( self, indices, columns ):
        """Write an array of integer indices to the file with increment."""
        format = self.integerFormat
        n      = 0
        for index in indices:
            self.file.write ( format.format ( index + 1 ) ) # . With increment.
            n += 1
            if n >= columns:
                self.file.write ( "\n" )
                n = 0
        if n != 0:
            self.file.write ( "\n" )

    def WriteTitle ( self ):
        """Write the title."""
        lines = []
        if self.system.label is not None: lines.append ( self.system.label )
        lines.append ( "PSF file created by pDynamo" )
        self.file.write ( self.sectionHeaderFormat1.format ( len ( lines ), self.__class__._SectionMapping["Title"] ) )
        for line in lines:
            self.file.write ( "* " + line + "\n" )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : CHARMMPSFFileWriter.PathFromSystem } , [ "psf", "psfx", "PSF", "PSFX" ], "Charmm PSF" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

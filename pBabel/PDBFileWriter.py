"""Write data to a PDB file."""

import itertools, time

from  pCore                 import TextFileWriter, __version__
from  pMolecule             import Sequence, System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Exporter

# . The occupancies and temperatureFactors arguments can be used to write PQR files if these contain the
# . charges and radii, respectively. Unfortunately, one can get precision problems as PDB uses 6.2f formatting
# . whereas PQR uses space delimiters.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Line formats.
_ATOMLINEFORMAT1    = "{:<6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4s}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:2s}{:2s}\n"
_ATOMLINEFORMAT2    = "{:<6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:5s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:<4s}{:2s}{:2s}\n"
_AUTHORLINEFORMAT   = "{:<6s}    {:s} {:s}\n"
_COMPOUNDLINEFORMAT = "{:<6s}    {:s}\n"
_ENDLINEFORMAT      = "{:<3s}\n"
_HEADERLINEFORMAT   = "{:<6s}    {:<40s}{:<9s}\n"
_MASTERLINEFORMAT   = "{:<6s}    {:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}{:5d}\n"
_TERLINEFORMAT1     = "{:<6s}{:5d}      {:3s} {:1s}{:4s}{:1s}\n"
_TERLINEFORMAT2     = "{:<6s}{:5d}      {:3s} {:1s}{:5s}\n"
_TITLELINEFORMAT    = "{:<6s}    {:s}\n"

# . Connection line data.
_CONECTLINEFORMAT1  = "{:6s}{:5d}"
_CONECTLINEFORMAT2  = "{:5d}"
_NCONECTS           = 4

# . Other data.
_MaximumTitleLength = 70

# . Standard residues.
_STANDARDRESIDUES = set ( [ "ALA", "ARG", "ASN", "ASP", "ASX", "CYS", "GLN", "GLU", "GLX", "GLY", \
                            "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", \
                            "TYR", "UNK", "VAL", "A"  , "+A" , "C"  , "+C" , "G"  , "+G" , "I"  , \
                            "+I", "T", "+T", "U", "+U", "UNK" ] )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBFileWriter ( TextFileWriter ):
    """PDBFileWriter is the class for writing PDB files."""

    _classLabel = "PDB File Writer"

    @classmethod
    def PathFromSystem ( selfClass                      ,
                         path                           ,
                         system                         ,
                         occupancies            = None  ,
                         QCONECT                = False ,
                         selection              = None  ,
                         suppressICodeField     = False ,
                         temperatureFactors     = None  ,
                         useSegmentEntityLabels = False ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSystem ( system                                          ,
                              occupancies            = occupancies            ,
                              QCONECT                = QCONECT                ,
                              selection              = selection              ,
                              suppressICodeField     = suppressICodeField     ,
                              temperatureFactors     = temperatureFactors     ,
                              useSegmentEntityLabels = useSegmentEntityLabels )

    # . Public methods.
    def WriteSystem ( self                           ,
                      system                         ,
                      data                   = None  ,
                      occupancies            = None  ,
                      QCONECT                = False ,
                      selection              = None  ,
                      suppressICodeField     = False ,
                      temperatureFactors     = None  ,
                      title                  = None  ,
                      useSegmentEntityLabels = False ):
        """Write out a system.

        |system|    is the system to be written.
        |data|      is the data to be written. The coordinates of system are written if |data| is absent.
        |QCONECT|   is the option for writing CONECT records. If |QCONECT| is True all connections are written (not just those to heteroatoms).
        |selection| gives the atoms to write. The default is to write all atoms.
        |title|     is an optional title.
        """
        # . Check |system|.
        if not isinstance ( system, System ): raise TypeError ( "Invalid |system| argument." )

        # . Get the sequence.
        sequence = getattr ( system, "sequence", None )
        if sequence is None: sequence = Sequence.FromAtoms ( system.atoms, componentLabel = "UNK.1" )

        # . Get the coordinate data.
        if isinstance ( data, Coordinates3 ): xyz = data
        else:                                 xyz = system.coordinates3
        if xyz is None: raise TypeError ( "Unable to obtain coordinate data from |system| or |data| arguments." )

        # . Check that the PDB model and coordinates are consistent.
        if len ( system.atoms ) != xyz.rows: raise TypeError ( "The PDB model and coordinate data are of different lengths." )

        # . Get the selection.
        if selection is None: towrite = range ( len ( system.atoms ) )
        else:                 towrite = selection

        # . Get polymer terminating atoms as the last atom of the terminating component.
        polymerTerminatingAtoms = set ( )
        for polymer in sequence.linearPolymers:
            polymerTerminatingAtoms.add ( polymer.rightTerminalComponent.children[-1] )

        # . Get the label.
        label = None
        if   isinstance ( title,        str ): label = title.upper        ( )
        elif isinstance ( system.label, str ): label = system.label.upper ( )
        if ( label is not None ) and ( len ( label ) > _MaximumTitleLength ): label = label[0:_MaximumTitleLength]

        # . Check for extra data.
        hasOccupancies        = ( occupancies        is not None ) and ( len ( occupancies        ) == len ( system.atoms ) )
        hasTemperatureFactors = ( temperatureFactors is not None ) and ( len ( temperatureFactors ) == len ( system.atoms ) )

        # . Start writing.
        self.Open ( )
        # . Header.
        self.file.write ( _HEADERLINEFORMAT.format ( "HEADER", "UNKNOWN", time.strftime ( "%d-%b-%y" ).upper ( ) ) )
        # . Title.
        if ( label is not None ): self.file.write ( _TITLELINEFORMAT.format ( "TITLE", label ) )
        # . Author.
        self.file.write ( _AUTHORLINEFORMAT.format ( "AUTHOR", "GENERATED BY PDYNAMO", __version__ ) )
        # . Compound.
        self.file.write ( _COMPOUNDLINEFORMAT.format ( "COMPND", "UNKNOWN" ) )
        # . Initialization.
        natm     = 0
        ncon     = 0
        nter     = 0
        cIndices = {}
        # . Loop over atom indices.
        for iatom in towrite:
            # . Get data.
            atom = system.atoms [iatom]
            # . Atom data.
            element = PeriodicTable.Symbol ( atom.atomicNumber ).upper ( )
            if element == "*": element = ""
            # . Component data.
            ( resName, resSeq, iCode ) = sequence.ParseLabel ( atom.parent.label, fields = 3 )
            # . Entity data.
            entityLabel = atom.parent.parent.label
            if useSegmentEntityLabels:
                chainID = ""
                segID   = entityLabel[0:4]
            else:
                chainID = entityLabel[0:1]
                segID   = ""
            # . Atom output.
            if resName in _STANDARDRESIDUES: recordname = "ATOM"
            else:                            recordname = "HETATM"
            # . Output label.
            label = atom.label
            if   len ( label ) >= 4    : outputlabel = label[0:4]
            elif label[0:1].isdigit ( ): outputlabel = label
            else:                        outputlabel = " " + label
            # . Extra data.
            if hasOccupancies:        occupancy         = occupancies[iatom]
            else:                     occupancy         = 0.0
            if hasTemperatureFactors: temperatureFactor = temperatureFactors[iatom]
            else:                     temperatureFactor = 0.0
            # . Output the line.
            if suppressICodeField:
                self.file.write ( _ATOMLINEFORMAT2.format ( recordname, natm + nter + 1, outputlabel, " ", resName[0:3], chainID, resSeq[:5],
                                                            xyz[iatom,0], xyz[iatom,1], xyz[iatom,2], occupancy, temperatureFactor, segID, element, " " ) )
            else:
                self.file.write ( _ATOMLINEFORMAT1.format ( recordname, natm + nter + 1, outputlabel, " ", resName[0:3], chainID, resSeq[:4], iCode,
                                                            xyz[iatom,0], xyz[iatom,1], xyz[iatom,2], occupancy, temperatureFactor, segID, element, " " ) )
            # . Save the atom index.
            cIndices[iatom] = ( natm + nter + 1 )
            # . Polymer termination output.
            if atom in polymerTerminatingAtoms:
                nter += 1
                if suppressICodeField: self.file.write ( _TERLINEFORMAT2.format ( "TER", natm + nter + 1, resName, chainID, resSeq        ) )
                else:                  self.file.write ( _TERLINEFORMAT1.format ( "TER", natm + nter + 1, resName, chainID, resSeq, iCode ) )
            # . Increment natm.
            natm += 1
        # . Connection data.
        if QCONECT and ( system.connectivity.bonds is not None ) and ( len ( system.connectivity.bonds ) > 0 ):
            system.connectivity.bonds.MakeConnections ( )
            for iatom in towrite:
                # . Get the connections.
                iconnections = system.connectivity.bonds.GetConnectedAtoms ( iatom ) # . Normal indices.
                # . Transform to pdb indices.
                iindex   = cIndices[iatom]
                jindices = []
                for j in iconnections:
                    jindex = cIndices.get ( j, -1 )
                    if jindex > 0: jindices.append ( jindex )
                # . Output.
                nindices = len ( jindices )
                if nindices > 0:
                    for i in range ( 0, nindices, _NCONECTS ):
                        self.file.write ( _CONECTLINEFORMAT1.format ( "CONECT", iindex ) )
                        for j in range ( i, min ( i + _NCONECTS, nindices ) ): self.file.write ( _CONECTLINEFORMAT2.format ( jindices[j] ) )
                        self.file.write ( "\n" )
                        ncon += 1
        # . Master.
        self.file.write ( _MASTERLINEFORMAT.format ( "MASTER", 0, 0, 0, 0, 0, 0, 0, 0, len ( system.atoms ), nter, ncon, 0 ) )
        # . End.
        self.file.write ( _ENDLINEFORMAT.format ( "END" ) )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : PDBFileWriter.PathFromSystem } , [ "ent", "ENT", "pdb", "PDB" ], "Protein Data Bank" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

"""A class for writing PQR files understandable by MEAD."""

from pBabel                import _Exporter
from pCore                 import TextFileWriter
from pMolecule             import Sequence       , \
                                  System
from pScientific.Geometry3 import Coordinates3

_ATOMLINEFORMAT = "{:<6s}{:5d} {:<4s}{:1s}{:3s} {:5d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:10.5f}{:8.3f}{:7d}{:7d}{:7d}{:7d}\n"

#===================================================================================================================================
# . PQR file writer class.
#===================================================================================================================================
class PQRFileWriter (TextFileWriter):
    """PQRFileWriter is a class for writing PQR files."""

    _classLabel = "PQR File Writer"

    @classmethod
    def PathFromSystem ( selfClass, path, system, charges = None, radii = None, selection = None ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSystem ( system, charges = charges, radii = radii, selection = selection )

    def WriteSystem (self, system, data=None, selection=None, charges=None, radii=None):
        """Write out a system.
        |system|    is the system to be written.
        |data|      is the coordinate data of the system.
        |charges|   is a sequence containing atomic charges. If not present, the charges are set to zero.
        |radii|     is a sequence containing atomic radii. If not present, the radii are set to zero.
        |selection| defines the atoms to write. The default is to write all atoms.
        """

        # Modified code from PDBFileWriter
        if not isinstance (system, System):
            raise TypeError ("Invalid |system| argument.")

        sequence = getattr (system, "sequence", None)
        if sequence is None:
            sequence = Sequence.FromAtomContainer (system.atoms, componentLabel = "UNK.1")

        if isinstance (data, Coordinates3):
            xyz = data
        else:
            xyz = system.coordinates3

        if xyz is None:
            raise TypeError ("Unable to obtain coordinate data from |system| or |data| arguments.")

        natoms = len (system.atoms)
        if natoms != xyz.rows:
            raise TypeError ("The PQR model and coordinate data are of different lengths.")

        if selection is None:
            towrite = xrange (natoms)
        else:
            towrite = selection

        if charges is None:
            charges = [0.] * natoms
        else:
            if natoms != len (charges):
                raise TypeError ("The PQR model and charge data are of different lengths.")

        if radii is None:
            radii = [0.] * natoms
        else:
            if natoms != len (radii):
                raise TypeError ("The PQR model and radii data are of different lengths.")

        # Shifting is necessary for multi-segment proteins because segments are not recognized by MEAD
        shifts   = {}
        segments = system.sequence.children
        for segmentIndex, segment in enumerate (segments):
            shifts[segment.label] = segmentIndex * 1000


        ParsePath = system.sequence.ParsePath
        self.Open ()

        for natm, iatom in enumerate (towrite):
            atom   = system.atoms [iatom]
            charge = charges      [iatom]
            radius = radii        [iatom]

            # FIXME
            segName, residue, atomName = ParsePath (atom.path)
            resName, resSerial = residue.split (".")

            resSerial = int (resSerial) + shifts[segName]
            atomIndex = atom.index

            label = atom.label
            if   len (label) >= 4:
                outputlabel = label[0:4]
            elif label[0:1].isdigit ():
                outputlabel = label
            else:
                outputlabel = " " + label
            x = xyz[iatom, 0]
            y = xyz[iatom, 1]
            z = xyz[iatom, 2]

            self.file.write ( _ATOMLINEFORMAT.format ( "ATOM", atomIndex, outputlabel, " ", resName[0:3], resSerial, "", x, y, z, charge, radius, 0, 0, 0, 0 ) )
        self.Close ()

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : PQRFileWriter.PathFromSystem } , [ "pqr", "PQR" ], "PQR" )

#===================================================================================================================================
# Testing
#===================================================================================================================================
if __name__ == "__main__": pass

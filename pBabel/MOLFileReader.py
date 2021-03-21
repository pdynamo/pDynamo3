#===================================================================================================================================
# . Classes and functions to read MOL and SDF files.
#===================================================================================================================================

from  pCore                 import logFile                  , \
                                   LogFileActive            , \
                                   TextFileReader
from  pMolecule             import Atom                     , \
                                   Bond                     , \
                                   BondType                 , \
                                   Connectivity             , \
                                   ConvertInputConnectivity , \
                                   System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Bond type definitions.
_MOLBondTypes = { 1 : ( BondType.Single   , False ) ,
                  2 : ( BondType.Double   , False ) ,
                  3 : ( BondType.Triple   , False ) ,
                  4 : ( BondType.Undefined, True  ) }

# . Charge definitions.
_ChargeCodes = { 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3 }

#===================================================================================================================================
# . MOL file reader class.
#===================================================================================================================================
class MOLFileReader ( TextFileReader ):
    """MOLFileReader is the class for MOL or SDF files that are to be read."""

    _classLabel = "MOL File Reader"

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            self.molRecords = []
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Parse all entries.
                notFinished = True
                while notFinished:
                    # . Header block.
                    # . Molecule name.
                    label = self.GetLine ( signalWarnings = False )
                    # . Skip lines 1 and 2 which are not parsed.
                    self.GetLine ( )
                    self.GetLine ( )
                    # . Connection table.
                    # . Counts line - only atom and bond numbers are parsed.
                    ( nAtoms, nBonds ) = self.GetFixedFormatTokens ( ( 0, 3, int, 0 ), ( 3, 6, int, 0 ) )
                    # . Initialize atom and coordinates data structures.
                    charges       = []
                    connectivity  = Connectivity ( )
                    coordinates3  = Coordinates3.WithExtent ( nAtoms )
                    mchg          = {}
                    miso          = {}
                    mrad          = {}
                    coordinates3.Set ( 0.0 )
                    # . Read the atom lines.
                    for n in range ( nAtoms ):
                        ( x, y, z, atomicNumber, charge ) = self.GetFixedFormatTokens ( ( 0, 10, float, 0.0 ), ( 10, 20, float, 0.0 ), ( 20, 30, float, 0.0 ), ( 31, 34, PeriodicTable.AtomicNumber, -1 ), ( 36, 39, int, 0 ) )
                        connectivity.AddNode ( Atom.WithOptions ( atomicNumber = atomicNumber ) )
                        charges.append ( _ChargeCodes.get ( charge, 0 ) )
                        coordinates3[n,0] = x
                        coordinates3[n,1] = y
                        coordinates3[n,2] = z
                    # . Read the bond lines.
                    for n in range ( nBonds ):
                        ( atom1, atom2, code ) = self.GetFixedFormatTokens ( ( 0, 3, int, 0 ), ( 3, 6, int, 0 ), ( 6, 9, int, 0 ) )
                        if ( atom1 <= 0 ) or ( atom1 > nAtoms ) or ( atom2 <= 0 ) or ( atom2 > nAtoms ):
                            self.Warning ( "Bond atom indices out of range: {:d}, {:d}.".format ( atom1, atom2 ), True )
                        ( bondType, bondIsAromatic ) = _MOLBondTypes.get ( code, ( BondType.Undefined, False ) )
                        connectivity.AddEdge ( Bond.WithNodes ( connectivity.atoms[atom1-1], connectivity.atoms[atom2-1], isAromatic = bondIsAromatic, type = bondType ) )
                    # . Properties lines.
                    while True:
                        line = self.GetLine ( )
                        if   line.startswith ( "M  CHG" ): self.ParsePropertiesLine ( line, mchg )
                        elif line.startswith ( "M  ISO" ): self.ParsePropertiesLine ( line, miso )
                        elif line.startswith ( "M  RAD" ): self.ParsePropertiesLine ( line, mrad )
                        elif line.startswith ( "M  END" ): break
                    # . Store the data.
                    self.molRecords.append ( ( label, connectivity, charges, coordinates3, mchg, miso, mrad ) )
                    # . Check for an SDF terminator.
                    notFinished = False
                    while True:
                        line = self.GetLine ( signalWarnings = False )
                        if line == "$$$$":
                            notFinished = True
                            break
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParsePropertiesLine ( self, line, mdict ):
        """Parse a CHG, ISO or RAD properties line."""
        try:    n = int ( line[8:9] )
        except: n = 0
        for i in range ( n ):
            o = i * 8 + 10
            try:
                a = int ( line[o  :o+3] ) - 1
                v = int ( line[o+4:o+7] )
                if a >= 0: mdict[a] = v
            except:
                pass

    @classmethod
    def PathToCoordinates3 ( selfClass, path, index = 0, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( index = index )

    @classmethod
    def PathToSystem ( selfClass, path, index = 0, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToSystem ( index = index )

    @classmethod
    def PathToSystems ( selfClass, path, log = logFile ):
        """Return all the systems from an SDF file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return [ inFile.ToSystem ( index = i ) for i in range ( len ( inFile.molRecords ) ) ]

    def ToAtoms ( self, connectivity, charges, mchg, miso, mrad ):
        """Set some extra atom attributes."""
        # . Only mchg for the moment.
        if self.isParsed:
            if len ( mchg ) > 0:
                for ( i, v ) in mchg.items ( ):
                    connectivity.atoms[i].formalCharge = v

    def ToCoordinates3 ( self, index = 0 ):
        """Return a coordinates3 object."""
        if self.isParsed:
            if index in range ( len ( self.molRecords ) ):
                data = self.molRecords[index]
                return data[3]
            else:
                raise IndexError ( "MOL record index {:d} not in range [0,{:d}].".format ( index, len ( self.molRecords ) - 1 ) )
        else:
            return None

    def ToSystem ( self, index = 0 ):
        """Return a system."""
        if self.isParsed:
            if index in range ( len ( self.molRecords ) ):
                # . Get data.
                ( label, connectivity, charges, coordinates3, mchg, miso, mrad ) = self.molRecords[index]
                # . Add extra atom attributes and finish connectivity.
                self.ToAtoms ( connectivity, charges, mchg, miso, mrad )
                ConvertInputConnectivity ( connectivity, {} )
                # . Make system.
                system              = System.FromConnectivity ( connectivity )
                system.label        = label
                system.coordinates3 = coordinates3
                return system
            else:
                raise IndexError ( "MOL record index {:d} not in range [0,{:d}].".format ( index, len ( self.molRecords ) - 1 ) )
        else:
            return None

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : MOLFileReader.PathToCoordinates3 ,
                         System       : MOLFileReader.PathToSystem       } ,
                       [ "mol", "MOL", "sdf", "SDF" ], "MDL MOL", defaultFunction = MOLFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

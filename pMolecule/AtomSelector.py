"""Atom selector."""

import sqlite3

from  pCore         import logFile            , \
                           LogFileActive      , \
                           Selection          , \
                           SummarizableObject
from .AtomSelection import AtomSelectionError
from .System        import System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomSelector ( SummarizableObject ):
    """An atom selector object."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Atom Selector"
    _attributable.update ( { "patterns"   : dict ,
                             "selections" : dict ,
                             "system"     : None } )
    _patterns = {}

    def __getattr__ ( self, name ): # or __getattribute__
        if   name in self.__dict__   : return self.__dict__  [name]
        elif name in self.selections : return self.selections[name]
        elif name in self.patterns   :
            selection = self.Where ( self.patterns[name] )
            self.selections[name] = selection
            return selection
        else: raise AttributeError ( "Unknown attribute: " + name + "." )

    def _CheckOptions ( self ):
        """Check options."""
        super ( AtomSelector, self )._CheckOptions ( )
        if not isinstance ( self.system, System ): raise AtomSelectionError ( "Invalid system passed to atom selector." )
        for ( key, value ) in self.__class__._patterns.items ( ): self.patterns[key] = value

    def AddPattern ( self, label, patternString ):
        """Add a predefined pattern."""
        self.patterns[label] = patternString
        self.selections.pop ( label, None )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Atoms"    , "{:d}".format ( len ( self.system.atoms ) ) ) ,
                 ( "Patterns" , "{:d}".format ( len ( self.patterns     ) ) ) ]

    def UpdateCoordinates3 ( self, coordinates3 = None ):
        """Update the coordinates in the selector."""
        # . Required if table other than in system.
        pass

    def Where ( self, selectionString ):
        """Return a selection."""
        pass

    @classmethod
    def WithSystem ( selfClass, system, **options ):
        """Constructor with system."""
        options = dict ( options )
        options["system"] = system
        return selfClass.WithOptions ( **options )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Deleting rows possible.
# . Merging done by row.

class SQLAtomSelector ( AtomSelector ):
    """An atom selector using SQL syntax."""

    _attributable = dict ( AtomSelector._attributable )
    _classLabel   = "SQL Atom Selector"
    _patterns     = dict ( AtomSelector._patterns     )
    _protein      = "ResNam IN ( 'ALA', 'ARG', 'ASN', 'ASP', 'CSE', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TPO', 'TRP', 'TYR', 'VAL' )"
    _attributable.update ( { "atomsColumns"    : None       ,
                             "connection"      : None       ,
                             "createFullTable" : True       ,
                             "path"            : ":memory:" } )
    _patterns.update ( { "aromatics"          : "IsAromatic = 1"                                   ,
                         "backbone"           : _protein + " and Label IN ( 'C', 'CA', 'N', 'O' )" ,
                         "boundaryAtoms"      : "IsBoundaryAtom = 1"                               ,
                         "counterions"        : "ResNam IN ( 'CL', 'K', 'NA' )"                    ,
                         "heavyAtoms"         : "( AtomicNumber > 1 ) and ( AtomicNumber < 110 )"  ,
                         "hydrogens"          : "AtomicNumber = 1"                                 ,
                         "linearPolymerAtoms" : "InLinearPolymer = 1"                              ,
                         "mmAtoms"            : "IsMMAtom = 1"                                     ,
                         "nucleicAcid"        : "ResNam IN ( 'A', 'C', 'G', 'I', 'T', 'U', 'DA', 'DC', 'DG', 'DI', 'DT', 'DU' )" ,
                         "protein"            : _protein                                           ,
                         "qcAtoms"            : "IsQCAtom = 1"                                     ,
                         "ringAtoms"          : "InRing = 1"                                       ,
                         "sidechain"          : _protein + " and Label NOT IN ( 'C', 'CA', 'N', 'O' ) and ( AtomicNumber != 1 )" ,
                         "water"              : "ResNam IN ( 'HOH', 'H2O', 'TIP', 'TIP3', 'WAT' )" } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( SQLAtomSelector, self )._CheckOptions ( )
        self.connection = sqlite3.connect ( self.path )
        if not self.CheckForExistingTable ( ): self.BuildAtomsTable ( )

    def BuildAtomsTable ( self ):
        """Build an atoms table."""
        # . Gather column data.
        columnData = self.DefineBasicAtomsColumns        ( ) + \
                     self.DefineConnectivityAtomsColumns ( ) + \
                     self.DefineMMAtomsColumns           ( ) + \
                     self.DefineSequenceAtomsColumns     ( )
        atomsColumns = []
        for item in columnData: atomsColumns.append ( item.split ( )[0] )
        atomsColumns.sort ( )
        self.atomsColumns = atomsColumns
        # . Create the table.
        cursor = self.connection.cursor ( )
        cursor.execute ( "CREATE TABLE Atoms ( " + ", ".join ( columnData ) + " )" )
        cursor.close ( )
        # . Fill the table.
        self.FillBasicAtomsColumns        ( )
        self.FillConnectivityAtomsColumns ( )
        self.FillMMAtomsColumns           ( )
        self.FillSequenceAtomsColumns     ( )
        self.UpdateCoordinates3           ( )

    def CheckForExistingTable ( self ):
        """Check for an existing table."""
        isOK    = False
        cursor  = self.connection.cursor ( )
        # . Check for table.
        cursor.execute ( "SELECT CASE WHEN tbl_name = 'Atoms' THEN 1 ELSE 0 END FROM sqlite_master WHERE tbl_name = 'Atoms' AND type = 'table'" )
#        print "TABLE DATA> ", cursor.fetchone ( )
        result = cursor.fetchone ( )
        if ( result is not None ) and ( result[0] == 1 ):
            # . Check number of rows.
            cursor.execute ( "SELECT COUNT(*) FROM Atoms" )
#            print "ROW DATA> ", cursor.fetchall ( )
            isOK = ( cursor.fetchone ( )[0] == len ( self.system.atoms ) )
            # . Get the column data.
            if isOK:
                cursor.execute ( "PRAGMA table_info ( Atoms )" )
                atomsColumns = []
                for item in cursor: atomsColumns.append ( item[1] )
#                print "COLUMN DATA> ", atomsColumns
                self.atomsColumns = atomsColumns
# . This also works but less information.
#                cursor.execute ( "SELECT * FROM Atoms" )
#                print "DESCRIPTION DATA> ", cursor.description
        # . Finish up.
        cursor.close ( )
        return isOK

    def Close ( self ): self.connection.close ( )

    def DefineBasicAtomsColumns ( self ):
        """Define basic columns in the atoms table."""
        return [ "RowID          INTEGER PRIMARY KEY"          ,
                 "AtomicNumber   INTEGER DEFAULT -1  NOT NULL" ,
                 "Label          TEXT    DEFAULT ''  NOT NULL" ,
                 "Mass           REAL    DEFAULT 0.0 NOT NULL" ,
                 "IsBoundaryAtom INTEGER DEFAULT 0   NOT NULL" ,
                 "IsFixed        INTEGER DEFAULT 1   NOT NULL" ,
                 "IsMMAtom       INTEGER DEFAULT 0   NOT NULL" ,
                 "IsQCAtom       INTEGER DEFAULT 0   NOT NULL" ,
                 "X              REAL    DEFAULT 0.0 NOT NULL" ,
                 "Y              REAL    DEFAULT 0.0 NOT NULL" ,
                 "Z              REAL    DEFAULT 0.0 NOT NULL" ]

    def DefineConnectivityAtomsColumns ( self ):
        """Define connectivity columns in the atoms table."""
        columns = []
        if ( len ( self.system.connectivity.bonds ) > 0 ) or self.createFullTable:
            columns = [ "Connections  INTEGER DEFAULT NULL" ,
                        "FormalCharge INTEGER DEFAULT NULL" ,
                        "Hydrogens    INTEGER DEFAULT NULL" ,
                        "InRing       INTEGER DEFAULT NULL" ,
                        "IsAromatic   INTEGER DEFAULT NULL" ,
                        "Isolate      INTEGER DEFAULT NULL" ,
                        "Valence      INTEGER DEFAULT NULL" ]
        return columns

    def DefineMMAtomsColumns ( self ):
        """Define MM columns in the atoms table."""
        columns = []
        if ( self.system.mmModel is not None ) or self.createFullTable:
            columns = [ "AtomType    TEXT DEFAULT NULL" ,
                        "Charge      REAL DEFAULT NULL" ,
                        "LJRadius    REAL DEFAULT NULL" ,
                        "LJWellDepth REAL DEFAULT NULL" ]
        return columns

    def DefineSequenceAtomsColumns ( self ):
        """Define sequence columns in the atoms table."""
        columns = []
        if ( self.system.sequence is not None ) or self.createFullTable:
            columns = [ "Component       TEXT    DEFAULT NULL" ,
                        "Entity          TEXT    DEFAULT NULL" ,
                        "InLinearPolymer INTEGER DEFAULT NULL" ,
                        "Path            TEXT    DEFAULT NULL" ,
                        "ICode           TEXT    DEFAULT NULL" ,
                        "ResNam          TEXT    DEFAULT NULL" ,
                        "ResSeq          INTEGER DEFAULT NULL" ]
        return columns

    def FillBasicAtomsColumns ( self ):
        """Fill the basic columns in an atoms table."""
        cursor = self.connection.cursor ( )
        # . Atom data.
        for atom in self.system.atoms:
            label   = atom.label.replace ( "'", "''" )
            cursor.execute ( "INSERT INTO Atoms ( RowID, AtomicNumber, Label, Mass ) VALUES ( {:d}, {:d}, '{:s}', {:f} )".format ( atom.index, atom.atomicNumber, label, atom.mass ) )
        # . Energy model.
        if self.system.mmModel is not None:
            for i in range ( len ( self.system.atoms ) ):
                cursor.execute ( "UPDATE Atoms SET IsMMAtom=1 WHERE RowID={:d}".format ( i ) )
        if self.system.qcModel is not None:
            for i in self.system.qcModel.pureQCAtoms:
                cursor.execute ( "UPDATE Atoms SET IsMMAtom=0, IsQCAtom=1 WHERE RowID={:d}".format ( i ) )
            for i in self.system.qcModel.boundaryAtoms:
                cursor.execute ( "UPDATE Atoms SET IsMMAtom=1, IsQCAtom=1 WHERE RowID={:d}".format ( i ) )
        # . Fixed atoms.
        if self.system.freeAtoms is None: freeAtoms = range ( len ( self.system.atoms ) )
        else:                             freeAtoms = self.system.freeAtoms
        for i in freeAtoms:
            cursor.execute ( "UPDATE Atoms SET IsFixed=0 WHERE RowID={:d}".format ( i ) )
        cursor.close ( )
        self.connection.commit ( )

    def FillConnectivityAtomsColumns ( self ):
        """Fill the connectivity columns in an atoms table."""
        if len ( self.system.connectivity.bonds ) > 0:
            isolates  = self.system.connectivity.isolates
            ringSets  = self.system.connectivity.ringSets
            inIsolate = { atom : n for ( n, isolate ) in enumerate ( isolates ) for atom in isolate }
            inRing    = { atom : 1 for      ringSet   in             ringSets   for atom in ringSet }
            cursor    = self.connection.cursor ( )
            for ( i, atom ) in enumerate ( self.system.atoms ):
                cursor.execute ( "UPDATE Atoms SET Connections={:d}, FormalCharge={:d}, Hydrogens={:d}, InRing={:d}, IsAromatic={:d}, Isolate={:d}, Valence={:d} WHERE RowID={:d}".format \
                                                  ( atom.connections, atom.formalCharge, atom.hydrogens, inRing.get ( atom, 0 ), atom.isAromatic, inIsolate[atom], atom.valence, i ) )
            cursor.close ( )
            self.connection.commit ( )

    def FillMMAtomsColumns ( self ):
        """Fill the MM columns in an atoms table."""
        if self.system.mmState is not None:
            # . Get the data.
            mmState  = self.system.mmState
            ljState  = mmState.ljParameters.__getstate__ ( )
            epsilons = ljState.get ( "epsilons", None )
            sigmas   = ljState.get ( "sigmas"  , None )
            if ( epsilons is None ) or ( sigmas is None ): raise AtomSelectionError ( "Unable to extract Lennard-Jones epsilon and sigma values." )
            # . Update the table.
            cursor  = self.connection.cursor ( )
            for atom in self.system.atoms:
                i        = atom.index
                atomType = mmState.atomTypes[mmState.atomTypeIndices[i]]
                charge   = mmState.charges[i]
                ljIndex  = mmState.ljTypeIndices[i]
                epsilon  = epsilons[ljIndex]
                sigma    = sigmas  [ljIndex]
                cursor.execute ( "UPDATE Atoms SET AtomType='{:s}', Charge={:f}, LJRadius={:f}, LJWellDepth={:f} WHERE RowID={:d}".format ( atomType, charge, sigma, epsilon, i ) )
            cursor.close ( )
            self.connection.commit ( )

    def FillSequenceAtomsColumns ( self ):
        """Fill the sequence columns in an atoms table."""
        if self.system.sequence is not None:
            # . Get the data.
            linearPolymerIndex = self.system.sequence.linearPolymerIndex
            # . Update the table.
            cursor = self.connection.cursor ( )
            for entity in self.system.sequence.children:
                entityLabel = entity.label.replace ( "'", "''" )
                for component in entity.children:
                    componentLabel = component.label.replace ( "'", "''" )
                    if component in linearPolymerIndex: inLinearPolymer = 1
                    else:                               inLinearPolymer = 0
                    ( resNam, resSeq, iCode ) = self.system.sequence.ParseLabel ( componentLabel, fields = 3 )
                    try   : resSeq = int ( resSeq )
                    except: resSeq = -1
                    for atom in component.children:
                        path = atom.path.replace ( "'", "''" )
                        cursor.execute ( "UPDATE Atoms SET Component='{:s}', Entity='{:s}', InLinearPolymer={:d}, Path='{:s}', ICode='{:s}', ResNam='{:s}', ResSeq='{:d}' WHERE RowID={:d}".format \
                                                                                    ( componentLabel, entityLabel, inLinearPolymer, path, iCode, resNam, resSeq, atom.index ) )
            cursor.close ( )
            self.connection.commit ( )

    def UpdateCoordinates3 ( self, coordinates3 = None, selection = None ):
        """Update the coordinates3 in the database."""
        # . Gather data - do more checks here?
        if selection    is None: selection = range ( len ( self.system.atoms ) )
        if coordinates3 is None: xyz       = self.system.coordinates3
        else:                    xyz       = coordinates3
        # . Update coordinates.
        if len ( selection ) > 0:
            cursor = self.connection.cursor ( )
            for i in selection:
                cursor.execute ( "UPDATE Atoms SET X={:f}, Y={:f}, Z={:f} WHERE RowID={:d}".format ( xyz[i,0], xyz[i,1], xyz[i,2], i ) )
            cursor.close ( )
            self.connection.commit ( )

    def Where ( self, selectionString ):
        """Get a selection."""
        cursor = self.connection.cursor ( )
        try   : cursor.execute ( "SELECT RowID FROM Atoms WHERE " + selectionString )
        except: raise AtomSelectionError ( "Invalid selection string: {:s}".format ( selectionString ) )
        items = set ( )
        for item in cursor: items.add ( item[0] )
        cursor.close ( )
        return Selection.FromIterable ( items )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

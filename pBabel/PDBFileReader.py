"""Read data from a PDB file and create a PDB model."""

import os.path

from  pCore                 import Clone                 , \
                                   logFile               , \
                                   LogFileActive         , \
                                   TextFileReader
from  pMolecule             import System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3          , \
                                   Transformation3
from .ExportImport          import _Importer
from .PDBModel              import PDBModelAtom          , \
                                   PDBModelAtomData      , \
                                   PDBModelComponent     , \
                                   PDBModelEntity        , \
                                   PDBModelLinearPolymer , \
                                   PDBModelLink          , \
                                   PDBModel

# . Helper functions for returning B-factors and occupancies?

#==================================================================================================================================
# . Parameters.
#==================================================================================================================================
# . Default link name.
_DefaultLinkName = "GenericSingle"

# . The multiplying factor for ANISOU records.
_PDBANISOUFactor    = 10000.0

# . Disulfide bridge data.
_PDBDisulfideBridge = "Disulfide_Bridge"
_PDBSSBondAtom      = "SG"

#==================================================================================================================================
# . Class.
#==================================================================================================================================
class PDBFileReader ( TextFileReader ):
    """PDBFileReader is the class for reading PDB files."""

    _attributable = dict ( TextFileReader._attributable )
    _classLabel   = "PDB File Reader"
    _withSections = True
    _attributable.update ( { "assemblies"             : None  ,
                             "frequencies"            : None  ,
                             "lastComponent"          : None  ,
                             "models"                 : None  ,
                             "seqResModel"            : None  ,
                             "suppressICodeField"     : False ,
                             "useSegmentEntityLabels" : False } )

    # . Public methods.
    def GetModel ( self, modelNumber = 0 ):
        """Return a model."""
        try:
            model            = self.models[modelNumber]
            model.assemblies = self.assemblies
            return model
        except:
            raise IndexError ( "Invalid model number: {:d}.".format ( modelNumber ) )

    def GetSeqResModel ( self ):
        """Return the SeqRes model."""
        return self.seqResModel

    def Parse ( self, log = logFile, suppressICodeField = False, useSegmentEntityLabels = False ):
        """Parse data from the file."""
        # . Check that parsing has not occurred.
        if not self.isParsed:
            # . Initialization.
            # . Options.
            self.suppressICodeField     = suppressICodeField
            self.useSegmentEntityLabels = useSegmentEntityLabels
            # . Local data.
            connections       = set ( )
            links             = []
            hets              = []
            self.currentModel = PDBModel.WithDefaults ( )
            self.lastLabel    = None
            self.serials      = {}
            self.title        = None
            if LogFileActive ( log ): self.log = log
            # . Global data.
            self.frequencies  = {}
            self.models       = [ self.currentModel ]
             # . Parse the file.
            self.Open ( )
            try:
                while True:
                    line = self.GetLine ( signalWarnings = False )
                    if line is None: break
                    if len ( line ) <= 0: continue
                    label = line[0:6].strip ( ).upper ( )
                    if   ( label == "ANISOU" ): self.__ParseAnisouRecord ( line              )
                    elif ( label == "ATOM"   ): self.__ParseAtomRecord   ( line, False       )
                    elif ( label == "CONECT" ): self.__ParseConectRecord ( line, connections )
                    elif ( label == "CRYST1" ): self.__ParseCryst1Record ( line              )
                    elif ( label == "END"    ): self.__ParseEndRecord    ( line              )
                    elif ( label == "ENDMDL" ): self.__ParseEndMdlRecord ( line              )
                    elif ( label == "HET"    ): self.__ParseHetRecord    ( line, hets        )
                    elif ( label == "HETATM" ): self.__ParseAtomRecord   ( line, True        )
                    elif ( label == "LINK"   ): self.__ParseLinkRecord   ( line, links       )
                    elif ( label == "MASTER" ): self.__ParseMasterRecord ( line              )
                    elif ( label == "MODEL"  ): self.__ParseModelRecord  ( line              )
                    elif ( label == "REMARK" ): self.__ParseRemarkRecord ( line              )
                    elif ( label == "SEQRES" ): self.__ParseSeqResRecord ( line              )
                    elif ( label == "SSBOND" ): self.__ParseSSBondRecord ( line, links       )
                    elif ( label == "TER"    ): self.__ParseTerRecord    ( line              )
                    elif ( label == "TITLE"  ): self.__ParseTitleRecord  ( line              )
                    self.lastLabel = label
                    if label in self.frequencies: self.frequencies[label] += 1
                    else:                         self.frequencies[label]  = 1
            except EOFError:
                pass
            # . Finish parsing.
            if len ( self.models ) > 1:
                if self.currentModel is not self.models[0]: self.Warning ( "Probable missing ENDMDL at end-of-file.", False )
                if self.models[0].numberOfAtoms > 0:    self.Warning ( "ATOM/HETATM records outside MODEL/ENDMDL records.", False )
            if self.lastLabel != "END": self.Warning ( "Missing END record at end-of-file.", False )
            # . Process the CONECT and HET record data.
            self.__ProcessConnections ( connections, links )
            self.__ProcessSeqResModel ( hets               )
            # . Verify the models.
            if self.seqResModel is not None:
                self.seqResModel.Verify ( )
                self.seqResModel.label = self.title
            for model in self.models:
                if model is not None:
                    model.Verify ( )
                    model.label = self.title
            # . Find the number of unknown atoms.
            numberOfUnknowns = 0
            for model in self.models: numberOfUnknowns += model.numberOfUnknowns
            if numberOfUnknowns > 0: self.Warning ( "There are {:d} atoms of unidentified element type.".format ( numberOfUnknowns ), False )
            # . Close the file and warning table.
            self.WarningStop ( )
            self.Close ( )
            # . Clean up.
            for attribute in ( "currentModel", "lastLabel", "serials", "title" ): delattr ( self, attribute )
            # . Everything is now parsed.
            self.log     = None
            self.isParsed = True

    @classmethod
    def PathToCoordinates3 ( selfClass                        ,
                             path                             ,
                             altLoc                 = None    ,
                             log                    = logFile ,
                             modelNumber            = 0       ,
                             suppressICodeField     = False   ,
                             useSegmentEntityLabels = False   ):
        """Return the coordinates from a PDB file."""
        # . Parse the file.
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log                    = log                    ,
                         suppressICodeField     = suppressICodeField     ,
                         useSegmentEntityLabels = useSegmentEntityLabels )
        inFile.Summary ( log = log )
        # . Get the model.
        if ( modelNumber == 0 ) and ( len ( inFile.models ) > 1 ): modelNumber = 1
        model = inFile.GetModel ( modelNumber = modelNumber )
        # . Finish up.
        return model.MakeCoordinates3 ( altLoc )

    @classmethod
    def PathToPDBModel ( selfClass                        ,
                         path                             ,
                         log                    = logFile ,
                         modelNumber            = 0       ,
                         suppressICodeField     = False   ,
                         useSegmentEntityLabels = False   ):
        """Return the PDB model from a PDB file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log                    = log                    ,
                         suppressICodeField     = suppressICodeField     ,
                         useSegmentEntityLabels = useSegmentEntityLabels )
        inFile.Summary ( log = log )
        # . Check the model number.
        if ( modelNumber == 0 ) and ( len ( inFile.models ) > 1 ): modelNumber = 1
        return inFile.GetModel ( modelNumber = modelNumber )

    @classmethod
    def PathToSystem ( selfClass                        ,
                       path                             ,
                       additionalData         = None    ,
                       altLoc                 = None    ,
                       embeddedHydrogens      = False   ,
                       libraryPaths           = None    ,
                       log                    = logFile ,
                       modelNumber            = 0       ,
                       suppressICodeField     = False   ,
                       useComponentLibrary    = False   ,
                       useSegmentEntityLabels = False   ):
        """Return a system defined by the atom data in a PDB file.

        This function will rarely work on an experimental PDB file. In these cases
        it will be necessary to pass by an explicit definition of the PDB model.
        """
        # . Parse the file.
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log                    = log                    ,
                         suppressICodeField     = suppressICodeField     ,
                         useSegmentEntityLabels = useSegmentEntityLabels )
        inFile.Summary ( log = log )
        # . Get the model.
        if ( modelNumber == 0 ) and ( len ( inFile.models ) > 1 ): modelNumber = 1
        model = inFile.GetModel ( modelNumber = modelNumber )
        # . Use the PDB component library.
        if useComponentLibrary:
            # . Build the model.
            rawModel = Clone ( model )
            model.ClearAtoms ( )
            model.MakeAtomicModelFromComponentLibrary ( libraryPaths = libraryPaths, log = log )
            model.ExtractAtomData ( rawModel, log = log )
            model.Summary ( log = log )
        # . Order hydrogens within each component.
        model.OrderHydrogens ( embedded = embeddedHydrogens )
        # . Make the system.
        system = model.MakeSystem ( altLoc = altLoc )
        # . Check if additional data is to be returned.
        if ( additionalData is not None ) and ( len ( additionalData ) > 0 ):
            results = [ system ]
            for attribute in additionalData:
                results.append ( model.MakeData ( attribute, altLoc = altLoc ) )
            results = tuple ( results )
        else:
            results = system
        # . Finish up.
        return results

    def SummaryItems ( self ):
        """Summary items."""
        # . Basic data.
        items = [ ( None       , False                              ) ,
                  ( "Records"  , "{:d}".format ( self.linesParsed ) ) ,
                  ( "Warnings" , "{:d}".format ( self.warnings    ) ) ]
        items.extend ( [ ( "{:<6s} Records".format ( key ), "{:d}".format ( self.frequencies[key] ) ) for key in sorted ( self.frequencies.keys ( ) ) ] )
        # . SeqRes data.
        if self.seqResModel is not None:
            items.append ( ( "SeqRes Sequence", True ) )
            items.extend ( self.seqResModel.SummaryItems ( ) )
        # . Model data.
        if len ( self.models ) > 0:
            for imodel in range ( len ( self.models ) ):
                model = self.models[imodel]
                if model.numberOfAtoms > 0:
                    title = "PDB Sequence Summary"
                    if imodel != 0: title += " for Model {:d}".format ( imodel )
                    items.append ( ( title, True ) )
                    items.extend ( model.SummaryItems ( ) )
        return items

    # . Private methods.
    def __GatherAtomData ( self, line ):
        """Gather data from an ATOM or HETATM record."""
        charge = 0.0
        radius = 0.0
        if self.suppressICodeField:
            iCode = ""
            ( serial            ,
              name              ,
              altLoc            ,
              resName           ,
              chainID           ,
              resSeq            ,
              x                 ,
              y                 ,
              z                 ,
              occupancy         ,
              temperatureFactor ,
              segID             ,
              atomicNumber      ,
              formalCharge      ) = self._ParseFixedFormatLine ( line                    ,
                                                                 (  6, 11, int  , None ) ,
                                                                 ( 12, 16, None , ""   ) ,
                                                                 ( 16, 17, None , ""   ) ,
                                                                 ( 17, 20, None , ""   ) ,
                                                                 ( 21, 22, None , ""   ) ,
                                                                 ( 22, 27, int  , None ) ,
                                                                 ( 30, 38, float, None ) ,
                                                                 ( 38, 46, float, None ) ,
                                                                 ( 46, 54, float, None ) ,
                                                                 ( 54, 60, float, 1.0  ) ,
                                                                 ( 60, 66, float, 0.0  ) ,
                                                                 ( 72, 76, None , ""   ) ,
                                                                 ( 76, 78, PeriodicTable.AtomicNumberFromSymbol, -1 ) ,
                                                                 ( 78, 80, None , ""   ) )
        else:
            ( serial            ,
              name              ,
              altLoc            ,
              resName           ,
              chainID           ,
              resSeq            ,
              iCode             ,
              x                 ,
              y                 ,
              z                 ,
              occupancy         ,
              temperatureFactor ,
              segID             ,
              atomicNumber      ,
              formalCharge      ) = self._ParseFixedFormatLine ( line                    ,
                                                                 (  6, 11, int  , None ) ,
                                                                 ( 12, 16, None , ""   ) ,
                                                                 ( 16, 17, None , ""   ) ,
                                                                 ( 17, 20, None , ""   ) ,
                                                                 ( 21, 22, None , ""   ) ,
                                                                 ( 22, 26, int  , None ) ,
                                                                 ( 26, 27, None , ""   ) ,
                                                                 ( 30, 38, float, None ) ,
                                                                 ( 38, 46, float, None ) ,
                                                                 ( 46, 54, float, None ) ,
                                                                 ( 54, 60, float, 1.0  ) ,
                                                                 ( 60, 66, float, 0.0  ) ,
                                                                 ( 72, 76, None , ""   ) ,
                                                                 ( 76, 78, PeriodicTable.AtomicNumberFromSymbol, -1 ) ,
                                                                 ( 78, 80, None , ""   ) )
        #. Entity label.
        if self.useSegmentEntityLabels: entityLabel = segID
        else:                           entityLabel = chainID
        # . Formal charge.
        try:    formalCharge = int ( formalCharge )
        except: formalCharge = 0
        return ( altLoc            ,
                 atomicNumber      ,
                 charge            ,
                 entityLabel       ,
                 formalCharge      ,
                 iCode             ,
                 name              ,
                 occupancy         ,
                 radius            ,
                 resName           ,
                 resSeq            ,
                 serial            ,
                 temperatureFactor ,
                 x                 ,
                 y                 ,
                 z                 )

    def __ParseAnisouRecord ( self, line ):
        """Process an ANISOU record."""
        # . Parsing.
        if self.suppressICodeField:
            iCode = ""
            ( serial       ,
              name         ,
              altLoc       ,
              resName      ,
              chainID      ,
              resSeq       ,
              u00          ,
              u11          ,
              u22          ,
              u01          ,
              u02          ,
              u12          ,
              segID        ,
              atomicNumber ,
              formalCharge ) = self._ParseFixedFormatLine ( line                    ,
                                                            (  6, 11, int  , None ) ,
                                                            ( 12, 16, None , ""   ) ,
                                                            ( 16, 17, None , ""   ) ,
                                                            ( 17, 20, None , ""   ) ,
                                                            ( 21, 22, None , ""   ) ,
                                                            ( 22, 27, int  , None ) ,
                                                            ( 28, 35, float, None ) ,
                                                            ( 35, 42, float, None ) ,
                                                            ( 42, 49, float, None ) ,
                                                            ( 49, 56, float, None ) ,
                                                            ( 56, 63, float, None ) ,
                                                            ( 63, 70, float, None ) ,
                                                            ( 72, 76, None , ""   ) ,
                                                            ( 76, 78, PeriodicTable.AtomicNumberFromSymbol, -1 ) ,
                                                            ( 78, 80, None , ""   ) )
        else:
            ( serial       ,
              name         ,
              altLoc       ,
              resName      ,
              chainID      ,
              resSeq       ,
              iCode        ,
              u00          ,
              u11          ,
              u22          ,
              u01          ,
              u02          ,
              u12          ,
              segID        ,
              atomicNumber ,
              formalCharge ) = self._ParseFixedFormatLine ( line                    ,
                                                            (  6, 11, int  , None ) ,
                                                            ( 12, 16, None , ""   ) ,
                                                            ( 16, 17, None , ""   ) ,
                                                            ( 17, 20, None , ""   ) ,
                                                            ( 21, 22, None , ""   ) ,
                                                            ( 22, 26, int  , None ) ,
                                                            ( 26, 27, None , ""   ) ,
                                                            ( 28, 35, float, None ) ,
                                                            ( 35, 42, float, None ) ,
                                                            ( 42, 49, float, None ) ,
                                                            ( 49, 56, float, None ) ,
                                                            ( 56, 63, float, None ) ,
                                                            ( 63, 70, float, None ) ,
                                                            ( 72, 76, None , ""   ) ,
                                                            ( 76, 78, PeriodicTable.AtomicNumberFromSymbol, -1 ) ,
                                                            ( 78, 80, None , ""   ) )
        if u00 is not None: u00 *= _PDBANISOUFactor
        if u11 is not None: u11 *= _PDBANISOUFactor
        if u22 is not None: u22 *= _PDBANISOUFactor
        if u01 is not None: u01 *= _PDBANISOUFactor
        if u02 is not None: u02 *= _PDBANISOUFactor
        if u12 is not None: u12 *= _PDBANISOUFactor
        # . Warnings.
        if ( resSeq is None ) or ( u00 is None ) or ( u11 is None ) or ( u22 is None ) or ( u01 is None ) or ( u02 is None ) or ( u12 is None ): self.Warning ( line, False )
        #. Entity label.
        if self.useSegmentEntityLabels: entityLabel = segID
        else:                           entityLabel = chainID
        # . Processing.
        # . Get data.
        data = self.currentModel.GetDescendantFromLabels ( entityLabel, self.currentModel.MakeLabel ( resName, "{:d}".format ( resSeq ), iCode ), name, altLoc )
        if data is None:
            path = self.currentModel.MakePath ( entityLabel, self.currentModel.MakeLabel ( resName, "{:d}".format ( resSeq ), iCode ), name, altLoc )
            self.Warning ( "The ATOM/HETATM record corresponding to the ANISOU record for \"" + path + "\" is missing.", False )
        else:
            data.u00 = u00
            data.u01 = u01
            data.u02 = u02
            data.u11 = u11
            data.u12 = u12
            data.u22 = u22

    def __ParseAtomRecord ( self, line, isHeteroatom ):
        """Process an ATOM or HETATM record."""
        # . Parsing.
        ( altLoc            ,
          atomicNumber      ,
          charge            ,
          entityLabel       ,
          formalCharge      ,
          iCode             ,
          name              ,
          occupancy         ,
          radius            ,
          resName           ,
          resSeq            ,
          serial            ,
          temperatureFactor ,
          x                 ,
          y                 ,
          z                 ) = self.__GatherAtomData ( line )
        # . Check for a warning.
        if ( resSeq is None ) or ( x is None ) or ( y is None ) or ( z is None ): self.Warning ( line, False )
        # . Processing.
        # . Get path items.
        componentLabel = self.currentModel.MakeLabel ( resName, "{:d}".format ( resSeq ), iCode )
        ( entity, component, atom ) = self.currentModel.GetDescendantsFromLabels ( entityLabel, componentLabel, name )
        if entity is None:
            entity = PDBModelEntity.WithOptions ( label = entityLabel )
            self.currentModel.AddChild ( entity )
        if component is None:
                   component = PDBModelComponent.WithOptions ( genericLabel = resName, label = componentLabel )
                   entity.AddChild ( component )
                   self.lastComponent = component
        if atom is None:
                    atom = PDBModelAtom.WithOptions ( atomicNumber = atomicNumber, formalCharge = formalCharge, label = name )
                    component.AddChild ( atom )
        # . Get data.
        data = atom.childIndex.get ( altLoc, None )
        if data is None:
                   data = PDBModelAtomData.WithOptions ( charge            = charge            ,
                                                         label             = altLoc            ,
                                                         occupancy         = occupancy         ,
                                                         radius            = radius            ,
                                                         temperatureFactor = temperatureFactor , x = x, y = y, z = z )
                   atom.AddChild ( data )
        else:
                   self.Warning ( "Duplicate ATOM/HETATM records for \"" + data.path + "\".", False )
        # . Serial number.
        if serial in self.serials:
            if len ( self.models ) == 1: self.Warning ( "There are duplicate atom serial numbers: {:d}.".format ( serial ), False )
        else:
            self.serials[serial] = atom

    def __ParseConectRecord ( self, line, connections ):
        """Process a CONECT record."""
        # . Parsing.
        tokens  = self._ParseFixedFormatLine ( line                  ,
                                               (  6, 11, int, None ) ,
                                               ( 11, 16, int, None ) ,
                                               ( 16, 21, int, None ) ,
                                               ( 21, 26, int, None ) ,
                                               ( 26, 31, int, None ) )
        serial1 = tokens.pop ( 0 )
        if ( serial1 is not None ):
            while len ( tokens ) > 0:
                serial2 = tokens.pop ( 0 )
                if ( serial2 is not None ) and ( serial1 != serial2 ):
                    maxSerial = max ( serial1, serial2 )
                    minSerial = min ( serial1, serial2 )
                    connections.add ( ( maxSerial, minSerial ) )

    def __ParseCryst1Record ( self, line ):
        """Process a CRYST1 record."""
        # . Parsing.
        ( a, b, c, alpha, beta, gamma, sGroup, zValue ) = self._ParseFixedFormatLine ( line,
                                                                                       (  6, 15, float, None ) ,
                                                                                       ( 15, 24, float, None ) ,
                                                                                       ( 24, 33, float, None ) ,
                                                                                       ( 33, 40, float, None ) ,
                                                                                       ( 40, 47, float, None ) ,
                                                                                       ( 47, 54, float, None ) ,
                                                                                       ( 55, 66, None , ""   ) ,
                                                                                       ( 66, 70, int  , None ) )
        # . Check syntax.
        if ( self.frequencies.get ( "CRYST1", 0 ) != 0 ): self.Warning ( "There are multiple CRYST1 records.", False )

    def __ParseEndRecord ( self, line ):
        """Process an END record."""
        # . Check syntax.
        if ( self.frequencies.get ( "MASTER", 0 ) != 1 ) or ( self.lastLabel != "MASTER" ): self.Warning ( "MASTER record missing before END record.", False )

    def __ParseEndMdlRecord ( self, line ):
        """Process an ENDMDL record."""
        # . Check syntax.
        if ( self.frequencies.get ( "MODEL", 0 ) != ( self.frequencies.get ( "ENDMDL", 0 ) + 1 ) ) or ( self.currentModel is self.models[0] ): self.Warning ( "Unpaired MODEL/ENDMDL records.", False )
        # . Processing.
        self.currentModel = self.models[0]

    # . Should move this to TextFileReader and rationalize with GetFixedFormatTokens.
    def _ParseFixedFormatLine ( self, *arguments ):
        """Get tokens in fixed format from a line."""
        # . Initialization.
        tokens = []
        length = 0
        for arg in arguments[1:]: length = max ( length, arg[1] )
        # . Get the line.
        line   = arguments[0].ljust ( length )
        # . Parse the line.
        for ( i, ( start, stop, converter, default ) ) in enumerate ( arguments[1:] ):
            word  = line[start:stop].strip ( )
            token = default
            if len ( word ) > 0:
                if converter is None:
                    token = word
                else:
                    try:    token = converter ( word )
                    except: self.Warning ( "Unable to convert token {:d}.".format ( i ), False )
            tokens.append ( token )
        return tokens

    def __ParseHetRecord ( self, line, hets ):
        """Process a HET record."""
        # . Parsing.
        ( resName, chainID, seqNum, iCode ) = self._ParseFixedFormatLine ( line                   ,
                                                                           (  7, 10, None, ""   ) ,
                                                                           ( 12, 13, None, ""   ) ,
                                                                           ( 13, 17, int , None ) ,
                                                                           ( 17, 18, None, ""   ) )
        # . Processing.
        hets.append ( ( chainID, resName, seqNum, iCode ) )

    def __ParseLinkRecord ( self, line, links ):
        """Process a LINK record."""
        # . Parsing.
        ( name1    ,
          altLoc1  ,
          resName1 ,
          chainID1 ,
          seqNum1  ,
          iCode1   ,
          name2    ,
          altLoc2  ,
          resName2 ,
          chainID2 ,
          seqNum2  ,
          iCode2   ,
          sym1     ,
          sym2     ) = self._ParseFixedFormatLine ( line                   ,
                                                    ( 12, 16, None, ""   ) ,
                                                    ( 16, 17, None, ""   ) ,
                                                    ( 17, 20, None, ""   ) ,
                                                    ( 21, 22, None, ""   ) ,
                                                    ( 22, 26, int , None ) ,
                                                    ( 26, 27, None, ""   ) ,
                                                    ( 42, 46, None, ""   ) ,
                                                    ( 46, 47, None, ""   ) ,
                                                    ( 47, 50, None, ""   ) ,
                                                    ( 51, 52, None, ""   ) ,
                                                    ( 52, 56, int , None ) ,
                                                    ( 56, 57, None, ""   ) ,
                                                    ( 59, 65, None, ""   ) ,
                                                    ( 66, 72, None, ""   ) )
        # . Warnings.
        if ( seqNum1 is None ) or ( seqNum2 is None ): self.Warning ( "Unidentified sequence numbers in LINK record.", False )
        if ( sym1 != sym2 ): self.Warning ( "Unable to process non-identical symmetry identifiers in LINK record.", False )
        # . Processing - altLocs are ignored.
        links.append ( ( _DefaultLinkName, ( chainID1, self.currentModel.MakeLabel ( resName1, "{:d}".format ( seqNum1 ), iCode1 ), name1 ),
                                           ( chainID2, self.currentModel.MakeLabel ( resName2, "{:d}".format ( seqNum2 ), iCode2 ), name2 ) ) )

    def __ParseMasterRecord ( self, line ):
        """Process a MASTER record."""
        # . Check syntax.
        if self.frequencies.get ( "MASTER", 0 ) > 0: self.Warning ( "There are multiple MASTER records.", False )

    def __ParseModelRecord ( self, line ):
        """Process a MODEL record."""
        # . If there is one previous model it must be of zero length.
        if len ( self.models ) == 1:
            if self.models[0].numberOfAtoms > 0: self.Warning ( "There are ATOM/HETATM records before the first MODEL record.", False )
        else:
            if ( self.lastLabel != "ENDMDL" ) or ( self.currentModel is not self.models[0] ): self.Warning ( "An ENDMDL record must precede this MODEL record.", False )
        # . Processing.
        self.currentModel  = PDBModel.WithDefaults ( )
        self.lastComponent = None
        self.models.append ( self.currentModel )

    def __ParseRemarkRecord ( self, line ):
        """Process a REMARK record."""
        # . Only 350 records.
        if line[7:10] == "350":
            if line.find ( "APPLY THE FOLLOWING TO CHAINS:" ) >= 0:
                chains = line.split ( ":", 1 )[-1].replace ( ",", " " ).split ( )
                chains.sort ( )
                if self.assemblies is None: self.assemblies = []
                self.assemblies.append ( ( tuple ( chains ), [] ) )
            elif line[13:18] == "BIOMT":
                ( i, n, ri0, ri1, ri2, ti ) = self.TokenizeLine ( line[18:], converters = [ int, int, float, float, float, float ] )
                if self.assemblies is not None:
                    transformations = self.assemblies[-1][-1]
                    if n > len ( transformations ): transformations.append ( Transformation3.Null ( ) )
                    transformation  = transformations[-1]
                    transformation.rotation[i-1,0]  = ri0
                    transformation.rotation[i-1,1]  = ri1
                    transformation.rotation[i-1,2]  = ri2
                    transformation.translation[i-1] = ti

    def __ParseSeqResRecord ( self, line ):
        """Process a SEQRES record."""
        if self.seqResModel is None: self.seqResModel = PDBModel.WithDefaults ( )
        # . Parsing.
        arguments = [ line, ( 8, 10, int, None ), ( 11, 12, None, "" ), ( 13, 17, int, None ) ]
        for i in range ( 19, 70, 4 ): arguments.append ( ( i, i+3, None, "" ) )
        tokens   = self._ParseFixedFormatLine ( *arguments )
        serNum   = tokens.pop ( 0 )
        chainID  = tokens.pop ( 0 )
        numRes   = tokens.pop ( 0 )
        resNames = []
        for token in tokens:
            if len ( token ) > 0: resNames.append ( token )
        # . Processing - all entities are chains.
        if len ( resNames ) > 0:
            entity = self.seqResModel.childIndex.get ( chainID, None )
            if entity is None:
                        entity = PDBModelEntity.WithOptions ( label = chainID )
                        self.seqResModel.AddChild ( entity )
            resSeq = len ( entity.children )
            for resName in resNames:
                resSeq += 1
                component = PDBModelComponent.WithOptions ( genericLabel = resName, label = self.seqResModel.MakeLabel ( resName, "{:d}".format ( resSeq ), "" ) )
                entity.AddChild ( component )

    def __ParseSSBondRecord ( self, line, links ):
        """Process a SSBOND record."""
        # . Parsing.
        ( serial   ,
          resName1 ,
          chainID1 ,
          seqNum1  ,
          iCode1   ,
          resName2 ,
          chainID2 ,
          seqNum2  ,
          iCode2   ,
          sym1     ,
          sym2     ) = self._ParseFixedFormatLine ( line                   ,
                                                    (  7, 10, int , None ) ,
                                                    ( 11, 14, None, ""   ) ,
                                                    ( 15, 16, None, ""   ) ,
                                                    ( 17, 21, int , None ) ,
                                                    ( 21, 22, None, ""   ) ,
                                                    ( 25, 28, None, ""   ) ,
                                                    ( 29, 30, None, ""   ) ,
                                                    ( 31, 35, int , None ) ,
                                                    ( 35, 36, None, ""   ) ,
                                                    ( 59, 65, None, ""   ) ,
                                                    ( 66, 72, None, ""   ) )
        # . Warnings.
        if ( seqNum1 is None ) or ( seqNum2 is None ): self.Warning ( line, False )
        if ( sym1 != sym2 ): self.Warning ( "Unable to process non-identical symmetry identifiers in LINK record.", False )
        # . Processing - altLocs are ignored.
        links.append ( ( _PDBDisulfideBridge, ( chainID1, self.currentModel.MakeLabel ( resName1, "{:d}".format ( seqNum1 ), iCode1 ), _PDBSSBondAtom ) ,
                                              ( chainID2, self.currentModel.MakeLabel ( resName2, "{:d}".format ( seqNum2 ), iCode2 ), _PDBSSBondAtom ) ) )

    def __ParseTerRecord ( self, line ):
        """Process a TER record."""
        # . Parsing.
        ( serial, resName, chainID, resSeq, iCode ) = self._ParseFixedFormatLine ( line                   ,
                                                                                   ( 6, 11, int  , None ) ,
                                                                                   ( 17, 20, None, ""   ) ,
                                                                                   ( 21, 22, None, ""   ) ,
                                                                                   ( 22, 26, int , None ) ,
                                                                                   ( 26, 27, None, ""   ) )
        # . Processing.
        # . Get the component.
        isOK = ( resSeq is not None )
        if isOK:
            componentLabel = self.currentModel.MakeLabel ( resName, "{:d}".format ( resSeq ), iCode )
            ( entity, component, atom ) = self.currentModel.GetDescendantsFromLabels ( chainID, componentLabel, "" )
            isOK = ( component is not None )
        if not isOK:
            self.Warning ( "Invalid component specification on TER record.", False )
            component = self.lastComponent
            if component is not None: entity = component.parent
        # . Add the chain to the model.
        if ( component is not None ):
            chain = PDBModelLinearPolymer.WithOptions ( leftTerminalComponent = entity.children[0], rightTerminalComponent = component )
            isOK  = self.currentModel.AddLinearPolymer ( chain )
            if not isOK: self.Warning ( "Error adding linear polymer to model.", False )
        # . Finish up.
        self.lastComponent = None

    def __ParseTitleRecord ( self, line ):
        """Process a TER record."""
        # . Parsing.
        ( continuation, title ) = self._ParseFixedFormatLine ( line, ( 8, 10, int, None ), ( 10, 70, None, "" ) )
        # . Processing.
        if len ( title ) > 0:
            if self.title is None: self.title  = title.title ( )
            else:                  self.title += title.title ( )

    def __ProcessConnections ( self, connections, links ):
        """Process the connections and links for the models."""
        # . Links.
        for model in self.models:
            if model is not None:
                for ( label, pathLabels1, pathLabels2 ) in links:
                    atom1 = model.GetDescendantFromLabels ( *pathLabels1 )
                    atom2 = model.GetDescendantFromLabels ( *pathLabels2 )
                    if ( atom1 is None ) or ( atom2 is None ):
                        self.Warning ( "Unknown atoms in link: {:s} and {:s}.".format ( model.MakePath ( *pathLabels1 ), model.MakePath ( *pathLabels2 ) ), False )
                    else:
                        model.AddLink ( PDBModelLink.WithOptions ( label = label, leftComponent = atom1.parent, rightComponent = atom2.parent ) )
        # . Connections.
        if len ( connections ) > 0:
            numberOfCross    = 0
            numberOfUnknowns = 0
            for ( serial1, serial2 ) in connections:
                atom1 = self.serials.get ( serial1, None )
                atom2 = self.serials.get ( serial2, None )
                if ( atom1 is None ) or ( atom2 is None ): numberOfUnknowns += 1
                else:
                    model1 = atom1.root
                    model2 = atom2.root
                    if model1 is model2: model1.AddLink ( PDBModelLink.WithOptions ( leftComponent = atom1.parent, rightComponent = atom2.parent ) )
                    else:                numberOfCross += 1
            if numberOfUnknowns > 0: self.Warning ( "There are {:d} connections with undefined atom serial numbers.".format ( numberOfUnknowns ), False )
            if numberOfCross    > 0: self.Warning ( "There are {:d} connections between different models."          .format ( numberOfCross    ), False )

    def __ProcessSeqResModel ( self, hets ):
        """Process the chains and het residues for the SeqRes model."""
        if len ( hets ) > 0:
            if self.seqResModel is None:
                self.seqResModel = PDBModel.WithDefaults ( )
            else:
                    # . Define chains.
                for entity in self.seqResModel.children:
                    chain = PDBModelLinearPolymer.WithOptions ( leftTerminalComponent = entity.children[0], rightTerminalComponent = entity.children[-1] )
                    self.seqResModel.AddLinearPolymer ( chain )
            # . Process hets.
            for ( chainID, resName, seqNum, iCode ) in hets:
                component = PDBModelComponent.WithOptions ( genericLabel = resName, label = self.seqResModel.MakeLabel ( resName, "{:d}".format ( seqNum ), iCode ) )
                entity    = self.seqResModel.childIndex.get ( chainID, None )
                # . Skip components that are already within chains.
                if ( entity is not None ) and ( component in entity.childIndex ): continue
                else:
                    if entity is None:
                        entity = PDBModelEntity.WithOptions ( label = chainID )
                        self.seqResModel.AddChild ( entity )
                    entity.AddChild ( component )

#==================================================================================================================================
# . Class.
#==================================================================================================================================
class PQRFileReader ( PDBFileReader ):
    """PQRFileReader is the class for reading PQR files."""

    # . The only difference with a PDB file is in the format of the ATOM/HETATM line.
    def __GatherAtomData ( self, line ):
        """Gather data from an ATOM or HETATM record."""
        altLoc            = ""
        atomicNumber      = -1
        formalCharge      =  0
        iCode             = ""
        occupancy         = 1.0
        temperatureFactor = 0.0
        ( recordName  ,
          serial      ,
          name        ,
          resName     ,
          entityLabel ,
          resSeq      ,
          x           ,
          y           ,
          z           ,
          charge      ,
          radius      ) = self.TokenizeLine ( line, converters = [ None, int, None, None, None, int, float, float, float, float, float ] )
        return ( altLoc            ,
                 atomicNumber      ,
                 charge            ,
                 entityLabel       ,
                 formalCharge      ,
                 iCode             ,
                 name              ,
                 occupancy         ,
                 radius            ,
                 resName           ,
                 resSeq            ,
                 serial            ,
                 temperatureFactor ,
                 x                 ,
                 y                 ,
                 z                 )

#==================================================================================================================================
# . Importer definitions.
#==================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : PDBFileReader.PathToCoordinates3 ,
                         PDBModel     : PDBFileReader.PathToPDBModel     , # . Not really used for the moment.
                         System       : PDBFileReader.PathToSystem       } ,
                       [ "ent", "ENT", "pdb", "PDB" ], "Protein Data Bank", defaultFunction = PDBFileReader.PathToSystem )
_Importer.AddHandler ( { Coordinates3 : PQRFileReader.PathToCoordinates3 ,
                         System       : PQRFileReader.PathToSystem       } ,
                       [ "pqr", "PQR" ], "PQR", defaultFunction = PQRFileReader.PathToSystem )

#==================================================================================================================================
# . Test.
#==================================================================================================================================
if __name__ == "__main__":
    pass

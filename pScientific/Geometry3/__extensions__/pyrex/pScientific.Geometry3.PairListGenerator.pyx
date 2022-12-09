"""Pairlist generation class and functions."""

from  pCore          import Align                , \
                            logFile              , \
                            LogFileActive        , \
                            RawObjectConstructor
from .Geometry3Error import Geometry3Error

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairListGenerator:
    """A pairlist generator."""

    def __copy__ ( self ):
        """Copying."""
        options = self.__getstate__ ( )
        new     = self.__class__ ( **options )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: PairListGenerator_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "sortIndices"          : self.cObject.sortIndices   == CTrue ,
                 "useGridByCell"        : self.cObject.useGridByCell == CTrue ,
                 "minimumCellExtent"    : self.cObject.minimumCellExtent      ,
                 "minimumPoints"        : self.cObject.minimumPoints          ,
                 "cellSize"             : self.cObject.cellSize               ,
                 "cutOff"               : self.cObject.cutOff                 ,
                 "cutOffCellSizeFactor" : self.cObject.cutOffCellSizeFactor   ,
                 "minimumCellSize"      : self.cObject.minimumCellSize        ,
                 "minimumExtentFactor"  : self.cObject.minimumExtentFactor    }

    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions ( **options )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def __str__ ( self ): return "Pair List Generator"

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = PairListGenerator_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def CrossPairListFromDoubleCoordinates3 ( self, Coordinates3         coordinates31 not None ,
                                                    Coordinates3         coordinates32 not None ,
                                                    RealArray1D          radii1                 ,
                                                    RealArray1D          radii2                 ,
                                                    Selection            andSelection1          ,
                                                    Selection            andSelection2          ,
                                                    Selection            orSelection1           ,
                                                    Selection            orSelection2           ,
                                                    PairList             exclusions             ,
                                                    RegularGrid          grid1                  ,
                                                    RegularGridOccupancy occupancy1             ):
        """Cross-pairlist from a double set of coordinates."""
        cdef CrossPairList          pairList
        cdef CPairList             *cExclusions    = NULL
        cdef CPairList             *cPairList      = NULL
        cdef CRealArray1D          *cRadii1        = NULL
        cdef CRealArray1D          *cRadii2        = NULL
        cdef CRegularGrid          *cGrid1         = NULL
        cdef CRegularGridOccupancy *cOccupancy1    = NULL
        cdef CSelection            *cAndSelection1 = NULL
        cdef CSelection            *cAndSelection2 = NULL
        cdef CSelection            *cOrSelection1  = NULL
        cdef CSelection            *cOrSelection2  = NULL
        cdef CStatus                status         = CStatus_OK
        if exclusions    is not None: cExclusions    = exclusions.cObject
        if radii1        is not None: cRadii1        = radii1.cObject
        if radii2        is not None: cRadii2        = radii2.cObject
        if grid1         is not None: cGrid1         = grid1.cObject
        if occupancy1    is not None: cOccupancy1    = occupancy1.cObject
        if andSelection1 is not None: cAndSelection1 = andSelection1.cObject
        if andSelection2 is not None: cAndSelection2 = andSelection2.cObject
        if orSelection1  is not None: cOrSelection1  = orSelection1.cObject
        if orSelection2  is not None: cOrSelection2  = orSelection2.cObject
        cPairList = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( self.cObject          ,
                                                                            coordinates31.cObject ,
                                                                            coordinates32.cObject ,
                                                                            cRadii1               ,
                                                                            cRadii2               ,
                                                                            cAndSelection1        ,
                                                                            cAndSelection2        ,
                                                                            cOrSelection1         ,
                                                                            cOrSelection2         ,
                                                                            cExclusions           ,
                                                                            cGrid1                ,
                                                                            cOccupancy1           ,
                                                                            &status               )
        if status != CStatus_OK: raise Geometry3Error ( "Error generating cross pairlist from a double set of coordinates." )
        pairList         = CrossPairList.Raw ( )
        pairList.cObject = cPairList
        pairList.isOwner = True
        return pairList

    def CrossPairListFromSingleCoordinates3 ( self, Coordinates3         coordinates3 not None ,
                                                    RealArray1D          radii                 ,
                                                    Selection            andSelection1         ,
                                                    Selection            andSelection2         ,
                                                    Selection            orSelection           ,
                                                    PairList             exclusions            ,
                                                                         excludeSelf           ,
                                                    RegularGrid          grid1                 ,
                                                    RegularGridOccupancy occupancy1            ):
        """Cross-pairlist from a single set of coordinates."""
        cdef CrossPairList          pairList
        cdef CBoolean               cExcludeSelf
        cdef CPairList             *cExclusions    = NULL
        cdef CPairList             *cPairList      = NULL
        cdef CRealArray1D          *cRadii         = NULL
        cdef CRegularGrid          *cGrid1         = NULL
        cdef CRegularGridOccupancy *cOccupancy1    = NULL
        cdef CSelection            *cAndSelection1 = NULL
        cdef CSelection            *cAndSelection2 = NULL
        cdef CSelection            *cOrSelection   = NULL
        cdef CStatus                status         = CStatus_OK
        if excludeSelf              : cExcludeSelf   = CTrue
        else:                         cExcludeSelf   = CFalse
        if exclusions    is not None: cExclusions    = exclusions.cObject
        if radii         is not None: cRadii         = radii.cObject
        if grid1         is not None: cGrid1         = grid1.cObject
        if occupancy1    is not None: cOccupancy1    = occupancy1.cObject
        if andSelection1 is not None: cAndSelection1 = andSelection1.cObject
        if andSelection2 is not None: cAndSelection2 = andSelection2.cObject
        if orSelection   is not None: cOrSelection   = orSelection.cObject
        cPairList = PairListGenerator_CrossPairListFromSingleCoordinates3 ( self.cObject         ,
                                                                            coordinates3.cObject ,
                                                                            cRadii               ,
                                                                            cAndSelection1       ,
                                                                            cAndSelection2       ,
                                                                            cOrSelection         ,
                                                                            cExclusions          ,
                                                                            cExcludeSelf         ,
                                                                            cGrid1               ,
                                                                            cOccupancy1          ,
                                                                            &status              )
        if status != CStatus_OK: raise Geometry3Error ( "Error generating cross pairlist from a single set of coordinates." )
        pairList         = CrossPairList.Raw ( )
        pairList.cObject = cPairList
        pairList.isOwner = True
        return pairList

    def OptionRecords ( self ):
        """Option records and subobjects that also have options."""
        return ( [ ( "cellSize"             , "Cell Size"               , "float", "{:.3f}".format ( self.cObject.cellSize               ) ) ,
                   ( "cutOff"               , "List CutOff"             , "float", "{:.3f}".format ( self.cObject.cutOff                 ) ) ,
                   ( "cutOffCellSizeFactor" , "CutOff Cell Size Factor" , "float", "{:.3f}".format ( self.cObject.cutOffCellSizeFactor   ) ) ,
                   ( "useGridByCell"        , "Grid Cell/Cell Method"   , "bool" , "{!r}"  .format ( self.cObject.useGridByCell == CTrue ) ) ,
                   ( "minimumCellExtent"    , "Minimum Cell Extent"     , "int"  , "{:d}"  .format ( self.cObject.minimumCellExtent      ) ) ,
                   ( "minimumCellSize"      , "Minimum Cell Size"       , "float", "{:.3f}".format ( self.cObject.minimumCellSize        ) ) ,
                   ( "minimumExtentFactor"  , "Minimum Extent Factor"   , "float", "{:.3f}".format ( self.cObject.minimumExtentFactor    ) ) ,
                   ( "minimumPoints"        , "Minimum Points"          , "int"  , "{:d}"  .format ( self.cObject.minimumPoints          ) ) ,
                   ( "sortIndices"          , "Sort Indices"            , "bool" , "{!r}"  .format ( self.cObject.sortIndices == CTrue   ) ) ], [] )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SelfPairListFromCoordinates3 ( self, Coordinates3         coordinates3 not None ,
                                             RealArray1D          radii                 ,
                                             Selection            andSelection          ,
                                             Selection            orSelection           ,
                                             PairList             exclusions            ,
                                             RegularGrid          grid                  ,
                                             RegularGridOccupancy occupancy             ):
        """Self-pairlist from a single set of coordinates."""
        cdef SelfPairList           pairList
        cdef CPairList             *cExclusions   = NULL
        cdef CPairList             *cPairList     = NULL
        cdef CRealArray1D          *cRadii        = NULL
        cdef CRegularGrid          *cGrid         = NULL
        cdef CRegularGridOccupancy *cOccupancy    = NULL
        cdef CSelection            *cAndSelection = NULL
        cdef CSelection            *cOrSelection  = NULL
        cdef CStatus                status        = CStatus_OK
        if exclusions   is not None: cExclusions   = exclusions.cObject
        if radii        is not None: cRadii        = radii.cObject
        if grid         is not None: cGrid         = grid.cObject
        if occupancy    is not None: cOccupancy    = occupancy.cObject
        if andSelection is not None: cAndSelection = andSelection.cObject
        if orSelection  is not None: cOrSelection  = orSelection.cObject
        cPairList = PairListGenerator_SelfPairListFromCoordinates3 ( self.cObject         ,
                                                                     coordinates3.cObject ,
                                                                     cRadii               ,
                                                                     cAndSelection        ,
                                                                     cOrSelection         ,
                                                                     cExclusions          ,
                                                                     cGrid                ,
                                                                     cOccupancy           ,
                                                                     &status              )
        if status != CStatus_OK: raise Geometry3Error ( "Error generating self pairlist from a set of coordinates." )
        pairList         = SelfPairList.Raw ( )
        pairList.cObject = cPairList
        pairList.isOwner = True
        return pairList

    def SetOptions ( self, **options ):
        """Set options for the model."""
        if "sortIndices"          in options:
            value = options.pop ( "sortIndices"   )
            if value: self.cObject.sortIndices   = CTrue
            else:     self.cObject.sortIndices   = CFalse
        if "useGridByCell"        in options:
            value = options.pop ( "useGridByCell" )
            if value: self.cObject.useGridByCell = CTrue
            else:     self.cObject.useGridByCell = CFalse
        if "minimumCellExtent"    in options: self.cObject.minimumCellExtent    = options.pop ( "minimumCellExtent"    )
        if "minimumPoints"        in options: self.cObject.minimumPoints        = options.pop ( "minimumPoints"        )
        if "cellSize"             in options: self.cObject.cellSize             = options.pop ( "cellSize"             )
        if "cutOff"               in options: self.cObject.cutOff               = options.pop ( "cutOff"               )
        if "cutOffCellSizeFactor" in options: self.cObject.cutOffCellSizeFactor = options.pop ( "cutOffCellSizeFactor" )
        if "minimumCellSize"      in options: self.cObject.minimumCellSize      = options.pop ( "minimumCellSize"      )
        if "minimumExtentFactor"  in options: self.cObject.minimumExtentFactor  = options.pop ( "minimumExtentFactor"  )
        if len ( options ) > 0: raise ValueError ( "Invalid options: " + ", ".join ( sorted ( options.keys ( ) ) ) + "." )
        # . Ensure the cell size is correct.
        self.cObject.cellSize = self.cObject.cutOffCellSizeFactor * self.cObject.cutOff

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Pairlist Generator Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Cell Size"               , "{:.3f}".format ( self.cObject.cellSize               ) ) ,
                 ( "List CutOff"             , "{:.3f}".format ( self.cObject.cutOff                 ) ) ,
                 ( "CutOff Cell Size Factor" , "{:.3f}".format ( self.cObject.cutOffCellSizeFactor   ) ) ,
                 ( "Grid Cell/Cell Method"   , "{!r}"  .format ( self.cObject.useGridByCell == CTrue ) ) ,
                 ( "Minimum Cell Extent"     , "{:d}"  .format ( self.cObject.minimumCellExtent      ) ) ,
                 ( "Minimum Cell Size"       , "{:.3f}".format ( self.cObject.minimumCellSize        ) ) ,
                 ( "Minimum Extent Factor"   , "{:.3f}".format ( self.cObject.minimumExtentFactor    ) ) ,
                 ( "Minimum Points"          , "{:d}"  .format ( self.cObject.minimumPoints          ) ) ,
                 ( "Sort Indices"            , "{!r}"  .format ( self.cObject.sortIndices == CTrue   ) ) ]

    def TableOfOptions ( self, log = logFile ):
        """Output a table with the default options."""
        if LogFileActive ( log ):
            ( records, subObjects ) = self.OptionRecords ( )
            alignments = [ Align.Left ] * ( len ( records[0] ) - 1 ) + [ Align.Right ]
            headers    = [ "Option", "Description", "Type", "Default" ]
            title      = "{:s} Option Table".format ( str ( self ) )
            log.TableOfRecords ( records, alignments = alignments, headers = headers, title = title )
            if len ( subObjects ) > 0:
                for subObject in subObjects:
                    subObject.TableOfOptions ( log = log )

    def UseGridSearch ( self, Coordinates3 coordinates3 not None ):
        """See if a grid search is to be performed for a set of coordinates."""
        cdef cUseGrid
        cUseGrid = PairListGenerator_DetermineMethod ( self.cObject, coordinates3.cObject, NULL )
        return ( cUseGrid == CTrue )

    @classmethod
    def WithOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    @property
    def cellSize ( self ): return self.cObject.cellSize

    @property
    def cutOff   ( self ): return self.cObject.cutOff

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def CrossPairList_FromDoubleCoordinates3 ( Coordinates3  coordinates31     ,
                                           Coordinates3  coordinates32     ,
                                           CrossPairList exclusions = None ,
                                           RealArray1D   radii1     = None , 
                                           RealArray1D   radii2     = None ,
                                           CReal         safety     = 0.0  ):
    """Create a cross-pairlist from two sets of coordinates."""
    cdef CrossPairList      pairList
    cdef CPairList         *cExclusions
    cdef CRealArray1D      *cRadii1
    cdef CRealArray1D      *cRadii2
    cdef PairListGenerator  self
    pairList = CrossPairList.Raw ( )
    if exclusions is None: cExclusions = NULL
    else:                  cExclusions = exclusions.cObject
    if radii1     is None: cRadii1     = NULL
    else:                  cRadii1     = radii1.cObject
    if radii2     is None: cRadii2     = NULL
    else:                  cRadii2     = radii2.cObject
    self = PairListGenerator.WithOptions ( cutOff = safety )
    pairList.cObject = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( self.cObject           ,
                                                                               coordinates31.cObject  ,
                                                                               coordinates32.cObject  ,
                                                                               cRadii1                ,
                                                                               cRadii2                ,
                                                                               NULL, NULL, NULL, NULL ,
                                                                               cExclusions            ,
                                                                               NULL, NULL, NULL       )
    pairList.isOwner = True
    return pairList

def CrossPairList_FromSingleCoordinates3 ( Coordinates3  coordinates3      ,
                                           CrossPairList exclusions = None ,
                                           RealArray1D   radii      = None ,
                                           Selection     selection1 = None ,
                                           Selection     selection2 = None ,
                                           CReal         safety     = 0.0  ):
    """Create a cross-pairlist from a single set of coordinates."""
    cdef CrossPairList  pairList
    cdef CPairList         *cExclusions
    cdef CSelection        *cSelection1
    cdef CSelection        *cSelection2
    cdef CRealArray1D      *cRadii
    cdef PairListGenerator  self
    pairList = CrossPairList.Raw ( )
    if exclusions is None: cExclusions = NULL
    else:                  cExclusions = exclusions.cObject
    if radii      is None: cRadii      = NULL
    else:                  cRadii      = radii.cObject
    if selection1 is None: cSelection1 = NULL
    else:                  cSelection1 = selection1.cObject
    if selection2 is None: cSelection2 = NULL
    else:                  cSelection2 = selection2.cObject
    self = PairListGenerator.WithOptions ( cutOff = safety )
    pairList.cObject = PairListGenerator_CrossPairListFromSingleCoordinates3 ( self.cObject            ,
                                                                               coordinates3.cObject    ,
                                                                               cRadii                  ,
                                                                               cSelection1             ,
                                                                               cSelection2             ,
                                                                               NULL                    ,
                                                                               cExclusions             ,
                                                                               CTrue, NULL, NULL, NULL )
    pairList.isOwner = True
    return pairList

def SelfPairList_FromCoordinates3 ( Coordinates3 coordinates3      ,
                                    SelfPairList exclusions = None ,
                                    RealArray1D  radii      = None ,
                                    CReal        safety     = 0.0  ):
    """Create a self-pairlist from a set of coordinates."""
    cdef CPairList         *cExclusions
    cdef CRealArray1D      *cRadii
    cdef PairListGenerator  self
    cdef SelfPairList       pairList
    pairList = SelfPairList.Raw ( )
    if exclusions is None: cExclusions = NULL
    else:                  cExclusions = exclusions.cObject
    if radii      is None: cRadii      = NULL
    else:                  cRadii      = radii.cObject
    self = PairListGenerator.WithOptions ( cutOff = safety )
    pairList.cObject = PairListGenerator_SelfPairListFromCoordinates3 ( self.cObject            ,
                                                                        coordinates3.cObject    , 
                                                                        cRadii                  ,
                                                                        NULL, NULL              ,
                                                                        cExclusions             ,
                                                                        NULL, NULL, NULL        )
    pairList.isOwner = True
    return pairList

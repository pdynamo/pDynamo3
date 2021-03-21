from pCore.CPrimitiveTypes                      cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.PairList                             cimport CPairList, CrossPairList, PairList, SelfPairList
from pCore.Selection                            cimport CSelection, Selection
from pCore.Status                               cimport CStatus, CStatus_OK
from pScientific.Arrays.RealArray1D             cimport CRealArray1D, RealArray1D
from pScientific.Arrays.RealArray2D             cimport CRealArray2D, RealArray2D
from pScientific.Geometry3.Coordinates3         cimport Coordinates3
from pScientific.Geometry3.RegularGrid          cimport CRegularGrid, RegularGrid
from pScientific.Geometry3.RegularGridOccupancy cimport CRegularGridOccupancy, RegularGridOccupancy

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairListGenerator.h":

    ctypedef struct CPairListGenerator "PairListGenerator":
        CBoolean  sortIndices
        CBoolean  useGridByCell
        CInteger  minimumCellExtent
        CInteger  minimumPoints
        CReal     cellSize
        CReal     cutOff
        CReal     cutOffCellSizeFactor
        CReal     minimumCellSize
        CReal     minimumExtentFactor

    cdef CPairListGenerator *PairListGenerator_Allocate                            ( )
    cdef CPairListGenerator *PairListGenerator_Clone                               ( CPairListGenerator    *self          )
    cdef CPairList          *PairListGenerator_CrossPairListFromDoubleCoordinates3 ( CPairListGenerator    *self          ,  
                                                                                     CRealArray2D          *coordinates31 ,  
                                                                                     CRealArray2D          *coordinates32 ,  
                                                                                     CRealArray1D          *radii1        ,  
                                                                                     CRealArray1D          *radii2        ,  
                                                                                     CSelection            *andSelection1 ,  
                                                                                     CSelection            *andSelection2 ,  
                                                                                     CSelection            *orSelection1  ,  
                                                                                     CSelection            *orSelection2  ,  
                                                                                     CPairList             *exclusions    ,  
                                                                                     CRegularGrid          *grid1         ,  
                                                                                     CRegularGridOccupancy *occupancy1    ,  
                                                                                     CStatus               *status        ) 
    cdef CPairList          *PairListGenerator_CrossPairListFromSingleCoordinates3 ( CPairListGenerator    *self          ,  
                                                                                     CRealArray2D          *coordinates3  ,  
                                                                                     CRealArray1D          *radii         ,  
                                                                                     CSelection            *andSelection1 ,  
                                                                                     CSelection            *andSelection2 ,  
                                                                                     CSelection            *orSelection   ,  
                                                                                     CPairList             *exclusions    ,  
                                                                                     CBoolean               excludeSelf   , 
                                                                                     CRegularGrid          *grid1         ,  
                                                                                     CRegularGridOccupancy *occupancy1    ,  
                                                                                     CStatus               *status        ) 
    cdef void                PairListGenerator_Deallocate                          ( CPairListGenerator   **self          )
    cdef CBoolean            PairListGenerator_DetermineMethod                     ( CPairListGenerator    *self          ,
                                                                                     CRealArray2D          *coordinates3  ,
                                                                                     CSelection            *andSelection  )
    cdef CPairList          *PairListGenerator_SelfPairListFromCoordinates3        ( CPairListGenerator    *self          ,  
                                                                                     CRealArray2D          *coordinates3  ,  
                                                                                     CRealArray1D          *radii         ,  
                                                                                     CSelection            *andSelection  ,  
                                                                                     CSelection            *orSelection   ,  
                                                                                     CPairList             *exclusions    ,  
                                                                                     CRegularGrid          *grid          ,  
                                                                                     CRegularGridOccupancy *occupancy     ,  
                                                                                     CStatus               *status        ) 

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairListGenerator ( object ):

    cdef CPairListGenerator *cObject
    cdef public object       isOwner

"""Coordinate arrays in 3-D."""

from .Geometry3Error import Geometry3Error

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Coordinates3 ( RealArray2D ):

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        clone = super ( Coordinates3, self ).__copy__ ( )
        if ( clone is not None ) and ( self.numberUndefined > 0 ): clone.undefined = set ( self.undefined )
        return clone

    def __deepcopy__ ( self, memo ):
        """Copying."""
        clone = super ( Coordinates3, self ).__deepcopy__ ( memo )
        if ( clone is not None ) and ( self.numberUndefined > 0 ): clone.undefined = set ( self.undefined )
        return clone

    def __getstate__ ( self ):
        """Return the state."""
        state = super ( Coordinates3, self ).__getstate__ ( )
        if self.numberUndefined > 0:
            indices = list ( self.undefined )
            indices.sort ( )
            state["undefined"] = indices
        return state

    def __setstate__ ( self, state ):
        """Set the state."""
        super ( Coordinates3, self ).__setstate__ ( state )
        if "undefined" in state: self.undefined = set ( state["undefined"] )

    def _Allocate ( self, rows, columns ):
        """Allocation."""
        if columns != 3: raise Geometry3Error ( "Invalid Coordinate3 shape." )
        super ( Coordinates3, self )._Allocate ( rows, columns )

    def _GetRawArray1D ( self ):
        """Get a raw 1-D array."""
        if ( MultiSlice_GetExtent ( self.cMultiSlice, 0, NULL ) == 1 ) and \
           ( MultiSlice_GetExtent ( self.cMultiSlice, 1, NULL ) == 3 ): return Vector3.RawWithCObject     ( )
        else:                                                           return RealArray1D.RawWithCObject ( )

    def _GetRawArray2D ( self ):
        """Get a raw 2-D array."""
        if MultiSlice_GetExtent ( self.cMultiSlice, 1, NULL ) == 3: return Coordinates3.RawWithCObject ( )
        else:                                                       return RealArray2D.RawWithCObject  ( )

    def _Initialize ( self ):
        """Initialization."""
        super ( Coordinates3, self )._Initialize ( )
        self.undefined = None

    def AddScaledRow ( self, CInteger i, CReal alpha, Vector3 v ):
        """Increment a row by a scaled vector."""
        self[i,0] = self[i,0] + alpha * v[0]
        self[i,1] = self[i,1] + alpha * v[1]
        self[i,2] = self[i,2] + alpha * v[2]

    def Angle ( self, CInteger i, CInteger j, CInteger k ):
        """Calculate the angle between three points."""
        return Coordinates3_Angle ( self.cObject, i, j, k )

    def BuildPointFromDistance ( self, CInteger i, CInteger j, CReal r, Vector3 direction ):
        """Build point |i| from point |j|, a distance |r| and a direction |direction|."""
        cdef CStatus status
        status = Coordinates3_BuildPointFromDistance ( self.cObject, i, j, r, direction.cObject )
        if status == CStatus_OK:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistanceAngle ( self, CInteger i, CInteger j, CInteger k, CReal r, CReal theta, Vector3 direction ):
        """Build point |i| from points |j| and |k|, a distance |r|, an angle |theta| and a direction |direction|."""
        cdef CStatus status
        status = Coordinates3_BuildPointFromDistanceAngle ( self.cObject, i, j, k, r, theta, direction.cObject )
        if status == CStatus_OK:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistanceAngleDihedral ( self, CInteger i, CInteger j, CInteger k, CInteger l, CReal r, CReal theta, CReal phi ):
        """Build point |i| from points |j|, |k| and |l|, a distance |r|, an angle |theta| and a dihedral |phi|."""
        cdef CStatus status
        status = Coordinates3_BuildPointFromDistanceAngleDihedral ( self.cObject, i, j, k, l, r, theta, phi )
        if status == CStatus_OK:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistancePlaneAngle ( self, CInteger i, CInteger j, CInteger k, CInteger l, CReal r, CReal planeAngle ):
        """Build point |i| from points |j|, |k| and |l|, a distance |r| and a plane angle |planeAngle|."""
        cdef CStatus status
        status = Coordinates3_BuildPointFromDistancePlaneAngle ( self.cObject, i, j, k, l, r, planeAngle )
        if status == CStatus_OK:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def BuildPointFromDistanceTetrahedralTripod ( self, CInteger i, CInteger j, CInteger k, CInteger l, CInteger m, CReal r ):
        """Build point |i| from points |j|, |k|, |l| and |m| and a distance |r| using a tetrahedral tripod."""
        cdef CStatus status
        status = Coordinates3_BuildPointFromDistanceTetrahedralTripod ( self.cObject, i, j, k, l, m, r )
        if status == CStatus_OK:
            self.FlagCoordinateAsDefined ( i )
            return True
        else:
            return False

    def Center ( self, Selection selection = None, RealArray1D weights = None ):
        """Determine the center of the matrix."""
        cdef Vector3       center
        cdef CRealArray1D *cWeights   = NULL
        cdef CSelection   *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        center = Vector3.Null ( )
        Coordinates3_Center ( self.cObject, cSelection, cWeights, &(center.cObject) )
        return center

    def Dihedral ( self, CInteger i, CInteger j, CInteger k, CInteger l ):
        """Calculate the dihedral between four points."""
        return Coordinates3_Dihedral ( self.cObject, i, j, k, l )

    def Displacement ( self, CInteger i, CInteger j, dr = None ):
        """Displacement between two points i and j."""
        cdef Vector3 vij
        if dr is None: vij = Vector3.Uninitialized ( )
        else:          vij = dr
        vij[0] = self[i,0] - self[j,0]
        vij[1] = self[i,1] - self[j,1]
        vij[2] = self[i,2] - self[j,2]
        return vij

    def Distance ( self, CInteger i, CInteger j ):
        """Calculate the distance between two points."""
        return Coordinates3_Distance ( self.cObject, i, j )

    def EnclosingOrthorhombicBox ( self, Selection selection = None, RealArray1D radii = None ):
        """Find the enclosing box around the matrix."""
        cdef Vector3       extents, origin
        cdef CRealArray1D *cRadii     = NULL
        cdef CSelection   *cSelection = NULL
        if radii     is not None: cRadii     = radii.cObject
        if selection is not None: cSelection = selection.cObject
        extents = Vector3.Uninitialized ( )
        origin  = Vector3.Uninitialized ( )
        Coordinates3_EnclosingOrthorhombicBox ( self.cObject, cSelection, cRadii, origin.cObject, extents.cObject )
        return ( origin, extents )

    def FlagCoordinateAsDefined ( self, i ):
        """Flag a coordinate as being defined."""
        if self.undefined is not None: self.undefined.discard ( i )

    def FlagCoordinateAsUndefined ( self, i ):
        """Flag a coordinate as being undefined."""
        if self.undefined is None: self.undefined = set ( )
        self.undefined.add ( i )

    @classmethod
    def FromGrid ( selfClass, RegularGrid grid, Selection selection = None ):
        """Constructor from grid."""
        cdef Coordinates3  self
        cdef CInteger      n
        cdef CSelection   *cSelection = NULL
        cdef CStatus       cStatus    = CStatus_OK
        n = grid.size
        if selection is not None:
            cSelection = selection.cObject
            n          = len ( selection )
        self = selfClass.WithExtent ( n )
        Coordinates3_FromRegularGrid ( self.cObject, grid.cObject, cSelection, &cStatus )
        if cStatus != CStatus_OK: raise Geometry3Error ( "Error making coordinates from grid." )
        return self

    def Gather ( self, Coordinates3 other, Selection selection = None ):
        """Gather items from a sparse other, indexed by selection, to a compact self."""
        if selection is None: Coordinates3_CopyTo ( other.cObject, self.cObject , NULL              )
        else:                 Coordinates3_Gather ( self.cObject , other.cObject, selection.cObject )

    def GatherAdd ( self, CReal alpha, Coordinates3 other, Selection selection = None ):
        """Gather items, with adding and scaling, from a sparse other, indexed by selection, to a compact self."""
        if selection is None: Coordinates3_Add       ( self.cObject, alpha, other.cObject, NULL              )
        else:                 Coordinates3_GatherAdd ( self.cObject, alpha, other.cObject, selection.cObject )

    def IdentifyOccupiedGridPoints ( self, RegularGrid grid, RealArray1D radii ):
        """Identify occupied grid points."""
        cdef Selection occupied
        occupied         = Selection.Raw ( )
        occupied.isOwner = True
        Coordinates3_IdentifyOccupiedGridPoints ( self.cObject, grid.cObject, radii.cObject, CFalse, &occupied.cObject )
        return occupied

    def IdentifyUnoccupiedGridPoints ( self, RegularGrid grid, RealArray1D radii ):
        """Identify unoccupied grid points."""
        occupied = self.IdentifyOccupiedGridPoints ( grid, radii )
        return occupied.Complement ( upperBound = grid.size )

    def InertiaMatrix ( self, Selection selection = None, RealArray1D weights = None  ):
        """Determine the inertia matrix."""
        cdef SymmetricMatrix  inertia
        cdef CRealArray1D    *cWeights   = NULL
        cdef CSelection      *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        inertia = SymmetricMatrix.WithExtent ( 3 )
        Coordinates3_InertiaMatrix ( self.cObject, cSelection, cWeights, inertia.cObject )
        return inertia

    def MakeGridAndOccupancy ( self, cellSize ):
        """Make a grid and its occupancy from the coordinates."""
        cdef RegularGrid            grid       = None
        cdef RegularGridOccupancy   occupancy  = None
        cdef CReal                  cCellSize
        cdef CRegularGrid          *cGrid      = NULL
        cdef CRegularGridOccupancy *cOccupancy = NULL
        cdef CStatus                status     = CStatus_OK
        cCellSize = float ( cellSize )
        Coordinates3_MakeGridAndOccupancy ( self.cObject, NULL, cCellSize, (&cGrid), (&cOccupancy), &status )
        if cGrid != NULL:
            grid         = RegularGrid.Raw ( )
            grid.cObject = cGrid
            grid.isOwner = True
        if cOccupancy != NULL:
            occupancy         = RegularGridOccupancy.Raw ( )
            occupancy.cObject = cOccupancy
            occupancy.isOwner = True
        if status != CStatus_OK: raise Geometry3Error ( "Error creating coordinates3 grid and occupancy." )
        return ( grid, occupancy )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging."""
        merged  = None
        lengths = [ item.rows for item in items if item is not None ]
        if len ( lengths ) == len ( items ):
            merged = selfClass.WithExtent ( sum ( lengths ) )
            m      = 0
            for ( n, item ) in zip ( lengths, items ):
                item.CopyTo ( merged[m:m+n] )
                m += n
        return merged
        # . Slow version!
        # . None not allowed.
#       cdef Coordinates3 item, new
#       cdef CInteger     i, j, n, n0
#       cdef CReal        value
#       isOK = True
#       new  = None
#       n0   = 0
#       for item in items:
#           if item is None:
#               isOK = False
#               break
#           else:
#               n0 += item.rows
#       if isOK:
#           new = selfClass.WithExtent ( n0 )
#           n0  = 0
#           for item in items:
#               n = item.rows
#               for i from 0 <= i < n:
#                   for j from 0 <= j < item.columns:
#                       value = Coordinates3_GetItem ( item.cObject, i, j, NULL )
#                       Coordinates3_SetItem ( new.cObject, i + n0, j, value, NULL )
#               n0 += n
#       return new

    def MomentsOfInertia ( self, Selection selection = None, RealArray1D weights = None  ):
        """Calculate the moments of inertia."""
        cdef Matrix33      axes
        cdef Vector3       moments
        cdef CRealArray1D *cWeights   = NULL
        cdef CSelection   *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        axes    = Matrix33.Null ( )
        moments = Vector3.Null  ( )
        Coordinates3_MomentsOfInertia ( self.cObject, cSelection, cWeights, moments.cObject, axes.cObject )
        return ( moments, axes )

    def Prune ( self, Selection selection not None, information = {} ):
        """Pruning."""
        pruned = self.__class__.WithExtent ( len ( selection ) )
        pruned.Gather ( self, selection )
        return pruned

    def RadiusOfGyration ( self, Selection selection = None, RealArray1D weights = None ):
        """Determine the radius of gyration."""
        cdef CReal         value
        cdef CRealArray1D *cWeights   = NULL
        cdef CSelection   *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        value = Coordinates3_RadiusOfGyration ( self.cObject, cSelection, cWeights )
        return value

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RootMeanSquareDeviation ( self, Coordinates3 other, Selection selection = None, RealArray1D weights = None ):
        """Determine the root mean square deviation between two matrices."""
        cdef CReal         rms
        cdef CRealArray1D *cWeights   = NULL
        cdef CSelection   *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        rms = Coordinates3_RootMeanSquareDeviation ( self.cObject, other.cObject, cSelection, cWeights )
        return rms

    def Rotate ( self, Matrix33 rotation, Selection selection = None ):
        """Rotation."""
        cdef CSelection *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        Coordinates3_Rotate ( self.cObject, rotation.cObject, cSelection )

    def RotationTranslationVectors ( self, qRotation, qTranslation, CInteger dimension = 0, RealArray1D weights = None ):
        """Construct a set of rotation and translation vectors."""
        cdef RealArray1D   reference
        cdef RealArray1D   scalars
        cdef RealArray2D   vectors
        cdef RealArray2D   vectorsN
        cdef CBoolean      qR, qT
        cdef CInteger      n
        cdef CRealArray1D *cWeights = NULL
        cdef CStatus       cStatus  = CStatus_OK
        n = 0
        if weights is not None: cWeights = weights.cObject
        if qRotation:
            qR = CTrue
            n += 3
        else:
            qR = CFalse
        if qTranslation:
            qT = CTrue
            n += 3
        else:
            qT = CFalse
        vectors   = RealArray2D.WithExtents ( max ( self.size, dimension ), n ) # . Dimension is required just in case there are more variables than coordinates.
        n         = Coordinates3_RotationTranslationVectors ( self.cObject, cWeights, qR, qR, qR, qT, qT, qT, vectors.cObject, &cStatus )
        reference = RealArray1D.WithExtent ( vectors.rows ) # . Flattened view possible?
        scalars   = RealArray1D.WithExtent ( n )
        vectorsN  = vectors[:,0:n]
        reference.Set ( 0.0 )
        self.iterator.CopyTo ( reference )
        RealArray2D_VectorMultiply ( CTrue, 1.0, vectorsN.cObject, reference.cObject, 0.0, scalars.cObject, &cStatus )
        if cStatus != CStatus_OK: raise Geometry3Error ( "Error making rotation/translation vectors." )
        return ( vectorsN, scalars )

    def Scatter ( self, Coordinates3 other, Selection selection = None ):
        """Scatter items from a compact self to a sparse other, indexed by selection."""
        if selection is None: Coordinates3_CopyTo  ( self.cObject, other.cObject, NULL              )
        else:                 Coordinates3_Scatter ( self.cObject, other.cObject, selection.cObject )

    def ScatterAdd ( self, CReal alpha, Coordinates3 other, Selection selection = None ):
        """Scatter items, with adding and scaling, from a compact self to a sparse other, indexed by selection."""
        if selection is None: Coordinates3_Add        ( other.cObject, alpha, self.cObject , NULL              )
        else:                 Coordinates3_ScatterAdd ( self.cObject , alpha, other.cObject, selection.cObject )

    def Superimpose ( self, Coordinates3 other, Selection selection = None, RealArray1D weights = None, returnTransformation = False ):
        """Superimpose the matrix onto a other matrix."""
        cdef Matrix33      rotation
        cdef Vector3       translation
        cdef CRealArray1D *cWeights     = NULL
        cdef CRealArray2D *cRotation    = NULL
        cdef CSelection   *cSelection   = NULL
        cdef CRealArray1D *cTranslation = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        if returnTransformation:
            rotation     = Matrix33.Null ( )
            translation  = Vector3.Null  ( )
            cRotation    = rotation.cObject
            cTranslation = translation.cObject
        else:
            cRotation    = NULL
            cTranslation = NULL
        Coordinates3_Superimpose ( self.cObject, other.cObject, cSelection, cWeights, cRotation, cTranslation )
        if returnTransformation: return ( rotation, translation )

    def ToPrincipalAxes ( self, Selection selection = None, RealArray1D weights = None ):
        """Do a principal axis transformation."""
        cdef CRealArray1D *cWeights   = NULL
        cdef CSelection   *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        Coordinates3_ToPrincipalAxes ( self.cObject, cSelection, cWeights )

    def Transform ( self, Transformation3 transformation3, Selection selection = None ):
        """Transformation."""
        cdef CSelection *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        Coordinates3_Transform ( self.cObject, transformation3.cObject, cSelection )

    def Translate ( self, Vector3 translation, Selection selection = None ):
        """Translation."""
        cdef CSelection *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        Coordinates3_Translate ( self.cObject, translation.cObject, cSelection )

    def TranslateToCenter ( self, Selection selection = None, RealArray1D weights = None ):
        """Translate the matrix to its center."""
        cdef CRealArray1D *cWeights   = NULL
        cdef CSelection   *cSelection = NULL
        if selection is not None: cSelection = selection.cObject
        if weights   is not None: cWeights   = weights.cObject
        Coordinates3_TranslateToCenter ( self.cObject, cSelection, cWeights )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent, 3 )

    # . Properties.
    @property
    def numberUndefined ( self ):
        if self.undefined is None: return 0
        else:                      return len ( self.undefined )

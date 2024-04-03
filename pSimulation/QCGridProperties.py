"""Classes and functions for calculating the properties of QC systems on a grid."""

import math, os.path

from pBabel                import ExportObjects      , \
                                  ExportSystem
from pCore                 import logFile            , \
                                  SummarizableObject
from pMolecule.QCModel     import QCModelError
from pScientific           import Units
from pScientific.Arrays    import Array              , \
                                  RealArrayND        , \
                                  Reshape
from pScientific.Geometry3 import Coordinates3       , \
                                  RegularGrid
from pScientific.Surfaces  import MarchingCubes_Isosurface3D

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Grid property defaults.
# . Grid generation.
_DefaultGridSpacing   = 0.2 # . Atomic units.
_DefaultRadiusFactor  = 3.0
# . Property tags.
_DefaultDensityTag    = "Grid Density"
_DefaultOrbitalTag    = "Grid Orbitals"
_DefaultPotentialTag  = "Grid Potential"
_DefaultIsosurfaceTag = "Isosurface"

#===================================================================================================================================
# . QC grid property classes.
#===================================================================================================================================
class QCBaseProperty ( SummarizableObject ):
    """Class for a QC grid property."""

    _attributable    = dict ( SummarizableObject._attributable )
    _exportAttribute = None
    _exportFormat    = None
    _exportFunction  = ExportObjects

    def Export ( self, path, format = format, log = logFile ):
        """Export the property."""
        # . Set the extension if the input path doesn't have one and no format is provided.
        if ( format is None ) and ( os.path.splitext ( path )[-1] == "" ):
            path = "{:s}.{:s}".format ( path, self.__class__._exportFormat )
        options = { "format" : format, "log" : log }
        options.update ( self.ExportOptions ( ) )
        self.__class__._exportFunction ( path, getattr ( self, self.__class__._exportAttribute ), **options )

    def ExportOptions ( self ):
        """Export options."""
        return {}

class QCGridProperty ( QCBaseProperty ):
    """A QC grid property."""

    _attributable    = dict ( QCBaseProperty._attributable )
    _exportAttribute = "system"
    _exportFormat    = "cube"
    _exportFunction  = ExportSystem
    _attributable.update ( { "grid"           : None ,
                             "gridPoints"     : None ,
                             "gridValues"     : None ,
                             "label"          : None ,
                             "orbitalIndices" : None ,
                             "system"         : None ,
                             "spinType"       : None } )

    def ExportOptions ( self ):
        """Export options."""
        return { "grid"           : self.grid           ,
                 "gridValues"     : self.gridValues     ,
                 "orbitalIndices" : self.orbitalIndices ,
                 "title"          : self.label          }

    def Isosurface ( self, isovalue, orbitalIndex = None ):
        """Generate an isosurface."""
        shape   = self.grid.shape
        sliceND = False
        if ( self.orbitalIndices is not None ) and ( len ( self.orbitalIndices ) > 1 ):
            shape   = list ( shape ) + [ len ( self.orbitalIndices ) ]
            sliceND = True
            try: index = self.orbitalIndices.index ( orbitalIndex )
            except: raise ValueError ( "Orbital index specified for isosurface generation not found." )
        dataND = Reshape ( self.gridValues, shape, resultClass = RealArrayND )
        if sliceND: dataND = dataND[:,:,:,index]
        surface = MarchingCubes_Isosurface3D ( self.grid, dataND, isovalue )
        return surface

class QCIsosurface ( QCBaseProperty ):
    """A QC isosurface."""

    _attributable    = dict ( QCBaseProperty._attributable )
    _exportAttribute = "isosurface"
    _exportFormat    = "stl"
    _attributable.update (  { "isosurface" : None ,
                              "isovalue"   : None } )

#===================================================================================================================================
# . QC grid property generator class.
#===================================================================================================================================
class QCGridPropertyGenerator ( SummarizableObject ):
    """Class for calculating properties of QC systems either on or derived from a grid."""

    _attributable = dict ( SummarizableObject._attributable )
    _attributable.update (  { "atomicNumbers"  : None ,
                              "grid"           : None ,
                              "gridPoints"     : None ,
                              "properties"     : dict ,
                              "qcCoordinates3" : None ,
                              "system"         : None } )

    def _CheckOptions ( self ):
        """Constructor."""
        super ( QCGridPropertyGenerator, self )._CheckOptions ( )
        try:
            self.atomicNumbers  = self.system.qcState.atomicNumbers
            self.qcCoordinates3 = self.system.scratch.qcCoordinates3AU
        except: raise TypeError ( "Invalid QC system." )

    def ClearGrid ( self ):
        """Clear the grid."""
        self.grid       = None
        self.gridPoints = None

    def DefineGrid ( self, gridSpacing = _DefaultGridSpacing, radiusFactor = _DefaultRadiusFactor ):
        """Define the grid."""
        # . Get radii.
        allRadii = Array.FromIterable ( [ atom.vdwRadius for atom in self.system.atoms ] )
        if len ( self.atomicNumbers ) == len ( allRadii ):
            radii = allRadii
        else:
            radii = Array.WithExtent ( len ( self.atomicNumbers ) )
            for ( i, s ) in enumerate ( self.qcAtoms ): radii[i] = allRadii[s]
        radii.Scale ( radiusFactor / Units.Length_Bohrs_To_Angstroms )
        # . Get grid axis data so that the grid encloses fully the QC system.
        gridAxes = []
        ( origin, extents ) = self.qcCoordinates3.EnclosingOrthorhombicBox ( radii = radii )
        for i in range ( len ( extents ) ):
            n = int ( math.ceil ( extents[i] / gridSpacing ) )
            l = float ( n ) * gridSpacing
            origin [i] -= 0.5 * ( l - extents[i] ) # . Readjust the origin and extents accordingly.
            extents[i]  = l
            gridAxes.append ( { "bins" : n, "binSize" : gridSpacing, "lower" : origin[i] } )
        # . Construct grid.
        self.grid       = RegularGrid.FromDimensionData ( gridAxes )
        self.gridPoints = Coordinates3.FromGrid ( self.grid )

    def ExportProperty ( self, path, tag, format = None, log = logFile ):
        """Export a property."""
        property = self.GetProperty ( tag )
        property.Export ( path, format = format, log = log )

    @classmethod
    def FromSystem ( selfClass, system ):
        """Constructor given system."""
        return selfClass.WithOptions ( system = system )

    def GetProperty ( self, tag ):
        """Return a property."""
        if tag in self.properties:
            return self.properties[tag]
        else:
            raise QCModelError ( "Property with tag \"{:s}\" not found.".format ( tag ) )

    def GridDensity ( self, spinType = None, tag = _DefaultDensityTag ):
        """The electronic density on the grid."""
        if self.gridPoints is None:
            raise QCModelError ( "Grid points are not defined." )
        else:
            qcModel = self.system.qcModel
            data    = qcModel.GridPointDensities ( self.system, self.gridPoints, spinType = spinType )
            self.properties[tag] = QCGridProperty.WithOptions ( grid       = self.grid       ,
                                                                gridPoints = self.gridPoints ,
                                                                gridValues = data            ,
                                                                label      = tag             ,
                                                                spinType   = spinType        ,
                                                                system     = self.system     )

    def GridOrbitals ( self, orbitalIndices, spinType = None, tag = _DefaultOrbitalTag ):
        """Electronic orbitals on the grid."""
        if self.gridPoints is None:
            raise QCModelError ( "Grid points are not defined." )
        else:
            qcModel = self.system.qcModel
            data    = qcModel.GridPointOrbitals ( self.system, self.gridPoints, orbitalIndices, spinType = spinType )
            self.properties[tag] = QCGridProperty.WithOptions  ( grid           = self.grid       ,
                                                                 gridPoints     = self.gridPoints ,
                                                                 gridValues     = data            ,
                                                                 label          = tag             ,
                                                                 orbitalIndices = orbitalIndices  ,
                                                                 spinType       = spinType        ,
                                                                 system         = self.system     )

    def GridPotential ( self, spinType = None, tag = _DefaultPotentialTag ):
        """The electrostatic potential on the grid."""
        if self.gridPoints is None:
            raise QCModelError ( "Grid points are not defined." )
        else:
            qcModel = self.system.qcModel
            data    = qcModel.GridPointPotentials ( self.system, self.gridPoints, spinType = spinType )
            self.properties[tag] = QCGridProperty.WithOptions  ( grid       = self.grid       ,
                                                                 gridPoints = self.gridPoints ,
                                                                 gridValues = data            ,
                                                                 label      = tag             ,
                                                                 spinType   = spinType        ,
                                                                 system     = self.system     )

    def Isosurface ( self, dataTag, isovalue, orbitalIndex = None, tag = _DefaultIsosurfaceTag ):
        """Generate an isosurface."""
        gridProperty = self.GetProperty ( dataTag )
        if hasattr ( gridProperty, "Isosurface" ):
            surface       = gridProperty.Isosurface ( isovalue, orbitalIndex = orbitalIndex )
            surface.label = "{:s} with isovalue {:.1f}".format ( tag.capitalize ( ), isovalue )
            self.properties[tag] = QCIsosurface.WithOptions ( isosurface = surface  ,                                               isovalue   = isovalue )
        else:
            raise QCModelError ( "Invalid QC property for isosurface generation." )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

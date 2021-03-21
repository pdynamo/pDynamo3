"""Basic class for periodic boundary conditions."""

from   pCore              import logFile                  , \
                                 LogFileActive            , \
                                 RawObjectConstructor     , \
                                 SummarizableObject
from  .CrystalSystem      import CrystalSystem            , \
                                 CrystalSystem_FromLabel
from  .SpaceGroup         import SpaceGroup_FromYAML
from  .SymmetryError      import SymmetryError
from  .SymmetryParameters import SymmetryParameters
from ..Geometry3          import Transformation3          , \
                                 Transformation3Container

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PeriodicBoundaryConditions ( SummarizableObject ):
    """Periodic boundary conditions."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Periodic Boundary Conditions"
    _attributable.update ( { "crystalSystem"   : None ,
                             "transformations" : None } )

    def CheckSymmetryParameters ( self, symmetryParameters ):
        """Check the validity of a set of symmetry parameters."""
        self.crystalSystem.AssignCrystalParameters ( symmetryParameters.a     ,
                                                     symmetryParameters.b     ,
                                                     symmetryParameters.c     ,
                                                     symmetryParameters.alpha ,
                                                     symmetryParameters.beta  ,
                                                     symmetryParameters.gamma )

    @classmethod
    def FromSpaceGroup ( selfClass, hermannMauguin, inPath = None ):
        """Constructor given a space group name."""
        spaceGroup      = SpaceGroup_FromYAML ( hermannMauguin, inPath = inPath )
        crystalSystem   = CrystalSystem_FromLabel ( spaceGroup.crystalSystem )
        tList           = [ Transformation3.FromSymmetryOperationString ( o ) for o in spaceGroup.transformations ]
        transformations = Transformation3Container.WithTransformations ( tList )
        return selfClass.WithCrystalSystem ( crystalSystem, transformations = transformations )

    def MakeSymmetryParameters ( self, a = None, b = None, c = None, alpha = None, beta = None, gamma = None ):
        """Make a set of symmetry parameters consistent with the crystal system."""
        ( a, b, c, alpha, beta, gamma ) = self.crystalSystem.AssignCrystalParameters ( a, b, c, alpha, beta, gamma )
        return SymmetryParameters.FromCrystalParameters ( a, b, c, alpha, beta, gamma )

    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( "Symmetry"      , True                     ) ,
                  ( "Crystal System", self.crystalSystem.label ) ]
        if self.transformations is not None: items.extend ( self.transformations.SummaryItems ( ) )
        return items

    @classmethod
    def WithCrystalSystem ( selfClass, crystalSystem, transformations = None ):
        """Constructor given a crystal system and an optional set of transformations."""
        if not ( isinstance ( crystalSystem, CrystalSystem ) and \
               ( ( transformations is None ) or \
               isinstance ( transformations, Transformation3Container ) ) ):
            raise SymmetryError ( "Invalid arguments to symmetry constructor." )
        if ( transformations is None ): transformations = Transformation3Container.Identity ( )
        return selfClass.WithOptions ( crystalSystem   = crystalSystem   ,
                                       transformations = transformations )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

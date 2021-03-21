"""Symmetry operations for point group identification."""

import math

from   pCore     import AttributableObject , \
                        Clone              , \
                        DataType
from ..Arrays    import Array
from ..Geometry3 import Matrix33           , \
                        Vector3

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class SymmetryElement ( AttributableObject ):
    """Base class for symmetry elements."""

    _attributable = dict ( AttributableObject._attributable )
    _label        = ""
    _attributable.update ( { "label"                : None ,
                             "mapping"              : None ,
                             "order"                :    1 ,
                             "selfMappings"         :    0 ,
                             "transformationMatrix" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( SymmetryElement, self )._CheckOptions ( )
        self.label = self.__class__._label
        self.MakeTransformationMatrix ( )

    def ApplyTo ( self, point ):
        """Apply the transformation to a point (in-place)."""
        self.transformationMatrix.ApplyTo ( point )

    def EstablishSymmetryRelatedPairs ( self, nodeGroups, coordinates3, tolerance ):
        """Establish the symmetry-related pairs for a transformation."""
        # . Initialization.
        numberNodes  = coordinates3.rows
        self.mapping = Array.WithExtent ( numberNodes, dataType = DataType.Integer )
        self.mapping.Set ( -1 )
        success      = True
        # . All nodes map to themselves.
        if self.order == 1:
            for i in range ( numberNodes ): self.mapping[i] = i
            self.selfMappings = numberNodes
        # . Build the pairs.
        else:
            xyz = Vector3.Null ( )
            for ( node, groups ) in nodeGroups.items ( ):
                for group in groups:
                    isUsed = set ( )
                    for iNode in group:
                        if self.mapping[iNode] < 0:
                            bestDistance = 2.0 * tolerance
                            bestJ        = -1
                            coordinates3[iNode].CopyTo ( xyz )
                            self.ApplyTo ( xyz )
                            for jNode in group:
                                if jNode not in isUsed:
                                    distance = math.sqrt ( ( xyz[0] - coordinates3[jNode,0] )**2 + \
                                                           ( xyz[1] - coordinates3[jNode,1] )**2 + \
                                                           ( xyz[2] - coordinates3[jNode,2] )**2 )
                                    if distance < bestDistance:
                                        bestJ        = jNode
                                        bestDistance = distance
                            if bestDistance > tolerance:
                                success = False
                                break
                            else:
                                self.mapping[iNode] = bestJ
                                isUsed.add ( bestJ )
            self.selfMappings = 0
            if success:
                for ( i, j ) in enumerate ( self.mapping ):
                    if ( i == j ): self.selfMappings += 1
        return success

    def IsEquivalent ( self, other ):
        """Check for equivalence of two transformations."""
        return True

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        self.transformationMatrix = Matrix33.Identity ( )

#===================================================================================================================================
class Identity ( SymmetryElement ):
    """Identity."""

    _label = "E"

    def ApplyTo ( self ): pass

#===================================================================================================================================
class Inversion ( SymmetryElement ):
    """Inversion."""

    _label = "i"

    def _CheckOptions ( self ):
        """Check options."""
        super ( Inversion, self )._CheckOptions ( )
        self.order = 2

    def ApplyTo ( self, point ):
        """Apply the transformation to a point (in-place)."""
        point.Scale ( -1.0 )

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( Inversion, self ).MakeTransformationMatrix ( )
        for i in range ( 3 ): self.transformationMatrix[i,i] = -1.0

#===================================================================================================================================
class ProperRotation ( SymmetryElement ):
    """Proper rotation."""

    _attributable = dict ( SymmetryElement._attributable )
    _label        = "C"
    _attributable.update ( { "angle" :  1.0  ,
                             "axis"  : None  } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( ProperRotation, self )._CheckOptions ( )
        if self.order == 1: tag = "inf"
        else:               tag = repr ( self.order )
        self.label += tag

    def IsEquivalent ( self, other, tolerance ):
        """Check for equivalence of two transformations."""
        self.axis.Normalize  ( )
        other.axis.Normalize ( )
        areSame = ( self.order == other.order ) and ( math.fabs ( math.fabs ( self.axis.Dot ( other.axis ) ) - 1.0 ) < tolerance )
        return areSame

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( ProperRotation, self ).MakeTransformationMatrix ( )
        if self.order > 1: self.angle = ( 2.0 * math.pi ) / float ( self.order )
        self.axis.Normalize ( )
        self.transformationMatrix.RotationAboutAxis ( self.angle, self.axis )

#===================================================================================================================================
class ImproperRotation ( ProperRotation ):
    """Improper rotation."""

    _label = "S"

    def _CheckOptions ( self ):
        """Check options."""
        super ( ImproperRotation, self )._CheckOptions ( )
        self.label += repr ( self.order )

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( ImproperRotation, self ).MakeTransformationMatrix ( )
        reflection = Matrix33.Null ( )
        reflection.Reflection ( self.axis )
        self.transformationMatrix.PreMultiplyBy ( reflection )

#===================================================================================================================================
class Reflection ( SymmetryElement ):
    """Reflection."""

    _attributable = dict ( SymmetryElement._attributable )
    _label        = "sigma"
    _attributable.update ( { "normal" : None } )

    def IsEquivalent ( self, other, tolerance ):
        """Check for equivalence of two transformations."""
        self.normal.Normalize  ( )
        other.normal.Normalize ( )
        areSame = ( self.order == other.order ) and ( math.fabs ( math.fabs ( self.normal.Dot ( other.normal ) ) - 1.0 ) < tolerance )
        return areSame

    def MakeTransformationMatrix ( self ):
        """Make the transformation matrix."""
        super ( Reflection, self ).MakeTransformationMatrix ( )
        self.normal.Normalize ( )
        self.transformationMatrix.Reflection ( self.normal )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

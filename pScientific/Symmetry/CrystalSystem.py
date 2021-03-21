"""Classes to represent crystal systems."""

from  pCore              import RawObjectConstructor
from .SymmetryError      import SymmetryError
from .SymmetryParameters import SymmetryParameters

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Base class.
class CrystalSystem:
    """Base class for crystal systems."""

    # . Class parameters.
    _alpha = 0.0
    _beta  = 0.0
    _gamma = 0.0
    _label = None

    def __getstate__ ( self ): return {}

    def __init__ ( self ):
        """Constructor."""
        self.label = self.__class__._label

    def __setstate__ ( self, state ): pass

    def __len__ ( self ):
        """Return the number of symmetry parameters needed to characterize the class."""
        return 0

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Form a complete set of crystal parameters."""
        return ( a, b, c, alpha, beta, gamma )

    def CreateSymmetryParameters ( self, a = None, b = None, c = None, alpha = None, beta = None, gamma = None ):
        """Create and fill a symmetry parameters object."""
        # . Get and check the parameter values.
        if alpha is None: alpha = self.__class__._alpha
        if beta  is None: beta  = self.__class__._beta
        if gamma is None: gamma = self.__class__._gamma
        ( a, b, c, alpha, beta, gamma ) = self.AssignCrystalParameters ( a, b, c, alpha, beta, gamma )
        self.VerifyCrystalParameterRanges ( a, b, c, alpha, beta, gamma )
        # . Construct the symmetry parameters.
        sp = SymmetryParameters ( )
        sp.SetCrystalParameters ( a, b, c, alpha, beta, gamma )
        return sp

    def IsOrthogonal ( self ):
        """Is this crystal class orthogonal?"""
        return False

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass ( )
        return self

    def ScaleVolume ( self, symmetryParameters, scale ):
        """Scale the box uniformly."""
        sp = symmetryParameters
        sp.SetCrystalParameters ( sp.a * scale, sp.b * scale, sp.c * scale, sp.alpha, sp.beta, sp.gamma )

    def VerifyCrystalParameterRanges ( self, a, b, c, alpha, beta, gamma ):
        """Basic checks on crystal parameter values."""
        # . Very basic checks.
        if ( a <= 0.0 ) or ( b <= 0.0 ) or ( c <= 0.0 ): raise SymmetryError ( "Unit cell lengths zero or negative."    )
        if ( alpha <= 0.0 ) or ( alpha >= 180.0 ) or \
           ( beta  <= 0.0 ) or ( beta  >= 180.0 ) or \
           ( gamma <= 0.0 ) or ( gamma >= 180.0 ):       raise SymmetryError ( "Unit cell angles not in range [0,180]." )

# . Derived classes.
class CrystalSystemCubic ( CrystalSystem ):
    """Cubic crystal class."""

    # . Class parameters.
    _alpha = 90.0
    _beta  = 90.0
    _gamma = 90.0
    _label = "Cubic"

    def __len__ ( self ): return 1

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        if b     is None: b     = a
        if c     is None: c     = a
        if alpha is None: alpha = self.__class__._alpha
        if beta  is None: beta  = self.__class__._beta
        if gamma is None: gamma = self.__class__._gamma
        QOK = ( a is not None ) and ( a == b ) and ( a == c ) and ( alpha == 90.0 ) and ( beta == 90.0 ) and ( gamma == 90.0 )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a = b = c ; alpha = beta = gamma = 90." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement] = symmetryParameterGradients.dEdA + symmetryParameterGradients.dEdB + symmetryParameterGradients.dEdC

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement] = symmetryParameters.a

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a = variables[vectorIncrement]
        symmetryParameters.SetCrystalParameters ( a, a, a, self.__class__._alpha, self.__class__._beta, self.__class__._gamma )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a }

    def GetVolumeDerivative ( self, symmetryParameters, symmetryParameterGradients ):
        """Get the symmetry parameter volume derivative."""
        return ( symmetryParameterGradients.dEdA + symmetryParameterGradients.dEdB + symmetryParameterGradients.dEdC ) / ( 3.0 * symmetryParameters.a**2 )

    def IsOrthogonal ( self ): return True

class CrystalSystemHexagonal ( CrystalSystem ):
    """Hexagonal crystal class."""

    # . Class parameters.
    _alpha =  90.0
    _beta  =  90.0
    _gamma = 120.0
    _label = "Hexagonal"

    def __len__ ( self ): return 2

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        if b     is None: b     = a
        if alpha is None: alpha = self.__class__._alpha
        if beta  is None: beta  = self.__class__._beta
        if gamma is None: gamma = self.__class__._gamma
        QOK = ( a is not None ) and ( a == b ) and ( c is not None ) and ( alpha == 90.0 ) and ( beta == 90.0 ) and ( gamma == 120.0 )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a = b != c ; alpha = beta = 90 ; gamma = 120." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement  ] = symmetryParameterGradients.dEdA + symmetryParameterGradients.dEdB
        gradients[vectorIncrement+1] = symmetryParameterGradients.dEdC

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement  ] = symmetryParameters.a
        variables[vectorIncrement+1] = symmetryParameters.c

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a = variables[vectorIncrement  ]
        c = variables[vectorIncrement+1]
        symmetryParameters.SetCrystalParameters ( a, a, c, self.__class__._alpha, self.__class__._beta, self.__class__._gamma )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a, "c" : symmetryParameters.c }

class CrystalSystemMonoclinic ( CrystalSystem ):
    """Monoclinic crystal class."""

    # . Class parameters.
    _alpha = 90.0
    _beta  =  0.0
    _gamma = 90.0
    _label = "Monoclinic"

    def __len__ ( self ): return 4

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        if alpha is None: alpha = self.__class__._alpha
        if gamma is None: gamma = self.__class__._gamma
        QOK = ( a is not None ) and ( b is not None ) and ( c is not None ) and ( alpha == 90.0 ) and ( gamma == 90.0 )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a != b != c ; alpha = gamma = 90 != beta." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement  ] = symmetryParameterGradients.dEdA
        gradients[vectorIncrement+1] = symmetryParameterGradients.dEdB
        gradients[vectorIncrement+2] = symmetryParameterGradients.dEdC
        gradients[vectorIncrement+3] = symmetryParameterGradients.dEdBeta

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement  ] = symmetryParameters.a
        variables[vectorIncrement+1] = symmetryParameters.b
        variables[vectorIncrement+2] = symmetryParameters.c
        variables[vectorIncrement+3] = symmetryParameters.beta

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a    = variables[vectorIncrement  ]
        b    = variables[vectorIncrement+1]
        c    = variables[vectorIncrement+2]
        beta = variables[vectorIncrement+3]
        symmetryParameters.SetCrystalParameters ( a, b, c, self.__class__._alpha, beta, self.__class__._gamma )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a, "b" : symmetryParameters.b, "c" : symmetryParameters.c, "beta" : symmetryParameters.beta }

class CrystalSystemOrthorhombic ( CrystalSystem ):
    """Orthorhombic crystal class."""

    # . Class parameters.
    _alpha = 90.0
    _beta  = 90.0
    _gamma = 90.0
    _label = "Orthorhombic"

    def __len__ ( self ): return 3

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        if alpha is None: alpha = self.__class__._alpha
        if beta  is None: beta  = self.__class__._beta
        if gamma is None: gamma = self.__class__._gamma
        QOK = ( a is not None ) and ( b is not None ) and ( c is not None ) and ( alpha == 90.0 ) and ( beta == 90.0 ) and ( gamma == 90.0 )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a != b != c ; alpha = beta = gamma = 90." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement  ] = symmetryParameterGradients.dEdA
        gradients[vectorIncrement+1] = symmetryParameterGradients.dEdB
        gradients[vectorIncrement+2] = symmetryParameterGradients.dEdC

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement  ] = symmetryParameters.a
        variables[vectorIncrement+1] = symmetryParameters.b
        variables[vectorIncrement+2] = symmetryParameters.c

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a = variables[vectorIncrement  ]
        b = variables[vectorIncrement+1]
        c = variables[vectorIncrement+2]
        symmetryParameters.SetCrystalParameters ( a, b, c, self.__class__._alpha, self.__class__._beta, self.__class__._gamma )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a, "b" : symmetryParameters.b, "c" : symmetryParameters.c }

    def IsOrthogonal ( self ): return True

class CrystalSystemTetragonal ( CrystalSystem ):
    """Tetragonal crystal class."""

    # . Class parameters.
    _alpha = 90.0
    _beta  = 90.0
    _gamma = 90.0
    _label = "Tetragonal"

    def __len__ ( self ): return 2

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        if b     is None: b     = a
        if alpha is None: alpha = self.__class__._alpha
        if beta  is None: beta  = self.__class__._beta
        if gamma is None: gamma = self.__class__._gamma
        QOK = ( a is not None ) and ( a == b ) and ( c is not None ) and ( alpha == 90.0 ) and ( beta == 90.0 ) and ( gamma == 90.0 )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a = b != c ; alpha = beta = gamma = 90." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement  ] = symmetryParameterGradients.dEdA + symmetryParameterGradients.dEdB
        gradients[vectorIncrement+1] = symmetryParameterGradients.dEdC

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement  ] = symmetryParameters.a
        variables[vectorIncrement+1] = symmetryParameters.c

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a = variables[vectorIncrement  ]
        c = variables[vectorIncrement+1]
        symmetryParameters.SetCrystalParameters ( a, a, c, self.__class__._alpha, self.__class__._beta, self.__class__._gamma )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a, "c" : symmetryParameters.c }

    def IsOrthogonal ( self ): return True

class CrystalSystemTriclinic ( CrystalSystem ):
    """Triclinic crystal class."""

    # . Class parameters.
    _alpha = 0.0
    _beta  = 0.0
    _gamma = 0.0
    _label = "Triclinic"

    def __len__ ( self ): return 6

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        QOK = ( a is not None ) and ( b is not None ) and ( c is not None )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a != b != c ; alpha != beta != gamma." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement  ] = symmetryParameterGradients.dEdA
        gradients[vectorIncrement+1] = symmetryParameterGradients.dEdB
        gradients[vectorIncrement+2] = symmetryParameterGradients.dEdC
        gradients[vectorIncrement+3] = symmetryParameterGradients.dEdAlpha
        gradients[vectorIncrement+4] = symmetryParameterGradients.dEdBeta
        gradients[vectorIncrement+5] = symmetryParameterGradients.dEdGamma

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement  ] = symmetryParameters.a
        variables[vectorIncrement+1] = symmetryParameters.b
        variables[vectorIncrement+2] = symmetryParameters.c
        variables[vectorIncrement+3] = symmetryParameters.alpha
        variables[vectorIncrement+4] = symmetryParameters.beta
        variables[vectorIncrement+5] = symmetryParameters.gamma

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a     = variables[vectorIncrement  ]
        b     = variables[vectorIncrement+1]
        c     = variables[vectorIncrement+2]
        alpha = variables[vectorIncrement+3]
        beta  = variables[vectorIncrement+4]
        gamma = variables[vectorIncrement+5]
        symmetryParameters.SetCrystalParameters ( a, b, c, alpha, beta, gamma )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a, "b" : symmetryParameters.b, "c" : symmetryParameters.c, "alpha" : symmetryParameters.alpha, "beta" : symmetryParameters.beta, "gamma" : symmetryParameters.gamma }

class CrystalSystemTrigonal ( CrystalSystem ):
    """Trigonal crystal class."""

    # . Class parameters.
    _alpha =  0.0
    _beta  = None
    _gamma = None
    _label = "Trigonal"

    def __len__ ( self ): return 2

    def AssignCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Verify the values of some crystal parameters."""
        if b     is None: b     = a
        if c     is None: c     = a
        if beta  is None: beta  = alpha
        if gamma is None: gamma = alpha
        QOK = ( a is not None ) and ( a == b ) and ( a == c ) and ( alpha == beta ) and ( alpha == gamma ) and ( alpha != 90.0 ) and ( alpha < 120.0 )
        if not QOK: raise SymmetryError ( "Invalid crystal parameters: a = b = c ; alpha = beta = gamma != 90 < 120." )
        return ( a, b, c, alpha, beta, gamma )

    def EmptyGradientsToVector ( self, gradients, symmetryParameterGradients, vectorIncrement ):
        """Put the symmetry parameter gradients into a vector."""
        gradients[vectorIncrement  ] = symmetryParameterGradients.dEdA     + symmetryParameterGradients.dEdB    + symmetryParameterGradients.dEdC
        gradients[vectorIncrement+1] = symmetryParameterGradients.dEdAlpha + symmetryParameterGradients.dEdBeta + symmetryParameterGradients.dEdGamma

    def EmptySymmetryParametersToVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Put the symmetry parameters into a vector."""
        variables[vectorIncrement  ] = symmetryParameters.a
        variables[vectorIncrement+1] = symmetryParameters.alpha

    def FillSymmetryParametersFromVector ( self, variables, symmetryParameters, vectorIncrement ):
        """Get the symmetry parameters from a vector."""
        a     = variables[vectorIncrement  ]
        alpha = variables[vectorIncrement+1]
        symmetryParameters.SetCrystalParameters ( a, a, a, alpha, alpha, alpha )

    def GetUniqueSymmetryParameters ( self, symmetryParameters ):
        """Return a dictionary with the unique symmetry parameters."""
        return { "a" : symmetryParameters.a, "alpha" : symmetryParameters.alpha }

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def CrystalSystem_FromLabel ( label ):
    """Return the crystal system given a label."""
    system = None
    if   label == "Cubic"        : return CrystalSystemCubic        ( )
    elif label == "Hexagonal"    : return CrystalSystemHexagonal    ( )
    elif label == "Monoclinic"   : return CrystalSystemMonoclinic   ( )
    elif label == "Orthorhombic" : return CrystalSystemOrthorhombic ( )
    elif label == "Tetragonal"   : return CrystalSystemTetragonal   ( )
    elif label == "Triclinic"    : return CrystalSystemTriclinic    ( )
    elif label == "Trigonal"     : return CrystalSystemTrigonal     ( )
    else: raise SymmetryError ( "Invalid crystal system label: {:s}.".format ( label ) )
    return system

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

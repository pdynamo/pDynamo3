"""Handle Gaussian basis sets."""

from  enum                      import Enum
from  pCore                     import RawObjectConstructor
from  pScientific.Arrays        import Array
from  pScientific.LinearAlgebra import OrthogonalizationMethod
from .GaussianBasisError        import GaussianBasisError
from .GaussianBasisUtilities    import AMLabelEncode    , \
                                       ShellLabelDecode , \
                                       ShellLabelEncode

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class GaussianBasisOperator ( Enum ):
    """Gaussian basis operators."""
    AntiCoulomb = 1 #( 1, "A", True  )
    Coulomb     = 2 #( 2, "C", True  )
    Dipole      = 3 #( 3, "D", False )
    Kinetic     = 4 #( 4, "K", True  )
    Overlap     = 5 #( 5, "O", True  )
    Poisson     = 6 #( 6, "P", True  )
    Quadrupole  = 7 #( 7, "Q", False )

# . Would be useful to have if can work out how to do it!
#    def __init__ ( self, order, code, isScalar ):
#        """Constructor."""
#        self.order    = order
#        self.code     = code
#        self.isScalar = isScalar

class GaussianBasisType ( Enum ):
    """The type of Gaussian basis."""
    Density = 1
    Orbital = 2

# . The YAML tag.
#_YAMLTag = "!GaussianBasis"

# . Basis sets of arbitary angular momentum are permitted.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasis:

    def __copy__ ( self ):
        """Copying."""
        cdef GaussianBasis new
        cdef CStatus       cStatus = CStatus_OK
        new         = self.__class__.Raw ( )
        new.cObject = GaussianBasis_Clone ( self.cObject, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error cloning Gaussian basis set." )
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            GaussianBasis_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        state = {}
        # . Basic data.
        state["Atomic Number"] = self.cObject.atomicNumber
        state["Basis Type"   ] = GaussianBasisType ( self.cObject.basisType ).name
        if self.cObject.isSpherical == CTrue: state["Is Spherical"] = True
        else:                                 state["Is Spherical"] = False
        if self.cObject.pNormalized == CTrue: state["Normalized Primitives"] = True
        else:                                 state["Normalized Primitives"] = False
        state["Primitive Fields"] = [ "Exponent", "Coefficients" ]
        # . Shell data.
        shells = []
        for i from 0 <= i < self.cObject.nShells:
            # . Primitive data.
            primitives = []
            for p from 0 <= p < self.cObject.shells[i].nPrimitives:
                coefficients = []
                for c from 0 <= c < ( self.cObject.shells[i].lHigh - self.cObject.shells[i].lLow + 1 ):
                    coefficients.append ( self.cObject.shells[i].primitives[p].coefficients0[c] )
                primitives.append ( [ self.cObject.shells[i].primitives[p].exponent0, coefficients ] )
            shells.append ( { "Primitives" : primitives ,
                              "Type"       : ShellLabelEncode ( self.cObject.shells[i].lLow, self.cObject.shells[i].lHigh ) } )
        state["Shells"] = shells
        return state

    def __init__ ( self                                    ,
                   atomicNumber                            ,
                   numberOfShells                          ,
                   basisType   = GaussianBasisType.Orbital ,
                   isSpherical = True                      ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfShells )
        self.cObject.atomicNumber = atomicNumber
        self.cObject.basisType    = basisType.value
        if isSpherical: self.cObject.isSpherical = CTrue
        else:           self.cObject.isSpherical = CFalse

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        # . Reallocate the object.
        shells = state["Shells"]
        self._Allocate ( len ( shells ) )
        # . Fill the object.
        # . Basic data.
        self.cObject.atomicNumber = state["Atomic Number"]
        if "Basis Type" in state: self.cObject.basisType = GaussianBasisType.__dict__[state["Basis Type"]].value
        if state["Is Spherical"         ]: self.cObject.isSpherical = CTrue
        else:                              self.cObject.isSpherical = CFalse
        if state["Normalized Primitives"]: self.cObject.pNormalized = CTrue
        else:                              self.cObject.pNormalized = CFalse
        # . Shell data.
        for ( i, shell ) in enumerate ( shells ):
            ( lLow, lHigh ) = ShellLabelDecode ( shell["Type"] )
            primitives = shell["Primitives"]
            self.AllocateShell ( i, lHigh, lLow, len ( primitives ) )
            # . Primitive data.
            for ( p, ( exponent, coefficients ) ) in enumerate ( primitives ):
                self.cObject.shells[i].primitives[p].exponent  = exponent
                self.cObject.shells[i].primitives[p].exponent0 = exponent
                # . Coefficient data.
                for ( c, coefficient ) in enumerate ( coefficients ):
                    self.cObject.shells[i].primitives[p].coefficients [c] = coefficient
                    self.cObject.shells[i].primitives[p].coefficients0[c] = coefficient
        # . Finish up.
        self.UnnormalizePrimitives ( )

    def _Allocate ( self, numberOfShells ):
        """Allocation."""
        self.cObject = GaussianBasis_Allocate ( numberOfShells )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def AllocateShell ( self, CInteger iShell, CInteger lHigh, CInteger lLow, CInteger nPrimitives ):
        """Start defining data for a shell."""
        if ( iShell >=                    0 ) and \
           ( iShell <  self.cObject.nShells ) and \
           ( lHigh  >=                    0 ) and \
           ( lHigh  >= lLow                 ) and \
           ( lLow   >=                    0 ):
            GaussianBasis_AllocateShell ( self.cObject, iShell, lHigh, lLow, nPrimitives )
        else: raise GaussianBasisError ( "Error allocating shell with index {:d}.".format ( iShell ) )

    @staticmethod
    def CartesianToSphericalTransformation ( CInteger lLow, CInteger lHigh ):
        """Return the C->S transformation for lLow to lHigh."""
        cdef RealArray2D   result
        cdef CRealArray2D *cResult
        cdef CView2D      *cView
        cdef CStatus       cStatus = CStatus_OK
        cResult = GaussianBasis_TransformationCartesianToSpherical ( lLow, lHigh, &cStatus )
        cView   = < CView2D * > cResult
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error generating C->S transformation." )
        result = Array.WithExtents ( View2D_GetRows ( cView ), View2D_GetColumns ( cView ) )
        RealArray2D_CopyTo     ( cResult, result.cObject, NULL )
        RealArray2D_Deallocate ( &cResult )
        return result

    def Finalize ( self ):
        """Finalize the basis ready for use."""
        cdef CStatus cStatus = CStatus_OK
        GaussianBasis_Finalize ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error finalizing basis set." )

    def MakeFunctionLabels ( self ):
        """Make an array of function labels."""
        # . Shell number followed by Cartesian or spherical harmonic identity.
        cdef CInteger i
        labels = []
        if self.cObject != NULL:
            isCartesian = ( self.cObject.isSpherical == CFalse )
            for i from 0 <= i < self.cObject.nShells:
                aLow  = self.cObject.shells[i].lLow
                aHigh = self.cObject.shells[i].lHigh
                if isCartesian:
                    for a in range ( aLow, aHigh + 1 ):
                        l = "{:d}".format ( i+1 )
                        if a == 0:
                            labels.append ( l + "s" )
                        else:
                            for z in range ( 0, a + 1 ):
                                for y in range ( 0, a - z + 1 ):
                                    x = a - y - z
                                    labels.append ( l + "".join ( x * [ "x" ] + y * [ "y" ] + z * [ "z" ] ) )
                else:
                    for a in range ( aLow, aHigh + 1 ):
                        l = "{:d}{:s}".format ( i+1, AMLabelEncode ( a ).lower ( ) )
                        if a == 0:
                            labels.append ( l )
                        else:
                            labels.append ( l + "0" )
                            for m in range ( 1, a+1 ):
                                labels.extend ( [ "{:s}{:d}c".format ( l, m ), "{:s}{:d}s".format ( l, m ) ] )
        return labels

    def Orthonormalize ( self                              ,
                         basisType                 = None  ,
                         orthogonalizationMethod   = None  ,
                         orthogonalizationOperator = None  ):
        """Orthonormalize the basis."""
        # . Principally for testing so results["Status"] must be checked on exit.
        cdef RealArray2D              M
        cdef RealArray2D              X
        cdef RealArray2D              Y
        cdef CInteger                 nIndependent = 0
        cdef CReal                    deviation    = 0.0
        cdef CGaussianBasisOperator   cOperator = GaussianBasisOperator_Overlap
        cdef COrthogonalizationMethod cMethod   = OrthogonalizationMethod_Symmetric
        cdef CStatus                  cStatus   = CStatus_OK
        if basisType                 is not None: self.cObject.basisType = basisType.value
        if orthogonalizationMethod   is not None: cMethod                = orthogonalizationMethod.value
        if orthogonalizationOperator is not None: cOperator              = orthogonalizationOperator.value
        d = self.numberOfFunctions
        M = RealArray2D.WithExtents ( d, d )
        X = RealArray2D.WithExtents ( d, d )
        Y = RealArray2D.WithExtents ( d, d )
        GaussianBasis_Orthonormalize ( self.cObject  ,
                                       cOperator     ,
                                       cMethod       ,
                                       &nIndependent ,
                                       &deviation    ,
                                       M.cObject     ,
                                       X.cObject     ,
                                       Y.cObject     ,
                                       &cStatus      )
        #if cStatus != CStatus_OK: raise GaussianBasisError ( "Error orthonormalizing Gaussian basis set." )
        return { "Basis Type"                 : self.cObject.basisType ,
                 "Independent Functions"      : nIndependent           ,
                 "Inverse Orthogonalizer"     : Y                      ,
                 "Maximum Deviation"          : deviation              ,
                 "Metric Matrix"              : M                      ,
                 "Orthogonalization Method"   : cMethod                ,
                 "Orthogonalization Operator" : cOperator              , 
                 "Orthogonalizer"             : X                      ,
                 "Status"                     : cStatus                ,
                 "Total Functions"            : d                      }

    @classmethod
    def Raw ( selfClass ):
        """Constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ScaleShellExponents ( self, CInteger index, CReal zeta ):
        """Scale the exponents of a shell."""
        cdef CStatus cStatus = CStatus_OK
        GaussianBasis_ScaleShellExponents ( self.cObject, index, zeta, &cStatus )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error scaling {:d} shell exponents.".format ( index ) )

    def SetShellCoefficients ( self, CInteger iShell, CInteger L, coefficients ):
        """Set the coefficients for a shell."""
        cdef CInteger p
        if ( iShell >= 0 ) and ( iShell < self.cObject.nShells ):
            if len ( coefficients ) == self.cObject.shells[iShell].nPrimitives:
                if ( L < self.cObject.shells[iShell].lLow ) or ( L > self.cObject.shells[iShell].lHigh ):
                    raise GaussianBasisError ( "Invalid angular momentum value, {:d}, for coefficients.".format ( L ) )
                L -= self.cObject.shells[iShell].lLow
                for p from 0 <= p < self.cObject.shells[iShell].nPrimitives:
                    self.cObject.shells[iShell].primitives[p].coefficients [L] = coefficients[p]
                    self.cObject.shells[iShell].primitives[p].coefficients0[L] = coefficients[p]

    def SetShellExponents ( self, CInteger iShell, exponents, CReal scale = 1.0 ):
        """Set the coefficients for a shell."""
        cdef CInteger p
        if ( iShell >= 0 ) and ( iShell < self.cObject.nShells ):
            if len ( exponents ) == self.cObject.shells[iShell].nPrimitives:
                for p from 0 <= p < self.cObject.shells[iShell].nPrimitives:
                    self.cObject.shells[iShell].primitives[p].exponent  = exponents[p] * scale * scale
                    self.cObject.shells[iShell].primitives[p].exponent0 = exponents[p] * scale * scale

    @staticmethod
    def SphericalToCartesianTransformation ( CInteger lLow, CInteger lHigh, RealArray2D c2s not None ):
        """Return the S->C transformation for lLow to lHigh."""
        cdef RealArray2D   result
        cdef CRealArray2D *cResult
        cdef CView2D      *cView
        cdef CStatus       cStatus = CStatus_OK
        cResult = GaussianBasis_TransformationSphericalToCartesian ( lLow, lHigh, c2s.cObject, &cStatus )
        cView   = < CView2D * > cResult
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error generating S->C transformation." )
        result = Array.WithExtents ( View2D_GetRows ( cView ), View2D_GetColumns ( cView ) )
        RealArray2D_CopyTo     ( cResult, result.cObject, NULL )
        RealArray2D_Deallocate ( &cResult )
        return result 

    @classmethod
    def Uninitialized ( selfClass, atomicNumber, numberOfShells, basisType = GaussianBasisType.Orbital, isSpherical = True ):
        """Constructor."""
        return selfClass ( atomicNumber, numberOfShells, basisType = basisType, isSpherical = isSpherical )

    def UnnormalizePrimitives ( self ):
        """Unnormalize the primitives of the basis."""
        GaussianBasis_UnnormalizePrimitives ( self.cObject )

    # . Properties.
    @property
    def atomicNumber ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.atomicNumber
    @atomicNumber.setter
    def atomicNumber ( self, CInteger value ):
        if self.cObject != NULL: self.cObject.atomicNumber = value

    @property
    def basisType ( self ):
        if self.cObject == NULL: return None
        else:                    return GaussianBasisType ( self.cObject.basisType )

    @property
    def isSpherical ( self ):
        if self.cObject == NULL: return False
        else:                    return ( self.cObject.isSpherical == CTrue )

    @property
    def maximumAngularMomentum ( self ):
        if self.cObject == NULL: return -1
        else:                    return self.cObject.lHigh

    @property
    def numberOfFunctions ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nBasis

    @property
    def numberOfShells ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nShells

    @property
    def numberOfWorkFunctions ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nCBF

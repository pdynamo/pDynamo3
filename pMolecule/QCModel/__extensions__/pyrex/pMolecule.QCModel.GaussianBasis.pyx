"""Handle Gaussian basis sets."""

from  enum         import Enum
from  pCore        import RawObjectConstructor
from .QCModelError import QCModelError

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class BasisType ( Enum ):
    """The type of basis."""
    Coulomb = 1
    Orbital = 2
    Poisson = 3

class NormalizationType ( Enum ):
    """The type of normalization."""
    Canonical = 1
    Diagonal  = 2
    Symmetric = 3

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The shell angular momenta.
_ShellAngularMomenta = [ "S", "P", "D", "F", "G" ]

# . The shell types.
_ShellTypes = [ "S", "P", "D", "F", "G", "SP", "SPD", "SPDF", "SPDFG" ]

# . Shell types accessed by low and high angular momenta.
_ShellTypeDictionary = { (0,0) : "S", (1,1) : "P", (2,2) : "D", (3,3) : "F", (4,4) : "G", (0,1) : "SP", (0,2) : "SPD", (0,3) : "SPDF", (0,4) : "SPDFG" }

# . The YAML tag.
#_YAMLTag = "!GaussianBasis"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasis:

    def __copy__ ( self ):
        """Copying."""
        cdef GaussianBasis new
        new         = self.__class__.Raw ( )
        new.cObject = GaussianBasis_Clone ( self.cObject )
        new.isOwner = True
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
        state["Atomic Number"     ] = self.cObject.atomicNumber
        state["Basis Type"        ] = BasisType         ( self.cObject.basisType         ).name
        state["Normalization Type"] = NormalizationType ( self.cObject.normalizationType ).name
        if self.cObject.QSPHERICAL            == CTrue: state["Is Spherical"] = True
        else:                                           state["Is Spherical"] = False
        if self.cObject.QNORMALIZEDPRIMITIVES == CTrue: state["Normalized Primitives"] = True
        else:                                           state["Normalized Primitives"] = False
        state["Primitive Fields"] = [ "Exponent", "Coefficients" ]
        # . Shell data.
        shells = []
        for i from 0 <= i < self.cObject.nshells:
            # . Primitive data.
            primitives = []
            for p from 0 <= p < self.cObject.shells[i].nprimitives:
                coefficients = []
                for c from 0 <= c < self.cObject.shells[i].type.angularmomentum_high - self.cObject.shells[i].type.angularmomentum_low + 1:
                    coefficients.append ( self.cObject.shells[i].primitives[p].coefficients0[c] )
                primitives.append ( [ self.cObject.shells[i].primitives[p].exponent0, coefficients ] )
            shells.append ( { "Primitives" : primitives , \
                              "Type"       : _ShellTypeDictionary[(self.cObject.shells[i].type.angularmomentum_low,self.cObject.shells[i].type.angularmomentum_high)] } )
        state["Shells"] = shells
        return state

    def __init__ ( self, atomicNumber, numberOfShells, basisType = BasisType.Orbital, isSpherical = True, normalizationType = NormalizationType.Symmetric ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfShells )
        self.cObject.atomicNumber = atomicNumber
        if isSpherical: self.cObject.QSPHERICAL = CTrue
        else:           self.cObject.QSPHERICAL = CFalse
        self.cObject.basisType         = basisType.value
        self.cObject.normalizationType = normalizationType.value

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
        if "Basis Type"         in state: self.cObject.basisType         = BasisType.__dict__[state["Basis Type"]].value
        if "Normalization Type" in state: self.cObject.normalizationType = NormalizationType.__dict__[state["Normalization Type"]].value
        if state["Is Spherical"         ]: self.cObject.QSPHERICAL            = CTrue
        else:                              self.cObject.QSPHERICAL            = CFalse
        if state["Normalized Primitives"]: self.cObject.QNORMALIZEDPRIMITIVES = CTrue
        else:                              self.cObject.QNORMALIZEDPRIMITIVES = CFalse
        # . Shell data.
        for ( i, shell ) in enumerate ( shells ):
            primitives = shell["Primitives"]
            self.AllocateShell ( i, shell["Type"], len ( primitives ) )
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
        # . Normalize here?
        #self.Normalize ( )

    def _Allocate ( self, numberOfShells ):
        """Allocation."""
        self.cObject = GaussianBasis_Allocate ( numberOfShells )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def AllocateShell ( self, CInteger ishell, object label, CInteger nprimitives ):
        """Start defining data for a shell."""
        if ( ishell >= 0 ) and ( ishell < self.cObject.nshells ):
            try:    ilabel = _ShellTypes.index ( label )
            except: raise QCModelError ( "Invalid shell type: " + label + "." )
            GaussianBasis_AllocateShell ( self.cObject, ishell, nprimitives, ilabel )

    def Normalize ( self, basisType = None, doReport = False, normalizationType = None ):
        """Normalize the basis."""
        cdef CBoolean cDoReport
        cdef CReal    deviation = 0.0
        cdef CStatus  status    = CStatus_OK
        if doReport: cDoReport = CTrue
        else:        cDoReport = CFalse
        if basisType         is not None: self.cObject.basisType         = basisType.value
        if normalizationType is not None: self.cObject.normalizationType = normalizationType.value
        deviation = GaussianBasis_Normalize ( self.cObject, cDoReport, &status )
        if doReport:
            report = { "Basis Type"         : self.cObject.basisType         ,
                       "Maximum Deviation"  : deviation                      ,
                       "Normalization Type" : self.cObject.normalizationType ,
                       "Status"             : status                         }
        else:
            if status != CStatus_OK: raise QCModelError ( "Error normalizing Gaussian basis set." )
            report = None
        return report

    @classmethod
    def Raw ( selfClass ):
        """Constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ScaleShellExponents ( self, CInteger  index, CReal zeta ):
        """Scale the exponents of a shell."""
        cdef CStatus cStatus = CStatus_OK
        GaussianBasis_ScaleShellExponents ( self.cObject, index, zeta, &cStatus )
        if cStatus != CStatus_OK: raise QCModelError ( "Error scaling {:d} shell exponents.".format ( index ) )

    def SetShellCoefficients ( self, CInteger ishell, object label, object coefficients ):
        """Set the coefficients for a shell."""
        if ( ishell >= 0 ) and ( ishell < self.cObject.nshells ):
            if len ( coefficients ) == self.cObject.shells[ishell].nprimitives:
                try:    am = _ShellAngularMomenta.index ( label )
                except: raise QCModelError ( "Invalid angular momentum type: " + label + "." )
                if ( am < self.cObject.shells[ishell].type.angularmomentum_low ) and ( am > self.cObject.shells[ishell].type.angularmomentum_high ):
                    raise QCModelError ( "Invalid angular momentum value for coefficients: " + label + "." )
                am = am - self.cObject.shells[ishell].type.angularmomentum_low
                for i from 0 <= i < self.cObject.shells[ishell].nprimitives:
                    self.cObject.shells[ishell].primitives[i].coefficients [am] = coefficients[i]
                    self.cObject.shells[ishell].primitives[i].coefficients0[am] = coefficients[i]

    def SetShellExponents ( self, CInteger ishell, object exponents, CReal scale = 1.0 ):
        """Set the coefficients for a shell."""
        if ( ishell >= 0 ) and ( ishell < self.cObject.nshells ):
            if len ( exponents ) == self.cObject.shells[ishell].nprimitives:
                for i from 0 <= i < self.cObject.shells[ishell].nprimitives:
                    self.cObject.shells[ishell].primitives[i].exponent  = exponents[i] * scale * scale
                    self.cObject.shells[ishell].primitives[i].exponent0 = exponents[i] * scale * scale

    @classmethod
    def Uninitialized ( selfClass, atomicNumber, numberOfShells, basisType = BasisType.Orbital, isSpherical = True ):
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
        else:                    return BasisType ( self.cObject.basisType )

    @property
    def isSpherical ( self ):
        if self.cObject == NULL: return False
        else:                    return ( self.cObject.QSPHERICAL == CTrue )

    @property
    def maximumAngularMomentum ( self ):
        if self.cObject == NULL: return -1
        else:                    return self.cObject.maximum_angularmomentum

    @property
    def normalizationType ( self ):
        if self.cObject == NULL: return None
        else:                    return NormalizationType ( self.cObject.normalizationType )

    @property
    def numberOfFunctions ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nbasis

    @property
    def numberOfShells ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nshells

    @property
    def numberOfWorkFunctions ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.nbasisw

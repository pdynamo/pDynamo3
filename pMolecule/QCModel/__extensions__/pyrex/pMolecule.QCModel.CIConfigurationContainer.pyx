"""A container for CI configurations."""

import os, os.path

from  pCore        import RawObjectConstructor
from .QCModelError import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CIConfigurationContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef CIConfigurationContainer new
        new         = self.__class__.Raw ( )
        new.cObject = CIConfigurationContainer_Clone ( self.cObject, NULL )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            CIConfigurationContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "Active Electrons" : self.numberOfActiveElectrons , 
                 "Active Orbitals"  : self.numberOfActiveOrbitals  , 
                 "MicroStates"      : self.microStates             }

    def __init__ ( self, activeOrbitals, configurations ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( activeOrbitals, configurations )

    def __len__ ( self ): return self.numberOfConfigurations

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._CObjectFromMicroStates ( state["MicroStates"], state["Active Electrons"], state["Active Orbitals"] )

    def _Allocate ( self, orbitals, configurations ):
        """Constructor."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = CIConfigurationContainer_Allocate ( orbitals, configurations, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error allocating CI configuration container." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject      = NULL
        self.isOwner      = False
        self._microStates = None

    def _CObjectFromMicroStates ( self, microStates not None, activeElectrons, activeOrbitals ):
        """C object constructor from a set of microstates."""
        cdef CInteger          i, j
        cdef CIntegerArray2D  *cMicroStates    = NULL
        cdef CStatus           cStatus         = CStatus_OK
        # . Convert input microStates to C.
        isOK = ( len ( microStates ) > 0 )
        if isOK:
            cMicroStates = IntegerArray2D_AllocateWithExtents ( len ( microStates ), 2 * activeOrbitals, &cStatus )
            isOK = ( cMicroStates != NULL )
            if isOK:
                IntegerArray2D_Set ( cMicroStates, 0 )
                for ( i, microState ) in enumerate ( microStates ):
                    indices = [ int ( c ) for c in microState ]
                    if ( len ( indices ) != ( 2 * activeOrbitals ) ) or ( sum ( indices ) != activeElectrons ):
                        isOK = False
                        break
                    for j from 0 <= j < ( 2 * activeOrbitals ): IntegerArray2D_SetItem ( cMicroStates, i, j, indices[j], NULL )
        if not isOK: raise QCModelError ( "Invalid microstate specification." )
        # . Construct the C object.
        self.cObject = CIConfigurationContainer_MakeUserSpecified ( cMicroStates, activeOrbitals, activeElectrons, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error making a CI configuration container from microstates." )
        IntegerArray2D_Deallocate ( &cMicroStates )

    def Characters ( self, RealArray2D orbitalTransformation not None ,
                           RealArray2D stateTransformation   not None ,
                           coreOrbitals = 0, includeCoreOrbitals = False ):
        """Determine characters for the CI configurations."""
        cdef CBoolean cIncludeCoreOrbitals
        cdef CStatus  cStatus = CStatus_OK
        if includeCoreOrbitals: cIncludeCoreOrbitals = CTrue
        else:                   cIncludeCoreOrbitals = CFalse
        CIConfigurationContainer_Characters ( self.cObject                  ,
                                              cIncludeCoreOrbitals          ,
                                              coreOrbitals                  ,
                                              orbitalTransformation.cObject ,
                                              stateTransformation.cObject   ,
                                              &cStatus                      )
        if cStatus != CStatus_OK: raise QCModelError ( "Error determining CI configuration characters." )

    def CIMatrixSparsity ( self ):
        """Get estimates of the sparsity of the CI matrix."""
        cdef CInteger  nonZero  = 0
        cdef CReal     sparsity = 0.0
        CIConfigurationContainer_GetCIMatrixSparsity ( self.cObject, &nonZero, &sparsity )
        return ( nonZero, sparsity )

    @classmethod
    def FromMicroStates ( selfClass, microStates not None, activeElectrons, activeOrbitals ):
        """Constructor from a set of microstates."""
        self = selfClass.Raw ( )
        self._CObjectFromMicroStates ( microStates, activeElectrons, activeOrbitals )
        return self

    def MakeCIMatrix ( self, SymmetricMatrix       fCoreMO not None ,
                             DoubleSymmetricMatrix moTEIs  not None ,
                             SymmetricMatrix       matrixFull       ,
                             SparseSymmetricMatrix matrixSparse     ):
        """Make the sparse and/or full CI matrices."""
        cdef CSparseSymmetricMatrix *cMatrixSparse = NULL
        cdef CSymmetricMatrix       *cMatrixFull   = NULL
        cdef CStatus                 cStatus       = CStatus_OK
        if matrixFull   is not None: cMatrixFull   = matrixFull.cObject
        if matrixSparse is not None: cMatrixSparse = matrixSparse.cObject
        CIConfigurationContainer_MakeCIMatrix ( self.cObject     ,
                                                fCoreMO.cObject  ,
                                                moTEIs.cObject   ,
                                                cMatrixFull      ,
                                                cMatrixSparse    ,
                                                &cStatus         )
        if cStatus != CStatus_OK: raise QCModelError ( "Error making CI matrix." )

    def MakeDensities ( self, RealArray1D           ciVector  not None ,
                              SymmetricMatrix       onePDMMOt not None ,
                              SymmetricMatrix       onePDMMOs not None ,
                              DoubleSymmetricMatrix twoPDM    not None ):
        """Make the total and spin CI densities in the MO basis."""
        cdef CStatus cStatus = CStatus_OK
        CIConfigurationContainer_MakeDensities ( self.cObject      ,
                                                 ciVector.cObject  ,
                                                 onePDMMOt.cObject ,
                                                 onePDMMOs.cObject ,
                                                 twoPDM.cObject    ,
                                                 &cStatus          )
        if cStatus != CStatus_OK: raise QCModelError ( "Error making CI densities." )

    @classmethod
    def MakeDoubles ( selfClass, alpha, beta, activeOrbitals ):
        """Constructor of doubles configurations given the number of active orbitals and of alpha and beta electrons."""
        cdef CIConfigurationContainer self
        cdef CStatus                  cStatus = CStatus_OK
        self         = selfClass.Raw ( )
        self.cObject = CIConfigurationContainer_MakeSinglesDoubles ( CFalse, CTrue, activeOrbitals, beta, ( alpha - beta ), &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error making a doubles CI configuration container." )
        return self

    @classmethod
    def MakeFull ( selfClass, alpha, beta, activeOrbitals ):
        """Constructor of all possible configurations given the number of active orbitals and of alpha and beta electrons."""
        cdef CIConfigurationContainer self
        cdef CStatus                  cStatus = CStatus_OK
        self         = selfClass.Raw ( )
        self.cObject = CIConfigurationContainer_MakeFull ( activeOrbitals, alpha, beta, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error making a full CI configuration container." )
        return self

    @classmethod
    def MakeSingles ( selfClass, alpha, beta, activeOrbitals ):
        """Constructor of singles configurations given the number of active orbitals and of alpha and beta electrons."""
        cdef CIConfigurationContainer self
        cdef CStatus                  cStatus = CStatus_OK
        self         = selfClass.Raw ( )
        self.cObject = CIConfigurationContainer_MakeSinglesDoubles ( CTrue, CFalse, activeOrbitals, beta, ( alpha - beta ), &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error making a singles and doubles CI configuration container." )
        return self

    @classmethod
    def MakeSinglesDoubles ( selfClass, alpha, beta, activeOrbitals ):
        """Constructor of singles and doubles configurations given the number of active orbitals and of alpha and beta electrons."""
        cdef CIConfigurationContainer self
        cdef CStatus                  cStatus = CStatus_OK
        self         = selfClass.Raw ( )
        self.cObject = CIConfigurationContainer_MakeSinglesDoubles ( CTrue, CTrue, activeOrbitals, beta, ( alpha - beta ), &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error making a singles and doubles CI configuration container." )
        return self

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def StateSpins ( self, RealArray2D vectors not None, RealArray1D spins not None ):
        """Find the spins of a set of CI state vectors."""
        cdef CStatus cStatus = CStatus_OK
        CIConfigurationContainer_StateSpins ( self.cObject, vectors.cObject, spins.cObject, &cStatus )
        if cStatus != CStatus_OK: raise QCModelError ( "Error finding CI vector state spins." )

    def TransitionDipoles ( self, SymmetricMatrix tdMOs not None, SymmetricMatrix tdMatrix not None ):
        """Make the transition dipoles."""
        CIConfigurationContainer_TransitionDipoles ( self.cObject, tdMOs.cObject, tdMatrix.cObject )

    @property
    def microStates ( self ):
        """Convert the container configurations to a list of 0/1 strings with alpha orbitals first."""
        cdef CInteger  a, b, c
        microStates = self._microStates
        if ( microStates is None ) and ( self.cObject != NULL ):
            microStates = []
            if ( self.cObject != NULL ) and ( self.cObject.nActive > 0 ) and ( self.cObject.nConfigurations > 0 ):
                for c from 0 <= c < self.cObject.nConfigurations:
                    state = []
                    for a from 0 <= a < self.cObject.nActive:
                        if ( IntegerArray1D_GetItem ( self.cObject.configurations[c].alphas, a, NULL ) == 0 ): state.append ( "0" )
                        else:                                                                                  state.append ( "1" )
                    for b from 0 <= b < self.cObject.nActive:
                        if ( IntegerArray1D_GetItem ( self.cObject.configurations[c].betas , b, NULL ) == 0 ): state.append ( "0" )
                        else:                                                                                  state.append ( "1" )
                    microStates.append ( "".join ( state ) )
            self._microStates = microStates
        return microStates

    @property
    def numberOfActiveElectrons ( self ):
        return CIConfigurationContainer_NumberOfActiveElectrons ( self.cObject )

    @property
    def numberOfActiveOrbitals ( self ):
        return CIConfigurationContainer_NumberOfActiveOrbitals ( self.cObject )

    @property
    def numberOfConfigurations ( self ):
        return CIConfigurationContainer_NumberOfConfigurations ( self.cObject )

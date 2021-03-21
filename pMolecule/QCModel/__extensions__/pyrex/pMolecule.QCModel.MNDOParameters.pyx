"""Handle the parameter data necessary for a MNDO NDDO-type calculation."""

from pCore import RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The YAML tag.
#_YAMLTag = "!MNDOParameters"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOParameters:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef MNDOParameters new
        new         = self.__class__.Raw ( )
        new.cObject = MNDOParameters_Clone ( self.cObject )
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            MNDOParameters_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getmodule__ ( self ):
        """Return the module name."""
        return "pMolecule.MNDOParameters"

    def __getstate__ ( self ):
        """Return the state."""
        state = {}
        # . Scalars.
        parameters = [ [ "atomicNumber" , self.cObject.atomicNumber , "none"          ] ,
                       [ "orbitals"     , self.cObject.norbitals    , "none"          ] ,
                       [ "iii"          , self.cObject.iii          , "none"          ] ,
                       [ "iiid"         , self.cObject.iiid         , "none"          ] ,
                       [ "ir016"        , self.cObject.ir016        , "none"          ] ,
                       [ "ir066"        , self.cObject.ir066        , "none"          ] ,
                       [ "ir244"        , self.cObject.ir244        , "none"          ] ,
                       [ "ir266"        , self.cObject.ir266        , "none"          ] ,
                       [ "ir466"        , self.cObject.ir466        , "none"          ] ,
                       [ "qnd"          , self.cObject.qnd          , "none"          ] ,
                       [ "qnp"          , self.cObject.qnp          , "none"          ] ,
                       [ "qns"          , self.cObject.qns          , "none"          ] ,
                       [ "ad"           , self.cObject.ad0          , "atomic"        ] ,  
                       [ "alp"          , self.cObject.alp0         , "A^-1"          ] ,
                       [ "am"           , self.cObject.am0          , "atomic"        ] ,
                       [ "aq"           , self.cObject.aq0          , "atomic"        ] ,
                       [ "betad"        , self.cObject.betad0       , "eV"            ] ,
                       [ "betap"        , self.cObject.betap0       , "eV"            ] ,
                       [ "betas"        , self.cObject.betas0       , "eV"            ] ,
                       [ "dd"           , self.cObject.dd0          , "atomic"        ] ,
                       [ "eheat"        , self.cObject.eheat0       , "kcal/mole"     ] ,
                       [ "eisol"        , self.cObject.eisol0       , "eV"            ] ,
                       [ "f0sd"         , self.cObject.f0sd0        , "eV"            ] ,
                       [ "gphot"        , self.cObject.gphot0       , "dimensionless" ] ,
                       [ "gpp"          , self.cObject.gpp0         , "eV"            ] ,
                       [ "gp2"          , self.cObject.gp20         , "eV"            ] ,
                       [ "gsp"          , self.cObject.gsp0         , "eV"            ] ,
                       [ "gss"          , self.cObject.gss0         , "eV"            ] ,
                       [ "g2sd"         , self.cObject.g2sd0        , "eV"            ] ,
                       [ "hsp"          , self.cObject.hsp0         , "eV"            ] ,
                       [ "pcore"        , self.cObject.pcore0       , "atomic"        ] ,
                       [ "qq"           , self.cObject.qq0          , "atomic"        ] ,
                       [ "udd"          , self.cObject.udd0         , "eV"            ] ,
                       [ "upp"          , self.cObject.upp0         , "eV"            ] ,
                       [ "uss"          , self.cObject.uss0         , "eV"            ] ,
                       [ "zcore"        , self.cObject.zcore0       , "atomic"        ] ,
                       [ "zetad"        , self.cObject.zetad0       , "atomic"        ] ,
                       [ "zetap"        , self.cObject.zetap0       , "atomic"        ] ,
                       [ "zetas"        , self.cObject.zetas0       , "atomic"        ] ,
                       [ "zdn"          , self.cObject.zdn0         , "atomic"        ] ,
                       [ "zpn"          , self.cObject.zpn0         , "atomic"        ] ,
                       [ "zsn"          , self.cObject.zsn0         , "atomic"        ] ]
        state["Scalar Parameter Fields"] = [ "Label", "Value", "Units" ]
        state["Scalar Parameters"      ] = parameters
        # . Arrays.
        if self.cObject.nam1pm3g > 0:
            parameters = []
            for i from 0 <= i < self.cObject.nam1pm3g:
                parameters.append ( [ self.cObject.fn10[i] ,
                                      self.cObject.fn20[i] ,
                                      self.cObject.fn30[i] ] )
            state["AM1/PM3 Gaussian Parameter Fields"] = [ "fn1"  , "fn2"  , "fn3" ]
            state["AM1/PM3 Gaussian Parameter Units" ] = [ "A eV" , "A^-2" , "A"   ]
            state["AM1/PM3 Gaussian Parameters"      ] = parameters
        if self.cObject.ndiatomic > 0:
            parameters = []
            for i from 0 <= i < self.cObject.ndiatomic:
                if self.cObject.QDIATOMICFLAGS[i] == CTrue:
                    parameters.append ( [ i, self.cObject.diatomicx0[i], self.cObject.diatomica0[i] ] )
            state["Diatomic Parameter Fields"] = [ "atomicNumber" , "coefficient"   , "exponent" ]
            state["Diatomic Parameter Units" ] = [ "none"         , "dimensionless" , "A^-1"     ]
            state["Diatomic Parameters"      ] = parameters
        if self.cObject.npddg > 0:
            parameters = []
            for i from 0 <= i < self.cObject.npddg:
                parameters.append ( [ self.cObject.pddgc0[i], self.cObject.pddge0[i] ] )
            state["PDDG Gaussian Parameter Fields"] = [ "coefficient", "exponentDistance" ]
            state["PDDG Gaussian Parameter Units" ] = [ "eV"         , "A"                ]
            state["PDDG Gaussian Parameters"      ] = parameters
        return state

    def __init__ ( self, atomicNumber, numberOfOrbitals ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( )
        self.cObject.atomicNumber = atomicNumber
        self.cObject.norbitals    = numberOfOrbitals

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        # . Allocate the object.
        self._Allocate ( )
        # . Fill the object.
        # . Scalars.
        parameters = state["Scalar Parameters"]
        for ( label, value, units ) in parameters:
            if   label == "atomicNumber" : self.cObject.atomicNumber = value
            elif label == "orbitals"     : self.cObject.norbitals    = value
            elif label == "iii"          : self.cObject.iii          = value
            elif label == "iiid"         : self.cObject.iiid         = value
            elif label == "ir016"        : self.cObject.ir016        = value
            elif label == "ir066"        : self.cObject.ir066        = value
            elif label == "ir244"        : self.cObject.ir244        = value
            elif label == "ir266"        : self.cObject.ir266        = value
            elif label == "ir466"        : self.cObject.ir466        = value
            elif label == "qnd"          : self.cObject.qnd          = value
            elif label == "qnp"          : self.cObject.qnp          = value
            elif label == "qns"          : self.cObject.qns          = value
            elif label == "ad"           : self.cObject.ad0          = value
            elif label == "alp"          : self.cObject.alp0         = value
            elif label == "am"           : self.cObject.am0          = value
            elif label == "aq"           : self.cObject.aq0          = value
            elif label == "betad"        : self.cObject.betad0       = value
            elif label == "betap"        : self.cObject.betap0       = value
            elif label == "betas"        : self.cObject.betas0       = value
            elif label == "dd"           : self.cObject.dd0          = value
            elif label == "eheat"        : self.cObject.eheat0       = value
            elif label == "eisol"        : self.cObject.eisol0       = value
            elif label == "f0sd"         : self.cObject.f0sd0        = value
            elif label == "gphot"        : self.cObject.gphot0       = value
            elif label == "gpp"          : self.cObject.gpp0         = value
            elif label == "gp2"          : self.cObject.gp20         = value
            elif label == "gsp"          : self.cObject.gsp0         = value
            elif label == "gss"          : self.cObject.gss0         = value
            elif label == "g2sd"         : self.cObject.g2sd0        = value
            elif label == "hsp"          : self.cObject.hsp0         = value
            elif label == "pcore"        : self.cObject.pcore0       = value
            elif label == "qq"           : self.cObject.qq0          = value
            elif label == "udd"          : self.cObject.udd0         = value
            elif label == "upp"          : self.cObject.upp0         = value
            elif label == "uss"          : self.cObject.uss0         = value
            elif label == "zcore"        : self.cObject.zcore0       = value
            elif label == "zetad"        : self.cObject.zetad0       = value
            elif label == "zetap"        : self.cObject.zetap0       = value
            elif label == "zetas"        : self.cObject.zetas0       = value
            elif label == "zdn"          : self.cObject.zdn0         = value
            elif label == "zpn"          : self.cObject.zpn0         = value
            elif label == "zsn"          : self.cObject.zsn0         = value
        # . Arrays.
        MNDOParameters_FillBetaUspd ( self.cObject )
        # . AM1/PM3.
        parameters = state.get ( "AM1/PM3 Gaussian Parameters", [] )
        size       = len ( parameters )
        if size > 0: self.FillFN123 ( size, parameters )
        # . Diatomics.
        parameters = state.get ( "Diatomic Parameters", [] )
        size       = len ( parameters )
        if size > 0: self.FillDiatomic ( size, parameters )
        # . PDDG.
        parameters = state.get ( "PDDG Gaussian Parameters", [] )
        size       = len ( parameters )
        if size > 0: self.FillPDDG ( size, parameters )
        # . Convert to atomic units.
        MNDOParameters_ToAtomicUnits ( self.cObject )
        # . Remaining data.
        MNDOParameters_CalculateOneCenterTEIs ( self.cObject )

    def _Allocate ( self ):
        """Allocation."""
        if self.cObject != NULL: MNDOParameters_Deallocate ( &self.cObject )
        self.cObject = MNDOParameters_Allocate ( )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def FillDiatomic ( self, CInteger  nterms, object data ):
        """Fill QDIATOMICFLAGS, diatomica, diatomicx."""
        cdef CInteger  j
        cdef CReal     a, x
        if nterms > 0:
            # . Find maximum atomic number.
            data.sort ( )
            maximumAtomicNumber = data[-1][0]
            # . Allocate and fill object.
            MNDOParameters_AllocateDiatomic ( self.cObject, maximumAtomicNumber + 1 )
            for ( j, x, a ) in data:
                self.cObject.QDIATOMICFLAGS[j] = CTrue
                self.cObject.diatomica0    [j] = a
                self.cObject.diatomicx0    [j] = x

    def FillFN123 ( self, CInteger  nterms, object data ):
        """Fill fn1, fn2 and fn3."""
        cdef CInteger  i
        if nterms > 0:
            MNDOParameters_AllocateFN123 ( self.cObject, nterms )
            for ( i, ( fn1, fn2, fn3 ) ) in enumerate ( data ):
                self.cObject.fn10[i] = fn1
                self.cObject.fn20[i] = fn2
                self.cObject.fn30[i] = fn3

    def FillPDDG ( self, CInteger  nterms, object data ):
        """Fill pddgc and pddge."""
        cdef CInteger  i
        if nterms > 0:
            MNDOParameters_AllocatePDDG ( self.cObject, nterms )
            for ( i, ( c, e ) ) in enumerate ( data ):
                self.cObject.pddgc0[i] = c
                self.cObject.pddge0[i] = e

    @classmethod
    def Raw ( selfClass ):
        """Constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def Uninitialized ( selfClass, atomicNumber, numberOfOrbitals ):
        """Constructor."""
        return selfClass ( atomicNumber, numberOfOrbitals )

    @property
    def chargeSeparations ( self ):
        """Charge separations in atomic units."""
        cdef CInteger  i
        data = []
        if self.cObject != NULL:
            for i from 0 <= i < 6: data.append ( self.cObject.ddp[i] )
        return data

    @property
    def diatomicTerms ( self ):
        """Diatomic terms."""
        cdef CInteger  i
        elements = None
        if ( self.cObject != NULL ) and ( self.cObject.ndiatomic > 0 ):
            elements = []
            for i from 0 <= i < self.cObject.ndiatomic:
                if self.cObject.QDIATOMICFLAGS[i] == CTrue:
                    elements.append ( i )
        return elements

    @property
    def energyBaseLine ( self ):
        """The energy base line."""
        if self.cObject == NULL: return 0.0
        else:                    return ( self.cObject.eheat - self.cObject.eisol )

    @property
    def nuclearCharge ( self ):
        """The nuclear charge."""
        if self.cObject == NULL: return 0.0
        else:                    return self.cObject.zcore

    @property
    def numberOfOrbitals ( self ):
        """The number of orbitals."""
        if self.cObject == NULL: return 0
        else:                    return self.cObject.norbitals

    @property
    def orbitalBasisLabel ( self ):
        """The orbital basis label."""
        if self.cObject == NULL:
            return None
        else:
            labels = []
            if self.cObject.qns != 0: labels.append ( "{:d}s".format ( self.cObject.qns ) )
            if self.cObject.qnp != 0: labels.append ( "{:d}p".format ( self.cObject.qnp ) )
            if self.cObject.qnd != 0: labels.append ( "{:d}d".format ( self.cObject.qnd ) )
            return ( "".join ( labels ).lower ( ) )

    @property
    def shellExponents ( self ):
        """The shell exponents."""
        exponents = []
        if self.cObject != NULL:
            if self.cObject.qns != 0: exponents.append ( self.cObject.zetas )
            if self.cObject.qnp != 0: exponents.append ( self.cObject.zetap )
            if self.cObject.qnd != 0: exponents.append ( self.cObject.zetad )
        return exponents


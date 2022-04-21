"""The System class is the essential class for simulating molecular systems.
It combines information about the atomic composition of the system along with
the ways in which its potential energy can be calculated."""

from  pCore                 import Clone                      , \
                                   logFile                    , \
                                   LogFileActive              , \
                                   Selection                  , \
                                   StorageNode                , \
                                   SummarizableObject
from  pScientific.Arrays    import Array
from  pScientific.Geometry3 import Coordinates3               , \
                                   Vector3
from  pScientific.Symmetry  import PeriodicBoundaryConditions , \
                                   SymmetryParameters         , \
                                   SymmetryParameterGradients
from .Atom                  import AtomContainer
from .Connectivity          import Connectivity
from .EnergyModel           import EnergyModel                , \
                                   EnergyModelPriority        , \
                                   EnergyModelState
from .MMModel               import MMModel
from .NBModel               import NBModel
from .QCModel               import ElectronicState            , \
                                   QCModel
from .RestraintModel        import RestraintModel
from .Sequence              import Sequence

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class System ( SummarizableObject ):
    """Define an atomic or molecular system and its energy models."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "System"
    _attributable.update ( { "_atoms"              : None        ,
                             "_connectivity"       : None        ,
                             "_coordinates3"       : None        ,
                             "_electronicState"    : None        ,
                             "_energyClosures"     : list        ,
                             "_energyModels"       : dict        ,
                             "_freeAtoms"          : None        ,
                             "_mmModel"            : None        ,
                             "_nbModel"            : None        ,
                             "_qcModel"            : None        ,
                             "_restraintModel"     : None        ,
                             "_scratch"            : StorageNode ,
                             "_sequence"           : None        ,
                             "_symmetry"           : None        ,
                             "_symmetryParameters" : None        } )
    _withSections = True

    def __getstate__ ( self ):
        """Set the state of the object."""
        # . Basic data.
        state = {}
        if self.sequence is not None: state["sequence"] = self.sequence.ToMapping ( )
        state.update ( self.atoms.ToMapping        ( ) )
        state.update ( self.connectivity.ToMapping ( ) )
        # . Other data.
        for ( key, value ) in self.__dict__.items ( ):
            if ( key not in ( "_atoms", "_connectivity", "_energyClosures", "_scratch", "_sequence" ) ) and ( value is not None ): state[key] = value
        return state

    def __setstate__ ( self, state ):
        """Set the state of the object."""
        # . Basic data.
        if "sequence" in state:
            sequence = Sequence.FromMapping ( state.pop ( "sequence" ) )
            atoms    = AtomContainer.FromIterable ( sequence.GatherAtoms ( ), attributes = state.pop ( "atoms", None ) )
        else:
            atoms    = AtomContainer.FromMapping ( state.pop ( "atoms", None ) )
            sequence = None
        connectivity = Connectivity.FromAtoms ( atoms, bonds = state.pop ( "bonds", None ) )
        self.__dict__["_atoms"       ] = atoms
        self.__dict__["_connectivity"] = connectivity
        self.__dict__["_sequence"    ] = sequence
        # . Other data.
        for ( key, value ) in state.items ( ):
            if value is not None:
                self.__dict__[key] = value
                if isinstance ( value, EnergyModelState ): value.target = self # . As target is unpicklable.
        self._UpdateEnergyClosures ( )

    def _AddEnergyModel ( self, key, value, **options ):
        """Add an energy model."""
        if value is not None:
            priority   = options.pop ( "priority"  , EnergyModelPriority.NullModel )
            valueClass = options.pop ( "valueClass", EnergyModel )
            if isinstance ( value, valueClass ):
                if hasattr ( value, "BuildModel" ): value.BuildModel ( self, **options )
                self._energyModels[key] = priority
                self.__dict__     [key] = value
            else: raise TypeError ( "Invalid \"{:s}\" argument.".format ( key ) )

    def _PopEnergyModel ( self, key, priority ):
        """Pop an energy model."""
        # . All incompatible models of higher priority are also removed.
        toDiscard = [ ( priority, key ) ]
        if priority > EnergyModelPriority.NullModel:
            toDiscard.extend ( [ ( priority0, key0 ) for ( key0, priority0 ) in self._energyModels.items ( ) if priority0 > priority ] )
        for ( priority0, key0 ) in sorted ( toDiscard, reverse = True ):
            self._energyModels.pop ( key0, None )
            if key0 in self.__dict__:
                value = self.__dict__[key0]
                if ( value is not None ) and hasattr ( value, "UnbuildModel" ): value.UnbuildModel ( self )
                self.__dict__[key0] = None
        self._energyClosures = []

    def _SetHandlerCoordinates3 ( self, value ):
        """Coordinates3 set handler."""
        if not ( ( value is None ) or ( isinstance ( value, Coordinates3 ) and ( value.rows == len ( self.atoms ) ) ) ):
            raise TypeError ( "Invalid coordinates3 argument." )
        self.__dict__["_coordinates3"] = value

    def _SetHandlerElectronicState ( self, value ):
        """Electronic state set handler."""
        if not ( ( value is not None ) and isinstance ( value, ElectronicState ) ):
            raise TypeError ( "Invalid electronic state argument." )
        self.__dict__["_electronicState"] = value
        if self.qcModel is not None: self.qcModel.ModifyElectronicState ( self )

    def _SetHandlerFreeAtoms ( self, value ):
        """Free atoms set handler."""
        if not ( ( value is None ) or ( isinstance ( value, Selection ) and ( max ( value ) < len ( self.atoms ) ) ) ):
            raise TypeError ( "Invalid free atoms argument." )
        self.scratch.Clear ( )
        for eKey in self._energyModels.keys ( ):
            model = self.__dict__[eKey]
            if isinstance ( model, EnergyModel ):
                label = model.__class__._stateName
                if label is None: state = None
                else:             state = getattr ( self, label, None )
                if ( state is not None ) and hasattr ( state, "UnfixAtoms" ):
                    state.UnfixAtoms ( )
                    if value is not None:
                        state.FixAtoms ( value )
        self.__dict__["_freeAtoms"] = value
        if value is not None: value.label = "Free Atoms"

    def _SetHandlerMMModel ( self, value, log = None ):
        """MM model set handler."""
        self.AddEnergyModel ( "_mmModel", value, log = log, priority = EnergyModelPriority.MMModel, valueClass = MMModel )

    def _SetHandlerNBModel ( self, value, assignQCMMModels = True ):
        """NB model set handler."""
        self.AddEnergyModel ( "_nbModel", value, assignQCMMModels = assignQCMMModels, priority = EnergyModelPriority.NBModel, valueClass = NBModel )

    def _SetHandlerQCModel ( self, value, qcSelection = None ):
        """QC model set handler."""
        self.AddEnergyModel ( "_qcModel", value, priority = EnergyModelPriority.QCModel, qcSelection = qcSelection, valueClass = QCModel )

    def _SetHandlerRestraintModel ( self, value ):
        """Restraint model set handler."""
        self.AddEnergyModel ( "_restraintModel", value, priority = EnergyModelPriority.Restraint, valueClass = RestraintModel )

    def _SetHandlerSymmetry ( self, value ):
        """Symmetry set handler."""
        if not ( ( value is None ) or isinstance ( value, PeriodicBoundaryConditions ) ): raise TypeError ( "Invalid symmetry argument." )
        if ( value is not None ) and ( self.symmetryParameters is not None ): value.CheckSymmetryParameters ( self.symmetryParameters )
        self.__dict__["_symmetry"          ] = value
        self.__dict__["_symmetryParameters"] = None
        self.scratch.Clear ( )

    def _SetHandlerSymmetryParameters ( self, value ):
        """Symmetry parameters set handler."""
        if not ( ( ( self.symmetry is     None ) and ( value is None ) ) or \
                 ( ( self.symmetry is not None ) and isinstance ( value, SymmetryParameters ) ) ):
            raise TypeError ( "Invalid symmetry parameters argument." )
        if value is not None: self.symmetry.CheckSymmetryParameters ( value )
        self.__dict__["_symmetryParameters"] = value

    def _UpdateEnergyClosures ( self ):
        """Update the list of energy closures."""
        closures = []
        for key in self._energyModels.keys ( ):
            closures.extend ( self.__dict__[key].EnergyClosures ( self ) )
        self._energyClosures = [ closure for ( priority, closure, label ) in sorted ( closures, key = lambda x: x[0] ) ]

    def AddEnergyModel ( self, key, value, **options ):
        """Add an energy model."""
        self._PopEnergyModel ( key, options.get ( "priority", EnergyModelPriority.NullModel ) )
        self._AddEnergyModel ( key, value, **options )
        self._UpdateEnergyClosures ( )

    def AtomicCharges ( self, qcChargeModel = None ):
        """Atomic charges of all atoms."""
        charges = None
        if self.mmModel is not None:
            charges   = self.mmModel.AtomicCharges ( self )
        if self.qcModel is not None:
            qcCharges = self.qcModel.AtomicCharges ( self, chargeModel = qcChargeModel )
            if charges is None:
                charges = qcCharges
            else:
                for ( s, q ) in zip ( self.qcState.qcAtoms, qcCharges ): charges[s] += q
        return charges

    def BondsFromCoordinates3 ( self, coordinates = None, radii = None, safety = 0.45 ):
        """Estimate bonds from coordinates using a distance search."""
        if coordinates is None:
            if hasattr ( self, "coordinates3" ): coordinates = self.coordinates3
        if radii is None: radii = Array.FromIterable ( [ atom.covalentRadius for atom in self.atoms ] )
        if self.connectivity is None: self.connectivity = Connectivity.FromAtoms ( self.atoms )
        self.connectivity.BondsFromCoordinates ( coordinates, radii, safety )

    # . Define methods only for extra options.
    def DefineMMModel ( self, value, log = logFile ):
        self._SetHandlerMMModel ( value, log = log )

    def DefineNBModel ( self, value, assignQCMMModels = True ):
        self._SetHandlerNBModel ( value, assignQCMMModels = assignQCMMModels )

    def DefineQCModel ( self, value, qcSelection = None ):
        self._SetHandlerQCModel ( value, qcSelection = qcSelection )

    def DefineRestraintModel ( self, value ):
        self._SetHandlerRestraintModel ( value )

    def DipoleMoment ( self, center = None ):
        """The dipole moment in Debyes."""
        dipole = None
        for attribute in ( "mmModel", "qcModel" ):
            model = getattr ( self, attribute, None )
            if model is not None: dipole = model.DipoleMoment ( self, center, dipole = dipole )
        return dipole

    def Energy ( self, doGradients = False, log = logFile ):
        """Calculate the energy and, optionally, the gradients for a system."""
        self.EnergyInitialize ( doGradients, log )
        for closure in self._energyClosures: closure ( )
        return self.EnergyFinalize ( )

    def EnergyFinalize ( self ):
        """Energy finalization."""
        doGradients     = self.scratch.doGradients
        log             = self.scratch.log
        energies        = self.scratch.energyTerms
        potentialEnergy = sum ( energies.values ( ) )
        if doGradients and ( self.freeAtoms is not None ): self.scratch.iFixedAtomGradients3.Set ( 0.0 )
        if log is not None:
            if doGradients: rmsGradient = "{:16.4f}".format ( self.scratch.gradients3.RootMeanSquare ( ) )
            else:           rmsGradient = "None"
            items = [ ( "Potential Energy", "{:16.4f}".format ( potentialEnergy ) ) ,
                       ( "RMS Gradient"    , rmsGradient                           ) ]
            items.extend ( [ ( key, "{:16.4f}".format ( energies[key] ) ) for key in sorted ( energies.keys ( ) ) ] )
            log.SummaryOfItems ( items, order = False, title = "Summary of Energy Terms" )
        energies["Potential Energy"] = potentialEnergy
        return potentialEnergy

    def EnergyInitialize ( self, doGradients, log ):
        """Energy initialization."""
        scratch = self.scratch
        if not LogFileActive ( log ): log = None
        scratch.doGradients = doGradients
        scratch.energyTerms = {}
        scratch.log         = log
        if doGradients:
            if not hasattr ( scratch, "gradients3" ):
                grd3               = Coordinates3.WithExtent ( len ( self.atoms ) )
                scratch.gradients3 = grd3
                if self.freeAtoms is not None:
                    fixedAtoms                   = self.freeAtoms.Complement ( upperBound = len ( self.atoms ) )
                    scratch.iFixedAtomGradients3 = grd3.RowIterator ( selection = fixedAtoms )
            scratch.gradients3.Set ( 0.0 )
            if self.symmetry is not None:
                if not hasattr ( scratch, "symmetryParameterGradients" ):
                    scratch.symmetryParameterGradients = SymmetryParameterGradients ( )
                scratch.symmetryParameterGradients.Clear ( )

    @classmethod
    def FromAtoms ( selfClass, atoms, bonds = None, withSequence = False ):
        """Constructor given a list of atoms."""
        self              = selfClass ( )
        self.__dict__["_atoms"       ] = AtomContainer.FromIterable ( atoms )
        self.__dict__["_connectivity"] = Connectivity.FromAtoms     ( self.atoms, bonds = bonds )
        if withSequence: self.__dict__["_sequence"] = Sequence.FromAtoms ( self.atoms )
        return self

    @classmethod
    def FromConnectivity ( selfClass, connectivity, withSequence = False ):
        """Constructor given a connectivity."""
        self = selfClass ( )
        self.__dict__["_atoms"       ] = AtomContainer.FromIterable ( connectivity.atoms )
        self.__dict__["_connectivity"] = connectivity
        if withSequence: self.__dict__["_sequence"] = Sequence.FromAtoms ( self.atoms )
        return self

    @classmethod
    def FromSequence ( selfClass, sequence, bonds = None ):
        """Constructor given a sequence."""
        atoms             = sequence.GatherAtoms ( )
        self              = selfClass ( )
        self.__dict__["_atoms"       ] = AtomContainer.FromIterable ( atoms )
        self.__dict__["_connectivity"] = Connectivity.FromAtoms     ( self.atoms, bonds = bonds )
        self.__dict__["_sequence"    ] = sequence
        return self

    #------------------------------------------------------------------------------------#
    # . Merging and pruning require Atoms to be done first as these store                #
    # . "Atom Container", "Atom Mapping" and "Index Increments" in information.          #
    # . Connectivity and Sequence use "Atom Container" and "Atom Mapping", respectively. #
    # . Many of the other entries use "Index Increments".                                #
    # . This is not explicitly done here as "Atoms" is alphabetically the first key.     #
    #------------------------------------------------------------------------------------#
    @classmethod
    def Merge ( selfClass, entries, information = {} ):
        """Merging."""
        # . Gather all keys to treat.
        keys = set ( )
        for entry in entries:
            for ( key, value ) in entry.__dict__.items ( ):
                if value is not None: keys.add ( key )
        # . Create object.
        new                   = selfClass ( )
        new.label             = "Merged {:s}".format ( new.__class__.__name__ )
        information["Target"] = new
        for key in sorted ( keys ):
            mergeEntries = [ entry.__dict__.get ( key, None ) for entry in entries ]
            notNone      = next ( ( entry for entry in mergeEntries if entry is not None ), None )
            if ( notNone is not None ) and hasattr ( notNone.__class__, "Merge" ):
                newEntry = notNone.__class__.Merge ( mergeEntries, information = information )
                if newEntry is not None: new.__dict__[key] = newEntry
        if ( new.freeAtoms is not None ) and ( new.mmState is not None ): new.mmState.FixAtoms ( new.freeAtoms )
        for ( key, priority ) in ( ( "_mmModel"       , EnergyModelPriority.MMModel   ) ,
                                   ( "_restraintModel", EnergyModelPriority.Restraint ) ):
            newEntry = new.__dict__.get ( key, None )
            if newEntry is not None: new._energyModels[key] = priority
        new._UpdateEnergyClosures ( )
        return new

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        new                   = self.__class__ ( )
        new.label             = "Pruned {:s}".format ( self.__class__.__name__ )
        information["Target"] = new
        for key in sorted ( self.__dict__.keys ( ) ):
            entry = self.__dict__.get ( key, None )
            if ( entry is not None ) and hasattr ( entry, "Prune" ):
                newEntry = entry.Prune ( selection, information = information )
                if ( newEntry is not None ): new.__dict__[key] = newEntry
        if ( new.freeAtoms is not None ) and ( new.mmState is not None ): new.mmState.FixAtoms ( new.freeAtoms )
        for ( key, priority ) in ( ( "_mmModel"       , EnergyModelPriority.MMModel   ) ,
                                   ( "_restraintModel", EnergyModelPriority.Restraint ) ):
            newEntry = new.__dict__.get ( key, None )
            if newEntry is not None: new._energyModels[key] = priority
        new._UpdateEnergyClosures ( )
        return new

    def PruneToQCRegion ( self, withAtomPaths = False ):
        """
        Return a minimal system consisting of the QC region, including boundary atoms.
        Nothing is done if no QC model is defined.
        """
        new = None
        if self.qcModel is not None:
            new              = self.__class__.FromAtoms ( self.qcState.atomicNumbers )
            new.coordinates3 = Coordinates3.WithExtent ( len ( new.atoms ) )
            self.qcModel.GatherQCCoordinates3 ( self.qcState, self.coordinates3, new.coordinates3 )
            if self.label is None: new.label = "QC Region"
            else:                  new.label = self.label + " - QC Region"
            if withAtomPaths:
                for ( i, s ) in enumerate ( self.qcState.qcAtoms ):
                    new.atoms[i].label = self.atoms[s].path
        return new

    def SummaryItems ( self ):
        """Summary sections."""
        # . Basic items then energy models and their states.
        items = []
        for name in ( "_atoms", "_freeAtoms", "_connectivity", "_electronicState", "_sequence", "_symmetry" ):
            item = getattr ( self, name, None )
            if item is not None: items.extend ( item.SummaryItems ( ) )
        for ( _, name ) in sorted ( [ ( priority, name ) for ( name, priority ) in self._energyModels.items ( ) ] ):
            item = getattr ( self, name, None )
            if item is not None:
                items.extend ( item.SummaryItems ( ) )
                if item.__class__._stateName is not None:
                    state = getattr ( self, item.__class__._stateName, None )
                    if state is not None: items.extend ( state.SummaryItems ( ) )
        # . Add code for remaining (sorted) attributes which are not explicitly included here?
        return items

    # . Properties.
    # . Getters.
    @property
    def atoms              ( self ): return self.__dict__["_atoms"             ]
    @property
    def connectivity       ( self ): return self.__dict__["_connectivity"      ]
    @property
    def coordinates3       ( self ): return self.__dict__["_coordinates3"      ]
    @property
    def electronicState    ( self ): return self.__dict__["_electronicState"   ]
    @property
    def freeAtoms          ( self ): return self.__dict__["_freeAtoms"         ]
    @property
    def mmModel            ( self ): return self.__dict__["_mmModel"           ]
    @property
    def nbModel            ( self ): return self.__dict__["_nbModel"           ]
    @property
    def qcModel            ( self ): return self.__dict__["_qcModel"           ]
    @property
    def restraintModel     ( self ): return self.__dict__["_restraintModel"    ]
    @property
    def scratch            ( self ): return self.__dict__["_scratch"           ]
    @property
    def sequence           ( self ): return self.__dict__["_sequence"          ]
    @property
    def symmetry           ( self ): return self.__dict__["_symmetry"          ]
    @property
    def symmetryParameters ( self ): return self.__dict__["_symmetryParameters"]

    # . Special getters.
    @property
    def energyModelLabel ( self ):
        """An energy model label."""
        items = []
        if self.qcModel is not None: items.append ( self.qcModel.__class__._classLabel )
        if self.mmModel is not None: items.append ( self.mmModel.__class__._classLabel )
        return "/".join ( items )

    # . Setters.
    @coordinates3.setter
    def coordinates3       ( self, value ): self._SetHandlerCoordinates3       ( value )
    @electronicState.setter
    def electronicState    ( self, value ): self._SetHandlerElectronicState    ( value )
    @freeAtoms.setter
    def freeAtoms          ( self, value ): self._SetHandlerFreeAtoms          ( value )
    @mmModel.setter
    def mmModel            ( self, value ): self._SetHandlerMMModel            ( value )
    @nbModel.setter
    def nbModel            ( self, value ): self._SetHandlerNBModel            ( value )
    @qcModel.setter
    def qcModel            ( self, value ): self._SetHandlerQCModel            ( value )
    @restraintModel.setter
    def restraintModel     ( self, value ): self._SetHandlerRestraintModel     ( value )
    @symmetry.setter
    def symmetry           ( self, value ): self._SetHandlerSymmetry           ( value )
    @symmetryParameters.setter
    def symmetryParameters ( self, value ): self._SetHandlerSymmetryParameters ( value )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

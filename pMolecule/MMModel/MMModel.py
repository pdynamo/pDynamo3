"""Defines MM energy model classes."""

import itertools, math, os, os.path

from   pCore                 import Clone              , \
                                    DataType           , \
                                    logFile            , \
                                    LogFileActive      , \
                                    Selection          , \
                                    SelectionContainer
from   pScientific           import Units
from   pScientific.Arrays    import Array
from   pScientific.Geometry3 import Vector3
from  .DYFFUtilities         import DYFFCosineAngleParameters      , \
                                    DYFFCosineDihedralParameters   , \
                                    DYFFCosineOutOfPlaneParameters , \
                                    DYFFHarmonicBondParameters
from  .LJParameterContainer  import LJForm
from  .MMAtomTyper           import MMAtomTyper
from  .MMModelError          import MMModelError
from  .MMParameterAssigner   import MMParameterAssigner
from ..ConnectivityUtilities import DetermineAtomGeometry          , \
                                    DetermineAtomOxidationState
from ..EnergyModel           import EnergyModel                    , \
                                    EnergyModelState

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The tolerance for ensuring that the total MM active charge is integral.
_DefaultChargeTolerance = 1.0e-4

# . Paths.
_ForceFieldPath = "forceFields"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelState ( EnergyModelState ):
    """An MM energy model state."""

    # . None values for mmAtoms and pureMMAtoms implies all atoms are MM. 
    _attributable = dict ( EnergyModelState._attributable )
    _attributable.update ( { "atomTypeIndices" : None ,
                             "atomTypes"       : None ,
                             "charges"         : None ,
                             "exclusions"      : None ,
                             "interactions14"  : None ,
                             "ljParameters"    : None ,
                             "ljParameters14"  : None ,
                             "ljTypeIndices"   : None , # . Add a duplicate for 1-4s?
                             "mmAtoms"         : None ,
                             "mmTerms"         : list ,
                             "pureMMAtoms"     : None } )

    def ActivateAllAtoms ( self ):
        """Activate all atoms."""
        self.mmAtoms     = None
        self.pureMMAtoms = None
        for mmTerm in self.mmTerms:
            if hasattr ( mmTerm, "ActivateTerms" ): mmTerm.ActivateTerms ( )

    def ActiveAtomTotalCharge ( self ):
        """The total charge of active atoms."""
        totalCharge = 0.0
        if self.charges is not None:
            if self.mmAtoms is None:
                totalCharge = sum ( self.charges )
            else:
                for i in self.mmAtoms:
                    totalCharge += self.charges[i]
        return totalCharge

    def CheckActiveAtomTotalCharge ( self, tolerance = _DefaultChargeTolerance ):
        """Throw an error if the total active charge is not integral."""
        totalCharge    = self.ActiveAtomTotalCharge ( )
        nearestInteger = round ( totalCharge )
        if math.fabs ( totalCharge - nearestInteger ) > tolerance:
            raise MMModelError ( "Total active MM charge is neither integral nor zero: {:.6f}.".format ( totalCharge ) )

    def DeactivateQCAtoms ( self, pureQCAtoms, baAtoms ):
        """Deactivate MM atoms and terms due to the presence of QC atoms."""
        if len ( pureQCAtoms ) > 0:
            # . To be improved.
            self.mmAtoms = Selection.FromIterable ( set ( range ( len ( self.charges ) ) ) - set ( pureQCAtoms ) )
            if len ( baAtoms ) > 0: self.pureMMAtoms = Selection.FromIterable ( set ( self.mmAtoms ) - set ( baAtoms ) )
            else:                   self.pureMMAtoms = self.mmAtoms
            # . Deactivate terms.
            for mmTerm in self.mmTerms:
                if hasattr ( mmTerm, "DeactivateQCAtomTerms" ):
                    mmTerm.DeactivateQCAtomTerms ( self.mmAtoms, self.pureMMAtoms )

    def FixAtoms ( self, freeAtoms ):
        """Fix atoms."""
        for mmTerm in self.mmTerms:
            if hasattr ( mmTerm, "DeactivateFixedAtomTerms" ):
                mmTerm.DeactivateFixedAtomTerms ( freeAtoms )

    def Get12Exclusions ( self ):
        """Get a 1-2 exclusion list from the MM terms."""
        exclusions = None
        indices    = []
        for mmTerm in self.mmTerms:
            if hasattr ( mmTerm, "Get12Indices" ): indices.extend ( mmTerm.Get12Indices ( ) )
        if len ( indices ) > 0: exclusions = SelfPairList.FromIndexPairs ( indices )
        return exclusions

    def IdentifyBoundaryAtoms ( self, selection ):
        """Identify atoms in |selection| that are covalently bound to atoms outside of |selection|."""
        results = {}
        for mmTerm in self.mmTerms:
            if hasattr ( mmTerm, "IdentifyBoundaryAtoms" ):
                mmTerm.IdentifyBoundaryAtoms ( selection, results )
        return results

    def MakeMMIsolates ( self ):
        """Make MM isolates."""
        # . QC atoms are kept together.
        # . Isolates with all free atoms are flagged.
        n          = len ( self.target.atoms )
        if self.exclusions is None: isolates = SelectionContainer.FromCapacity ( n )
        else:                       isolates = self.exclusions.GetConnectedComponents ( upperBound = n )
        if self.target.qcModel is not None:
            toFuse = isolates.MakeMembershipFlags ( self.target.qcState.qcAtoms )
            isolates.FuseItems ( toFuse )
        if self.target.freeAtoms is not None:
            freeIsolates = isolates.MakeMembershipFlags ( self.target.freeAtoms, andTest = True )
            self.target.scratch.freeMMIsolates = freeIsolates
        self.target.scratch.mmIsolates = isolates

    # . Merging and pruning are done piecewise due to the diversity of entries.
    # . The new models are delivered with all atoms activated.
    # . Attribute names are handled implicitly.
    @classmethod
    def Merge ( selfClass, entries, information = {} ):
        """Merging."""
        new = None
        if ( None not in entries ):
            targets = [ entry.target for entry in entries ]
            if any ( x is None for x in targets ): raise MMModelError ( "Incompatible MM states for merging." )
            # . Create object.
            new        = selfClass ( )
            new.target = information.get ( "Target", None )
            # . LJ parameter merging produces "LJ Type Increments" and "LJ Type Mapping" (same for both sets of parameters).
            for attribute in ( "charges", "exclusions", "interactions14", "ljParameters", "ljParameters14", "ljTypeIndices" ):
                mergeEntries = [ getattr ( entry, attribute, None ) for entry in entries ]
                notNone      = next ( ( entry for entry in mergeEntries if entry is not None ), None )
                if ( notNone is not None ) and hasattr ( notNone.__class__, "Merge" ):
                    newEntry = notNone.__class__.Merge ( mergeEntries, information = information )
                    if newEntry is not None: setattr ( new, attribute, newEntry )
            # . MM terms.
            mergeEntries = [ getattr ( entry, "mmTerms", None ) for entry in entries ]
            states       = [ {} if entry is None else { value.label : value for value in entry } for entry in mergeEntries ]
            keys         = set ( itertools.chain.from_iterable ( [ state.keys ( ) for state in states ] ) )
            if len ( keys ) > 0:
                for key in sorted ( keys ):
                    toMerge = [ state.get ( key, None ) for state in states ]
                    notNone = next ( ( entry for entry in toMerge if entry is not None ), None )
                    if ( notNone is not None ) and hasattr ( notNone.__class__, "Merge" ):
                        newEntry = notNone.__class__.Merge ( toMerge, information = information )
                        if newEntry is not None: new.mmTerms.append ( newEntry )
            # . Atom types.
            oldTypes = []
            for entry in entries: oldTypes.extend ( entry.atomTypes )
            newTypes      = sorted ( set ( oldTypes ) )
            new.atomTypes = newTypes
            # . Atom type indices.
            oldToNew = [ newTypes.index ( o ) for o in oldTypes ]
            indices  = Array.WithExtent ( len ( new.charges ), dataType = DataType.Integer )
            n        = 0
            for entry in entries:
                value = entry.atomTypeIndices
                for ( i, v ) in enumerate ( value ): indices[i+n] = oldToNew[v]
                n += len ( value )
            new.atomTypeIndices = indices
            # . LJ type indices - modify existing array if extra information found.
            increments = information.get ( "LJ Type Increments" , None )
            mapping    = information.get ( "LJ Type Mapping"    , None )
            if ( increments is not None ) and ( mapping is not None ):
                indices = new.ljTypeIndices
                n       = 0
                for ( increment, entry ) in zip ( increments, entries ):
                        value = entry.ljTypeIndices
                        for ( i, v ) in enumerate ( value ):
                            indices[i+n] = mapping[v + increment]
                        n += len ( value )
            # . Finish up.
            new.ActivateAllAtoms ( )
        return new

    # . Atom types and LJ parameters remain unchanged implying the atom and LJ type indices do too.
    def Prune ( self, selection, information = {} ):
        """Pruning."""
        new        = self.__class__ ( )
        new.target = information.get ( "Target", None )
        entry      = getattr ( self, "atomTypes", None )
        if entry is not None: new.atomTypes = Clone ( entry )
        for attribute in ( "atomTypeIndices" ,
                           "charges"         ,
                           "exclusions"      ,
                           "interactions14"  ,
                           "ljParameters"    ,
                           "ljParameters14"  ,
                           "ljTypeIndices"   ):
            entry =  getattr ( self, attribute, None )
            if entry is not None: setattr ( new, attribute, entry.Prune ( selection, information = information ) )
        for entry in self.mmTerms:
            if entry is not None:
                newEntry = entry.Prune ( selection, information = information )
                if newEntry is not None: new.mmTerms.append ( newEntry )
        new.ActivateAllAtoms ( )
        return new

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.target is not None:
            nT = len ( self.charges )
            if self.mmAtoms     is None: nM = nT
            else:                        nM = len ( self.mmAtoms     )
            if self.pureMMAtoms is None: nP = nT
            else:                        nP = len ( self.pureMMAtoms )
            items.extend ( [ ( "Number of MM Atoms"     , "{:d}".format ( nT                     ) ) ,
                             ( "Number of MM Atom Types", "{:d}".format ( len ( self.atomTypes ) ) ) ] )
            if nM < nT:
                items.extend ( [ ( "Number of Pure MM Atoms" , "{:d}".format ( nP      ) ) ,
                                 ( "Number of Other MM Atoms", "{:d}".format ( nM - nP ) ) ] )
            items.append ( ( "Total MM Charge", "{:.2f}".format ( self.ActiveAtomTotalCharge ( ) ) ) )
            for attribute in ( "exclusions", "interactions14", "ljParameters", "ljParameters14" ):
                entry = getattr ( self, attribute, None )
                if entry is not None: items.extend ( entry.SummaryItems ( ) )
            for mmTerm in self.mmTerms: items.extend ( mmTerm.SummaryItems ( ) )
        return items

    def UnfixAtoms ( self ):
        """Unfix atoms."""
        for mmTerm in self.mmTerms:
            if hasattr ( mmTerm, "ActivateTerms" ): mmTerm.ActivateTerms ( )
        if ( self.mmAtoms is not None ) and ( self.pureMMAtoms is not None ):
            for mmTerm in self.mmTerms:
                if hasattr ( mmTerm, "DeactivateQCAtomTerms" ):
                    mmTerm.DeactivateQCAtomTerms ( self.mmAtoms, self.pureMMAtoms )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModel ( EnergyModel ):
    """An MM energy model."""

    _attributable         = dict ( EnergyModel._attributable )
    _classLabel           = "MM Model"
    _electrostaticScale14 = 1.0
    _forceField           = None
    _lennardJonesScale14  = 1.0
    _lennardJonesStyle    = LJForm.OPLS
    _stateName            = "mmState"
    _stateObject          = MMModelState
    _summarizable         = dict ( EnergyModel._summarizable )
    _attributable.update ( { "electrostaticScale14" : None ,
                             "forceField"           : None ,
                             "lennardJonesScale14"  : None ,
                             "lennardJonesStyle"    : None ,
                             "parameterSet"         : None ,
                             "parameterSetPath"     : None } )
    _summarizable.update ( { "electrostaticScale14" : ( "Electrostatic 1-4 Scaling", "{:.3f}" ) ,
                             "forceField          " :   "Force Field"                           ,
                             "lennardJonesScale14 " : ( "Lennard-Jones 1-4 Scaling", "{:.3f}" ) ,
                             "lennardJonesStyle   " :   "Lennard-Jones Style"                   ,
                             "parameterSet"         :   "Parameter Set"                       } )


    def _Initialize ( self ):
        """Initialization."""
        super ( MMModel, self )._Initialize ( )
        self.electrostaticScale14 = self.__class__._electrostaticScale14
        self.forceField           = self.__class__._forceField
        self.lennardJonesScale14  = self.__class__._lennardJonesScale14
        self.lennardJonesStyle    = self.__class__._lennardJonesStyle

    def AtomicCharges ( self, target ):
        """The atomic charges of active atoms."""
        charges = Array.WithExtent ( len ( target.atoms ) )
        state   = getattr ( target, self.__class__._stateName )
        if state.mmAtoms is None:
            state.charges.CopyTo ( charges )
        else:
            charges.Set ( 0.0 )
            for s in state.mmAtoms: charges[s] = state.charges[s]
        return charges

    def ActiveAtomTotalCharge ( self ):
        """The total charge of active atoms."""
        totalCharge = 0.0
        if self.charges is not None:
            if self.mmAtoms is None:
                totalCharge = sum ( self.charges )
            else:
                for i in self.mmAtoms:
                    totalCharge += self.charges[i]
        return totalCharge


    def BuildModel ( self, target, log = logFile ):
        """Build the model given the parameter set and system information."""
        state = super ( MMModel, self ).BuildModel ( target )
        if self.parameterSetPath is not None:
            # . Atom typing.
            typer = MMAtomTyper.FromPath ( self.parameterSetPath )
            ( atomTypes, atomCharges ) = typer.TypeAtoms ( target.connectivity, target.sequence, log )
            # . Parameter assignment.
            assigner = MMParameterAssigner.FromPath ( self.parameterSetPath                                               ,
                                                      lennardJonesScale14 = self.lennardJonesScale14                      ,
                                                      parameterFactories  = self.ParameterFactories ( typer.mmAtomTypes ) )
            assigner.AssignParameters ( target.connectivity, atomTypes, atomCharges, state, log )
            if target.freeAtoms is not None: state.FixAtoms ( target.freeAtoms )
            state.CheckActiveAtomTotalCharge ( )

    def CompleteConnectivity ( self, target ):
        """Complete the connectivity of target as much as possible."""
        connectivity = getattr ( target, "connectivity", None )
        if connectivity is not None:
            state = getattr ( target, self.__class__._stateName )
            bonds = [] # . Initially a list of indices but eventually could include bond types if this information is known.
            for mmTerm in state.mmTerms:
                if hasattr ( mmTerm, "GetBonds" ): bonds.extend ( mmTerm.GetBonds ( ) )
            connectivity.BondsFromIterable ( bonds ) # . Empty lists are assumed to correspond to no bonds at all.
            connectivity.CompleteConnectivity ( )

    def DipoleMoment ( self, target, center, dipole = None ):
        """The dipole moment in Debyes."""
        state = getattr ( target, self.__class__._stateName )
        if dipole is None: dipole = Vector3.Null ( )
        if state.charges is not None:
            coordinates3 = target.coordinates3
            units        = ( Units.Dipole_Atomic_Units_To_Debyes * Units.Length_Angstroms_To_Bohrs )
            if center is None:
                xC = yC = zC = 0.0
            else:
                xC = center[0] ; yC = center[1] ; zC = center[2]
            if state.mmAtoms is None: indices = range ( len ( state.charges ) )
            else:                     indices = state.mmAtoms
            for i in indices:
                q = ( state.charges[i] * units )
                dipole[0] += q * ( coordinates3[i,0] - xC )
                dipole[1] += q * ( coordinates3[i,1] - yC )
                dipole[2] += q * ( coordinates3[i,2] - zC )
        return dipole

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        closures = []
        state    = getattr ( target, self.__class__._stateName )
        for mmTerm in state.mmTerms:
            closures.extend ( mmTerm.EnergyClosures ( target ) )
        return closures

    @classmethod
    def Merge ( selfClass, entries, information = {} ):
        """Merging."""
        # . Parameter set and path only need to be identical if the model has not yet been built ...
        new = None
        if ( None not in entries ):
            isOK = all ( x == entries[0].electrostaticScale14 for x in [ entry.electrostaticScale14 for entry in entries ] ) and \
                   all ( x == entries[0].forceField           for x in [ entry.forceField           for entry in entries ] ) and \
                   all ( x == entries[0].lennardJonesScale14  for x in [ entry.lennardJonesScale14  for entry in entries ] ) and \
                   all ( x == entries[0].lennardJonesStyle    for x in [ entry.lennardJonesStyle    for entry in entries ] )
# . It is assumed if the force field is the same, then everything else is OK.
#                   all ( x == entries[0].parameterSet         for x in [ entry.parameterSet         for entry in entries ] ) and \
#                   all ( x == entries[0].parameterSetPath     for x in [ entry.parameterSetPath     for entry in entries ] )
            if not isOK: raise MMModelError ( "Incompatible MM models for merging." )
            new = Clone ( entries[0] )
        return new

    def ParameterFactories ( self, atomTypes ):
        """Return a dictionary of parameter factories for the force field."""
        return {}

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        return Clone ( self )

    @classmethod
    def WithParameterSet ( selfClass, parameterSet, path = None ):
        """Constructor with a parameter set and optional path where the parameter set is to be found."""
        self = selfClass ( )
        if path is None:
            path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), _ForceFieldPath, self.forceField.lower () )
        fullPath = os.path.join ( path, parameterSet )
        if not ( os.path.exists ( fullPath ) and os.path.isdir ( fullPath ) ):
            raise MMModelError ( "Missing or invalid force field parameter set." )
        else:
            self.parameterSet     = parameterSet
            self.parameterSetPath = fullPath
        return self

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelAMBER ( MMModel ):
    """Define an AMBER MM energy model."""

    _classLabel           = "AMBER MM Model"
    _electrostaticScale14 = ( 1.0 / 1.2 )
    _forceField           = "AMBER"
    _lennardJonesScale14  = ( 1.0 / 2.0 )
    _lennardJonesStyle    = LJForm.Amber

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelCHARMM ( MMModel ):
    """Define a CHARMM MM energy model."""

    _classLabel        = "CHARMM MM Model"
    _forceField        = "CHARMM"
    _lennardJonesStyle = LJForm.Amber

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelDYFF ( MMModel ):
    """Define an DYFF MM energy model."""

    _classLabel           = "DYFF MM Model"
    _electrostaticScale14 = 0.5
    _forceField           = "DYFF"
    _lennardJonesScale14  = 0.5

    def BuildModel ( self, target, log = logFile ):
        """Build the model given the parameter set and system information."""
        # . Make sure extra connectivity data exists for matching.
        DetermineAtomGeometry       ( target.connectivity )
        DetermineAtomOxidationState ( target.connectivity )
        return super ( MMModelDYFF, self ).BuildModel ( target, log = log )

    def ParameterFactories ( self, atomTypes ):
        """DYFF parameter factories."""
        return { "cosineAngleParameters"      : DYFFCosineAngleParameters      ( atomTypes ) ,
                 "cosineDihedralParameters"   : DYFFCosineDihedralParameters   ( atomTypes ) ,
                 "cosineOutOfPlaneParameters" : DYFFCosineOutOfPlaneParameters ( atomTypes ) ,
                 "harmonicBondParameters"     : DYFFHarmonicBondParameters     ( atomTypes ) }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMModelOPLS ( MMModel ):
    """Define an OPLS MM energy model."""

    _classLabel           = "OPLS MM Model"
    _electrostaticScale14 = 0.5
    _forceField           = "OPLS"
    _lennardJonesScale14  = 0.5

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""Classes for a charge and spin restraint energy model."""

from   pCore                      import AttributableObject    , \
                                         Clone                 , \
                                         DataType              , \
                                         logFile               , \
                                         LogFileActive         , \
                                         RawObjectConstructor
from   pScientific                import Units
from   pScientific.Arrays         import Array                 , \
                                         IntegerArray1D        , \
                                         RealArray1D           , \
                                         StorageType
from  .LoewdinMultipoleEvaluator  import LoewdinMultipoleEvaluator
from  .MNDOMultipoleEvaluator     import MNDOMultipoleEvaluator
from  .MullikenMultipoleEvaluator import MullikenMultipoleEvaluator
from  .QCDefinitions              import ChargeModel           , \
                                         FockClosurePriority
from  .QCModelError               import QCModelError
from ..EnergyModel                import EnergyClosurePriority , \
                                         EnergyModel           , \
                                         EnergyModelState
from ..Restraint                  import RestraintEnergyModel

#===================================================================================================================================
# . Notes:
#===================================================================================================================================
#
# . Both charge and spin restraints are of the form:
#
#     R_k = Tr ( P * W_k ) + N_k - Q_k                        (1)
#
#   where P is either the total or spin density, W_k is the symmetric weight matrix appropriate for the charge model,
#   N_k is the sum of nuclear charges, which is non-zero only for charge restraints, and Q_k is the target value.
#   For charge restraints W is negative (as electrons are negative) whereas for spin restraints W is positive.
#
#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Allowed charge models and the overall default for models with more than one choice.
_AllowedChargeModels = { "DFT QC Model"  : ( ChargeModel.Loewdin, ChargeModel.Mulliken ) ,
                         "MNDO QC Model" : ( ChargeModel.MNDO   ,                      ) }
_DefaultChargeModel  = ChargeModel.Loewdin #Mulliken
_DefaultChargeModels = { "MNDO QC Model" : ChargeModel.MNDO } # . Models with specific defaults.

# . Evaluators.
_Evaluators = { ChargeModel.Loewdin  : LoewdinMultipoleEvaluator  ,
                ChargeModel.MNDO     : MNDOMultipoleEvaluator     ,
                ChargeModel.Mulliken : MullikenMultipoleEvaluator }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ChargeRestraint ( AttributableObject ):
    """A charge restraint."""
    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "energyModel" : None  , # . Atomic units!
                             "indices"     : None  ,
                             "isSpin"      : False ,
                             "target"      : 0.0   ,
                             "weights"     : None  } )

    def _CheckOptions ( self ):
        """Check options."""
        n = len ( self.indices )
        # . It is easier to use IntegerArray1D here than Selection as indices and weights are in pairs
        # . and because the indices do not act like selections for merging and pruning.
        # . If weights is absent then create it full of ones.
        if self.weights is None:
            self.weights = Array.WithExtent ( n, dataType = DataType.Real )
            self.weights.Set ( 1.0 )
        # . General checks.
        isOK = isinstance ( self.energyModel, RestraintEnergyModel ) and \
               isinstance ( self.indices    , IntegerArray1D       ) and \
               isinstance ( self.weights    , RealArray1D          ) and \
               ( len (       self.weights   ) == n ) and \
               ( len ( set ( self.indices ) ) == n ) and \
               ( min ( self.indices ) >= 0 )
        if not isOK: raise ValueError ( "Invalid arguments to charge restraint constructor." )

    def Merge ( self, increment ):
        """Merging."""
        # . Just increment indices.
        new = Clone ( self )
        new.indices.Add ( increment )
        return new

    def Prune ( self, selection ):
        """Pruning."""
        # . Only can prune if all indices of the charge restraint are present in selection.
        new        = None
        newIndices = [ selection.Position ( oldIndex ) for oldIndex in self.indices ]
        if all ( [ newIndex >= 0 for newIndex in newIndices ] ):
            new = Clone ( self )
            for ( i, newIndex ) in enumerate ( newIndices ): new.indices[i] = newIndex
        return new

    def Value ( self, charges, spins ):
        """Calculate the value of the restraint given charge and spin arrays."""
        if self.isSpin: data = spins
        else:           data = charges
        value = - self.target
        for ( index, weight ) in zip ( self.indices, self.weights ):
            value += weight * data[index]
        return value

    @property
    def upperBound ( self ):
        return ( max ( self.indices ) + 1 ) # . + 1 because it refers to the size of the set from which the indices are drawn.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ChargeRestraintModelState ( EnergyModelState ):
    """A charge restraint model state."""

    _attributable = dict ( EnergyModelState._attributable )
    _attributable.update ( { "chargeModel" : None ,
                             "evaluator"   : None } )

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( ChargeRestraintModelState, self ).SummaryItems ( )
        if self.chargeModel is not None: items.append ( ( "Charge Model", self.chargeModel.name ) )
        return items

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ChargeRestraintModel ( EnergyModel ):
    """A class defining a charge restraint energy model."""

    _attributable = dict ( EnergyModel._attributable )
    _classLabel   = "Charge Restraint Model"
    _stateName    = "chargeRestraintModelState"
    _stateObject  = ChargeRestraintModelState
    _summarizable = dict ( EnergyModel._summarizable )
    _attributable.update ( { "chargeModel" : None ,
                             "restraints"  : dict } )
    _summarizable.update ( { "chargeModel" : "Charge Model" } )

    def __contains__ ( self, key ):
        """Membership."""
        return ( key in self.restraints )

    def __delitem__  ( self, key ):
        """Delete key."""
        self.restraints.pop ( key, None )

    def __getitem__  ( self, key ):
        """Get an item."""
        return self.restraints.get ( key, None )

    def __len__ ( self ):
        """The number of restraints."""
        return len ( self.restraints )

    def __setitem__ ( self, key, value ):
        """Set an item."""
        if isinstance ( key, str ) and isinstance ( value, ChargeRestraint ):
            self.restraints.pop ( key, None )
            self.restraints[key] = value
        else: raise TypeError ( "Invalid restraint key or value type." )

    def BuildModel ( self, target, **options ):
        """Build the model."""
        state = super ( ChargeRestraintModel, self ).BuildModel ( target )
        if not hasattr ( target, "qcModel" ):
            raise QCModelError ( "Cannot add charge restraints to a system without a QC model." )
        target.qcState.AddFockModel ( self.__class__._classLabel, self )
        # . Check the charge model and set the evaluator.
        qcName = str ( target.qcModel )
        # . Set a default charge model in state.
        if self.chargeModel is None:
            chargeModel       = _DefaultChargeModels.get ( qcName, _DefaultChargeModel )
            state.chargeModel = chargeModel
        # . Check that the charge model is compatible with the QC model.
        else:
            chargeModel = self.chargeModel
            allowed     = _AllowedChargeModels.get ( qcName, None )
            if ( allowed is None ) or ( chargeModel not in allowed ):
                raise QCModelError ( "A {:s} charge model is incompatible with the system's {:s}".format ( chargeModel.name, qcName ) )
        state.evaluator = _Evaluators[chargeModel] ( )
        return state

    def ClearScratch ( self, scratch ):
        """Clear scratch."""
        for attribute in ( "crDLambdas"  ,
                           "crLambdas"   ,
                           "crWMatrices" ): 
            scratch.Delete ( attribute )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.SetUpMakeFock ( target )
        closures = [ ( EnergyClosurePriority.QCPreEnergy, a, self.__class__._classLabel ) ]
        state = getattr ( target, self.__class__._stateName )
        if hasattr ( state.evaluator, "ChargeRestraintWeightedDensity" ):
            def b ( ): self.WeightedDensity ( target )
            closures.append ( ( EnergyClosurePriority.QCPreGradients, b, self.__class__._classLabel ) )
        return closures

    def FockClosures ( self, target ):
        """Fock closures."""
        def a ( ): return self.MakeFock ( target )
        return [ ( FockClosurePriority.Medium, a ) ]

    def MakeFock ( self, target ):
        """Charge restraint contribution to the Fock matrices and the electronic energy."""
        eCR       = 0.0
        scratch   = target.scratch
        onePDMP   = scratch.Get ( "onePDMP", None )
        onePDMQ   = scratch.Get ( "onePDMQ", None )
        lambdas   = scratch.crLambdas
        wMatrices = scratch.crWMatrices
        for ( key, restraint ) in self.restraints.items ( ):
            ( W, N ) = wMatrices[key]
            if restraint.isSpin:
                P = onePDMQ.density
                F = onePDMQ.fock
            else:
                P = onePDMP.density
                F = onePDMP.fock
            value = P.TraceOfProduct ( W ) + N - restraint.target # . The restraint value.
            ( e, g ) = restraint.energyModel.Energy ( value )
            F.Add ( W, scale = g )
            eCR += e
            lambdas[key] = ( e, g )
        scratch.energyTerms   ["Charge Restraint Energy"]  = ( eCR * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        scratch.qcEnergyReport["Charge Restraint Energy"]  = eCR
        return eCR

    # . No checking is done for duplicate keys.
    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging."""
        new         = None
        new._target = information.get ( "Target"          , None )
        increments  = information.get ( "Index Increments", None )
        if ( increments is not None ) and ( len ( increments ) == len ( items ) ):
            newItems = {}
            for ( item, increment ) in zip ( items, increments ):
                if item is not None:
                    for ( key, value ) in item.restraints.items ( ):
                        if hasattr ( value, "Merge" ):
                            newItems[key] = value.Merge ( increment )
            if len ( newItems ) > 0:
                new = selfClass ( )
                new.restraints.update ( newItems )
        return new

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        pruned = None
        items  = {}
        for ( key, value ) in self.restraints.items ( ):
            if hasattr ( value, "Prune" ):
                item = value.Prune ( selection )
                if item is not None: items[key] = item
        if len ( items ) > 0:
            pruned = self.__class__ ( )
            pruned.restraints.update ( items )
        return pruned

    def SetUpMakeFock ( self, target ):
        """Set up the charge restraint energy and Fock calculation."""
        restraints = self.restraints
        scratch    = target.scratch
        state      = getattr ( target, self.__class__._stateName )
        onePDMP    = scratch.Get ( "onePDMP", None )
        onePDMQ    = scratch.Get ( "onePDMQ", None )
        if ( onePDMP is None ) and ( self.chargeRestraints > 0 ):
            raise QCModelError ( "Missing total density matrix for charge restraint calculation." )
        if ( onePDMQ is None ) and ( self.spinRestraints > 0 ):
            raise QCModelError ( "Missing spin density matrix for spin restraint calculation." )
        scratch.Set ( "crLambdas"  , { key : None                                                   for ( key, item ) in self.restraints.items ( ) } )
        scratch.Set ( "crWMatrices", { key : state.evaluator.ChargeRestraintMatrix ( target, item ) for ( key, item ) in self.restraints.items ( ) } )

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( ChargeRestraintModel, self ).SummaryItems ( )
        items.extend ( [ ( "Charge Restraints", "{:d}".format ( self.chargeRestraints ) ) ,
                         ( "Spin Restraints"  , "{:d}".format ( self.spinRestraints   ) ) ] )
        return items

    def UnbuildModel ( self, target ):
        """Unbuild the model."""
        self.ClearScratch ( target.scratch )
        target.qcState.AddFockModel ( self.__class__._classLabel, None )

    def Values ( self, charges, spins ):
        """Calculate the values of the restraints given charge and spin arrays."""
        return { key : item.Value ( charges, spins ) for ( key, item ) in self.restraints.items ( ) }

    def WeightedDensity ( self, target ): 
        """Get the weighted density."""
        scratch = target.scratch
        if scratch.doGradients:
            state = getattr ( target, self.__class__._stateName )
            state.evaluator.ChargeRestraintWeightedDensity ( target                  ,
                                                             self.restraints         ,
                                                             scratch.crLambdas       , # . Assume latest values OK!
                                                             scratch.weightedDensity )

    @property
    def chargeRestraints ( self ):
        n = self.__dict__.get ( "_chargeRestraints", None )
        if n is None:
            n = [ item.isSpin for item in self.restraints.values ( ) ].count ( False )
            self.__dict__["_chargeRestraints"] = n
        return n

    @property
    def spinRestraints ( self ):
        n = self.__dict__.get ( "_spinRestraints", None )
        if n is None:
            n = [ item.isSpin for item in self.restraints.values ( ) ].count ( True )
            self.__dict__["_spinRestraints"] = n
        return n

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

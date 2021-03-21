"""Class for a restraint energy model."""

from  collections import defaultdict
from  pCore       import logFile               , \
                         LogFileActive         , \
                         RawObjectConstructor
from .EnergyModel import EnergyClosurePriority , \
                         EnergyModel
from .Restraint   import Restraint

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RestraintModel ( EnergyModel ):
    """A class defining a restraint energy model."""

    _attributable = dict ( EnergyModel._attributable )
    _classLabel   = "Restraint Model"
    _attributable.update ( { "restraints" : dict } )

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
        if isinstance ( key, str ) and isinstance ( value, Restraint ):
            self.restraints.pop ( key, None )
            self.restraints[key] = value
        else: raise TypeError ( "Invalid key or value type." )

    def Energy ( self, coordinates3, gradients3 ):
        """Energy and gradients."""
        energy  = 0.0
        scTerms = { key : value.Energy ( coordinates3, gradients3 ) for ( key, value ) in self.restraints.items ( ) }
        for ( e, d ) in scTerms.values ( ): energy += e
        return ( ( "Restraints", energy ), scTerms )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def f ( ):
            coordinates3 = target.coordinates3
            gradients3   = target.scratch.Get ( "gradients3", None )
            ( ( label, energy ), state )      = self.Energy ( coordinates3, gradients3 )
            target.scratch.energyTerms[label] = energy
            target.scratch.restraintTerms     = state
        return [ ( EnergyClosurePriority.IndependentEnergyTerm, f, self.__class__._classLabel ) ]

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

    def SummaryItems ( self ):
        """Summary items."""
        frequencies = defaultdict ( int )
        for ( key, value ) in self.restraints.items ( ): frequencies[value.__class__._classLabel] += 1
        return [ ( key, "{:d}".format ( value ) ) for ( key, value ) in frequencies.items ( ) ]

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

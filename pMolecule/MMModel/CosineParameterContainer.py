"""Cosine parameter classes."""

import math

from  pScientific          import Units
from .CosineTermContainer  import CosineAngleContainer      , \
                                  CosineDihedralContainer   , \
                                  CosineOutOfPlaneContainer
from .MMParameterContainer import MMParameterContainer
from .MMModelError         import MMModelError

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
_PeriodForceConstantLabel = "Period / Force Constant Pairs"

#===================================================================================================================================
# . Base class.
#===================================================================================================================================
class CosineParameterContainer ( MMParameterContainer ):
    """A container for cosine parameters."""

    # . Defaults.
    _attributable    = dict ( MMParameterContainer._attributable )
    _containerClass  = None
    _numberOfIndices = 0
    _termLabel       = "Cosine Generic"
    _attributable.update ( { "containerClass"  : None ,
                             "itemsWild"       : None ,
                             "numberOfIndices" : None } )

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        mapping["Parameter Fields"] = [ "Atom Type {:d}".format ( i+1 ) for i in range ( self.numberOfIndices ) ] + [ _PeriodForceConstantLabel ]
        # . Rows.
        rows = []
        for key in sorted ( self.rawItems.keys ( ) ):
            rows.append ( list ( key ) + self.rawItems[key] )
        mapping["Parameter Values"] = rows
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        for key in ( "Analytic Form", "Units", "Wild Card" ):
            if key in self.properties:
                mapping[key] = self.properties[key]
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self._Initialize ( )
            columns         = mapping.pop ( "Parameter Fields" )
            rows            = mapping.pop ( "Parameter Values" )
            self.label      = mapping.pop ( "Label", None )
            self.properties = dict ( mapping )
            indices         = [ columns.index ( "Atom Type {:d}".format ( i+1 ) ) for i in range ( self.numberOfIndices ) ]
            pairs           = columns.index ( _PeriodForceConstantLabel )
            # . Create the raw items.
            self.rawItems = {}
            if rows is not None:
                for row in rows:
                    key = self.MakeKey ( *[ row[i] for i in indices ] )
                    self.rawItems[key] = pairs
                if len ( self.rawItems ) != len ( rows ): raise
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
        except Exception as e:
            print ( e.args )
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def _Initialize ( self ):
        """Initialization."""
        super ( CosineParameterContainer, self )._Initialize ( )
        self.containerClass  = self.__class__._containerClass
        self.numberOfIndices = self.__class__._numberOfIndices

    @staticmethod
    def ConvertItemUnits ( items, fCFactor ):
        """Unit conversion."""
        newItems = {}
        for ( key, pairs ) in items.items ( ):
            newItems[key] = [ ( p, fCFactor * fC ) for ( p, fC ) in pairs ]
        return newItems

    @staticmethod
    def GetTermsFromIndices ( indices ):
        """Get all possible terms from a set of indices."""
        return [ indices ]

    def MakeMMTerms ( self, atomTypeLabels, termIndices, connectivity ):
        """Make the appropriate MM terms given a list of term indices."""
        missingParameters = set ( )
        mmTerms           = []
        if len ( termIndices ) > 0:
            # . Wild cards.
            hasWildCard = ( len ( self.itemsWild ) > 0 )
            # . Initialization.
            parameterKeys = {}
            parameters    = []
            terms         = []
            # . Generate both term and parameter data.
            for indices in termIndices:
                typeLabels = [ atomTypeLabels[i] for i in indices ]
                key        = self.MakeKey ( *typeLabels )
                pairs      = self.items.get ( key, None )
                if ( pairs is None ) and hasWildCard:
                    numberFound = 0
                    for newKey in self.MakeWildCardKeys ( typeLabels ):
                        newPairs = self.itemsWild.get ( newKey, None )
                        if newPairs is not None:
                            if numberFound == 0: pairs = newPairs
                            numberFound += 1
                    # . Error - temporary.
                    if numberFound > 1: pairs = None
                    if pairs is not None: key = newKey
                # . Automatic generation.
                if ( pairs is None ) and ( self.parameterFactory is not None ):
                    ( pairs, key ) = self.parameterFactory ( typeLabels, indices, connectivity, key )
                # . Still no parameters.
                if pairs is None:
                    if self.useStrictAssignment:
                        missingParameters.add ( ( self.termLabel, key ) )
                # . Parameters found.
                # . Zero-length parameters allowed but ignored.
                elif len ( pairs ) > 0:
                    p = parameterKeys.get ( key, -1 )
                    # . Add the parameters if necessary.
                    # . The parameters are already adjusted for term multiplicity.
                    if p == -1:
                        p = len ( parameterKeys )
                        parameterKeys[key] = p
                        parameters.append ( pairs )
                    # . Add the terms.
                    for newIndices in self.GetTermsFromIndices ( indices ):
                        terms.append ( [ newIndices, p, True ] )
            # . Construct the parameter keys.
            newKeys = [ None for i in range ( len ( parameterKeys ) ) ]
            for ( key, value ) in parameterKeys.items ( ):
                newKeys[value] = self.keySeparator.join ( key )
            # . Construct the container.
            numberOfParameters = len ( parameters )
            numberOfTerms      = len ( terms      ) 
            if ( numberOfParameters > 0 ) and ( numberOfTerms > 0 ) and ( len ( missingParameters ) <= 0 ):
                state = { "indices"       : self.numberOfIndices ,
                          "label"         : self.termLabel       ,
                          "parameterKeys" : newKeys              ,
                          "parameters"    : parameters           ,
                          "terms"         : sorted ( terms )     } # . Simple sorting.
                #self.PrintMMTerms ( atomTypeLabels, state )
                mm = self.containerClass.Raw ( )
                mm.__setstate__ ( state )
                mmTerms.append ( mm )
        return ( mmTerms, missingParameters )

    def MakeMMTermsFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make the appropriate MM terms given a connectivity."""
        termIndices = self.GetTermIndices ( connectivity )
        return self.MakeMMTerms ( atomTypeLabels, termIndices, connectivity )

    def PrintMMTerms ( self, atomTypeLabels, state ):
        """Print MM terms - for debugging only."""
        print ( "\nMM Terms for {:s}:\n".format ( state["label"] ) )
        for ( i, ( indices, p, _ ) ) in enumerate ( state["terms"] ):
            key   = state["parameterKeys"][p]
            pairs = state["parameters"][p]
            print ( "{:30s} {:30s} {:s}".format ( repr ( indices ), key, repr ( pairs ) ) )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        # . Get options.
        units    = self.properties.get ( "Units"    , None )
        wildCard = self.properties.get ( "Wild Card", None )
        # . Do nothing.
        if ( units is None ) and ( wildCard is None ):
            self.items     = self.rawItems
            self.itemsWild = {}
        # . Process entries.
        else:
            # . Units.
            convert  = False
            fCFactor = 1.0
            if units is not None:
                if units.get ( "Force Constant", "" ).find ( "kcal" ) > -1:
                    convert  = True
                    fCFactor = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
            # . Wild card.
            if wildCard is None:
                items     = self.rawItems
                itemsWild = {}
            else:
                items     = {}
                itemsWild = {}
                for ( key, value ) in self.rawItems.items ( ):
                    if wildCard in key: itemsWild[key] = value
                    else:               items    [key] = value
            # . Convert.
            if convert:
                self.items     = self.ConvertItemUnits ( items    , fCFactor )
                self.itemsWild = self.ConvertItemUnits ( itemsWild, fCFactor )
            else:
                self.items     = items
                self.itemsWild = itemsWild

#===================================================================================================================================
# . Specific classes.
#===================================================================================================================================
class CosineAngleParameterContainer ( CosineParameterContainer ):
    """A container for cosine angle parameters."""

    _containerClass  = CosineAngleContainer
    _numberOfIndices = 3
    _termLabel       = "Cosine Angle"

    @staticmethod
    def GetTermsFromIndices ( indices ):
        """Get all possible terms from a set of indices."""
        ( i, j, k ) = indices
        if i > k: return [ [ i, j, k ] ]
        else:     return [ [ k, j, i ] ]

    def GetTermIndices ( self, connectivity ):
        """Get the term indices."""
        return connectivity.angleIndices # . i-j-k.

    @staticmethod
    def MakeKey ( *labels ):
        """Make a key."""
        ( label1, label2, label3 ) = labels
        if ( label1 >= label3 ): return ( label1, label2, label3 )
        else:                    return ( label3, label2, label1 )

    def MakeWildCardKeys ( self, typeLabels ):
        """Make wild card keys."""
        # . X-j-X.
        keys     = []
        wildCard = self.properties.get ( "Wild Card", None )
        if wildCard is not None:
            keys.append ( self.MakeKey ( wildCard, typeLabels[1], wildCard ) )
        return keys

class CosineDihedralParameterContainer ( CosineParameterContainer ):
    """A container for cosine dihedral parameters."""

    _containerClass  = CosineDihedralContainer
    _numberOfIndices = 4
    _termLabel       = "Cosine Dihedral"

    @staticmethod
    def GetTermsFromIndices ( indices ):
        """Get all possible terms from a set of indices."""
        ( i, j, k, l ) = indices
        if   ( j >  k ): return [ ( i, j, k, l ) ]
        elif ( j <  k ): return [ ( l, k, j, i ) ]
        elif ( j == k ):
            if ( i > l ): return [ ( i, j, k, l ) ]
            else:         return [ ( l, k, j, i ) ]

    def GetTermIndices ( self, connectivity ):
        """Get the term indices."""
        return connectivity.dihedralIndices # . i-j-k-l.

    @staticmethod
    def MakeKey ( *labels ):
        """Make a key."""
        ( label1, label2, label3, label4 ) = labels
        if   ( label2 >  label3 ): return ( label1, label2, label3, label4 )
        elif ( label2 <  label3 ): return ( label4, label3, label2, label1 )
        elif ( label2 == label3 ):
            if ( label1 > label4 ): return ( label1, label2, label3, label4 )
            else:                   return ( label4, label3, label2, label1 )

    def MakeWildCardKeys ( self, typeLabels ):
        """Make wild card keys."""
        # . X-j-k-X.
        keys     = []
        wildCard = self.properties.get ( "Wild Card", None )
        if wildCard is not None:
            keys.append ( self.MakeKey ( wildCard, typeLabels[1], typeLabels[2], wildCard ) )
        return keys

class CosineOutOfPlaneParameterContainer ( CosineParameterContainer ):
    """A container for cosine out-of-plane parameters."""

    _containerClass  = CosineOutOfPlaneContainer
    _numberOfIndices = 4
    _termLabel       = "Cosine Out-Of-Plane"

    @staticmethod
    def GetTermsFromIndices ( indices ):
        """Get all possible terms from a set of indices."""
        # . All possible terms around j which is the central atom.
        ( j, i, k, l ) = indices
        newIndices     = []
        for ( i1, i2, i3 ) in ( ( i, k, l ), ( k, l, i ), ( l, i, k ) ):
            if i2 > i3: newIndices.append ( ( i1, j, i2, i3 ) )
            else:       newIndices.append ( ( i1, j, i3, i2 ) )
        return newIndices

    def GetTermIndices ( self, connectivity ):
        """Get the term indices."""
        return connectivity._IndicesNeighbor3 ( ) # . All 3-coordinate atoms such that i-(j,k,l). 

    @staticmethod
    def MakeKey ( *labels ):
        """Make a key."""
        return tuple ( [ labels[0] ] + sorted ( labels[1:] ) )

    def MakeWildCardKeys ( self, typeLabels ):
        """Make wild card keys."""
        # . i-j-X-X.
        keys     = []
        wildCard = self.properties.get ( "Wild Card", None )
        if wildCard is not None:
            for label in typeLabels[1:]:
                keys.append ( self.MakeKey ( typeLabels[0], label, wildCard, wildCard ) )
        return keys

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

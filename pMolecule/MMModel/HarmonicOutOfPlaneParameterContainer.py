"""Harmonic improper dihedral out-of-plane parameter classes."""

import math

from  pScientific               import Units
from .HarmonicImproperContainer import HarmonicImproperContainer
from .MMModelError              import MMModelError
from .MMParameterContainer      import MMParameterContainer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class HarmonicOutOfPlaneParameterContainer ( MMParameterContainer ):
    """A container for harmonic improper dihedral out-of-plane parameters."""

    _attributable = dict ( MMParameterContainer._attributable )
    _termLabel    = "Harmonic Out-Of-Plane"
    _attributable.update ( { "itemsWild" : None } )

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        mapping["Parameter Fields"] = [ "Central Atom Type", "Peripheral Atom Type 1", "Peripheral Atom Type 2", "Peripheral Atom Type 3", "Force Constant" ]
        # . Rows.
        keys = list ( self.rawItems.keys ( ) )
        keys.sort ( )
        rows = []
        for key in keys:
            rows.append ( [ key[0], key[1], key[2], key[3], self.rawItems[key] ] )
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
            # . Get token indices.
            l1 = columns.index ( "Central Atom Type"      )
            l2 = columns.index ( "Peripheral Atom Type 1" )
            l3 = columns.index ( "Peripheral Atom Type 2" )
            l4 = columns.index ( "Peripheral Atom Type 3" )
            fC = columns.index ( "Force Constant"               )
            # . Create the raw items.
            self.rawItems = {}
            if rows is not None:
                for row in rows:
                    key = self.MakeKey ( row[l1], row[l2], row[l3], row[l4] )
                    self.rawItems[key] = float ( row[fC] )
                if len ( self.rawItems ) != len ( rows ): raise
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
        except Exception as e:
            print ( e[0] )
#        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @staticmethod
    def ConvertItemUnits ( items, fCFactor ):
        """Unit conversion."""
        newItems = {}
        for ( key, fC ) in items.items ( ):
            newItems[key] = fCFactor * fC
        return newItems

    @staticmethod
    def MakeKey ( label1, label2, label3, label4 ):
        """Make a key."""
        peripherals = [ label2, label3, label4 ]
        peripherals.sort ( )
        return tuple ( [ label1 ] + peripherals )

    def MakeMMTerms ( self, atomTypeLabels, termIndices, connectivity ):
        """Make the appropriate MM terms given a list of term indices."""
        missingParameters = set ( )
        mmTerms           = []
        if len ( termIndices ) > 0:
            # . Wild cards.
            hasWildCard =  ( len ( self.itemsWild ) > 0 )
            if hasWildCard: wildCard = self.properties.get ( "Wild Card", None )
            # . Initialization.
            parameterKeys = {}
            parameters    = []
            terms         = []
            # . Generate both term and parameter data.
            for ( i, j, k, l ) in termIndices:
                ti  = atomTypeLabels[i]
                tj  = atomTypeLabels[j]
                tk  = atomTypeLabels[k]
                tl  = atomTypeLabels[l]
                key = self.MakeKey ( ti, tj, tk, tl )
                fC  = self.items.get ( key, None )
                if ( fC is None ) and hasWildCard:
                    # . Do all three possible permutations.
                    numberFound = 0
                    toTry       = set ( [ tj, tk, tl ] )
                    for to in toTry:
                        newKey = self.MakeKey ( ti, to, wildCard, wildCard )
                        newFC  = self.itemsWild.get ( newKey, None )
                        if newFC is not None:
                            if numberFound == 0: fC = newFC
                            numberFound += 1
                    # . Error - temporary.
                    if numberFound > 1: fC = None
                    if fC is not None: key = newKey
                if fC is None:
                    if self.useStrictAssignment:
                        missingParameters.add ( ( self.termLabel, key ) )
                elif fC != 0.0:
                    p = parameterKeys.get ( key, -1 )
                    if p == -1:
                        p = len ( parameterKeys )
                        parameterKeys[key] = p
                        parameters.append ( ( 0.0, fC / 3.0 ) )
                    # . Add three symmetrical terms.
                    for ( i1, i2, i3 ) in ( ( j, k, l ), ( k, l, j ), ( l, j, k ) ):
                        terms.append ( ( i, i1, i2, i3, p, True ) )
            # . Construct the parameter keys.
            newKeys = [ None for i in range ( len ( parameterKeys ) ) ]
            for ( key, value ) in parameterKeys.items ( ):
                newKeys[value] = self.keySeparator.join ( key )
            # . Construct the container.
            numberOfParameters = len ( parameters )
            numberOfTerms      = len ( terms      ) 
            if ( numberOfParameters > 0 ) and ( numberOfTerms > 0 ) and ( len ( missingParameters ) <= 0 ):
                state = { "label"         : self.termLabel ,
                          "parameterKeys" : newKeys        ,
                          "parameters"    : parameters     ,
                          "terms"         : terms          }
                mm = HarmonicImproperContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
                mmTerms.append ( mm )
        return ( mmTerms, missingParameters )

    def MakeMMTermsFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make the appropriate MM terms given a connectivity."""
        termIndices = connectivity._IndicesNeighbor3 ( )
        return self.MakeMMTerms ( atomTypeLabels, termIndices, connectivity )

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
                if units.get ( "Force Constant", "" ).find ( "kcal"    ) > -1:
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
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

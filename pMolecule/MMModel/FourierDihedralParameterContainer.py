"""Fourier dihedral parameter classes."""

from  pScientific              import Units
from .FourierDihedralContainer import FourierDihedralContainer
from .MMParameterContainer     import MMParameterContainer
from .MMModelError             import MMModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class FourierDihedralParameterContainer ( MMParameterContainer ):
    """A container for fourier dihedral parameters."""

    _attributable = dict ( MMParameterContainer._attributable )
    _termLabel    = "Fourier Dihedral"
    _attributable.update ( { "itemsWild" : None } )

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        mapping["Parameter Fields"] = [ "Atom Type 1", "Atom Type 2", "Atom Type 3", "Atom Type 4", "Force Constant", "Periodicity", "Phase" ]
        # . Rows.
        keys = list ( self.rawItems.keys ( ) )
        keys.sort ( )
        rows = []
        for key in keys:
            for ( fC, period, phase ) in self.rawItems[key]:
                rows.append ( [ key[0], key[1], key[2], key[3], fC, period, phase ] )
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
            l1     = columns.index ( "Atom Type 1" )
            l2     = columns.index ( "Atom Type 2" )
            l3     = columns.index ( "Atom Type 3" )
            l4     = columns.index ( "Atom Type 4" )
            fC     = columns.index ( "Force Constant"    )
            period = columns.index ( "Periodicity"       )
            phase  = columns.index ( "Phase"             )
            # . Create the raw items.
            numberOfDuplicates = 0
            self.rawItems      = {}
            if rows is not None:
                for row in rows:
                    key   = self.MakeKey ( row[l1], row[l2], row[l3], row[l4] )
                    terms = self.rawItems.get ( key, [] )
                    n     = int ( row[period] )
                    foundDuplicate = False
                    for term in terms:
                        if n == term[1]:
                            foundDuplicate = True
                            numberOfDuplicates += 1
                            break
                    if not foundDuplicate:
                        terms.append ( ( float ( row[fC] ), n, float ( row[phase] ) ) )
                        self.rawItems[key] = terms
                if numberOfDuplicates > 0: raise
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @staticmethod
    def ConvertItemUnits ( items, fCFactor, phFactor ):
        """Unit conversion."""
        newItems = {}
        for ( key, terms ) in items.items ( ):
            newTerms = []
            for ( fC, period, phase ) in terms:
                newTerms.append ( ( fCFactor * fC, period, phFactor * phase ) )
            newItems[key] = newTerms
        return newItems

    @staticmethod
    def MakeKey ( label1, label2, label3, label4 ):
        """Make a key."""
        if   ( label2 >  label3 ): return ( label1, label2, label3, label4 )
        elif ( label2 <  label3 ): return ( label4, label3, label2, label1 )
        elif ( label2 == label3 ):
            if ( label1 > label4 ): return ( label1, label2, label3, label4 )
            else:                   return ( label4, label3, label2, label1 )

    def MakeMMTerms ( self, atomTypeLabels, termIndices, connectivity ):
        """Make the appropriate MM terms given a list of term indices."""
        missingParameters = set ( )
        mmTerms           = []
        if len ( termIndices ) > 0:
            # . Wild cards.
            hasWildCard       =  ( len ( self.itemsWild ) > 0 )
            if hasWildCard: wildCard = self.properties.get ( "Wild Card", None )
            # . Initialization.
            parameterKeys = {}
            parameters    = []
            terms         = []
            # . Generate both term and parameter data.
            for ( i, j, k, l ) in termIndices:
                ti        = atomTypeLabels[i]
                tj        = atomTypeLabels[j]
                tk        = atomTypeLabels[k]
                tl        = atomTypeLabels[l]
                key       = self.MakeKey ( ti, tj, tk, tl )
                parameter = self.items.get ( key, None )
                if ( parameter is None ) and hasWildCard:
                    newKey    = self.MakeKey ( wildCard, tj, tk, wildCard )
                    parameter = self.itemsWild.get ( newKey, None )
                    if parameter is not None: key = newKey
                if parameter is None:
                    if self.useStrictAssignment:
                        missingParameters.add ( ( self.termLabel, key ) )
                else:
                    pIndices = parameterKeys.get ( key, None )
                    if pIndices is None:
                        pIndex   = len ( parameters )
                        pIndices = []
                        for ( fC, period, phase ) in parameter:
                            if fC != 0.0:
                                parameters.append ( ( fC, period, phase ) )
                                pIndices.append ( pIndex )
                                pIndex += 1
                        parameterKeys[key] = pIndices
                    for p in pIndices:
                        terms.append ( ( i, j, k, l, p, True ) )
            # . Construct the parameter keys.
            newKeys = [ None for i in range ( len ( parameters ) ) ]
            for ( key, pIndices ) in parameterKeys.items ( ):
                for pIndex in pIndices:
                    newKeys[pIndex] = self.keySeparator.join ( list ( key ) + [ repr ( parameters[pIndex][1] ) ] )
            # . Construct the container.
            numberOfParameters = len ( parameters )
            numberOfTerms      = len ( terms      ) 
            if ( numberOfParameters > 0 ) and ( numberOfTerms > 0 ) and ( len ( missingParameters ) <= 0 ):
                state = { "label"         : self.termLabel ,
                          "parameterKeys" : newKeys        ,
                          "parameters"    : parameters     ,
                          "terms"         : terms          }
                mm = FourierDihedralContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
                mmTerms.append ( mm )
        return ( mmTerms, missingParameters )

    def MakeMMTermsFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make the appropriate MM terms given a connectivity."""
        termIndices = connectivity.dihedralIndices
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
            phFactor = 1.0
            if units is not None:
                if units.get ( "Force Constant", "" ).find ( "kcal"    ) > -1:
                    convert  = True
                    fCFactor = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
                if units.get ( "Phase"         , "" ).find ( "degrees" ) > -1:
                    convert  = True
                    phFactor = Units.Angle_Degrees_To_Radians
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
                self.items     = self.ConvertItemUnits ( items    , fCFactor, phFactor )
                self.itemsWild = self.ConvertItemUnits ( itemsWild, fCFactor, phFactor )
            else:
                self.items     = items
                self.itemsWild = itemsWild

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

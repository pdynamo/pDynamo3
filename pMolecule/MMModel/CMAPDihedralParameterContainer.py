"""CMAP dihedral parameter classes."""

from  pScientific           import Units
from .CMAPDihedralContainer import CMAPDihedralContainer
from .MMParameterContainer  import MMParameterContainer
from .MMModelError          import MMModelError

# . At the moment it is assumed that there are always three common atoms in each of the dihedrals.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CMAPDihedralParameterContainer ( MMParameterContainer ):
    """A container for CMAP dihedral parameters."""


    _attributable     = dict ( MMParameterContainer._attributable )
    _strictAssignment = False
    _termLabel        = "CMAP Dihedral"
    _attributable.update ( { "rawXIntervals" : None ,
                             "rawYIntervals" : None ,
                             "xIntervals"    : None ,
                             "yIntervals"    : None } )

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Intervals.
        mapping["X Intervals"] = self.rawXIntervals
        mapping["Y Intervals"] = self.rawYIntervals
        # . Parameters.
        parameters = []
        keys       = list ( self.rawItems.keys ( ) )
        keys.sort ( )
        for key in keys:
            values      = self.rawItems[key]
            listKey     = list ( key )
            commonTypes = set ( listKey[0:4] ).intersection ( set ( listKey[4:] ) )
            parameters.append ( { "Atom Types" : listKey, "Common Atom Types" : commonTypes, "Values" : values } )
        mapping["Parameters"] = parameters
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        for key in ( "Units", ):
            if key in self.properties:
                mapping[key] = self.properties[key]
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self._Initialize ( )
            parameters         = mapping.pop ( "Parameters"  )
            self.rawXIntervals = mapping.pop ( "X Intervals" )
            self.rawYIntervals = mapping.pop ( "Y Intervals" )
            self.label         = mapping.pop ( "Label", None )
            self.properties    = dict ( mapping )
            # . Process parameters.
            nXY   = len ( self.rawXIntervals ) * len ( self.rawYIntervals )
            items = {}
            if parameters is not None:
                for data in parameters:
                    key    = self.MakeKey ( data["Atom Types"] )
                    values = data["Values"]
                    if ( key in items ) or  ( len ( values ) != nXY ): raise
                    items[key] = values
            self.rawItems = items
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @staticmethod
    def ConvertIntervalUnits ( intervals, factor ):
        """Unit conversion."""
        newIntervals = []
        for v in intervals: newIntervals.append ( v * factor )
        return newIntervals

    @staticmethod
    def MakeKey ( labels ):
        """Make a key."""
        ( l1, l2, l3, l4, l5, l6, l7, l8 ) = labels
        if   l2 >  l3: key1 = [ l1, l2, l3, l4 ]
        elif l2 == l3: key1 = [ max ( l1, l4 ), l2, l3, min ( l1, l4 ) ]
        else:          key1 = [ l4, l3, l2, l1 ]
        if   l6 >  l7: key2 = [ l5, l6, l7, l8 ]
        elif l6 == l7: key2 = [ max ( l5, l8 ), l6, l7, min ( l5, l8 ) ]
        else:          key2 = [ l8, l7, l6, l5 ]
        if key1 < key2:
            t = key1 ; key1 = key2 ; key2 = t
        return tuple ( key1 + key2 )

    def MakeMMTerms ( self, atomTypeLabels, termIndices, connectivity ):
        """Make the appropriate MM terms given a list of term indices."""
        # . Assume are 1-5 terms for moment.
        missingParameters = set ( )
        mmTerms           = []
        if len ( termIndices ) > 0:
            # . Initialization.
            parameterKeys = {}
            parameters    = []
            terms         = []
            # . Generate both term and parameter data.
            for localIndices in termIndices:
                labels    = []
                for i in localIndices: labels.append ( atomTypeLabels[i] )
                key       = self.MakeKey ( labels )
                parameter = self.items.get ( key, None )
                if parameter is None:
                    if self.useStrictAssignment:
                        missingParameters.add ( ( self.termLabel, key ) )
                else:
                    p = parameterKeys.get ( key, -1 )
                    if p == -1:
                        p = len ( parameterKeys )
                        parameterKeys[key] = p
                        parameters.append ( ( self.xIntervals, self.yIntervals, parameter ) )
                    terms.append ( list ( localIndices ) + [ p, True ] )
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
                mm = CMAPDihedralContainer.Raw ( )
                mm.__setstate__ ( state )
                mm.Sort ( )
                mmTerms.append ( mm )
        return ( mmTerms, missingParameters )

    def MakeMMTermsFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make the appropriate MM terms given a connectivity."""
        termIndices = self.MakeTermIndicesFromConnectivity ( atomTypeLabels, connectivity )
        return self.MakeMMTerms ( atomTypeLabels, termIndices, connectivity )

    def MakeTermIndicesFromConnectivity ( self, atomTypeLabels, connectivity ):
        """Make term indices from connectivity."""
        # . Unique atom types.
        uniqueAtomTypes = set ( )
        for key in self.items.keys ( ): uniqueAtomTypes.update ( key )
        # . For the moment assume there are three atoms common to all CMAP terms
        # . which means calculating all 1-2-3-4-5 terms.
        atoms           = connectivity.atoms
        dihedralIndices = connectivity.dihedralIndices
        neighbors       = connectivity.adjacentNodes
        nodeToIndex     = connectivity.nodeIndices
        termIndices     = []
        for ( i, j, k, l ) in dihedralIndices:
            types = set ( [ atomTypeLabels[i], atomTypeLabels[j], atomTypeLabels[k], atomTypeLabels[l] ] )
            if types.issubset ( uniqueAtomTypes ):
                iBonds = set ( [ nodeToIndex[n] for n in neighbors[atoms[i]] ] ).difference ( set ( [ j, k, l ] ) )
                lBonds = set ( [ nodeToIndex[n] for n in neighbors[atoms[l]] ] ).difference ( set ( [ i, j, k ] ) ).difference ( iBonds )
                for x in iBonds:
                    if atomTypeLabels[x] in uniqueAtomTypes:
                        termIndices.append ( ( x, i, j, k, i, j, k, l ) )
                for x in lBonds:
                    if atomTypeLabels[x] in uniqueAtomTypes:
                        termIndices.append ( ( i, j, k, l, j, k, l, x ) )
        return termIndices

    def ProcessRawItems ( self ):
        """Process the raw items."""
        units = self.properties.get ( "Units", None )
        if units is not None:
            # . Intervals.
            if units.get ( "Intervals", "" ).find ( "degrees" ) > -1:
                factor = Units.Angle_Degrees_To_Radians
                self.xIntervals = self.ConvertIntervalUnits ( self.rawXIntervals, factor )
                self.yIntervals = self.ConvertIntervalUnits ( self.rawYIntervals, factor )
            else:
                self.xIntervals = self.rawXIntervals
                self.yIntervals = self.rawYIntervals
            # . Parameters.
            if units.get ( "Energies", "" ).find ( "kcal"    ) > -1:
                factor = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole
                items  = {}
                for ( key, values ) in self.rawItems.items ( ):
                    newValues = []
                    for v in values: newValues.append ( factor * v )
                    items[key] = newValues
                self.items = items
            else:
                self.items = self.rawItems

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

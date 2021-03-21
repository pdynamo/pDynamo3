"""Lennard-Jones parameter classes."""

import re

from  pCore                import AttributableObject
from  pScientific          import Units
from .LJParameterContainer import LJForm               , \
                                  LJParameterContainer
from .MMModelError         import MMModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LennardJonesParameterContainer ( AttributableObject ):
    """A container for Lennard-Jones parameters."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "items"                 :            None ,
                             "properties"            :            None ,
                             "label"                 :            None ,
                             "rawItems"              :            None ,
                             "regularExpressions"    :            None ,
                             "style"                 :            None ,
                             "termLabel"             : "Lennard-Jones" ,
                             "useRegularExpressions" :          False  ,
                             "useStrictAssignment"   :           True  } )

    #yaml_tag = "!LennardJonesParameterContainer"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        mapping["Parameter Fields"] = [ "Atom Type", "Epsilon", "Sigma" ]
        # . Rows.
        keys = list ( self.rawItems.keys ( ) )
        keys.sort ( )
        rows = []
        for key in keys:
            ( epsilon, sigma ) = self.rawItems[key]
            rows.append ( [ key, epsilon, sigma ] )
        mapping["Parameter Values"] = rows
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        if self.style is not None: mapping["Style"] = self.style.name.upper ( )
        mapping["Regular Expressions"] = self.useRegularExpressions
        mapping["Strict Assignment"  ] = self.useStrictAssignment
        for key in ( "Analytic Form", "Units" ):
            if key in self.properties:
                mapping[key] = self.properties[key]
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self._Initialize ( )
            self.CheckStyle            ( mapping.pop ( "Style"               , None  ) )
            columns                    = mapping.pop ( "Parameter Fields" )
            rows                       = mapping.pop ( "Parameter Values" )
            self.label                 = mapping.pop ( "Label"               , None  )
            self.useRegularExpressions = mapping.pop ( "Regular Expressions" , False )
            self.useStrictAssignment   = mapping.pop ( "Strict Assignment"   , True  )
            self.properties            = dict ( mapping )
            # . Get token indices.
            l = columns.index ( "Atom Type" )
            e = columns.index ( "Epsilon"   )
            s = columns.index ( "Sigma"     )
            # . Create the raw items.
            self.rawItems = {}
            if rows is not None:
                for row in rows:
                    key = row[l]
                    self.rawItems[key] = ( float ( row[e] ), float ( row[s] ) )
                if len ( self.rawItems ) != len ( rows ): raise
            # . Create the items that are to be used.
            self.ProcessRawItems           ( )
            self.ProcessRegularExpressions ( )
#        except Exception as e:
#            print e[0]
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def CheckStyle ( self, style ):
        """Set the style."""
        if isinstance ( style, str ):
            if   style.upper ( ) == "AMBER": self.style = LJForm.Amber
            elif style.upper ( ) == "OPLS" : self.style = LJForm.OPLS
        if self.style is None:
            raise MMModelError ( "Invalid Lennard-Jones parameter style: {:s}.".format ( str ( style ) ) )

    @staticmethod
    def MakeKey ( label1, label2 ):
        """Make a key."""
        if ( label1 >= label2 ): return ( label1, label2 )
        else:                    return ( label2, label1 )

    def MakeParameterContainer ( self, atomTypeLabels ):
        """Make a LJ parameter container."""
        missingParameters = set ( )
        mm                = None
        if len ( atomTypeLabels ) > 0:
            epsilons      = []
            parameterKeys = []
            sigmas        = []
            for key in atomTypeLabels:
                parameters = self.items.get ( key, None )
                if ( parameters is None ) and self.useRegularExpressions:
                    matches = [ rKey for ( rKey, rValue ) in self.regularExpressions.items ( ) if re.fullmatch ( rValue, key ) ]
                    if len ( matches ) == 1: parameters = self.items[matches[0]]
                if parameters is None:
                    if self.useStrictAssignment:
                        missingParameters.add ( ( self.termLabel, tuple ( [ key ] ) ) )
                else:
                    ( epsilon, sigma ) = parameters
                    epsilons.append      ( epsilon )
                    parameterKeys.append ( key     )
                    sigmas.append        ( sigma   )
            numberOfParameters = len ( epsilons )
            if ( numberOfParameters > 0 ) and ( len ( missingParameters ) <= 0 ):
                state = { "label"         : self.termLabel    ,
                          "analyticForm"  : self.style        ,
                          "epsilons"      : epsilons          ,
                          "numberOfTypes" : numberOfParameters,
                          "parameterKeys" : parameterKeys     ,
                          "sigmas"        : sigmas            }
                mm = LJParameterContainer.Raw ( )
                mm.__setstate__ ( state )
        return ( mm, missingParameters )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        # . Create the items that are to be used and convert to the correct units if necessary.
        units = self.properties.get ( "Units", None )
        if units is not None:
            if units.get ( "Epsilon", "" ).find ( "kcal" ) > -1:
                self.items = {}
                for ( key, ( epsilon, sigma ) ) in self.rawItems.items ( ):
                    self.items[key] = ( Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole * epsilon, sigma )
            else:
                self.items = self.rawItems

    def ProcessRegularExpressions ( self ):
        """Process regular expressions."""
        if self.useRegularExpressions:
            self.regularExpressions = { key : re.compile ( key ) for key in self.items.keys ( ) }

    def ScaleEnergies ( self, scale ):
        """Scale the energy parameters."""
        if scale != 1.0:
            for ( key, ( epsilon, sigma ) ) in self.items.items ( ):
                self.items[key] = ( scale * epsilon, sigma )

    def UpdateParameters ( self, other ):
        """Update the parameters with additional values from another container."""
        if self.items is None:
            self.items = {}
        if other.items is not None:
            self.items.update ( other.items )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

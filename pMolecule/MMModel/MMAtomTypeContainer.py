"""MM atom types."""

from  pCore        import AttributableObject
from .MMModelError import MMModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMAtomType ( AttributableObject ):
    """An MM atom type."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "atomicNumber"  :   -1 ,
                             "charge"        :  0.0 ,
                             "hydrogenType"  : None ,
                             "label"         : None ,
                             "properties"    : None } )

    @classmethod
    def FromKeysValues ( selfClass, keys, values ):
        """Constructor from keys and values."""
        self = selfClass ( )
        self.properties = {}
        for ( key, value ) in zip ( keys, values ): self.properties[key] = value
        self.ProcessProperties ( )
        return self

    def GetProperty ( self, label ):
        """Get a property."""
        value = None
        if self.properties is not None:
            value = self.properties.get ( label, None )
        return value

    def ProcessProperties ( self ):
        """Process properties for needed attributes."""
        if self.properties is not None:
            self.atomicNumber = self.properties.get ( "Atomic Number" ,   -1 )
            self.charge       = self.properties.get ( "Charge"        ,  0.0 )
            self.label        = self.properties.get ( "Label"         , None )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMAtomTypeContainer ( AttributableObject ):
    """A container for MM atom types."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "items"           :            None ,
                             "parameterFields" :            None ,
                             "properties"      :            None ,
                             "label"           :            None ,
                             "rawItems"        :            None ,
                             "termLabel"       : "MM Atom Type"  } )

    #yaml_tag = "!MMAtomTypeContainer"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Columns.
        columnKeys                  = self.parameterFields
        mapping["Parameter Fields"] = columnKeys
        # . Rows.
        rowKeys = list ( self.rawItems.keys ( ) )
        rowKeys.sort ( )
        rows = []
        for rowKey in rowKeys:
            row = []
            for columnKey in columnKeys:
                row.append ( self.rawItems[rowKey].properties[columnKey] )
            rows.append ( row )
        mapping["Parameter Values"] = rows
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Basic construction.
            self.__init__ ( )
            columns              = mapping.pop ( "Parameter Fields" )
            rows                 = mapping.pop ( "Parameter Values" )
            self.label           = mapping.pop ( "Label", None )
            self.parameterFields = columns
            self.properties      = dict ( mapping )
            # . Get token indices.
            l = columns.index ( "Label" )
            # . Create the raw items.
            self.rawItems = {}
            for row in rows:
                key = row[l]
                self.rawItems[key] = MMAtomType.FromKeysValues ( columns, row )
            if len ( self.rawItems ) != len ( rows ): raise
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
#            for ( key, value ) in self.items.items ( ):
#                print ( key, value.__dict__ )
        except:
            raise MMModelError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def GetItem ( self, label ):
        """Get an item given its label."""
        return self.items.get ( label, None )

    def ProcessRawItems ( self ):
        """Process the raw items."""
        self.items = self.rawItems
        if self.items is not None:
            if "Hydrogen Type" in self.parameterFields:
                for item in self.items.values ( ):
                    hydrogenLabel     = item.properties["Hydrogen Type"]
                    item.hydrogenType = self.GetItem ( hydrogenLabel )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

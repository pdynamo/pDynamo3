"""PDB Components."""

import os, os.path

from pCore             import AttributableObject  , \
                              logFile             , \
                              LogFileActive
from pMolecule         import BondType            , \
                              Sequence            , \
                              System
from pMolecule.MMModel import MMSequenceAtom      , \
                              MMSequenceComponent , \
                              MMSequenceLink      , \
                              MMSequenceVariant

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Bond label/object mappings.
_BondLabelObjectMapping = { "Aromatic"  : ( BondType.Single    , True  ) ,
                            "Double"    : ( BondType.Double    , False ) ,
                            "Single"    : ( BondType.Single    , False ) ,
                            "Triple"    : ( BondType.Triple    , False ) ,
                            "Quadruple" : ( BondType.Quadruple , False ) ,
                            "Undefined" : ( BondType.Undefined , False ) }
_BondObjectLabelMapping = { ( BondType.Single    , True  ) : "Aromatic"  ,
                            ( BondType.Double    , False ) : "Double"    ,
                            ( BondType.Single    , False ) : "Single"    ,
                            ( BondType.Quadruple , False ) : "Quadruple" , 
                            ( BondType.Triple    , False ) : "Triple"    ,
                            ( BondType.Undefined , False ) : "Undefined" }

# . Key separators.
_KeyTokenSeparator = "_"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentError ( Exception ):
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponent ( AttributableObject ):
    """A PDB component."""

    _attributable = dict ( AttributableObject._attributable )
    _mappable     = { "Component Class"      : "componentClass"     ,
                      "Formal Charge"        : "formalCharge"       ,
                      "Label"                : "label"              ,
                      "Left Link"            : "leftLink"           ,
                      "Left Termination"     : "leftTermination"    ,
                      "Name"                 : "name"               ,
                      "PDB Class"            : "pdbClass"           ,
                      "Is In Chain"          : "isInChain"          ,
                      "Is Chain Terminating" : "isChainTerminating" ,
                      "Is Heteroatom"        : "isHeteroatom"       ,
                      "Right Link"           : "rightLink"          ,
                      "Right Termination"    : "rightTermination"   ,
                      "Variants"             : "variants"           }
    _attributable.update ( { "atomFields"           : None ,
                             "atoms"                : None ,
                             "bondFields"           : None ,
                             "bonds"                : None ,
                             "componentClass"       : None ,
                             "formalCharge"         : None ,
                             "label"                : None ,
                             "leftLink"             : None ,
                             "leftTermination"      : None ,
                             "name"                 : None ,
                             "pdbClass"             : None ,
                             "isInChain"            : None ,
                             "isChainTerminating"   : None ,
                             "isHeteroatom"         : None ,
                             "rightLink"            : None ,
                             "rightTermination"     : None ,
                             "variants"             : None } )

    #yaml_tag = "!PDBComponent"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Atom and bonds.
        for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atoms", "atomFields", "Atoms", "Atom Fields", PDBComponentAtom ) ,
                                                              ( "bonds", "bondFields", "Bonds", "Bond Fields", PDBComponentBond ) ):
            items = getattr ( self, tag, None )
            if items is not None:
                listItems = []
                keys      = getattr ( self, keyTag, None )
                if keys is None: keys = itemObject._mappable0
                for item in items:
                    listItems.append ( item.ToList ( keys ) )
                mapping[label   ] = listItems
                mapping[keyLabel] = keys
        # . Other attributes.
        for ( fullKey, key ) in self.__class__._mappable.items ( ):
            attribute = getattr ( self, key, None )
            if attribute is not None: mapping[fullKey] = attribute
        # . Finish up.
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Atom and bonds.
            for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atoms", "atomFields", "Atoms", "Atom Fields", PDBComponentAtom ),
                                                                  ( "bonds", "bondFields", "Bonds", "Bond Fields", PDBComponentBond ) ):
                listItems = mapping.get ( label, None )
                if listItems is not None:
                    items = []
                    keys  = mapping[keyLabel]
                    for listItem in listItems:
                        items.append ( itemObject.FromList ( listItem, keys ) )
                    setattr ( self, tag   , items )
                    setattr ( self, keyTag, keys  )
            # . Other attributes.
            for ( fullKey, key ) in self.__class__._mappable.items ( ):
                attribute = mapping.pop ( fullKey, None )
                if attribute is not None: setattr ( self, key, attribute )
#        except Exception as e:
#            print e[0]
        except:
            raise PDBComponentError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @classmethod
    def FromSystem ( selfClass, system, label = None ):
        """Generate a component from a system."""
        # . Atoms.
        atoms        = []
        formalCharge = 0
        labelIndices = {}
        for ( i, atom ) in enumerate ( system.atoms ):
            formalCharge += atom.formalCharge
            atomLabel     = atom.label
            if len ( atomLabel ) > 3: pdbAlign = 0
            else:                     pdbAlign = 1
            labelIndices[atom] = atomLabel
            atoms.append ( PDBComponentAtom.WithOptions ( atomicNumber = atom.atomicNumber, formalCharge = atom.formalCharge, label = atomLabel, pdbAlign = pdbAlign ) )
        # . Bonds.
        bonds = []
        if system.connectivity is not None:
            for bond in system.connectivity.bonds:
                bonds.append ( PDBComponentBond.WithOptions ( atomLabel1 = labelIndices[bond.node1], atomLabel2 = labelIndices[bond.node2], bondType = bond.type ) )
        # . Create the component.
        self = selfClass.WithOptions ( atoms = atoms, bonds = bonds, formalCharge = formalCharge, label = label )
        return self

    @staticmethod
    def MakeKey ( label ):
        """Make a key."""
        tokens = label.lower ( ).split ( )
        return _KeyTokenSeparator.join ( tokens )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( [ ( "Atoms"        , "{:d}".format ( len ( self.atoms ) ) ) ,
                                   ( "Bonds"        , "{:d}".format ( len ( self.bonds ) ) ) ,
                                   ( "Formal Charge", "{:d}".format ( self.formalCharge  ) ) ] ,
                                   title = "PDB Component Summary" )

    def ToMMSequenceObject ( self ):
        """Make an MM sequence object."""
        mmAtoms = []
        for atom in self.atoms:
            mmAtoms.append ( MMSequenceAtom.WithOptions ( charge = 0.0, label = atom.label, typeLabel = atom.label ) ) 
        mmItem = MMSequenceComponent.WithOptions ( atoms = mmAtoms, label = self.label )
        return mmItem

    def ToSystem ( self ):
        """Generate a system from the component."""
        # . Bonds.
        labelIndices = {}
        for ( i, atom ) in enumerate ( self.atoms ): labelIndices[atom.label] = i
        bonds = []
        for bond in self.bonds: bonds.append ( ( labelIndices[bond.atomLabel1], labelIndices[bond.atomLabel2], bond.bondType ) )
        # . System from atoms and bonds.
        system = System.FromAtoms ( self.atoms, bonds = bonds, withSequence = True )
        # . Set additional atom data.
        for ( pAtom, sAtom ) in zip ( self.atoms, system.atoms ):
            sAtom.formalCharge = pAtom.formalCharge
            sAtom.label        = pAtom.label
        # . Label.
        system.label = "PDB Component"
        if self.label is not None: system.label += " " + self.label
        # . Define component labels.
        for entity in system.sequence.children:
            for component in entity.children:
                component.label = self.label
        # . Finish up.
        return system

    # . Key property.
    @property
    def key ( self ):
        """Get a key."""
        if self.__dict__.get ( "key", None ) is None:
            self.__dict__["key"] = self.MakeKey ( self.label )
        return self.__dict__["key"]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentAtom ( AttributableObject ):
    """Atom in a PDB component."""

    # . To follow is position in component. Only useful for links and variants which can involve changes in order.
    # . Also cannot be used to put atom at beginning of a residue unless use special notation - e.g. toFollow = "start"?

    _attributable = dict ( AttributableObject._attributable )
    _mappable     = { "Atomic Number" : "atomicNumber" ,
                      "Formal Charge" : "formalCharge" ,
                      "Label"         : "label"        ,
                      "PDB Align"     : "pdbAlign"     ,
                      "To Follow"     : "toFollow"     }
    _mappable0    = ( "Atomic Number", "Formal Charge", "Label" )
    _attributable.update ( { "atomicNumber"  :   -1 ,
                             "formalCharge"  :    0 ,
                             "label"         : None ,
                             "pdbAlign"      :    0 ,
                             "toFollow"      : None } )

    @classmethod
    def FromList ( selfClass, attributeList, attributes ):
        """Constructor from list."""
        options = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value is not None: options[selfClass._mappable[key]] = value
        self = selfClass.WithOptions ( **options )
        return self

    def ToList ( self, attributes ):
        """Get a list representation of the object."""
        selfList = []
        for attribute in attributes:
            selfList.append ( getattr ( self, self.__class__._mappable[attribute], None ) )
        return selfList

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentBond ( AttributableObject ):
    """Bond in a PDB component."""

    _attributable = dict ( AttributableObject._attributable )
    _mappable     = { "Atom 1"        : "atomLabel1"    ,
                      "Atom 2"        : "atomLabel2"    ,
                      "Type"          : "bondTypeLabel" }
    _mappable0    = ( "Atom 1", "Atom 2", "Type" )
    _attributable.update ( { "atomLabel1"    : None  ,
                             "atomLabel2"    : None  ,
                             "bondType"      : None  ,
                             "bondTypeLabel" : None  ,
                             "isAromatic"    : False } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( PDBComponentBond, self )._CheckOptions ( )
        self.CheckBond ( )

    def CheckBond ( self ):
        """Check the bond attributes."""
        if   self.bondType      is not None: self.bondTypeLabel = _BondObjectLabelMapping[( self.bondType, self.isAromatic )]
        elif self.bondTypeLabel is not None: ( self.bondType, self.isAromatic ) = _BondLabelObjectMapping[self.bondTypeLabel]

    @classmethod
    def FromList ( selfClass, attributeList, attributes ):
        """Constructor from list."""
        options = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value is not None: options[selfClass._mappable[key]] = value
        self = selfClass.WithOptions ( **options )
        return self

    def ToList ( self, attributes ):
        """Get a list representation of the object."""
        selfList = []
        for attribute in attributes:
            selfList.append ( getattr ( self, self.__class__._mappable[attribute], None ) )
        return selfList

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentLink ( AttributableObject ):
    """Non-directional link between two PDB components."""

    _attributable = dict ( AttributableObject._attributable )
    _mappable     = { "Bond Type"      : "bondTypeLabel"  ,
                      "Is Aromatic"    : "isAromatic"     ,
                      "Label"          : "label"          ,
                      "Left Atom"      : "leftAtomLabel"  ,
                      "Left Variant"   : "leftVariant"    ,
                      "Right Atom"     : "rightAtomLabel" ,
                      "Right Variant"  : "rightVariant"   }
    _attributable.update ( { "bondType"       : None  ,
                             "bondTypeLabel"  : None  ,
                             "isAromatic"     : False ,
                             "label"          : None  ,
                             "leftAtomLabel"  : None  ,
                             "leftVariant"    : None  ,
                             "rightAtomLabel" : None  ,
                             "rightVariant"   : None  } )


    #yaml_tag = "!PDBComponentLink"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping         = {}
        previousState   = None
        previousVariant = None
        for ( fullKey, key ) in self.__class__._mappable.items ( ):
            attribute = getattr ( self, key, None )
            if attribute is not None:
                if key.endswith ( "Variant" ):
                    if attribute is previousVariant:
                        state           = previousState
                    else:
                        state           = attribute.__getstate__ ( )
                        previousState   = state
                        previousVariant = attribute
                    attribute = state
                mapping[fullKey] = attribute
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Build the object.
            previousState   = None
            previousVariant = None
            for ( fullKey, key ) in self.__class__._mappable.items ( ):
                attribute = mapping.pop ( fullKey, None )
                if attribute is not None:
                    if key.endswith ( "Variant" ):
                        if attribute is previousState:
                            variant = previousVariant
                        else:
                            variant = PDBComponentVariant.WithDefaults ( )
                            variant.__setstate__ ( attribute )
                            previousState   = attribute
                            previousVariant = variant
                        attribute = variant
                    setattr ( self, key, attribute )
            # . Checks.
            self.CheckBond ( )
#        except Exception as e:
#            print e[0]
        except:
            raise PDBComponentError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def _CheckOptions ( self ):
        """Check options."""
        super ( PDBComponentLink, self )._CheckOptions ( )
        self.CheckBond ( )

    def CheckBond ( self ):
        """Check the bond attributes."""
        if   self.bondType      is not None: self.bondTypeLabel = _BondObjectLabelMapping[( self.bondType, self.isAromatic )]
        elif self.bondTypeLabel is not None: ( self.bondType, self.isAromatic ) = _BondLabelObjectMapping[self.bondTypeLabel]

    @staticmethod
    def MakeKeys ( linkLabel, leftComponentLabel, rightComponentLabel ):
        """Make keys."""
        # . Get the labels.
        label1 = _KeyTokenSeparator.join ( leftComponentLabel.lower  ( ).split ( ) )
        label2 = _KeyTokenSeparator.join ( rightComponentLabel.lower ( ).split ( ) )
        # . Get the keys.
        key00 = _KeyTokenSeparator.join ( linkLabel.lower ( ).split ( ) )
        keyX0 = _KeyTokenSeparator.join ( [ key00, label1        ] )
        key0Y = _KeyTokenSeparator.join ( [ key00, ""   , label2 ] )
        keyXY = _KeyTokenSeparator.join ( [ keyX0,        label2 ] )
        # . Most specific to most general.
        return ( keyXY, keyX0, key0Y, key00 )

    def ToMMSequenceObject ( self ):
        """Make an MM sequence object."""
        leftVariant = self.leftVariant.ToMMSequenceObject ( )
        if self.leftVariant is self.rightVariant: rightVariant = leftVariant
        else: rightVariant = self.rightVariant.ToMMSequenceObject ( )
        mmItem = MMSequenceLink.WithOptions ( label = self.label, leftVariant = leftVariant, rightVariant = rightVariant )
        return mmItem

    # . Key property.
    @property
    def key ( self ):
        """Get a key."""
        if self.__dict__.get ( "key", None ) is None:
            tokens = self.label.lower ( ).split ( )
            label1 = getattr ( self.leftVariant , "componentLabel", None )
            label2 = getattr ( self.rightVariant, "componentLabel", None )
            if label1 is not None: tokens1 = label1.lower ( ).split ( )
            if label2 is not None: tokens2 = label2.lower ( ).split ( )
            if label1 is None:
                if label2 is not None: tokens.extend ( [ "" ] + tokens2 )
            else:
                tokens.extend ( tokens1 )
                if label2 is not None: tokens.extend ( tokens2 )
            self.__dict__["key"] = _KeyTokenSeparator.join ( tokens )
        return self.__dict__["key"]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentVariant ( AttributableObject ):
    """Variant of a PDB component."""

    _attributable = dict ( AttributableObject._attributable )
    _mappable     = { "Atoms To Delete"      : "atomsToDelete"      ,
                      "Bonds To Delete"      : "bondsToDelete"      ,
                      "Component"            : "componentLabel"     ,
                      "Formal Charges"       : "formalCharges"      ,
                      "Is Chain Terminating" : "isChainTerminating" ,
                      "Label"                : "label"              }
    _attributable.update ( { "atomsToAddFields"     : None ,
                             "atomsToAdd"           : None , # . List of PDBComponentAtom objects.
                             "atomsToDelete"        : None , # . List of atom labels.
                             "bondsToAddFields"     : None ,
                             "bondsToAdd"           : None , # . List of PDBComponentBond objects.
                             "bondsToDelete"        : None , # . List of lists of atom label pairs.
                             "bondTypesFields"      : None ,
                             "bondTypes"            : None , # . List of PDBComponentBond objects.
                             "componentLabel"       : None ,
                             "formalCharges"        : None , # . Mapping of formal charges to atom labels
                             "isChainTerminating"   : None ,
                             "label"                : None } )

    #yaml_tag = "!PDBComponentVariant"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Atom and bonds.
        for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atomsToAdd", "atomsToAddFields", "Atoms To Add", "Atom To Add Fields", PDBComponentAtom ),
                                                              ( "bondsToAdd", "bondsToAddFields", "Bonds To Add", "Bond To Add Fields", PDBComponentBond ),
                                                              ( "bondTypes" , "bondTypesFields" , "Bond Types"  , "Bond Type Fields"  , PDBComponentBond ) ):
            items = getattr ( self, tag, None )
            if items is not None:
                listItems = []
                keys      = getattr ( self, keyTag, None )
                if keys is None: keys = itemObject._mappable0
                for item in items:
                    listItems.append ( item.ToList ( keys ) )
                mapping[label   ] = listItems
                mapping[keyLabel] = keys
        # . Other attributes.
        for ( fullKey, key ) in self.__class__._mappable.items ( ):
            attribute = getattr ( self, key, None )
            if attribute is not None: mapping[fullKey] = attribute
        # . Finish up.
        return mapping

    def __setstate__ ( self, mapping ):
        """Set state from a mapping."""
        # . There are assumed to be no errors!
        try:
            # . Atom and bonds.
            for ( tag, keyTag, label, keyLabel, itemObject ) in ( ( "atomsToAdd", "atomsToAddFields", "Atoms To Add", "Atom To Add Fields", PDBComponentAtom ) ,
                                                                  ( "bondsToAdd", "bondsToAddFields", "Bonds To Add", "Bond To Add Fields", PDBComponentBond ) ,
                                                                  ( "bondTypes" , "bondTypesFields" , "Bond Types"  , "Bond Type Fields"  , PDBComponentBond ) ):
                listItems = mapping.get ( label, None )
                if listItems is not None:
                    items = []
                    keys  = mapping[keyLabel]
                    for listItem in listItems:
                        items.append ( itemObject.FromList ( listItem, keys ) )
                    setattr ( self, tag   , items )
                    setattr ( self, keyTag, keys  )
            # . Other attributes.
            for ( fullKey, key ) in self.__class__._mappable.items ( ):
                attribute = mapping.pop ( fullKey, None )
                if attribute is not None: setattr ( self, key, attribute )
#        except Exception as e:
#            print e[0]
        except:
            raise PDBComponentError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    @staticmethod
    def MakeKeys ( variantLabel, componentLabel ):
        """Make keys."""
        tokens0 = variantLabel.lower   ( ).split ( )
        tokens  = tokens0 + componentLabel.lower ( ).split ( )
        key0    = _KeyTokenSeparator.join ( tokens0 )
        key     = _KeyTokenSeparator.join ( tokens  )
        # . Most specific to most general.
        return ( key, key0 )

    def ToMMSequenceObject ( self ):
        """Make an MM sequence object."""
        mmAtoms = []
        if self.atomsToAdd is not None:
            for atom in self.atomsToAdd:
                mmAtoms.append ( MMSequenceAtom.WithOptions ( charge = 0.0, label = atom.label, typeLabel = atom.label ) ) 
        mmItem = MMSequenceVariant.WithOptions ( atoms = mmAtoms, componentLabel = self.componentLabel, label = self.label )
        return mmItem

    # . Key property.
    @property
    def key ( self ):
        """Get a key."""
        if self.__dict__.get ( "key", None ) is None:
            tokens = self.label.lower ( ).split ( )
            label  = getattr ( self, "componentLabel", None )
            if label is not None: tokens.extend ( label.lower ( ).split ( ) )
            self.__dict__["key"] = _KeyTokenSeparator.join ( tokens )
        return self.__dict__["key"]

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

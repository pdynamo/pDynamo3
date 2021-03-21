"""Classes and functions for handling a PDB component library."""

import os, os.path

from  pCore                     import AttributableObject         , \
                                       logFile                    , \
                                       LogFileActive              , \
                                       YAMLMappingFile_FromObject , \
                                       YAMLMappingFile_ToObject   , \
                                       YAMLPickleFileExtension
from  pMolecule                 import BondType
from .PDBComponent              import PDBComponent               , \
                                       PDBComponentAtom           , \
                                       PDBComponentBond           , \
                                       PDBComponentLink           , \
                                       PDBComponentVariant
from .PDBComponentCIFFileReader import PDBComponentCIFFileReader

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Amino acids - note UNK is treated as an amino acid in the PDB standard!
_AminoAcids = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "UNK", "VAL" ]

# . Default variants for the amino acids.
# . All amino acids appear to be fully protonated in the cif file. This means ARG and LYS do not need defaults.
_DefaultVariants = { "ASP" : "Deprotonated", "GLU" : "Deprotonated", "HIS" : "Delta Protonated" }

# . Reduced set of components to save in the default distribution - amino acids, water and some counterions.
_ReducedComponentSet = set ( _AminoAcids + [ "HOH", "CL", "K", "NA" ] )

# . Paths.
_ComponentPath      = "components"
_DefaultCIFFileName = "components.cif"
_LibraryPath        = "pdbComponents"
_LinkPath           = "links"
_VariantPath        = "variants"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PDBComponentLibrary ( AttributableObject ):
    """Library of PDB components, links and variants."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "componentPaths" :  None ,
                             "components"     :  None ,
                             "isSetup"        : False ,
                             "libraryPaths"   :  None ,
                             "linkPaths"      :  None ,
                             "links"          :  None ,
                             "paths"          :  None ,
                             "readOnly"       :  True ,
                             "variantPaths"   :  None ,
                             "variants"       :  None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( PDBComponentLibrary, self )._CheckOptions ( )
        self.DefineLibraryPaths ( self.paths )
        self.SetupPaths ( makePaths = not self.readOnly )
        self.ResetCache ( )

    def AddItems ( self, items, libraryPath = None ):
        """Add items to the library."""
        # . Find library path.
        if libraryPath is None:
            libraryPath = self.libraryPaths[0]
        # . Loop over items.
        for item in items:
            # . Identify item type and path.
            if   isinstance ( item, PDBComponent        ): path = os.path.join ( libraryPath, _ComponentPath )
            elif isinstance ( item, PDBComponentLink    ): path = os.path.join ( libraryPath, _LinkPath      )
            elif isinstance ( item, PDBComponentVariant ): path = os.path.join ( libraryPath, _VariantPath   )
            else: raise TypeError ( "Invalid PDB component, link or variant." )
            # . Save the item.
            YAMLMappingFile_FromObject ( os.path.join ( path, item.key + YAMLPickleFileExtension ), "!" + item.__class__.__name__, item )

    def DefineLibraryPaths ( self, paths ):
        """Define the library paths."""
        if paths is None: self.libraryPaths = [ self.MakeStandardLibraryPath ( ) ]
        else:             self.libraryPaths = paths

    def GetComponent ( self, componentLabel, missingItems = None ):
        """Get a component."""
        key = PDBComponent.MakeKey ( componentLabel )
        if key not in self.components:
            for componentPath in self.componentPaths:
                fileName = os.path.join ( componentPath, key + YAMLPickleFileExtension )
                if os.path.exists ( fileName ):
                    self.components[key] = YAMLMappingFile_ToObject ( fileName, PDBComponent )
                    break
        component = self.components.get ( key, None )
        if ( component is None ) and ( missingItems is not None ): missingItems.add ( ( "Component", key ) )
        return component

    def GetLink ( self, linkLabel, leftComponentLabel, rightComponentLabel, missingItems = None ):
        """Get a link."""
        keys = PDBComponentLink.MakeKeys ( linkLabel, leftComponentLabel, rightComponentLabel )
        key  = keys[0]
        if key not in self.links:
            found = False
            for tag in keys:
                if found: break
                for linkPath in self.linkPaths:
                    fileName = os.path.join ( linkPath, tag + YAMLPickleFileExtension )
                    if os.path.exists ( fileName ):
                        self.links[key] = YAMLMappingFile_ToObject ( fileName, PDBComponentLink )
                        found = True
                        break
        link = self.links.get ( key, None )
        if ( link is None ) and ( missingItems is not None ): missingItems.add ( ( "Link", key ) )
        return link

    def GetVariant ( self, variantLabel, componentLabel, missingItems = None ):
        """Get a variant."""
        keys = PDBComponentVariant.MakeKeys ( variantLabel, componentLabel )
        key  = keys[0]
        if key not in self.variants:
            found = False
            for tag in keys:
                if found: break
                for variantPath in self.variantPaths:
                    fileName = os.path.join ( variantPath, tag + YAMLPickleFileExtension )
                    if os.path.exists ( fileName ):
                        self.variants[key] = YAMLMappingFile_ToObject ( fileName, PDBComponentVariant )
                        found = True
                        break
        variant = self.variants.get ( key, None )
        if ( variant is None ) and ( missingItems is not None ): missingItems.add ( ( "Variant", key ) )
        return variant

    def MakeStandardLibraryPath ( self ):
        """Make the standard library path."""
        return os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), _LibraryPath )

    def ResetCache ( self ):
        """Clear stored components, links and variants."""
        self.components = {}
        self.links      = {}
        self.variants   = {}

    def SetupPaths ( self, makePaths = False ):
        """Set up all paths."""
        # . Define paths.
        paths = list ( self.libraryPaths )
        for ( attribute, tailPath ) in ( ( "componentPaths", _ComponentPath ) ,
                                         ( "linkPaths"     , _LinkPath      ) ,
                                         ( "variantPaths"  , _VariantPath   ) ):
            localPaths = []
            for rootPath in self.libraryPaths:
                path = os.path.join ( rootPath, tailPath )
                localPaths.append ( path )
            paths.extend ( localPaths )
            setattr ( self, attribute, localPaths )
        # . Make paths.
        if makePaths:
            for path in paths:
                if not os.path.exists ( path ): os.mkdir ( path )

#===================================================================================================================================
# . Functions for making the default library.
#===================================================================================================================================
def _MakeDefaultLinks ( ):
    """Make links for the default distribution."""
    # . Initialization.
    links = []
    # . Disulfide bridge.
    variant = PDBComponentVariant.WithOptions ( atomsToDelete = [ "HG" ], componentLabel = "CYS" )
    links.append ( PDBComponentLink.WithOptions ( label = "Disulfide Bridge", leftAtomLabel = "SG", leftVariant = variant, rightAtomLabel = "SG", rightVariant = variant, bondType = BondType.Single ) )
    # . Peptide bonds.
    # . The left variant is always the same.
    leftVariant = PDBComponentVariant.WithOptions ( atomsToDelete = [ "HXT", "OXT" ] )
    # . General case, PRO and UNK.
    rightVariant = PDBComponentVariant.WithOptions ( atomsToDelete = [ "H2" ] )
    links.append ( PDBComponentLink.WithOptions ( label = "Peptide", leftAtomLabel = "C", leftVariant = leftVariant, rightAtomLabel = "N", rightVariant = rightVariant, bondType = BondType.Single ) )
    rightVariant = PDBComponentVariant.WithOptions ( atomsToDelete = [ "H"  ], componentLabel = "PRO" )
    links.append ( PDBComponentLink.WithOptions ( label = "Peptide", leftAtomLabel = "C", leftVariant = leftVariant, rightAtomLabel = "N", rightVariant = rightVariant, bondType = BondType.Single ) )
    rightVariant = PDBComponentVariant.WithOptions ( atomsToDelete = [ "H2" ], componentLabel = "UNK" )
    links.append ( PDBComponentLink.WithOptions ( label = "Peptide", leftAtomLabel = "C", leftVariant = leftVariant, rightAtomLabel = "N", rightVariant = rightVariant, bondType = BondType.Single ) )
    # . Finish up.
    return links

def _MakeDefaultVariants ( ):
    """Make variants for the default distribution."""
    # . Initialization.
    variants = []
    # . CTerminal.
    variants.append ( PDBComponentVariant.WithOptions ( label = "C Terminal", atomsToDelete = [ "HXT" ], formalCharges = { "OXT" : -1 } ) )
    # . NTerminal - general case and UNK.
    atomsToAdd = [ PDBComponentAtom.WithOptions ( atomicNumber = 1, label = "H1", pdbAlign = 0, toFollow = "N"  ),
                   PDBComponentAtom.WithOptions ( atomicNumber = 1, label = "H2", pdbAlign = 0, toFollow = "H1" ),
                   PDBComponentAtom.WithOptions ( atomicNumber = 1, label = "H3", pdbAlign = 0, toFollow = "H2" )  ]
    bondsToAdd = [ PDBComponentBond.WithOptions ( atomLabel1 = "H1", atomLabel2 = "N", bondType = BondType.Single ),
                   PDBComponentBond.WithOptions ( atomLabel1 = "H2", atomLabel2 = "N", bondType = BondType.Single ),
                   PDBComponentBond.WithOptions ( atomLabel1 = "H3", atomLabel2 = "N", bondType = BondType.Single )  ]
    variants.append ( PDBComponentVariant.WithOptions (                         label = "N Terminal", atomsToAdd = atomsToAdd, atomsToDelete = [ "H", "H2" ], formalCharges = { "N" : +1 }, bondsToAdd = bondsToAdd ) )
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "UNK", label = "N Terminal", atomsToAdd = atomsToAdd, atomsToDelete = [ "H", "H2" ], formalCharges = { "N" : +1 }, bondsToAdd = bondsToAdd ) )
    # . NTerminal - PRO.
    atomsToAdd = [ PDBComponentAtom.WithOptions ( atomicNumber = 1, label = "H2", pdbAlign = 0, toFollow = "N"  ),
                   PDBComponentAtom.WithOptions ( atomicNumber = 1, label = "H3", pdbAlign = 0, toFollow = "H2" )  ]
    bondsToAdd = [ PDBComponentBond.WithOptions ( atomLabel1 = "H2", atomLabel2 = "N", bondType = BondType.Single ),
                   PDBComponentBond.WithOptions ( atomLabel1 = "H3", atomLabel2 = "N", bondType = BondType.Single )  ]
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "PRO", label = "N Terminal", atomsToAdd = atomsToAdd, atomsToDelete = [ "H" ], formalCharges = { "N" : +1 }, bondsToAdd = bondsToAdd ) )
    # . Deprotonated ASP and GLU.
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "ASP", label = "Deprotonated", atomsToDelete = [ "HD2" ], formalCharges = { "OD2" : -1 } ) )
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "GLU", label = "Deprotonated", atomsToDelete = [ "HE2" ], formalCharges = { "OE2" : -1 } ) )
    # . Histidine variants.
    # . Doubly protonated - leave as is but shift charge to Nepsilon.
    bondTypes = [ PDBComponentBond.WithOptions ( atomLabel1 = "CE1", atomLabel2 = "ND1", bondType = BondType.Single ),
                  PDBComponentBond.WithOptions ( atomLabel1 = "CE1", atomLabel2 = "NE2", bondType = BondType.Double ) ]
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "HIS", label = "Doubly Protonated",                             formalCharges = { "ND1" : 0, "NE2" : +1 }, bondTypes = bondTypes ) )
    # . Delta protonated.
    bondTypes = [ PDBComponentBond.WithOptions ( atomLabel1 = "CE1", atomLabel2 = "ND1", bondType = BondType.Single ),
                  PDBComponentBond.WithOptions ( atomLabel1 = "CE1", atomLabel2 = "NE2", bondType = BondType.Double ) ]
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "HIS", label = "Delta Protonated",   atomsToDelete = [ "HE2" ], formalCharges = { "ND1" : 0, "NE2" :  0 }, bondTypes = bondTypes ) )
    # . Epsilon protonated.
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "HIS", label = "Epsilon Protonated", atomsToDelete = [ "HD1" ], formalCharges = { "ND1" : 0, "NE2" :  0 } ) )
    # . Fully deprotonated.
    variants.append ( PDBComponentVariant.WithOptions ( componentLabel = "HIS", label = "Fully Deprotonated", atomsToDelete = [ "HD1", "HE2" ], formalCharges = { "ND1" : 0, "NE2" : -1 } ) )
    # . Finish up.
    return variants

def _ModifyAminoAcidComponents ( components ):
    """Modify the amino acid components for the default distribution."""
    for label in _AminoAcids:
        component = components[label]
        component.leftAtom         = "N"
        component.leftLink         = "Peptide"
        component.leftTermination  = "N Terminal"
        component.isInChain        = True
        component.isHeteroatom     = False
        component.rightAtom        = "C"
        component.rightLink        = "Peptide"
        component.rightTermination = "C Terminal"
        defaultvariant = _DefaultVariants.get ( component.label, None )
        if defaultvariant is not None: component.variants = [ defaultvariant ]

def MakeDefaultPDBComponentLibrary ( cifPath = None, fullLibrary = False, log = logFile, libraryPaths = None, outPath = None ):
    """Make the default library that comes with the pDynamo distribution."""
    # . Get the cifPath.
    if cifPath is None:
        cifPath = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), _LibraryPath, _DefaultCIFFileName )
    # . Get the components.
    components = PDBComponentCIFFileReader.PathToComponents ( cifPath, asDictionary = True, log = log )
    # . Process the amino acid components.
    _ModifyAminoAcidComponents ( components )
    # . Create the default links and variants.
    links    = _MakeDefaultLinks    ( )
    variants = _MakeDefaultVariants ( )
    # . Get the items to put in the library.
    # . All items.
    if fullLibrary:
        items = list ( components.values ( ) )
        items.extend ( links    )
        items.extend ( variants )
    # . Items in the reduced set only.
    else:
        items = []
        # . Components.
        for label in _ReducedComponentSet: items.append ( components[label] )
        # . Links.
        for item in links:
            if ( ( item.leftVariant.componentLabel  is None ) or ( item.leftVariant.componentLabel  in _ReducedComponentSet ) ) and \
               ( ( item.rightVariant.componentLabel is None ) or ( item.rightVariant.componentLabel in _ReducedComponentSet ) ): items.append ( item )
        # . Variants.
        for item in variants:
            if ( item.componentLabel is None ) or ( item.componentLabel in _ReducedComponentSet ): items.append ( item )
    # . Create the library.
    library = PDBComponentLibrary.WithOptions ( paths = libraryPaths, readOnly = False )
    library.AddItems ( items, libraryPath = outPath )
    return library

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":

    # . Make the default library.
    MakeDefaultPDBComponentLibrary ( fullLibrary = False )

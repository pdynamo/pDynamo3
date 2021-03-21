"""Classes and functions for writing SMILES strings."""

import string

from  pCore             import AttributableObject       , \
                               logFile
from  pMolecule         import BondType                 , \
                               CheckForKekuleOutput
from  pScientific       import PeriodicTable
from .SMILESDefinitions import _BondTypeToSymbol        , \
                               _ChiralityDefaultClasses , \
                               _ElementsOrganic         , \
                               _NullBondToken           , \
                               _ValenciesOrganic

#===================================================================================================================================
# . SMILES writer class.
#===================================================================================================================================
class SMILESWriter ( AttributableObject ):
    """SMILESWriter is the class for writing SMILES strings."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "removeExplicitHydrogens" : True  ,
                             "useKekule"               : False } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( SMILESWriter, self )._CheckOptions ( )
        self.Reset ( )

    def AtomToken ( self, atom, bonds ):
        """Return a SMILES token for an atom."""
        # . Gather all data.
        atomicNumber      = getattr ( atom, "atomicNumber"     ,    -1 )
        chiralityClass    = getattr ( atom, "chiralityClass"   ,  None )
        chiralityNumber   = getattr ( atom, "chiralityNumber"  ,     0 )
        connections       = getattr ( atom, "connections"      ,     0 )
        formalCharge      = getattr ( atom, "formalCharge"     ,     0 )
        implicitHydrogens = getattr ( atom, "implicitHydrogens",     0 )
        isAromatic        = getattr ( atom, "isAromatic"       , False ) and ( not self.forceKekule )
        isotope           = getattr ( atom, "isotope"          ,     0 )
        # . Remove explicit hydrogens.
        if self.removeExplicitHydrogens:
            currentBonds = []
            for bond in bonds:
                other = bond.Opposite ( atom )
                if ( getattr ( other, "atomicNumber"     , -1 ) == 1 ) and \
                   ( getattr ( other, "connections"      ,  0 ) == 1 ) and \
                   ( getattr ( other, "formalCharge"     ,  0 ) == 0 ) and \
                   ( getattr ( other, "isotope"          ,  0 ) == 0 ):
                    implicitHydrogens += 1
                    self.qAtoms[other] = True
                else:
                    currentBonds.append ( bond )
        # . Keep explicit hydrogens - in other words do nothing.
        else:
            currentBonds = bonds
        # . Check for a reduced atom representation.
        # . Basic checks.
        isReduced = ( atomicNumber in _ElementsOrganic ) and ( formalCharge == 0 ) and ( isotope <= 0 )
        # . Check for a reduced chirality representation.
        if chiralityClass is not None:
            isChiralityReduced = ( _ChiralityDefaultClasses.get ( connections, None ) == chiralityClass ) and ( chiralityNumber <= 4 )
            isReduced          = isReduced and isChiralityReduced
            # . Generate the string.
            if isChiralityReduced: cToken = chiralityNumber * "@"
            else:                  cToken = "@" + chiralityClass + "{:d}".format ( chiralityNumber )
        else: cToken = ""
        # . Encode the string.
        if isAromatic: aToken = PeriodicTable.Symbol ( atomicNumber ).lower ( )
        else:          aToken = PeriodicTable.Symbol ( atomicNumber )
        tokens = [ aToken, cToken ]
        if not isReduced:
            tokens[0:0] = "["
            if isotope != 0: tokens[1:1] = "{:d}".format ( isotope )
            if   implicitHydrogens == 1: tokens.append ( "H" )
            elif implicitHydrogens  > 1: tokens.append ( "H" + "{:d}".format ( implicitHydrogens ) )
            if   formalCharge       < 0: tokens.append (       "{:d}".format ( formalCharge      ) )
            elif formalCharge       > 0: tokens.append ( "+" + "{:d}".format ( formalCharge      ) )
            tokens.append ( "]" )
        return ( "".join ( tokens ), currentBonds )

    def CreateBranch ( self, current, previous, bonds ):
        """Return tokens for a branch."""
        self.qAtoms[current]           = True
        ( currentToken, currentBonds ) = self.AtomToken ( current, bonds[current] )
        tokens                         = [ currentToken ]
        branches                       = []
        for bond in currentBonds:
            other     = bond.Opposite ( current )
            bondOrder = bond.type.bondOrder
            bondToken = None
            # . Bond token cases.
            # . Single and aromatic bonds are never represented except for cases where there
            # . is a single non-aromatic bond between two aromatic atoms (e.g. biphenyl).
            if ( self.forceKekule and ( bondOrder > 1 ) ) or \
               ( ( not bond.isAromatic ) and \
               ( ( current.isAromatic and other.isAromatic and ( bondOrder > 0 ) ) or ( bondOrder > 1 ) ) ):
                bondToken = _BondTypeToSymbol[bond.type]
            if self.qAtoms.get ( other, False ):
                if other is not previous:
                    if bondToken is not None: tokens.append ( bondToken )
                    crossLink = self.qBonds.pop ( bond, self.crossLinks + 1 )
                    if crossLink > self.crossLinks:
                        self.crossLinks += 1
                        self.qBonds[bond] = self.crossLinks
                    tokens.append ( "%" + "{:d}".format ( crossLink ) )
            else: # . Need to explicitly put single bond if two aromatic atoms but only single bond (e.g. biphenyl).
                branch = []
                if bondToken is not None: branch.append ( bondToken )
                branch.extend ( self.CreateBranch ( other, current, bonds ) )
                branches.append ( branch )
        if len ( branches ) > 0:
            for branch in branches[0:-1]: tokens.extend ( [ "(" ] + branch + [ ")" ] )
            tokens.extend ( branches[-1] )
        return tokens

    def FromConnectivity ( self, connectivity ):
        """Generate a SMILES from a connectivity."""
        # . Initialization.
        self.Reset ( )
        tokens = []
        # . Loop over isolates.
        self.forceKekule = CheckForKekuleOutput ( connectivity, self.useKekule )
        for current in connectivity.atoms:
            if not self.qAtoms.get ( current, False ):
                isolate = self.CreateBranch ( current, None, connectivity.adjacentEdges )
                if len ( tokens ) > 0: tokens.append ( _NullBondToken )
                tokens.extend ( isolate )
        self.Reset ( )
        return "".join ( tokens )

    def Reset ( self ):
        """Reset all attributes."""
        self.qAtoms      = {}
        self.qBonds      = {}
        self.crossLinks  = 0
        self.forceKekule = False

    @classmethod
    def StringFromSystem ( selfClass, system, **options ):
        """Generate a SMILES string for a system."""
        writer = selfClass.WithOptions ( **options )
        return writer.FromConnectivity ( system.connectivity )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

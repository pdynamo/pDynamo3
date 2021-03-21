"""Miscellaneous SMILES parameters and functions."""

from pCore       import logFile       , \
                        LogFileActive
from pMolecule   import Atom          , \
                        BondType      , \
                        System
from pScientific import PeriodicTable

#===================================================================================================================================
# . SMILES parameters.
#===================================================================================================================================
# . Bond type symbols.
# . The flag is for aromaticity.
_BondTypeFromSymbol = { "."  : ( BondType.Null      , False ) ,
                        "-"  : ( BondType.Single    , False ) ,
                        "/"  : ( BondType.Single    , False ) ,
                        "\\" : ( BondType.Single    , False ) ,
                        ":"  : ( BondType.Undefined , True  ) ,
                        "="  : ( BondType.Double    , False ) ,
                        "#"  : ( BondType.Triple    , False ) ,
                        "$"  : ( BondType.Quadruple , False ) }
_BondTypeToSymbol   = { BondType.Null      : "."  ,
                        BondType.Single    : "-"  ,
                        BondType.Double    : "="  ,
                        BondType.Triple    : "#"  ,
                        BondType.Quadruple : "$"  }
_AromaticBondToken  = ":"
_NullBondToken      = "."

# . Charge.
_ChargeTokens = { "+" :  1, "++" :  2, "+++" :  3, "++++" :  4, "-" : -1, "--" : -2, "---" : -3, "----" : -4 }

# . Chirality data.
_ChiralityFull    = { "AL" : "AL", "OH" : "OH", "SP" : "SP", "TB" : "TB", "TH" : "TH" }
_ChiralityReduced = { "@" : ( "??", 1 ), "@@" : ( "??", 2 ), "@@@" : ( "??", 3 ), "@@@@" : ( "??", 4 ) }

# . The allowed number of connections for each class.
_ChiralityClassConnections = { "AL" : 2, "OH" : 6, "SP" : 4, "TB" : 5, "TH" : 4 }

# . The default class for each connectivity.
_ChiralityDefaultClasses = { 2 : "AL", 4 : "TH", 5 : "TB", 6 : "OH" }

# . The maximum number permitted for each class.
_ChiralityMaximumNumbers = { "AL" : 2, "OH" : 30, "SP" : 3, "TB" : 20, "TH" : 2 }

# . Element data.
# . Atomic numbers of aromatic elements.
_ElementsAromatic = ( 5, 6, 7, 8, 15, 16, 33, 34, 51, 52 )

# . Organic elements.
_ElementsOrganic  = ( 5, 6, 7, 8, 9, 15, 16, 17, 35, 53 )

# . Element tokens.
_ElementTokens        = { PeriodicTable.Symbol ( n )           : ( n, False ) for n in PeriodicTable.atomicNumbers }
_ElementTokens.update ( { PeriodicTable.Symbol ( n ).lower ( ) : ( n, True  ) for n in _ElementsAromatic           } )

# . Hydrogen count.
_HCountTokens = { "H" : 1 }

# . Reduced element tokens (outside of square brackets).
_ReducedElementTokens        = { PeriodicTable.Symbol ( n )           : ( n, False ) for n in       _ElementsOrganic                               }
_ReducedElementTokens.update ( { PeriodicTable.Symbol ( n ).lower ( ) : ( n, True  ) for n in set ( _ElementsOrganic ) & set ( _ElementsAromatic ) } )

# . The maximum possible valencies for aromatic atoms.
_ValenciesAromatic = { 5 : 3, 6 : 4, 7 : 5, 8 : 2, 15 : 5, 16 : 6, 33 : 5, 34 : 6, 51 : 5, 52 : 6 }

# . The possible valencies for atoms in the organic subset when connected hydrogens are specified implicitly.
# . Note that the minimum number of hydrogens are added so as to be consistent with explicitly-specified bonds.
_ValenciesOrganic = { 5: ( 3, ), 6: ( 4, ), 7: ( 3, 5 ), 8: ( 2, ), 9: ( 1, ), 15: ( 3, 5 ), 16: ( 2, 4, 6 ), 17: ( 1, ), 35: ( 1, ), 53: ( 1, ) }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""DYFF utilities."""

import math

from  enum           import Enum
from  pMolecule.Bond import BondType
from  pScientific    import PeriodicTable  , \
                            Units
from .MMModelError   import MMModelError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Hybridizations.
_Octahedral   = "Oct"
_Resonant     = "Res"
_SquarePlanar = "SPl"
_Tetrahedral  = "Tet"
_Trigonal     = "Tri"

# . Bond orders.
_AmideBondOrder    = 1.41   # . N-C=O.
_AromaticBondOrder = 0.50   # . To be added or subtracted from non-aromatic bond order ( 1 -> 1.5, 2 -> 1.5, 3 -> 2.5 ).
_BondOrderScaling  = 0.1332 # . The scaling factor for the bond-order term in the DYFF expression for bond equilibrium values.

# . Angles.
_AnglePiTolerance  = 0.01   # . The tolerance in radians for determining whether an angle has a value of pi.

# . Dihedrals.
# . All values from original paper.
_Dihedral_Group16sp3       = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole *  6.8
_Dihedral_Osp3             = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole *  2.0
_Dihedral_PlanarPlanar     = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole *  5.0
_Dihedral_ResonantResonant = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole * 25.0
_Dihedral_sp2_A            = 5.0
_Dihedral_sp2_B            = 4.18
_Dihedral_sp3sp2           = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole *  1.0
_Dihedral_sp3sp2sp2        = Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole *  2.0

# . Out-of-planes.
#
# . These terms have been reformulated due to numerical imprecisions in the formulation of the original paper.
#
#   Each term is now E = c0 + c2 * cos 2g + c4 * cos 4g where g is the angle between the ij bond and the normal to the jkl plane.
#   The sum is designed so that E is periodic in the range [0,pi/2] and identical for both + and - angles. I.e. the direction of
#   ij wrt the jk and jl bonds is irrelevant. Unchemical conformations will be prevented by the presence of other terms such as
#   bond angles.
#
#   In the general case parametrization is done using three conditions:
#
#   E(g0)       = c0 + c2 cos 2g0 + c4 cos 4g0 = 0        the minimum energy value at g=g0
#   E(pi/2)     = c0 - c2 + c4 = Eb                       the inversion barrier Eb at g=pi/2
#   dE(g0)/dg   = 0 => c2 sin 2g0 + 2 c4 sin 4g0 = 0      the minimum condition at g=g0
#
#   When g0 = pi/2, c0 = c2 = 0.5 * Eb and c4 = 0.
#
# . All force constants are divided by three because there are three inversion terms per center.
#
# . Force constants in kJ/mol and periodicities.
_OOP_CNsp2     = ( (  4.184, 0 ), (   4.184, 2 )               ) # . C/N    with w0    = pi/2  and Eb =  6 kcal/mol at w = 0.
_OOP_CNsp2Osp2 = ( ( 34.867, 0 ), (  34.867, 2 )               ) # . C/N=O  with w0    = pi/2  and Eb = 50 kcal/mol at w = 0.
_OOP_Nsp3      = ( (  2.792, 0 ), (  -3.803, 2 ), ( 1.772, 4 ) ) # . N  sp3 with theta = 106.7 and Eb =  6 kcal/mol at w = pi/2.
_OOP_Psp3      = ( ( 11.445, 0 ), ( -15.340, 2 ), ( 3.897, 4 ) ) # . P  sp3 with theta =  93.5 and Eb = 22 kcal/mol at w = pi/2.
_OOP_Assp3     = ( ( 11.490, 0 ), ( -15.341, 2 ), ( 3.851, 4 ) ) # . As sp3 with theta =  91.8 and Eb = 22 kcal/mol at w = pi/2.
_OOP_Sbsp3     = ( ( 11.492, 0 ), ( -15.341, 2 ), ( 3.849, 4 ) ) # . Sb sp3 with theta =  91.7 and Eb = 22 kcal/mol at w = pi/2.
_OOP_Bisp3     = ( ( 11.505, 0 ), ( -15.341, 2 ), ( 3.837, 4 ) ) # . Bi sp3 with theta =  90.5 and Eb = 22 kcal/mol at w = pi/2.
_OOP_Mcsp3     = ( ( 11.505, 0 ), ( -15.341, 2 ), ( 3.837, 4 ) ) # . Mc sp3 with theta =  90.5 and Eb = 22 kcal/mol at w = pi/2.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Unused at present.
class CoordinationType ( Enum ):
    """Coordination types."""
    CubicAntiprism         = ( "CubicAntiprism"         , 8 , "CAp" )
    Linear                 = ( "Linear"                 , 2 , "Lin" )
    Octahedral             = ( "Octahedral"             , 6 , "Oct" )
    PentagonalBipyramidal  = ( "PentagonalBipyramidal"  , 7 , "PBp" )
    Resonant               = ( "Resonant"               , 3 , "Res" )
    SquarePlanar           = ( "SquarePlanar"           , 4 , "SPl" )
    SquarePyramidal        = ( "SquarePyramidal"        , 5 , "SPy" )
    Tetrahedral            = ( "Tetrahedral"            , 4 , "Tet" ) 
    TricappedTrigonalPrism = ( "TricappedTrigonalPrism" , 9 , "TTP" )
    Trigonal               = ( "Trigonal"               , 3 , "Tri" )
    TrigonalBipyramidal    = ( "TrigonalBipyramidal"    , 5 , "TBp" )
    Undefined              = ( "Undefined"              , 0 , ""    )

    def __init__ ( self, label, neighbors, tag ):
        """Constructor."""
        self.label     = label
        self.neighbors = neighbors
        self.tag       = tag

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def DYFFBondOrder ( i, j, connectivity ):
    """The DYFF bond order."""
    # . Get an initial guess for the bond order (correct most of the time).
    iAtom    = connectivity.nodes[i]
    jAtom    = connectivity.nodes[j]
    bond     = connectivity.GetEdge ( iAtom, jAtom )
    if bond is None: raise MMModelError ( "Bond not found between atoms {:d} and {:d}.".format ( i, j ) )
    bondType = bond.type
    bO       = float ( bondType.bondOrder )
    # . Special cases.
    # . Aromatic - assumes one electron extra in each bond (which is clearly not always the case).
    if bond.isAromatic:
        if bO == 1.0: bO += _AromaticBondOrder
        else:         bO -= _AromaticBondOrder
    # . Possible amide - N-C=O.
    else:
        nI = iAtom.atomicNumber
        nJ = jAtom.atomicNumber
        if ( bondType is BondType.Single ) and ( ( ( nI == 6 ) and ( nJ == 7 ) ) or ( ( nI == 7 ) and ( nJ == 6 ) ) ):
            if nI == 6: cAtom = iAtom
            else:       cAtom = jAtom
            for bond in connectivity.adjacentEdges[cAtom]:
                other = bond.Opposite ( cAtom )
                if ( other.atomicNumber == 8 ) and ( bondType is BondType.Double ):
                    bO = _AmideBondOrder
                    break
    return bO

def DYFFHybridization ( label ):
    """The hybridization from an atom type label."""
    tokens = label.split ( ":" )
    if len ( tokens ) > 1: return tokens[1]
    else:                  return ""

def DYFFIsResonant ( label ):
    """The resonance of the atom type."""
    return ( DYFFHybridization ( label ) == _Resonant )

def DYFFNaturalBondLength ( dyffAtoms, iT, jT, bondOrder ):
    """The natural bond length between two elements."""
    # . Atom type data.
    iAtom = dyffAtoms.GetItem ( iT )
    jAtom = dyffAtoms.GetItem ( jT )
    rI    = iAtom.GetProperty ( "Bond Radius"       )
    rJ    = jAtom.GetProperty ( "Bond Radius"       )
    xI    = iAtom.GetProperty ( "Electronegativity" )
    xJ    = jAtom.GetProperty ( "Electronegativity" )
    # . Factors.
    rIJ0 = rI + rJ
    rBO  = - _BondOrderScaling * rIJ0 * math.log ( bondOrder )
    rEN  = ( rI * rJ ) * ( math.sqrt ( xI ) - math.sqrt ( xJ ) )**2 / ( xI * rI + xJ * rJ )
    return ( rIJ0 + rBO - rEN ) # . Note minus sign.

#===================================================================================================================================
# . Factory functions.
#===================================================================================================================================
def DYFFCosineAngleParameters ( dyffAtoms ):
    """Factory function for cosine angle parameters."""
    def f ( types, indices, connectivity ):
        """The cosine angle parameters for three elements."""
        # . Bond orders and connections.
        ( i, j, k )   = indices
        ijBondOrder   = DYFFBondOrder ( i, j, connectivity )
        jkBondOrder   = DYFFBondOrder ( j, k, connectivity )
        nJConnections = len ( connectivity.adjacentNodes[connectivity.nodes[j]] )
        # . Atom type data.
        ( iT, jT, kT ) = types
        jAtom   = dyffAtoms.GetItem ( jT )
        jHybrid = DYFFHybridization ( jT )
        theta   = math.radians ( jAtom.GetProperty ( "Bond Angle" ) ) # . Input angle in degrees.
        zI      = dyffAtoms.GetItem ( iT ).GetProperty ( "Effective Charge" )
        zK      = dyffAtoms.GetItem ( kT ).GetProperty ( "Effective Charge" )
        # . Pair data.
        rIJ     = DYFFNaturalBondLength ( dyffAtoms, iT, jT, ijBondOrder )
        rJK     = DYFFNaturalBondLength ( dyffAtoms, jT, kT, jkBondOrder )
        # . Factors.
        results = []
        sinT    = math.sin ( theta )
        cosT    = math.cos ( theta )
        rIK     = math.sqrt ( rIJ**2 + rJK**2 - 2.0 * rIJ * rJK * cosT )
        fC0     = 2.0 * Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zI * zK * ( 3.0 * rIJ * rJK * ( 1.0 - cosT**2 ) - rIK**2 * cosT ) / rIK**5
        if fC0 != 0.0:
            # . Linear central atom.
            if math.fabs ( theta - math.pi ) < _AnglePiTolerance:
                results.append ( ( fC0, 0 ) )
                results.append ( ( fC0, 1 ) )
            # . Trigonal-planar or resonant central atom.
            elif ( ( jHybrid == _Trigonal ) or DYFFIsResonant ( jT ) ) and ( nJConnections == 3 ):
                results.append ( (   fC0 / 9.0, 0 ) )
                results.append ( ( - fC0 / 9.0, 3 ) )
            # . Square-planar or octahedral central atom.
            elif ( ( jHybrid == _SquarePlanar ) and ( nJConnections == 4 ) ) or ( ( jHybrid == _Octahedral ) and ( nJConnections == 6 ) ):
                results.append ( (   fC0 / 16.0, 0 ) )
                results.append ( ( - fC0 / 16.0, 4 ) )
            # . General case.
            else:
                fC0 /= ( 4.0 * sinT**2 )
                results.append ( (   fC0 * ( 1.0 + 2.0 * cosT**2 ), 0 ) )
                results.append ( ( - fC0 * 4.0 * cosT             , 1 ) )
                results.append ( (   fC0                          , 2 ) )
        return results
    return f

def DYFFCosineDihedralParameters ( dyffAtoms ):
    """Factory function for cosine dihedral parameters."""
    def f ( types, indices, connectivity ):
        """The cosine dihedral parameters for four elements."""
        # . Bond orders and connections.
        ( i, j, k, l ) = indices
        jkBondOrder    = DYFFBondOrder ( j, k, connectivity )
        nJConnections  = len ( connectivity.adjacentNodes[connectivity.nodes[j]] )
        nKConnections  = len ( connectivity.adjacentNodes[connectivity.nodes[k]] )
        # . Get the dihedral multiplicity (the factor of 2 is due to the 1/2 outside the UFF force constant definition).
        weight = 2.0 * float ( nJConnections - 1 ) * float ( nKConnections - 1 )
        # . Atom type data.
        ( iT, jT, kT, lT ) = types
        jAtom     = dyffAtoms.GetItem ( jT )
        jHybrid   = DYFFHybridization ( jT )
        jIsPlanar = ( ( jHybrid == _Resonant ) or ( jHybrid == _Trigonal ) ) 
        jTorsion2 = jAtom.GetProperty ( "Torsion 2" )
        jTorsion3 = jAtom.GetProperty ( "Torsion 3" )
        nJ        = jAtom.atomicNumber
        kAtom     = dyffAtoms.GetItem ( kT )
        kHybrid   = DYFFHybridization ( kT )
        kIsPlanar = ( ( kHybrid == _Resonant ) or ( kHybrid == _Trigonal ) ) 
        kTorsion2 = kAtom.GetProperty ( "Torsion 2" )
        kTorsion3 = kAtom.GetProperty ( "Torsion 3" )
        nK        = kAtom.atomicNumber
        # . Get flags to indicate group 16 elements (O, S, Se, Te, Po, Lv).
        isJGroup16 = ( PeriodicTable[nJ].group == 16 )
        isKGroup16 = ( PeriodicTable[nK].group == 16 )
        # . sp3-sp3 case.
        if ( jHybrid == _Tetrahedral ) and ( kHybrid == _Tetrahedral ):
            # . Two group 16 elements.
            if isJGroup16 and isKGroup16:
                if nJ == 8: jTorsion3 = _Dihedral_Osp3
                else:       jTorsion3 = _Dihedral_Group16sp3
                if nK == 8: kTorsion3 = _Dihedral_Osp3
                else:       kTorsion3 = _Dihedral_Group16sp3
                period = 2
                phi0   = 0.5 * math.pi
            # . Other combinations.
            else:
                period = 3
                phi0   = math.pi
            fC = math.sqrt ( jTorsion3 * kTorsion3 )
            #print ( "DIHSP3SP3> {:5d} {:5d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format ( nJConnections, nKConnections, jTorsion3, kTorsion3, fC, weight ) )
        # . sp3-sp2/resonant case.
        elif ( ( jHybrid == _Tetrahedral ) and kIsPlanar ) or ( ( kHybrid == _Tetrahedral ) and jIsPlanar ):
            # . j or k are group 16 sp3.
            if ( isJGroup16 and ( jHybrid == _Tetrahedral ) ) or ( isKGroup16 and ( kHybrid == _Tetrahedral ) ):
                fC     = _Dihedral_sp2_A * math.sqrt ( jTorsion2 * kTorsion2 ) * ( 1.0 + _Dihedral_sp2_B * math.log ( jkBondOrder ) )
                period = 2
                phi0   = math.pi
            # . j is sp3.
            elif jHybrid == _Tetrahedral:
                factor = float ( nKConnections - 1 )
                # . l is sp2/resonant.
                if DYFFHybridization ( lT ) in ( _Resonant, _Trigonal ):
                    fC     = factor * _Dihedral_sp3sp2sp2
                    period = 3
                    phi0   = math.pi
                else:
                    fC     = factor * _Dihedral_sp3sp2
                    period = 6
                    phi0   = 0.0
            # . k is sp3.
            elif kHybrid == _Tetrahedral:
                factor = float ( nJConnections - 1 )
                # . i is sp2/resonant.
                if DYFFHybridization ( iT ) in ( _Resonant, _Trigonal ):
                    fC     = factor * _Dihedral_sp3sp2sp2
                    period = 3
                    phi0   = math.pi
                else:
                    fC     = factor * _Dihedral_sp3sp2
                    period = 6
                    phi0   = 0.0
            #print ( "DIHSP3SP2> {:5d} {:5d} {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format ( nJConnections, nKConnections, jTorsion2, kTorsion2, fC, weight ) )
        # . Resonant atom-resonant atom case.
        elif ( jHybrid == _Resonant ) and ( kHybrid == _Resonant ):
            fC     = _Dihedral_ResonantResonant
            period = 2
            phi0   = math.pi
        # . sp2-sp2 case.
        elif ( jHybrid == _Trigonal ) and ( kHybrid == _Trigonal ):
            fC     = _Dihedral_sp2_A * math.sqrt ( jTorsion2 * kTorsion2 ) * ( 1.0 + _Dihedral_sp2_B * math.log ( jkBondOrder ) )
            period = 2
            phi0   = math.pi
        # . sp2/resonant-sp2/resonant atom case.
        elif ( jIsPlanar and kIsPlanar ):
            fC     = _Dihedral_PlanarPlanar
            period = 2
            phi0   = math.pi
        # . Other cases.
        else:
            fC = 0.0
        if fC == 0.0:
            results = []
        else:
            fC /= weight
            results = [ ( fC, 0 ), ( - fC * math.cos ( float ( period ) * phi0 ), period ) ]
        return results
    return f

def DYFFCosineOutOfPlaneParameters ( dyffAtoms ):
    """Factory function for cosine out-of-plane parameters."""
    def f ( types, indices, connectivity ):
        """The cosine out-of-plane parameters for four elements."""
        # . Atom data.
        ( jT, iT, kT, lT ) = types # . jT is central atom.
        jAtom              = dyffAtoms.GetItem ( jT )
        jHybrid            = DYFFHybridization ( jT )
        nJ                 = jAtom.atomicNumber
        # . sp2 and resonant carbon and nitrogen.
        if ( ( nJ == 6 ) or ( nJ == 7 ) ) and ( ( jHybrid == _Resonant ) or ( jHybrid == _Trigonal ) ):
            # . One of the other atoms is an sp2 oxygen.
            if any ( [ ( ( dyffAtoms.GetItem ( t ).atomicNumber == 8 ) and ( DYFFHybridization ( t ) ) ) for t in ( iT, kT, lT ) ] ):
                results = _OOP_CNsp2Osp2
            else:
                results = _OOP_CNsp2
        # . Group 15 tetrahedral.
        elif ( jHybrid == _Tetrahedral ) and ( PeriodicTable[nJ].group == 15 ):
            if   nJ ==   7: results = _OOP_Nsp3
            elif nJ ==  15: results = _OOP_Psp3
            elif nJ ==  33: results = _OOP_Assp3
            elif nJ ==  51: results = _OOP_Sbsp3
            elif nJ ==  83: results = _OOP_Bisp3
            elif nJ == 115: results = _OOP_Mcsp3
        # . All other elements.
        else: results = []
        return results
    return f

def DYFFHarmonicBondParameters ( dyffAtoms ):
    """Factory function for harmonic bond parameters."""
    def f ( types, indices, connectivity ):
        """The harmonic bond parameters for two elements."""
        # . Bond order.
        bondOrder = DYFFBondOrder ( indices[0], indices[1], connectivity )
        # . Bond length.
        ( iT, jT ) = types
        eq         = DYFFNaturalBondLength ( dyffAtoms, iT, jT, bondOrder )
        # . Force constant.
        zI = dyffAtoms.GetItem ( iT ).GetProperty ( "Effective Charge" )
        zJ = dyffAtoms.GetItem ( jT ).GetProperty ( "Effective Charge" )
        fC = Units.Energy_E2Angstroms_To_Kilojoules_Per_Mole * zI * zJ / eq**3
        if fC == 0.0: return None
        else:         return ( eq, fC )
    return f

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""Classes for dealing with electronic states."""

import math

from   enum               import Enum                  , \
                                 IntEnum
from   pCore              import AttributableObject    , \
                                 logFile               , \
                                 LogFileActive
from   pScientific        import Units
from   pScientific.Arrays import Array
from  .QCDefinitions      import FockClosurePriority
from  .QCModelError       import QCModelError
from ..EnergyModel        import EnergyClosurePriority

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_CardinalOccupancyTolerance =  1.0e-10
_FermiBisectionTolerance    =  1.0e-12
_FermiChargeTolerance       =  1.0e-10
_FermiEnergyLower           = -1.0e+03
_FermiEnergyUpper           =  1.0e+03
_FermiMaximumExponent       =  500.0
_FermiMaximumIterations     =  500
_TotalChargeTolerance       =  1.0e-10

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class MOMMetric ( Enum ):
    """The MOM metric to use."""
    AbsoluteMaximum = 1
    Norm2           = 2
    PMOM            = 3

class OccupancyType ( Enum ):
    """The type of occupancy."""
    Cardinal           = 1 # . Cardinal occupancy.
    FractionalFixed    = 2 # . Certain specified orbitals have non-cardinal occupancy.
    FractionalVariable = 3 # . Occupancies can be fractional and are determined on-the-fly. 
    MOM                = 4 # . Maximum overlap method - cardinal occupancy.

class SpinMultiplicity ( IntEnum ):
    """Spin multiplicities."""
    Singlet     =  1
    Doublet     =  2
    Triplet     =  3
    Quartet     =  4
    Quintet     =  5
    Sextet      =  6
    Septet      =  7
    Octet       =  8
    Nonet       =  9
    Decuplet    = 10
    Undecuplet  = 11
    Duodecuplet = 12
    Tredecuplet = 13

class SpinType ( Enum ):
    """The type of spin."""
    Alpha = 1
    Beta  = 2
    Spin  = 3
    Total = 4

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ElectronicState ( AttributableObject ):
    """A QC electronic state."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "charge"                     : 0                      ,
                             "fermiBroadening"            : 1000.0                 ,
                             "isSpinRestricted"           : True                   ,
                             "multiplicity"               : 1                      ,
                             "numberFractionalAlphaHOOs"  : 0                      ,
                             "numberFractionalAlphaLUOs"  : 0                      ,
                             "numberFractionalBetaHOOs"   : 0                      ,
                             "numberFractionalBetaLUOs"   : 0                      ,
                             "numberFractionalHOOs"       : 0                      ,
                             "numberFractionalLUOs"       : 0                      ,
                             "occupancyType"              : OccupancyType.Cardinal ,
                             "permitRestrictedNonSinglet" : False                  } )

    def _OccupancyEnergyClosure ( self, target ):
        """Set up the MOM occupancy handlers."""
        scratch = target.scratch
        for attribute in ( "orbitalsP", "orbitalsQ" ):
            orbitals = scratch.Get ( attribute, None )
            if orbitals is not None:
                orbitals.occupancyHandler.SetUpProcessOrbitals ( orbitals, scratch )

    def _OccupancyFockClosure ( self, target ):
        """The occupancy energy."""
        energy  = 0.0
        scratch = target.scratch
        for attribute in ( "orbitalsP", "orbitalsQ" ):
            orbitals = scratch.Get ( attribute, None )
            if orbitals is not None: energy += orbitals.occupancyHandler.OccupancyEnergy ( orbitals.occupancies )
        scratch.energyTerms["QC Orbital Occupancy"] = ( energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        return energy

    def EnergyClosures ( self, target ):
        """Fock closures."""
        closures = []
        if self.occupancyType is OccupancyType.MOM:
            def a ( ):
                return self._OccupancyEnergyClosure ( target )
            closures.append ( ( EnergyClosurePriority.QCPreEnergy, a, "QC MOM Set Up" ) )
        return closures

    def FockClosures ( self, target ):
        """Fock closures."""
        closures = []
        if self.occupancyType is OccupancyType.FractionalVariable:
            def b ( ):
                return self._OccupancyFockClosure ( target )
            closures.append ( ( FockClosurePriority.Medium, b ) )
        return closures

    def FractionalOrbitals ( self, spinType ):
        """Get the fractional orbital counters for a given spin type."""
        if   spinType is SpinType.Total: return ( self.numberFractionalHOOs      , self.numberFractionalLUOs      )
        elif spinType is SpinType.Alpha: return ( self.numberFractionalAlphaHOOs , self.numberFractionalAlphaLUOs )
        elif spinType is SpinType.Beta : return ( self.numberFractionalBetaHOOs  , self.numberFractionalBetaLUOs  )
        else: return None

    def OccupancyHandlers ( self, alphaCharge, betaCharge ):
        """Get the occupancy handlers for the state."""
        handlers = None
        if self.occupancyType is not None:
            if self.isSpinRestricted: options = [ { "spinType" : SpinType.Total, "totalCharge" : alphaCharge + betaCharge } ]
            else:                     options = [ { "spinType" : SpinType.Alpha, "totalCharge" : alphaCharge              } ,
                                                  { "spinType" : SpinType.Beta , "totalCharge" :               betaCharge } ]
            if self.occupancyType is OccupancyType.Cardinal:
                handlerClass = OccupancyHandlerCardinal
            elif self.occupancyType is OccupancyType.FractionalFixed:
                handlerClass = OccupancyHandlerFractionalFixed
                for k in options: ( k["numberFractionalHOOs"], k["numberFractionalLUOs"] ) = self.FractionalOrbitals ( k["spinType"] )
            elif self.occupancyType is OccupancyType.FractionalVariable:
                handlerClass = OccupancyHandlerFractionalVariable
                for k in options: k["fermiBroadening"] = self.fermiBroadening
            elif self.occupancyType is OccupancyType.MOM:
                handlerClass = OccupancyHandlerMOM
            handlers = [ handlerClass.WithOptions ( **k ) for k in options ] # . Old way using FromCharge better as no unnecessary attributes?
        return handlers

    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( "Electronic State"   , True                                               ) ,
                  ( "Electronic Charge"  , "{:d}".format ( self.charge                      ) ) ,
                  ( "Spin Multiplicity"  , "{:d}".format ( self.multiplicity                ) ) ,
                  ( "Is Spin Restricted" , "{:s}".format ( repr ( self.isSpinRestricted   ) ) ) ,
                  ( "Occupancy Type"     , "{:s}".format (        self.occupancyType.name   ) ) ]
        if self.occupancyType is OccupancyType.FractionalFixed:
            if self.isSpinRestricted:
                items.extend ( [ ( "Number Fractional HOOs" , "{:d}".format ( self.numberFractionalHOOs ) ) ,
                                 ( "Number Fractional LUOs" , "{:d}".format ( self.numberFractionalLUOs ) ) ] )
            else:
                items.extend ( [ ( "Number Fractional Alpha HOOs" , "{:d}".format ( self.numberFractionalAlphaHOOs ) ) ,
                                 ( "Number Fractional Alpha LUOs" , "{:d}".format ( self.numberFractionalAlphaLUOs ) ) ,
                                 ( "Number Fractional Beta HOOs"  , "{:d}".format ( self.numberFractionalBetaHOOs  ) ) ,
                                 ( "Number Fractional Beta LUOs"  , "{:d}".format ( self.numberFractionalBetaLUOs  ) ) ] )
        if ( self.occupancyType != OccupancyType.Cardinal ) and ( self.multiplicity > 1 ) and self.isSpinRestricted:
            items.append ( ( "Permit Restricted Non-Singlet", "{:s}".format ( repr ( self.permitRestrictedNonSinglet ) ) ) )
        return items

    def Verify ( self, nuclearCharge ):
        """Verify the state versus the electronic charge."""
        numberOfElectrons = round ( nuclearCharge )
        electronicCharge  = float ( numberOfElectrons ) - float ( self.charge )
        alphaCharge       = ( electronicCharge + ( float ( self.multiplicity ) - 1.0 ) ) / 2.0
        betaCharge        = ( electronicCharge - ( float ( self.multiplicity ) - 1.0 ) ) / 2.0
        if self.isSpinRestricted and ( self.multiplicity > 1 ):
            if self.permitRestrictedNonSinglet:
                alphaCharge = electronicCharge / 2.0
                betaCharge  = electronicCharge / 2.0
            else:
                raise QCModelError ( "A spin-unrestricted calculation is required for non-singlet states." )
        return ( alphaCharge, betaCharge )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class OccupancyHandler ( AttributableObject ):
    """Base class for orbital occupancy handlers."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "numberOccupied"  : 0    ,
                             "occupancyFactor" : 1.0  ,
                             "spinType"        : None ,
                             "totalCharge"     : 0.0  } )

    def CheckOccupancies ( self, occupancies ):
        """Check the occupancies."""
        dQ = sum ( occupancies ) - self.totalCharge
        if math.fabs ( dQ ) > _TotalChargeTolerance:
            raise QCModelError ( "Orbital occupancy charge mismatch ({:.5f}).".format ( dQ ) )

    def ProcessOrbitals ( self, energies, occupancies, orbitals, scratch ):
        """Process orbitals after construction."""
        # . By default only the Fermi energy is returned.
        fermiEnergy = 0.0
        if self.numberOccupied > 0:
            homo = self.numberOccupied - 1
            lumo = self.numberOccupied
            if lumo < len ( energies ): fermiEnergy = 0.5 * ( energies[homo] + energies[lumo] )
            else:                       fermiEnergy = energies[homo]
        return fermiEnergy

    def SetOccupancies ( self, occupancies ):
        """Set occupancies."""
        occupancies.Set ( 0.0 )

#-----------------------------------------------------------------------------------------------------------------------------------
class OccupancyHandlerCardinal ( OccupancyHandler ):
    """Cardinal occupancies."""

    def _CheckOptions ( self ):
        """Check options."""
        charge = self.totalCharge
        if ( math.fabs ( round ( charge ) - charge ) > _CardinalOccupancyTolerance ):
            raise QCModelError ( "Cardinal occupancy non-integral charge {:.1f}.".format ( charge ) )
        if self.spinType is SpinType.Total: self.occupancyFactor = 2.0
        self.numberOccupied = int ( math.ceil ( charge / self.occupancyFactor ) )

    def SetOccupancies ( self, occupancies ):
        """Set occupancies."""
        occupancies.Set ( 0.0 )
        for i in range ( self.numberOccupied ): occupancies[i] = self.occupancyFactor
        self.CheckOccupancies ( occupancies )

#-----------------------------------------------------------------------------------------------------------------------------------
class OccupancyHandlerFractionalFixed ( OccupancyHandler ):
    """Fixed fractional occupancies."""

    _attributable = dict ( OccupancyHandler._attributable )
    _attributable.update ( { "fractionalCharge"     : 0.0 ,
                             "numberCardinal"       : 0   ,
                             "numberFractional"     : 0   ,
                             "numberFractionalHOOs" : 0   ,
                             "numberFractionalLUOs" : 0   } )

    def _CheckOptions ( self ):
        """Check options."""
        charge = self.totalCharge
        if self.spinType is SpinType.Total: self.occupancyFactor = 2.0
        numberOccupied        = int ( math.ceil ( charge / self.occupancyFactor ) )
        self.numberCardinal   =      numberOccupied       - self.numberFractionalHOOs      
        self.numberFractional = self.numberFractionalHOOs + self.numberFractionalLUOs
        self.numberOccupied   = self.numberCardinal + self.numberFractional
        self.fractionalCharge = ( ( charge / self.occupancyFactor ) - float ( self.numberCardinal ) ) / float ( max ( self.numberFractional, 1 ) )

    def SetOccupancies ( self, occupancies ):
        """Set occupancies."""
        occupancies.Set ( 0.0 )
        for i in range ( self.numberCardinal   ): occupancies[i                    ] =                         self.occupancyFactor
        for i in range ( self.numberFractional ): occupancies[i+self.numberCardinal] = self.fractionalCharge * self.occupancyFactor
        self.CheckOccupancies ( occupancies )

#-----------------------------------------------------------------------------------------------------------------------------------
class OccupancyHandlerFractionalVariable ( OccupancyHandler ):
    """Variable fractional occupancies."""

    _attributable = dict ( OccupancyHandler._attributable )
    _attributable.update ( { "fermiBroadening" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        charge = self.totalCharge
        if self.spinType is SpinType.Total: self.occupancyFactor = 2.0
        self.numberOccupied = int ( math.ceil ( charge / self.occupancyFactor ) )

    def _NumberOccupied ( self, occupancies ):
        """Determine the number of occupied orbitals."""
        n = 0
        for o in occupancies:
            if o <= _FermiChargeTolerance: break
            n+= 1
        return n

    def _Occupancies ( self, fermiEnergy, energies, occupancies ):
        """Determine the occupancies from the orbital energies."""
        for ( i, e ) in enumerate ( energies ):
            e = - self.fermiBroadening * ( fermiEnergy - energies[i] )
            if   e >   _FermiMaximumExponent: occupancies[i] = 0.0
            elif e < - _FermiMaximumExponent: occupancies[i] = 1.0
            else: occupancies[i] = 1.0 / ( 1.0 + math.exp ( e ) )
        return sum ( occupancies )

    def OccupancyEnergy ( self, occupancies ):
        """Calculate the occupancy energy."""
        e = 0.0
        for p in occupancies:
            p /= self.occupancyFactor
            q  = 1.0 - p
            if p > _FermiChargeTolerance: e += p * math.log ( p )
            if q > _FermiChargeTolerance: e += q * math.log ( q )
        return ( e * self.occupancyFactor / self.fermiBroadening )

    def ProcessOrbitals ( self, energies, occupancies, orbitals, scratch ):
        """Set occupancies given the orbital energies."""
        fermiEnergy    = 0.0
        oTarget        = self.totalCharge / self.occupancyFactor
        occupancies.Set ( 0.0 )
        if oTarget > 0.0:
            efLower = _FermiEnergyLower
            efUpper = _FermiEnergyUpper
            oLower  = self._Occupancies ( efLower, energies, occupancies ) - oTarget
            oUpper  = self._Occupancies ( efUpper, energies, occupancies ) - oTarget
            if ( oUpper < 0.0 ) or ( oLower > 0.0 ):
                raise QCModelError ( "Extreme orbital energies do not bracket target charge." )
            dE = ( efUpper - efLower )
            for i in range ( _FermiMaximumIterations ):
                dE *= 0.5
                fermiEnergy = efLower + dE
                oMiddle     = self._Occupancies ( fermiEnergy, energies, occupancies ) - oTarget
                if oMiddle <= 0.0: efLower = fermiEnergy
                if math.fabs ( oMiddle ) <= _FermiBisectionTolerance: break
            if math.fabs ( oMiddle ) > _FermiBisectionTolerance:
                raise QCModelError ( "Unable to find Fermi energy to required accuracy." )
            self.numberOccupied = self._NumberOccupied ( occupancies )
            if self.occupancyFactor != 1.0: occupancies.Scale ( self.occupancyFactor )
        return fermiEnergy

#-----------------------------------------------------------------------------------------------------------------------------------
# .  MOM method: ATB Gilbert, NA Besley, PMW Gill. JPCA (2008) 112, 13164-13171.
#   PMOM method: HH Corzo, AA Taka, A Pribram-Jones, HP Hratchian. JCC (2022) 43, 382-390.
#   The last paper also details the IMOM method and the absolute maximum and norm2 metrics.
#-----------------------------------------------------------------------------------------------------------------------------------
class OccupancyHandlerMOM ( OccupancyHandlerCardinal ):
    """Maximum overlap cardinal occupancies."""

    _attributable = dict ( OccupancyHandlerCardinal._attributable )
    _attributable.update ( { "metric"  : MOMMetric.PMOM ,
                             "useIMOM" : True           } )

    def ProcessOrbitals ( self, energies, occupancies, orbitals, scratch ):
        """Reorder the orbitals and determine the Fermi energy."""
        # . Create new orbital projections from the overlap between old occupied and new orbitals.
        if self.spinType in ( SpinType.Alpha, SpinType.Total ): momTag = "P"
        else:                                                   momTag = "Q" 
        m       = self.numberOccupied
        overlap = scratch.Get ( "overlapMatrix"                         , None )
        CSC     = scratch.Get ( "momOverlapMatrix{:s}".format ( momTag ), None )
        SCOld   = scratch.Get ( "momOrbitals{:s}".format      ( momTag ), None )
        CSC.MatrixMultiply ( orbitals, SCOld, xTranspose = True )
        if self.metric is MOMMetric.AbsoluteMaximum:
            projections0 = [ ( CSC[i,:].AbsoluteMaximum ( ), i ) for i in range ( CSC.rows ) ]
        elif self.metric is MOMMetric.Norm2:
            projections0 = [ ( CSC[i,:].Norm2 ( ), i ) for i in range ( CSC.rows ) ]
        elif self.metric is MOMMetric.PMOM:
            CSC2 = scratch.Get ( "momPMOMMatrix{:s}".format ( momTag ), None )
            CSC2.MatrixMultiply ( CSC, CSC, yTranspose = True )
            projections0 = [ ( CSC2[i,:].Sum ( ), i ) for i in range ( CSC2.rows ) ]
        projections0.sort ( reverse = True )
        pOccupied    = [ ( i, p ) for ( p, i ) in projections0[0:m] ] ; pOccupied.sort ( )
        pVirtual     = [ ( i, p ) for ( p, i ) in projections0[m: ] ] ; pVirtual.sort  ( )
        projections  = pOccupied + pVirtual
        oldPositions = { i : i for i in range ( len ( projections ) ) }
        # . Order orbitals and their energies by maximum projection - in place!
        for ( new, ( old, _ ) ) in enumerate ( projections ):
            if new != old:
                index = oldPositions[old]
                ( energies[new], energies[index] ) = ( energies[index], energies[new] )
                orbitals[:,new].Swap ( orbitals[:,index] )
                oldPositions[oldPositions[new]] = index
        # . Update the old occupied orbitals.
        if not self.useIMOM:
            if overlap is None: orbitals[:,0:m].CopyTo ( SCOld )
            else: overlap.PostMatrixMultiply ( orbitals[:,0:m], SCOld )
        # . Finish up.
        return super ( OccupancyHandlerMOM, self ).ProcessOrbitals ( energies, occupancies, orbitals, scratch )

    def SetUpProcessOrbitals ( self, orbitals, scratch ):
        """Set up the reordering calculation."""
        if self.spinType in ( SpinType.Alpha, SpinType.Total ): momTag = "P"
        else:                                                   momTag = "Q"
        orbitals = orbitals.orbitals
        m        = self.numberOccupied
        n        = orbitals.columns
        overlap  = scratch.Get ( "overlapMatrix", None ) # . The overlap matrix must have been calculated!
        CSC      = scratch.Get ( "momOverlapMatrix{:s}".format ( momTag ), None )
        if ( CSC is None ) or ( CSC.rows != n ) or ( CSC.columns != m ):
            CSC = Array.WithExtents ( n, m )
            scratch.Set ( "momOverlapMatrix{:s}".format ( momTag ), CSC )
        SCOld = scratch.Get ( "momOrbitals{:s}".format ( momTag ), None )
        if ( SCOld is None ) or ( SCOld.rows != n ) or ( SCOld.columns != m ):
            SCOld = Array.WithExtents ( n, m )
            scratch.Set ( "momOrbitals{:s}".format ( momTag ), SCOld )
        if self.metric is MOMMetric.PMOM:
            CSC2 = scratch.Get ( "momPMOMMatrix{:s}".format ( momTag ), None )
            if ( CSC2 is None ) or ( CSC2.rows != n ) or ( CSC2.columns != n ):
                CSC2 = Array.WithExtents ( n, n )
                scratch.Set ( "momPMOMMatrix{:s}".format ( momTag ), CSC2 )
        if overlap is None: orbitals[:,0:m].CopyTo ( SCOld )
        else:               overlap.PostMatrixMultiply ( orbitals[:,0:m], SCOld )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

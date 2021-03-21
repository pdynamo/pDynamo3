"""Classes for dealing with electronic states."""

import math

from  enum          import Enum                , \
                           IntEnum
from  pCore         import AttributableObject  , \
                           logFile             , \
                           LogFileActive
from  pScientific   import Units
from .QCDefinitions import FockClosurePriority
from .QCModelError  import QCModelError

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
class OccupancyType ( Enum ):
    """The type of occupancy."""
    Cardinal           = 1
    FractionalFixed    = 2
    FractionalVariable = 3

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
        """The occupancy energy."""
        energy  = 0.0
        scratch = target.scratch
        for attribute in ( "orbitalsP", "orbitalsQ" ):
            orbitals = scratch.Get ( attribute, None )
            if orbitals is not None: energy += orbitals.occupancyHandler.OccupancyEnergy ( orbitals.occupancies )
        scratch.energyTerms["QC Orbital Occupancy"] = ( energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        return energy

    def FockClosures ( self, target ):
        """Fock closures."""
        def a ( ):
            return self._OccupancyEnergyClosure ( target )
        if self.occupancyType == OccupancyType.FractionalVariable: return [ ( FockClosurePriority.Medium, a ) ]
        else: return []

    def FractionalOrbitals ( self, spinType ):
        """Get the fractional orbital counters for a given spin type."""
        if   spinType == SpinType.Total: return ( self.numberFractionalHOOs      , self.numberFractionalLUOs      )
        elif spinType == SpinType.Alpha: return ( self.numberFractionalAlphaHOOs , self.numberFractionalAlphaLUOs )
        elif spinType == SpinType.Beta : return ( self.numberFractionalBetaHOOs  , self.numberFractionalBetaLUOs  )
        else: return None

    def OccupancyHandlers ( self, alphaCharge, betaCharge ):
        """Get the occupancy handlers for the state."""
        handlers = None
        if self.occupancyType is not None:
            if self.isSpinRestricted: arguments = [ [ ( SpinType.Total, alphaCharge + betaCharge ), {} ] ]
            else:                     arguments = [ [ ( SpinType.Alpha, alphaCharge ), {} ], [ ( SpinType.Beta, betaCharge ), {} ] ]
            if self.occupancyType == OccupancyType.Cardinal:
                handlerClass = OccupancyHandlerCardinal
            elif self.occupancyType == OccupancyType.FractionalFixed:
                handlerClass = OccupancyHandlerFractionalFixed
                for ( a, k ) in arguments: ( k["numberFractionalHOOs"], k["numberFractionalLUOs"] ) = self.FractionalOrbitals ( a[0] )
            elif self.occupancyType == OccupancyType.FractionalVariable:
                handlerClass = OccupancyHandlerFractionalVariable
                for ( a, k ) in arguments: k["fermiBroadening"] = self.fermiBroadening
            handlers = [ handlerClass.FromCharge ( *a, **k ) for ( a, k ) in arguments ]
        return handlers

    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( "Electronic State"   , True                                               ) ,
                  ( "Electronic Charge"  , "{:d}".format ( self.charge                      ) ) ,
                  ( "Spin Multiplicity"  , "{:d}".format ( self.multiplicity                ) ) ,
                  ( "Is Spin Restricted" , "{:s}".format ( repr ( self.isSpinRestricted   ) ) ) ,
                  ( "Occupancy Type"     , "{:s}".format (        self.occupancyType.name   ) ) ]
        if self.occupancyType == OccupancyType.FractionalFixed:
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
    _attributable.update ( { "numberOccupied"  : 0   ,
                             "occupancyFactor" : 1.0 ,
                             "totalCharge"     : 0.0 } )

    def CheckOccupancies ( self, occupancies ):
        """Check the occupancies."""
        dQ = sum ( occupancies ) - self.totalCharge
        if math.fabs ( dQ ) > _TotalChargeTolerance:
            raise QCModelError ( "Orbital occupancy charge mismatch ({:.5f}).".format ( dQ ) )

    def SetOccupancies ( self, occupancies ):
        """Set occupancies."""
        occupancies.Set ( 0.0 )

    def SetOccupanciesFromEnergies ( self, energies, occupancies ):
        """Set occupancies given the orbital energies."""
        # . By default only the Fermi energy is returned.
        fermiEnergy = 0.0
        if self.numberOccupied > 0:
            homo = self.numberOccupied - 1
            lumo = self.numberOccupied
            if lumo < len ( energies ): fermiEnergy = 0.5 * ( energies[homo] + energies[lumo] )
            else:                       fermiEnergy = energies[homo]
        return fermiEnergy

class OccupancyHandlerCardinal ( OccupancyHandler ):
    """Cardinal occupancies."""

    @classmethod
    def FromCharge ( selfClass, spinType, charge ):
        """Constructor given a charge and spin type."""
        self = selfClass ( )
        if ( math.fabs ( round ( charge ) - charge ) > _CardinalOccupancyTolerance ):
            raise QCModelError ( "Cardinal occupancy non-integral charge {:.1f}.".format ( charge ) )
        if spinType == SpinType.Total: self.occupancyFactor = 2.0
        self.numberOccupied = int ( math.ceil ( charge / self.occupancyFactor ) )
        self.totalCharge    = charge
        return self

    def SetOccupancies ( self, occupancies ):
        """Set occupancies."""
        occupancies.Set ( 0.0 )
        for i in range ( self.numberOccupied ): occupancies[i] = self.occupancyFactor
        self.CheckOccupancies ( occupancies )

class OccupancyHandlerFractionalFixed ( OccupancyHandler ):
    """Fixed fractional occupancies."""

    _attributable = dict ( OccupancyHandler._attributable )
    _attributable.update ( { "fractionalCharge" : 0.0 ,
                             "numberCardinal"   : 0   ,
                             "numberFractional" : 0   } )

    @classmethod
    def FromCharge ( selfClass, spinType, charge, numberFractionalHOOs = 0, numberFractionalLUOs = 0 ):
        """Constructor given a charge and spin type."""
        self = selfClass ( )
        if spinType == SpinType.Total: self.occupancyFactor = 2.0
        numberOccupied        = int ( math.ceil ( charge / self.occupancyFactor ) )
        self.numberCardinal   = numberOccupied       - numberFractionalHOOs      
        self.numberFractional = numberFractionalHOOs + numberFractionalLUOs
        self.numberOccupied   = self.numberCardinal + self.numberFractional
        self.fractionalCharge = ( ( charge / self.occupancyFactor ) - float ( self.numberCardinal ) ) / float ( max ( self.numberFractional, 1 ) )
        self.totalCharge      = charge
        return self

    def SetOccupancies ( self, occupancies ):
        """Set occupancies."""
        occupancies.Set ( 0.0 )
        for i in range ( self.numberCardinal   ): occupancies[i                    ] =                         self.occupancyFactor
        for i in range ( self.numberFractional ): occupancies[i+self.numberCardinal] = self.fractionalCharge * self.occupancyFactor
        self.CheckOccupancies ( occupancies )

class OccupancyHandlerFractionalVariable ( OccupancyHandler ):
    """Variable fractional occupancies."""

    _attributable = dict ( OccupancyHandler._attributable )
    _attributable.update ( { "fermiBroadening" : None } )

    @classmethod
    def FromCharge ( selfClass, spinType, charge, fermiBroadening = 1000.0 ):
        """Constructor given a charge and spin type."""
        self = selfClass ( )
        if spinType == SpinType.Total: self.occupancyFactor = 2.0
        self.fermiBroadening = fermiBroadening
        self.numberOccupied  = int ( math.ceil ( charge / self.occupancyFactor ) )
        self.totalCharge     = charge
        return self

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

    def SetOccupanciesFromEnergies ( self, energies, occupancies ):
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

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

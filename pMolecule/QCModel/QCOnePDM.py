"""A class for handling one-particle density matrices."""

from  pScientific.Arrays import Array       , \
                                StorageType
from .ElectronicState    import SpinType
from .QCModelError       import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCOnePDM:
    """The one-particle density matrix."""

    _attributable  = { "density"        : None                          ,
                       "fock"           : False                         ,
                       "isValid"        : False                         ,
                       "label"          : "One-Particle Density Matrix" ,
                       "numberOrbitals" : 0                             ,
                       "spinType"       : SpinType.Total                ,
                       "totalCharge"    : None                          }

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__._attributable.items ( ): setattr ( self, key, value )

    # . These next two methods are primarily for internal use.
    def AlphaBetaFromTotalSpin ( self, other ):
        """Convert total and spin densities to alpha and beta densities."""
        if other is not None:
            qA = 0.5 * ( self.totalCharge + other.totalCharge )
            qB = 0.5 * ( self.totalCharge - other.totalCharge )
            self.density.Add  ( other.density ) ; self.density.Scale ( 0.5 )
            self.spinType     = SpinType.Alpha
            self.totalCharge  = qA
            other.density.Add ( self.density, scale = -1.0 ) ; other.density.Scale ( -1.0 )
            other.spinType    = SpinType.Beta 
            other.totalCharge = qB

    def AlphaBetaToTotalSpin ( self, other ):
        """Convert alpha and beta densities to total and spin densities."""
        if other is not None:
            qT = self.totalCharge + other.totalCharge
            qS = self.totalCharge - other.totalCharge
            self.density.Add  ( other.density )
            self.spinType     = SpinType.Total
            self.totalCharge  = qT
            other.density.Scale ( -2.0 ) ; other.density.Add ( self.density )
            other.spinType    = SpinType.Spin 
            other.totalCharge = qS

    @classmethod
    def FromDiagonalGuess ( selfClass      ,
                            numberOrbitals ,
                            spinType       ,
                            totalCharge    ):
        """Constructor from a diagonal guess."""
        n = float ( numberOrbitals )
        if spinType == SpinType.Total: maximumCharge = 2.0 * n
        else:                          maximumCharge =       n
        if totalCharge > maximumCharge: raise QCModelError ( "Density matrix charge overflow {:.1f}/{:.1f}.".format ( totalCharge, maximumCharge ) )
        self             = selfClass.WithExtent ( numberOrbitals )
        self.spinType    = spinType
        self.totalCharge = totalCharge
        self.density.Set ( 0.0 )
        population = ( totalCharge / n )
        for i in range ( numberOrbitals ): self.density[i,i] = population
        return self

    # . This requires total and spin densities in their orthogonal representation.
    # . Otherwise the trace is of PSPS - QSQS.
    def SpinExpectationValues ( self, other, overlap = None ):
        """Spin expectation values."""
        if ( self.spinType == SpinType.Total ) and ( other.spinType == SpinType.Spin ):
            qTotal = self.totalCharge  ; qSpin = other.totalCharge ; sign =  1.0
        elif ( self.spinType == SpinType.Spin ) and ( other.spinType == SpinType.Total ):
            qTotal = other.totalCharge ; qSpin = self.totalCharge  ; sign = -1.0
        else:
            raise QCModelError ( "Spin expectation values require total and spin densities." )
        if overlap is None:
            trace = self.density.TraceOfProduct ( self.density ) - other.density.TraceOfProduct ( other.density )
        else:
            SDS = Array.WithShape ( overlap.shape, storageType = StorageType.Symmetric )
            self.density.SymmetricTransform  ( overlap, SDS ) ; trace  = self.density.TraceOfProduct  ( SDS )
            other.density.SymmetricTransform ( overlap, SDS ) ; trace -= other.density.TraceOfProduct ( SDS )
        Sz =           0.5 * qSpin
        S2 = Sz * Sz + 0.5 * qTotal - 0.25 * sign * trace
        return ( S2, Sz )

    @classmethod
    def WithExtent ( selfClass, numberOrbitals ):
        """Constructor with a given number of orbitals."""
        self = selfClass ( )
        if numberOrbitals > 0:
            self.numberOrbitals = numberOrbitals
            self.density        = Array.WithExtent ( numberOrbitals, storageType = StorageType.Symmetric )
            self.fock           = Array.WithExtent ( numberOrbitals, storageType = StorageType.Symmetric )
        return self

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

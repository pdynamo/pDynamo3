"""A class for handling orbital sets."""

from  pScientific.Arrays        import Array       , \
                                       StorageType
from  pScientific.LinearAlgebra import EigenPairs
from .ElectronicState           import SpinType
from .QCModelError              import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCOrbitals:
    """A set of QC orbitals."""

    _attributable = { "energies"         : None       ,
                      "fermiEnergy"      : 0.0        ,
                      "label"            : "Orbitals" ,
                      "numberOrbitals"   : 0          ,
                      "occupancies"      : None       ,
                      "occupancyHandler" : None       ,
                      "orbitals"         : None       }

    def __init__ ( self ):
        """Constructor."""
        for ( key, value ) in self.__class__._attributable.items ( ): setattr ( self, key, value )

    def MakeFromFock ( self, fock, orthogonalizer = None, preserveInput = True ):
        """Make the orbitals from a Fock matrix."""
        eigenValues = self.energies
        if orthogonalizer is None:
            eigenVectors = self.orbitals
            fockLocal    = fock
        else:
            n            = orthogonalizer.columns
            eigenVectors = Array.WithExtents ( n, n )
            fockLocal    = Array.WithExtent  ( n, storageType = StorageType.Symmetric )
            fock.Transform ( orthogonalizer, fockLocal, useTranspose = False )
        EigenPairs ( fockLocal, eigenValues, eigenVectors, preserveInput = preserveInput )
        if orthogonalizer is not None:
            self.orbitals[:,0:eigenVectors.shape[1]].MatrixMultiply ( orthogonalizer, eigenVectors )
        self.fermiEnergy = self.occupancyHandler.SetOccupanciesFromEnergies ( self.energies, self.occupancies )

    def MakeWeightedDensity ( self, wDM ):
        """Make the weighted density matrix."""
        # . wDM must be initialized on entry.
        n       = self.occupancyHandler.numberOccupied
        scratch = Array.WithExtent ( n )
        for i in range ( n ): scratch[i] = -2.0 * self.energies[i] * self.occupancies[i]
        wDM.MakeFromEigenSystem ( n, scratch, self.orbitals, initialize = False )

    @classmethod
    def WithExtents ( selfClass, numberBasisFunctions, numberOrbitals, occupancyHandler ):
        """Constructor given the number of basis functions and orbitals, and an occupancy handler."""
        self = selfClass ( )
        if numberOrbitals > 0:
            self.numberOrbitals   = numberOrbitals
            self.occupancyHandler = occupancyHandler
            self.energies         = Array.WithExtent  ( numberOrbitals )
            self.occupancies      = Array.WithExtent  ( numberOrbitals )
            self.orbitals         = Array.WithExtents ( numberBasisFunctions, numberOrbitals )
            self.occupancyHandler.SetOccupancies ( self.occupancies )
        return self

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""QC definitions."""

from enum  import Enum    , \
                  IntEnum

# . These are defined here so as to avoid cyclic references (ElectronicState <-> QCModelBase, etc.).

# . Some definitions are general to all QC models, some are relevant to built-in QC models only.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class BasisRepresentation ( Enum ):
    """The basis represenation."""
    Actual     = 10 # . The orthogonalized atom-centered basis.
    Orthogonal = 20 # . The molecular basis resulting from orthogonalization (during the SCF for example).
    Work       = 30 # . The primitive atom-centered work basis.

class ChargeModel ( Enum ):
    """Charge types."""
    CHelpG   =  0
    Loewdin  = 10
    MNDO     = 20
    Mulliken = 30

class FockClosurePriority ( IntEnum ):
    """Fock closure priorities."""
    VeryHigh =  0 # . Includes initialization.
    High     = 10
    Medium   = 20
    Low      = 30
    VeryLow  = 40

class OrthogonalizationType ( Enum ):
    """The type of orthogonalization."""
    Canonical = 1
    Symmetric = 2

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

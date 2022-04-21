"""QC definitions."""

from enum import Enum, IntEnum

# . These are defined here so as to avoid cyclic references (ElectronicState <-> QCModelBase, etc.).

# . Some definitions are general to all QC models, some are relevant to built-in QC models only.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
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

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

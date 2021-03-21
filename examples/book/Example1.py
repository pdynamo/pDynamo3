"""Example 1."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Create a water molecule.
water       = System.FromAtoms ( [ 8, 1, 1 ] )
water.label = "Water"
water.Summary ( )

# . Create a water dimer.
waterDimer       = MergeByAtom ( water, water )
waterDimer.label = "Water Dimer"
waterDimer.Summary ( )

# . Create a hydroxyl.
oh = Selection.FromIterable ( [ 0, 1 ] )
hydroxyl       = PruneByAtom ( water, oh )
hydroxyl.label = "Hydroxyl"
hydroxyl.Summary ( )

# . Footer.
logFile.Footer ( )

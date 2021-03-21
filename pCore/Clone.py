"""Clone functions."""

# . Clone function.
import copy

# . Generic.
Clone = copy.deepcopy

# . Specific.
ShallowClone = copy.copy
DeepClone    = copy.deepcopy

del copy

"""A package for numerical integration."""

from .LebedevLaikovGrid    import LebedevLaikovGrid_GetGridPoints
from .NumericalIntegration import AdaptiveSimpsonsRule , \
                                  TrapezoidalRule

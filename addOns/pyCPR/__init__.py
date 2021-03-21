"""Conjugate-peak refinement."""

#===================================================================================================================================
# . Modules contributed by Florian Gisdon, Martin Culka and Matthias Ullmann.
#
#   Source, documentation and tests from:
#
#   http://www.bisb.uni-bayreuth.de/index.php?page=data/PyCPR/PyCPR
#
#   Reference:
#
#   Florian J. Gisdon, Martin Culka, G. Matthias Ullmann (2016)
#   PyCPR - A Python-based Implementation of the Conjugate Peak Refinement (CPR) Algorithm for Finding Transition State Structures.
#   J Mol Model 22, 242
#   DOI: 10.1007/s00894-016-3116-8
#
#===================================================================================================================================

from .ConjugatePeakRefinement      import CPRObjectiveFunction                , \
                                          CPRRefinement                       , \
                                          ConjugatePeakRefinementOptimizePath
from .CPRSaddlePointRefinement     import CPRSaddlePointRefinement            , \
                                          CPRSaddlePointRefinementState       , \
                                          LineSearchObjectiveFunction
from .MoreThuenteLineSearchWithMax import MoreThuenteLineSearcherWithMax      , \
                                          MoreThuenteLineSearcherWithMaxState

"""A sub-package for dense and sparse linear algebra."""

from .CGLinearEquationSolver        import CGLinearEquationSolver             , \
                                           CGLinearEquationSolverState
from .DenseLinearAlgebra            import Determinant                        , \
                                           EigenPairs                         , \
                                           EigenValues                        , \
                                           LinearEquations                    , \
                                           LinearLeastSquaresBySVD            , \
                                           MatrixPower                        , \
                                           MatrixPowerInverse                 , \
                                           MatrixPseudoInverse
from .LinearAlgebraError            import LinearAlgebraError
from .MachineConstants              import MachineConstants
from .OrthogonalizingTransformation import OrthogonalizationMethod            , \
                                           OrthogonalizingTransformation_Make

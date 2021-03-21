"""A conjugate-gradient linear equation solver."""

# . Basic no-frills implementation.

from pCore              import AttributableObject , \
                               SummarizableObject
from pScientific.Arrays import Array

#==================================================================================================================================
# . Solver state.
#=================================================================================================================================*/
class CGLinearEquationSolverState ( AttributableObject ):
    """State for a conjugate-gradient linear equation solver."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "b"                 : None  ,
                             "doPreconditioning" : False ,
                             "h"                 : None  ,
                             "r"                 : None  ,
                             "x"                 : None  ,
                             "target"            : None  } )

    @classmethod
    def FromTarget ( selfClass, target, rhs, solution, doPreconditioning = False ):
        """Constructor from target."""
        self                   = selfClass ( )
        self.doPreconditioning = doPreconditioning
        self.target            = target
        # . Aliases.
        self.b                 = rhs
        self.x                 = solution
        # . Allocate space.
        n                      = len ( self.b )
        self.h                 = Array.WithExtent ( n ) 
        self.r                 = Array.WithExtent ( n ) 
        return self

#==================================================================================================================================
# . Solver.
#=================================================================================================================================*/
class CGLinearEquationSolver ( SummarizableObject ):
    """A conjugate-gradient linear equation solver."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Conjugate Gradient Linear Equation Solver"
    _summarizable = dict ( SummarizableObject._summarizable )
    _attributable.update ( { "convergenceMode"   :    1                 , 
                             "errorTolerance"    :    1.0e-10           ,
                             "maximumIterations" :  500                 } )
    _summarizable.update ( { "convergenceMode"   : "Convergence Mode"   ,
                             "errorTolerance"    : "Error Tolerance"    ,
                             "maximumIterations" : "Maximum Iterations" } )

    def IsConverged ( self, rNorm2, hNorm2, bNorm2 ):
        """Check for convergence."""
        if   self.convergenceMode == 1: stepTest = rNorm2
        elif self.convergenceMode == 2: stepTest = rNorm2 / bNorm2
        elif self.convergenceMode == 3: stepTest = hNorm2
        elif self.convergenceMode == 4: stepTest = hNorm2 / bNorm2
        return ( stepTest <= self.errorTolerance )

    def Solve ( self, state ):
        """Solver for SPD matrices with or without preconditioning given a state."""
        # . H can be R if there is no preconditioning.
        # . Initialization.
        doPreconditioning = state.doPreconditioning
        target            = state.target
        b                 = state.b
        h                 = state.h
        r                 = state.r
        x                 = state.x
        # . RHS norm2 and the denominator for convergence checks.
        bNorm2            = b.Norm2 ( )
        denominator       = 1.0
        if   self.convergenceMode == 2: denominator = bNorm2
        elif self.convergenceMode == 4:
            if doPreconditioning:
                target.ApplyPreconditioner ( b, r )
                denominator = r.Norm2 ( )
            else:
                denominator = bNorm2
        # . Compute initial residual r = b - A*x.
        target.ApplyMatrix ( x, r )
        r.Add   ( b, scale = -1.0 )
        r.Scale ( -1.0    )
        r0Norm2 = r.Norm2 ( )
        rNorm2  = r0Norm2
        # . Preconditioning h = p * r (in b).
        if doPreconditioning:
            target.ApplyPreconditioner ( r, b )
            h0Norm2 = b.Norm2 ( )
        else:
            r.CopyTo ( b )
            h0Norm2 = r0Norm2
        # . Initial convergence check.
        isConverged = self.IsConverged ( r0Norm2, h0Norm2, denominator )
        iterations  = 0
        if not isConverged:
            rDotH = r.Dot ( b )
            for iterations in range ( 1, self.maximumIterations + 1 ):
                # . New x.
                target.ApplyMatrix ( b, h )
                alpha = rDotH / h.Dot ( b )
                x.Add ( b, scale = alpha )
                # . New r.
                r.Add ( h, scale = -alpha )
                rNorm2 = r.Norm2 ( )
                # . New h.
                if doPreconditioning:
                    target.ApplyPreconditioner ( r, h )
                    hNorm2 = h.Norm2 ( )
                else:
                    r.CopyTo ( h )
                    hNorm2 = rNorm2
                # . Check for termination.
                isConverged = self.IsConverged ( rNorm2, hNorm2, denominator )
                if isConverged: break
                # . New p.
                oldRH = rDotH
                rDotH = r.Dot ( h )
                beta  = rDotH / oldRH
                b.Scale ( beta )
                b.Add ( h )
        # . Finish up.
        report = { "Final Residual"   : rNorm2      ,
                   "Initial Residual" : r0Norm2     ,
                   "Is Converged"     : isConverged ,
                   "Iterations"       : iterations  ,
                   "RHS Norm2"        : bNorm2      }
        return report

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

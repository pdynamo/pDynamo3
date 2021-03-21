"""Helper functions for geometry optimizations."""

from pCore                                  import logFile
from pMolecule                              import SystemGeometryObjectiveFunction
from pScientific.ObjectiveFunctionIterators import BakerOptimizer             , \
                                                   ConjugateGradientMinimizer , \
                                                   FIREMinimizer              , \
                                                   LBFGSMinimizer             , \
                                                   QuasiNewtonMinimizer       , \
                                                   SteepestDescentMinimizer

#===================================================================================================================================
# . Helper function for setting up an optimization.
#===================================================================================================================================
def _SetUpOptimization ( system, defaultOptions, defaultsToChange, inputOptions, removeRotationTranslation ):
    """Generic function to set up an optimization."""
    # . Get default options - overridden with more suitable values if necessary.
    options = dict ( defaultOptions )
    options.update ( { "logFrequency" : 1, "maximumIterations" : 50, "rmsGradientTolerance" : 1.5 } )
    options.update ( defaultsToChange )
    # . Update with the input options.
    options.update ( inputOptions )
    # . Get some non-optimizer options.
    log                        = options.pop ( "log"                       , logFile )
    optimizeSymmetryParameters = options.pop ( "optimizeSymmetryParameters", False   )
    trajectories               = options.pop ( "trajectories"              , []      )
    # . Make sure label is not None.
    if options.get ( "label", None ) is None: options.pop ( "label", None )
    # . Set up the objective function.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    if removeRotationTranslation : of.RemoveRotationTranslation ( )
    if optimizeSymmetryParameters: of.IncludeSymmetryParameters ( )
    for ( trajectory, saveFrequency ) in trajectories:
        of.DefineTrajectory ( trajectory, saveFrequency )
    # . Finish up.
    return ( of, options, log )

#===================================================================================================================================
# . Baker saddle optimization.
#===================================================================================================================================
# . Keyword arguments are those from BakerOptimizer._attributable along with "findMinimum", "log", "optimizeSymmetryParameters" and "trajectories".
def BakerSaddleOptimize_SystemGeometry ( system, **options ):
    """Baker stationary point geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, BakerOptimizer._attributable, {}, options, True )
    # . Minimization.
    if options.pop ( "findMinimum", False ):
        options["followMode"           ] = -1
        options["numberOfNegativeModes"] =  0
        if "hessianUpdatingOption" not in options: options["hessianUpdatingOption"] = "BFGS"
    # . Saddle point.
    else:
        options["followMode"           ] = max ( options["followMode"], 0 )
        options["numberOfNegativeModes"] = 1
        if "hessianUpdatingOption" not in options: options["hessianUpdatingOption"] = "BOFILL"
    optimizer = BakerOptimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Conjugate gradient minimization.
#===================================================================================================================================
# . Keyword arguments are those from ConjugateGradientMinimizer._attributable along with "log", "optimizeSymmetryParameters" and "trajectories".
# . It seems to be necessary to limit the initialStep here so as to prevent problems with SCF convergence when there are QC atoms.
def ConjugateGradientMinimize_SystemGeometry ( system, **options ):
    """Conjugate gradient geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, ConjugateGradientMinimizer._attributable, { "initialStep" : 0.1 }, options, True )
    optimizer = ConjugateGradientMinimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . FIRE minimization.
#===================================================================================================================================
# . Keyword arguments are those from FIREMinimizer._attributable along with "log", "optimizeSymmetryParameters" and "trajectories".
def FIREMinimize_SystemGeometry ( system, **options ):
    """FIRE geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, FIREMinimizer._attributable, { "maximumTimeStep" : 0.01, "timeStep" : 0.001 }, options, True )
    optimizer = FIREMinimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . LBFGS minimization.
#===================================================================================================================================
# . Keyword arguments are those from LBFGSMinimizer._attributable along with "log", "optimizeSymmetryParameters" and "trajectories".
def LBFGSMinimize_SystemGeometry ( system, **options ):
    """LBFGS geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, LBFGSMinimizer._attributable, {}, options, True )
    optimizer = LBFGSMinimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Quasi-Newton minimization.
#===================================================================================================================================
# . Keyword arguments are those from QuasiNewtonMinimizer._attributable along with "hessian", "log", "optimizeSymmetryParameters" and "trajectories".
def QuasiNewtonMinimize_SystemGeometry ( system, **options ):
    """Quasi-Newton geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, QuasiNewtonMinimizer._attributable, {}, options, True )
    hessian   = options.pop ( "hessian", None )
    optimizer = QuasiNewtonMinimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    if hessian is None:
        variables = of.VariablesAllocate ( )
        of.VariablesGet ( variables )
        hessian   = of.StartingHessian ( variables )
    of.SetStartingHessian ( hessian )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Steepest descent minimization.
#===================================================================================================================================
# . Keyword arguments are those from SteepestDescentMinimizer._attributable along with "log", "optimizeSymmetryParameters" and "trajectories".
def SteepestDescentMinimize_SystemGeometry ( system, **options ):
    """Steepest descent geometry optimization."""
    ( of, options, log ) = _SetUpOptimization ( system, SteepestDescentMinimizer._attributable, {}, options, False )
    optimizer = SteepestDescentMinimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

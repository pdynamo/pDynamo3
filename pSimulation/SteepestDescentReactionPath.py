"""Steepest-descent method for finding reaction paths."""

from pCore                                  import logFile
from pMolecule                              import SystemGeometryObjectiveFunction
from pScientific.ObjectiveFunctionIterators import SteepestDescentPathFinder

#===================================================================================================================================
# . Steepest descent path.
#===================================================================================================================================
# . Currently, structures are saved on the trajectory in the order they are produced without regard to which saddle point
# . branch they are in. This needs to be changed. At the moment, the best solution would be to store values on the trajectory
# . according to their identifier (in this case, pathPosition) if this were possible. Otherwise, the trajectory should be
# . reversed after the first branch has been computed. However, not all input trajectories may be reversible in this fashion.
def SteepestDescentPath_SystemGeometry ( system                      ,
                                         fromSaddle        = True    ,
                                         functionStep      =  0.2    ,
                                         logFrequency      =  1      ,
                                         maximumIterations = 50      ,
                                         pathStep          =  0.025  ,
                                         useMassWeighting  = False   ,
                                         saveFrequency     =  0      ,
                                         trajectory        = None    ,
                                         log               = logFile ):
    """Determine a steepest descent path."""
    # . Create an object function.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    if useMassWeighting: of.DefineWeights ( )
    of.RemoveRotationTranslation ( )
    # . Define trajectories.
    if ( saveFrequency > 0 ) and ( saveFrequency < maximumIterations ) and ( trajectory is not None ): of.DefineTrajectory ( trajectory, saveFrequency )
    # . Set up the iterator.
    optimizer = SteepestDescentPathFinder.WithOptions ( fromSaddle        = fromSaddle        ,
                                                        functionStep      = functionStep      ,
                                                        logFrequency      = logFrequency      ,
                                                        maximumIterations = maximumIterations ,
                                                        pathStep          = pathStep          )
    optimizer.Summary ( log = log )
    # . Iterate.
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

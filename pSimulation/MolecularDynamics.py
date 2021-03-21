"""Helper functions for performing molecular dynamics simulations."""

from pCore                                  import logFile
from pMolecule                              import SystemGeometryObjectiveFunction
from pScientific.ObjectiveFunctionIterators import LangevinVelocityVerletIntegrator , \
                                                   LeapFrogIntegrator               , \
                                                   VelocityVerletIntegrator
from pScientific.Symmetry                   import CrystalSystemCubic

#===================================================================================================================================
# . Helper function for setting up a molecular dynamics simulation.
#===================================================================================================================================
def _SetUpSimulation ( system, defaultOptions, defaultsToChange, inputOptions ):
    """Generic function to set up an optimization."""
    # . Get default options - overridden with more suitable values if necessary.
    options = dict ( defaultOptions )
    options.update ( { "logFrequency" : 1, "steps" : 1000, "timeStep" : 0.001 } )
    options.update ( defaultsToChange )
    # . Update with the input options.
    options.update ( inputOptions )
    # . Change steps to maximumIterations.
    options["maximumIterations"] = options.pop ( "steps" )
    # . Get some non-optimizer options.
    log                       = options.pop ( "log"                       , logFile )
    removeRotationTranslation = options.pop ( "removeRotationTranslation" , True    )
    trajectories              = options.pop ( "trajectories"              , []      )
    # . Make sure label is not None.
    if options.get ( "label", None ) is None: options.pop ( "label", None )
    # . Set up the objective function.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    of.DefineWeights ( )
    if removeRotationTranslation: of.RemoveRotationTranslation ( )
    for ( trajectory, saveFrequency ) in trajectories:
        of.DefineTrajectory ( trajectory, saveFrequency )
    # . Finish up.
    return ( of, options, log )

#===================================================================================================================================
# . Langevin molecular dynamics.
#===================================================================================================================================
# . Keyword arguments are those from LangevinVelocityVerletIntegrator._attributable along with "log", "removeRotationTranslation" and "trajectories".
def LangevinDynamics_SystemGeometry ( system, **options ):
    """Molecular dynamics using the velocity Verlet algorithm."""
    ( of, options, log ) = _SetUpSimulation ( system, LangevinVelocityVerletIntegrator._attributable, {}, options )
    of.VelocitiesAssign ( options["temperature"], normalDeviateGenerator = options.get ( "normalDeviateGenerator", None ) )
    optimizer = LangevinVelocityVerletIntegrator.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Molecular dynamics using the leapfrog algorithm.
#===================================================================================================================================
# . Keyword arguments are those from LeapFrogIntegrator._attributable along with "log", "normalDeviateGenerator", "removeRotationTranslation" and "trajectories".
def LeapFrogDynamics_SystemGeometry ( system, **options ):
    """Molecular dynamics using the leap-frog algorithm."""
    ( of, options, log ) = _SetUpSimulation ( system, LeapFrogIntegrator._attributable, {}, options )
    # . The system has symmetry.
    if hasattr ( system, "symmetry" ) and ( system.symmetry is not None ):
        # . For the moment restrict calculation to systems with cubic symmetry.
        if options.get ( "pressureControl", False ) and not isinstance ( system.symmetry.crystalSystem, CrystalSystemCubic ):
            raise ValueError ( "Pressure coupling only works currently for systems with cubic symmetry." )
    # . Turn off pressure coupling for systems with no symmetry.
    else:
        options["pressureControl"] = False
    of.VelocitiesAssign ( options["temperature"], normalDeviateGenerator = options.pop ( "normalDeviateGenerator", None ) )
    optimizer = LeapFrogIntegrator.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Molecular dynamics using the velocity Verlet algorithm.
#===================================================================================================================================
# . Keyword arguments are those from VelocityVerletIntegrator._attributable along with "log", "normalDeviateGenerator", "removeRotationTranslation" and "trajectories".
def VelocityVerletDynamics_SystemGeometry ( system, **options ):
    """Molecular dynamics using the velocity Verlet algorithm."""
    ( of, options, log ) = _SetUpSimulation ( system, VelocityVerletIntegrator._attributable, {}, options )
    of.VelocitiesAssign ( options["temperatureStart"], normalDeviateGenerator = options.pop ( "normalDeviateGenerator", None ) )
    optimizer = VelocityVerletIntegrator.WithOptions ( **options )
    optimizer.Summary ( log = log )
    return optimizer.Iterate ( of, log = log )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""Define classes for Langevin velocity Verlet dynamics."""

import math

from  .MultiDimensionalDynamics       import MultiDimensionalDynamics      , \
                                             MultiDimensionalDynamicsState
from  .ObjectiveFunctionIteratorError import ObjectiveFunctionIteratorError
from ..Arrays                         import Array
from ..RandomNumbers                  import NormalDeviateGenerator        , \
                                             RandomNumberGenerator

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LangevinVelocityVerletIntegratorState ( MultiDimensionalDynamicsState ):
    """State for the Langevin velocity Verlet integrator."""

    _attributable = dict ( MultiDimensionalDynamicsState._attributable )
    _attributable.update ( { "sdR"  :  0.0 ,
                             "sdV1" :  0.0 ,
                             "sdV2" :  0.0 ,
                             "w"    : None } )

    def SetUp ( self ):
        """Set up the state."""
        super ( LangevinVelocityVerletIntegratorState, self ).SetUp ( )
        n = self.numberOfVariables
        if self.w is None: self.w = Array.WithExtent ( n ) ; self.w.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LangevinVelocityVerletIntegrator ( MultiDimensionalDynamics ):
    """Class for Langevin velocity Verlet dynamics."""

    _attributable = dict ( MultiDimensionalDynamics._attributable )
    _classLabel   = "Langevin Velocity Verlet Integrator"
    _stateObject  = LangevinVelocityVerletIntegratorState
    _summarizable = dict ( MultiDimensionalDynamics._summarizable )
    _attributable.update ( { "collisionFrequency"     :  25.0                 ,   
                             "facR1"                  :   0.0                 ,   
                             "facR2"                  :   0.0                 ,   
                             "facV1"                  :   0.0                 ,   
                             "facV2"                  :   0.0                 ,   
                             "facV3"                  :   0.0                 ,   
                             "normalDeviateGenerator" :  None                 ,   
                             "temperature"            : 300.0                 } ) 
    _summarizable.update ( { "collisionFrequency"     : "Collision Frequency" ,
                             "temperature"            : "Temperature"         ,
                             "timeStep"               : "Time Step"           } )

    def CalculateIntegrationConstants ( self, state ):
        """Calculate some quantities necessary for the integration."""
        # . Initialization.
        fact = self.collisionFrequency * self.timeStep
        # . Use exponential formulae.
        if fact > 0.009:
            # . Local factors.
            c0   = math.exp ( - fact )
            c1   = ( 1.0 - c0 ) / fact
            c2   = ( 1.0 - c1 ) / fact
            # . Random force constants.
            sdR  = math.sqrt ( self.timeStep**2 * ( 2.0 - ( 3.0 - 4.0 * c0 + c0 * c0 ) / fact ) / fact )
            sdV  = math.sqrt ( ( 1.0 - c0 * c0 ) )
            cRV1 = self.timeStep * ( 1.0 - c0 )**2 / ( fact * sdR * sdV )
            cRV2 = math.sqrt ( 1.0 - cRV1 * cRV1 )
        # . Use series expansions (valid to fact^5).
        # . The first is gamma * sdR = f^(3/2) * ( ). So bring out fact and divide by gamma to avoid problems when gamma small.
        # . Note that there was a slight error in cRV1 in initial implementations of the power series.
        else:
            # . Local factors.
            c0 = 1.0 - fact       + fact**2 /  2.0 - fact**3 /   6.0 + fact**4 /  24.0 - fact**5 /  120.0
            c1 = 1.0 - fact / 2.0 + fact**2 /  6.0 - fact**3 /  24.0 + fact**4 / 120.0 - fact**5 /  720.0
            c2 = 0.5 - fact / 6.0 + fact**2 / 24.0 - fact**3 / 120.0 + fact**4 / 720.0 - fact**5 / 5040.0
            # . Random force constants.
            sdR  = (     2.0                        - \
                         3.0 * fact    /        4.0 + \
                        67.0 * fact**2 /      320.0 - \
                       119.0 * fact**3 /     2560.0   ) * math.sqrt ( fact / 6.0 ) * self.timeStep
            sdV  = (     2.0 - fact                 + \
                         5.0 * fact**2 /       12.0 - \
                               fact**3 /        8.0 + \
                        79.0 * fact**4 /     2880.0   ) * math.sqrt ( fact / 2.0 )
            cRV1 = (     3.0           /        2.0 - \
                         3.0 * fact    /       16.0 - \
                        51.0 * fact**2 /     1280.0 + \
                        17.0 * fact**3 /     2048.0 + \
                     40967.0 * fact**4 / 11468800.0 - \
                     57203.0 * fact**5 / 91750400.0   ) / math.sqrt ( 3.0 )
            cRV2 = (     0.5                       +   
                         3.0 * fact    /       16.0 - \
                         9.0 * fact**2 /     1280.0 - \
                       109.0 * fact**3 /    10240.0 + \
                     10077.0 * fact**4 / 11468800.0 + \
                     14887.0 * fact**5 / 18350080.0   )
        # . Multiply sdR and sdV by the appropriate conversion factors.
        fact = math.sqrt ( state.objectiveFunction.TemperatureConversionFactor ( ) * self.temperature )
        state.sdR  =        sdR * fact
        state.sdV1 = cRV1 * sdV * fact
        state.sdV2 = cRV2 * sdV * fact
        # . Calculate some integration constants.
        self.facR1 = c1 * self.timeStep
        self.facR2 = c2 * self.timeStep**2
        self.facV1 = c0
        self.facV2 = ( c1 - c2 ) * self.timeStep
        self.facV3 = c2 * self.timeStep

    def Iteration ( self, state ):
        """Perform one dynamics step."""
        try:
            state.numberOfIterations += 1
            state.time               += self.timeStep
            # . Calculate the new variables.
            state.x.Add ( state.v, scale = self.facR1 )
            state.x.Add ( state.a, scale = self.facR2 )
            # . Calculate the half-velocities.
            state.v.Scale  ( self.facV1 )
            state.v.Add ( state.a, scale = self.facV2 )
            # . Add in random forces.
            self.RandomForces ( state )
            # . Calculate the function and accelerations at the new point.
            state.f = state.objectiveFunction.Accelerations ( state.x, state.a )
            # . Complete the velocity calculation.
            state.v.Add ( state.a, scale = self.facV3 )
            # . Calculate the temperature, etc.
            ( state.kineticEnergy, state.temperature ) = state.objectiveFunction.Temperature ( state.v )
            state.totalEnergy = state.f + state.kineticEnergy
        except Exception as error:
            state.HandleError ( error )

    def RandomForces ( self, state ):
        """Add in random forces."""
        # . First vector.
        self.normalDeviateGenerator.NextStandardDeviates ( state.w )
        state.objectiveFunction.ApplyConstraintsToVector ( state.w )
        state.x.Add ( state.w, scale = state.sdR  )
        state.v.Add ( state.w, scale = state.sdV1 )
        # . Second vector.
        self.normalDeviateGenerator.NextStandardDeviates ( state.w )
        state.objectiveFunction.ApplyConstraintsToVector ( state.w )
        state.v.Add ( state.w, scale = state.sdV2 )

    def TemperatureHandlingOptions ( self ):
        """Set up the temperature-handling options."""
        # . Check the temperature handling options.
        if ( self.collisionFrequency <= 0.0 ) or ( self.temperature < 0.0 ): raise ObjectiveFunctionIteratorError ( "Invalid temperature handling options." )
        # . Check for a valid normal deviate generator.
        if self.normalDeviateGenerator is None:
            self.normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithRandomSeed ( ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

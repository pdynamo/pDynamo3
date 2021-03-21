"""Define classes for velocity Verlet dynamics."""

import math

from .MultiDimensionalDynamics import MultiDimensionalDynamics      , \
                                      MultiDimensionalDynamicsState

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Temperature-handling options.
_TemperatureHandlingOptions = ( "Constant", "Exponential", "Linear" )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class VelocityVerletIntegratorState ( MultiDimensionalDynamicsState ):
    """State for the velocity Verlet integrator."""
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class VelocityVerletIntegrator ( MultiDimensionalDynamics ):
    """Class for velocity Verlet dynamics."""

    _attributable = dict ( MultiDimensionalDynamics._attributable )
    _classLabel   = "Velocity Verlet Integrator"
    _stateObject  = VelocityVerletIntegratorState
    _summarizable = dict ( MultiDimensionalDynamics._summarizable )
    _attributable.update ( { "facX1"                     :  0.0                          ,
                             "facX2"                     :  0.0                          ,
                             "facV"                      :  0.0                          ,
                             "temperatureScale"          :  0.0                          ,
                             "temperatureScaleFrequency" :    0                          ,
                             "temperatureScaleOption"    : None                          ,
                             "temperatureStart"          : None                          ,
                             "temperatureStop"           : None                          } )
    _summarizable.update ( { "temperatureScaleOption"    : "Temperature Handling"        ,
                             "temperatureScale"          : "Temperature Scale"           ,
                             "temperatureScaleFrequency" : "Temperature Scale Frequency" ,
                             "temperatureStart"          : "Temperature Start"           ,
                             "temperatureStop"           : "Temperature Stop"            } )

    def CalculateIntegrationConstants ( self, state ):
        """Calculate constants for the integration."""
        self.facX1 =       self.timeStep
        self.facX2 = 0.5 * self.timeStep**2
        self.facV  = 0.5 * self.timeStep

    def Iteration ( self, state ):
        """Perform one dynamics step."""
        state.numberOfIterations += 1
        state.time               += self.timeStep
        # . Calculate the new variables.
        state.x.Add ( state.v, scale = self.facX1 )
        state.x.Add ( state.a, scale = self.facX2 )
        # . Calculate the half-velocities.
        state.v.Add ( state.a, scale = self.facV )
        # . Calculate the function and accelerations at the new point.
        state.f = state.objectiveFunction.Accelerations ( state.x, state.a )
        # . Complete the velocity calculation.
        state.v.Add ( state.a, scale = self.facV )
        # . Calculate other quantities and perform temperature scaling if necessary.
        ( state.kineticEnergy, state.temperature ) = state.objectiveFunction.Temperature ( state.v )
        if ( self.temperatureScaleOption is not None ) and ( state.numberOfIterations % self.temperatureScaleFrequency == 0 ):
            scale = self.TargetTemperature ( state.time ) / state.temperature
            state.v.Scale ( math.sqrt ( scale ) )
            state.kineticEnergy *= scale
            state.temperature   *= scale
        # . Finish up.
        state.totalEnergy = state.f + state.kineticEnergy

    def TargetTemperature ( self, time ):
        """Get the target temperature for temperature scaling."""
        if   self.temperatureScaleOption == "Constant"    : tNeeded = self.temperatureStart
        elif self.temperatureScaleOption == "Exponential" : tNeeded = self.temperatureStart * math.exp ( self.temperatureScale * time )
        elif self.temperatureScaleOption == "Linear"      : tNeeded = self.temperatureStart +            self.temperatureScale * time
        return max ( 0.0, tNeeded )

    def TemperatureHandlingOptions ( self ):
        """Set up the temperature-handling options."""
        # . Convert the option to capitalized form.
        if isinstance ( self.temperatureScaleOption, str ): self.temperatureScaleOption = self.temperatureScaleOption.capitalize ( )
        # . No handling.
        if ( self.temperatureScaleFrequency <= 0 ) or ( self.temperatureScaleFrequency >= self.maximumIterations ) or \
           ( self.temperatureScaleOption is None ) or ( self.temperatureScaleOption not in _TemperatureHandlingOptions ):
            self.temperatureScaleOption = None
        # . Constant temperature.
        elif self.temperatureScaleOption == "Constant":
            self.temperatureStop = self.temperatureStart
        # . Exponential changes.
        elif self.temperatureScaleOption == "Exponential":
            self.temperatureScale = math.log ( self.temperatureStop / self.temperatureStart ) / ( float ( self.maximumIterations ) * self.timeStep )
        # . Linear changes.
        elif self.temperatureScaleOption == "Linear":
            self.temperatureScale =          ( self.temperatureStop - self.temperatureStart ) / ( float ( self.maximumIterations ) * self.timeStep )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

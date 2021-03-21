"""Define classes for leapfrog Verlet dynamics."""

import math

from .MultiDimensionalDynamics import MultiDimensionalDynamics      , \
                                      MultiDimensionalDynamicsState

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LeapFrogIntegratorState ( MultiDimensionalDynamicsState ):
    """State for the leap frog integrator."""

    _attributable = dict ( MultiDimensionalDynamicsState._attributable )
    _attributable.update ( { "pressure" : 0.0 ,
                             "volume"   : 0.0 } )

    def AddAccumulables ( self ):
        """Define extra quantities to accumulate."""
        self.accumulableAttributes.extend ( [ "pressure", "volume" ] )
        self.accumulableLabels.extend     ( [ "Pressure", "Volume" ] )
        self.averages.extend              ( [ 0.0, 0.0 ] )
        self.rmsDeviations.extend         ( [ 0.0, 0.0 ] )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LeapFrogIntegrator ( MultiDimensionalDynamics ):
    """Class for leap frog dynamics."""

    _attributable = dict ( MultiDimensionalDynamics._attributable )
    _classLabel   = "Leapfrog Verlet Integrator"
    _stateObject  = LeapFrogIntegratorState
    _summarizable = dict ( MultiDimensionalDynamics._summarizable )
    _attributable.update ( { "pressure"            :     1.0                ,
                             "pressureControl"     :   False                ,
                             "pressureCoupling"    :  2000.0                ,
                             "temperature"         :   300.0                ,
                             "temperatureControl"  :   False                ,
                             "temperatureCoupling" :     0.1                } )
    _summarizable.update ( { "pressure"            : "Pressure"             ,
                             "pressureControl"     : "Pressure Control"     ,
                             "pressureCoupling"    : "Pressure Coupling"    ,
                             "temperature"         : "Temperature"          ,
                             "temperatureControl"  : "Temperature Control"  ,
                             "temperatureCoupling" : "Temperature Coupling" } )

    def _CheckOptions ( self ):
        """Check the attributes."""
        super ( LeapFrogIntegrator, self )._CheckOptions ( )
        self.PressureHandlingOptions ( )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        if self.pressureControl: state.objectiveFunction.DefinePressure ( )
        super ( LeapFrogIntegrator, self ).Initialize ( state )
        ( state.pressure, state.volume ) = state.objectiveFunction.Pressure ( state.kineticEnergy )

    def Iteration ( self, state ):
        """Perform one dynamics step."""
        state.numberOfIterations += 1
        state.time               += self.timeStep
        # . Calculate the new velocities.
        state.v.Add ( state.a, scale = self.timeStep )
        # . Temperature scaling.
        self.TemperatureScale ( state )
        # . Calculate the new variables.
        state.x.Add ( state.v, scale = self.timeStep )
        # . Pressure scaling.
        self.PressureScale  ( state )
        # . Calculate the function and accelerations at the new point.
        state.f = state.objectiveFunction.Accelerations ( state.x, state.a )
        # . Calculate the temperature, pressure, etc.
        ( state.kineticEnergy, state.temperature ) = state.objectiveFunction.Temperature ( state.v             )
        ( state.pressure     , state.volume      ) = state.objectiveFunction.Pressure    ( state.kineticEnergy )
        state.totalEnergy = state.f + state.kineticEnergy

    def PressureHandlingOptions ( self ):
        """Set up the pressure-handling options."""
        self.pressureControl = self.pressureControl and ( self.pressure > 0.0 ) and ( self.pressureCoupling > 0.0 )

    def PressureScale ( self, state ):
        """Do pressure scaling."""
        if self.pressureControl:
            zetaP = 1.0 - ( self.timeStep / ( 3.0 * self.pressureCoupling ) ) * ( self.pressure - state.pressure ) # . Use first-order formula.
            state.objectiveFunction.VolumeScale ( zetaP )
            state.x.Scale                       ( zetaP )

    def StateFromObjectiveFunction ( self, objectiveFunction ):
        """Set up the state."""
        state = super ( LeapFrogIntegrator, self ).StateFromObjectiveFunction ( objectiveFunction )
        if self.pressureControl: state.AddAccumulables ( )
        return state

    def TemperatureHandlingOptions ( self ):
        """Set up the temperature-handling options."""
        self.temperatureControl = self.temperatureControl and ( self.temperature > 0.0 ) and ( self.temperatureCoupling > 0.0 )

    def TemperatureScale ( self, state ):
        """Do temperature scaling."""
        if self.temperatureControl:
            zetaT = math.sqrt ( 1.0 + ( self.timeStep / self.temperatureCoupling ) * ( ( self.temperature / state.temperature ) - 1.0 ) )
            state.v.Scale ( zetaT )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""Define classes for multidimensional dynamics algorithms."""

import math

from   pCore                     import Align                          , \
                                        logFile                        , \
                                        LogFileActive
from  .ObjectiveFunctionIterator import ObjectiveFunctionIterator      , \
                                        ObjectiveFunctionIteratorState
from ..Arrays                    import Array

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Output options.
_TableColumnWidth = 20

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MultiDimensionalDynamicsState ( ObjectiveFunctionIteratorState ):
    """State for multidimensional dynamics integrators."""

    _attributable = dict ( ObjectiveFunctionIteratorState._attributable )
    _attributable.update ( { "a"                     : None ,
                             "accumulableAttributes" : None ,
                             "accumulableLabels"     : None ,
                             "averages"              : None ,
                             "kineticEnergy"         :  0.0 ,
                             "numberOfSamples"       :  0   ,
                             "rmsDeviations"         : None ,
                             "temperature"           :  0.0 ,
                             "time"                  :  0.0 ,
                             "totalEnergy"           :  0.0 ,
                             "v"                     : None } )

    def DefineAccumulables ( self, **options ):
        """Define the quantities to accumulate."""
        self.accumulableAttributes = [ "totalEnergy"  , "kineticEnergy"  , "f"                , "temperature" ]
        self.accumulableLabels     = [ "Total Energy" , "Kinetic Energy" , "Potential Energy" , "Temperature" ]

    def SetUp ( self ):
        """Set up the state."""
        # . Arrays.
        n = self.numberOfVariables
        if self.a is None: self.a = Array.WithExtent ( n ) ; self.a.Set ( 0.0 )
        # . Accumulables.
        self.DefineAccumulables ( )
        n = len ( self.accumulableAttributes )
        self.averages      = [ 0.0 for i in range ( n ) ]
        self.rmsDeviations = [ 0.0 for i in range ( n ) ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MultiDimensionalDynamics ( ObjectiveFunctionIterator ):
    """Class for multidimensional dynamics."""

    _attributable = dict ( ObjectiveFunctionIterator._attributable )
    _classLabel   = "Multi-Dimensional Dynamics"
    _stateObject  = MultiDimensionalDynamicsState
    _summarizable = dict ( ObjectiveFunctionIterator._summarizable )
    _attributable.update ( { "timeStep" : 0.0         } )
    _summarizable.update ( { "timeStep" : "Time Step" } )

    def _CheckOptions ( self ):
        """Check the attributes."""
        self.TemperatureHandlingOptions ( )

    def CalculateIntegrationConstants ( self, state ): pass

    def Initialize ( self, state ):
        """Initialization before iteration."""
        # . Calculate integration constants.
        self.CalculateIntegrationConstants ( state )
        # . Get velocities.
        state.v = state.objectiveFunction.VelocitiesAllocate ( )
        # . First function evaluation.
        state.objectiveFunction.VariablesGet ( state.x )
        state.f = state.objectiveFunction.Accelerations ( state.x, state.a )
        # . Initialize the data at zero time.
        state.time = 0.0
        ( state.kineticEnergy, state.temperature  ) = state.objectiveFunction.Temperature ( state.v )
        state.totalEnergy = state.f + state.kineticEnergy

    def LogIteration ( self, state ):
        """Log an iteration."""
        state.objectiveFunction.LogIteration ( state.numberOfIterations, identifier = state.time )
        if state.table is not None:
            # . Statistics.
            values = []
            for ( i, attribute ) in enumerate ( state.accumulableAttributes ):
                value = getattr ( state, attribute )
                values.append ( value )
                state.averages     [i] += value
                state.rmsDeviations[i] += value**2
            state.numberOfSamples += 1
            # . Output.
            if state.numberOfIterations % self.logFrequency == 0:
                state.table.Entry ( "{:20.8f}".format ( state.time ) )
                for value in values: state.table.Entry ( "{:20.8f}".format ( value ) )

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        state.objectiveFunction.LogStart ( )
        if ( self.logFrequency > 0 ) and LogFileActive ( log ):
            state.log   = log
            state.table = log.GetTable ( columns = ( len ( state.accumulableAttributes ) + 1 ) * [ _TableColumnWidth ] )
            state.table.Start ( )
            state.table.Title ( self.__class__._classLabel + " Results" )
            state.table.Heading ( "Time" )
            for tag in state.accumulableLabels: state.table.Heading ( tag )

    def LogStop ( self, state ):
        """Stop logging."""
        state.objectiveFunction.LogStop ( )
        if state.log is not None:
            if state.numberOfSamples > 0:
                # . Statistics.
                fact = float ( state.numberOfSamples )
                for i in range ( len ( state.accumulableAttributes ) ):
                    average      =                   state.averages     [i] / fact
                    rmsDeviation = math.sqrt ( max ( state.rmsDeviations[i] / fact - average**2, 0.0 ) )
                    state.averages     [i] = average
                    state.rmsDeviations[i] = rmsDeviation
                # . Output.
                state.table.Entry ( "Averages",       align = Align.Left )
                for value in state.averages     : state.table.Entry ( "{:20.8f}".format ( value ) )
                state.table.Entry ( "RMS Deviations", align = Align.Left )
                for value in state.rmsDeviations: state.table.Entry ( "{:20.8f}".format ( value ) )
            state.table.Stop ( )

    def TemperatureHandlingOptions ( self ): pass

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

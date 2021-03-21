"""System class with energy timings."""

from  pCore       import CPUTime       , \
                         logFile       , \
                         LogFileActive , \
                         Timings
from .EnergyModel import EnergyModelState
from .System      import System

# . For testing only.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemWithTimings ( System ):
    """A system with a timing analysis of its energy calculation."""

    _attributable = dict ( System._attributable )
    _classLabel   = "System With Timings"
    _attributable.update ( { "cpuTimer"            : None ,
                             "numberOfEnergyCalls" :    0 ,
                             "timings"             : None } )

    def _UpdateEnergyClosures ( self ):
        """Update the list of energy closures."""
        closures = []
        for key in self._energyModels.keys ( ):
            closures.extend ( self.__dict__[key].EnergyClosures ( self ) )
        self._energyClosures = [ ( closure, label ) for ( priority, closure, label ) in sorted ( closures, key = lambda x: x[0] ) ]

    def Energy ( self, doGradients = False, log = logFile ):
        """Calculate the energy and, optionally, the gradients for a system."""
        tGlobal = self.cpuTimer.Current ( )
        self.EnergyInitialize ( doGradients, log )
        for ( closure, label ) in self._energyClosures:
            tLocal = self.cpuTimer.Current ( )
            closure ( )
            self.timings[label] += ( self.cpuTimer.Current ( ) - tLocal )
        energy = self.EnergyFinalize ( )
        self.timings["Energy"]   += ( self.cpuTimer.Current ( ) - tGlobal )
        self.numberOfEnergyCalls += 1
        return energy

    def EnergyFinalize ( self ):
        """Energy finalization."""
        tLocal = self.cpuTimer.Current ( )
        energy = super ( SystemWithTimings, self ).EnergyFinalize ( )
        self.timings["Finalization"] += ( self.cpuTimer.Current ( ) - tLocal )
        return energy

    def EnergyInitialize ( self, doGradients, log ):
        """Energy initialization."""
        tLocal = self.cpuTimer.Current ( )
        super ( SystemWithTimings, self ).EnergyInitialize ( doGradients, log )
        self.timings["Initialization"] += ( self.cpuTimer.Current ( ) - tLocal )

    @classmethod
    def FromSystem ( selfClass, system ):
        """Constructor from system."""
        # . The input system is unusable after this.
        self = selfClass ( )
        self.__dict__.update ( system.__dict__ )
        # . Reassign states and update energy closures.
        for value in self.__dict__.values ( ):
            if ( value is not None ) and isinstance ( value, EnergyModelState ): value.target = self
        self._UpdateEnergyClosures ( )
        return self

    def TimingStart ( self ):
        """Start timing."""
        self.cpuTimer            = CPUTime ( )
        self.numberOfEnergyCalls = 0
        self.timings             = Timings ( )
        self.timings.Clear ( )

    def TimingStop ( self ):
        """Stop timing."""
        self.timings["Total"] += self.cpuTimer.Current ( )

    def TimingSummary ( self, log = logFile, orderByMagnitude = False ):
        """Timing summary."""
        self.TimingStop ( )
        self.timings.Summary ( log = log, orderByMagnitude = orderByMagnitude )
        if LogFileActive ( log ): log.Paragraph ( "Number of Energy Calls = {:d}".format ( self.numberOfEnergyCalls ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

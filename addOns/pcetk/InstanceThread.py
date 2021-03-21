import threading, time

from pCore import logFile, LogFileActive

class InstanceThread (threading.Thread):
    """A class for the parallel calculation of electrostatic energy terms."""

    def __init__ (self, instance, log=logFile):
        """Constructor."""
        threading.Thread.__init__ (self)
        self.instance = instance
        self.log      = log
        self.time     = 0


    def run (self):
        """A single thread."""
        time0     = time.time ()
        instance  = self.instance

        instance.CalculateModelCompound (log=self.log)
        instance.CalculateProtein       (log=self.log)
        instance.CalculateGintr         (log=self.log)

        # . Calculate the time of execution
        self.time = (time.time () - time0)


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

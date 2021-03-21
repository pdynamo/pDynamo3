"""SGOF process classes and functions for various parallelization modes."""

from  pScientific.Arrays     import Array
from .ParallelizationOptions import HasFutures, HasMPI2, HasMultiprocessing

# . It does not seem that @staticmethod works for process functions.

#===================================================================================================================================
# . Properties.
#===================================================================================================================================
# . Process exit codes.
ProcessExitCode_Failure = -1
ProcessExitCode_Success =  0

# . Process task codes.
ProcessTaskCode_Exit              = -1
ProcessTaskCode_FunctionGradients =  0

# . For ORCA (and similar external programs), must set scratch/job combination per process. Maybe also invokation command.
# . Have a separate task to do this and also in futures function?

#===================================================================================================================================
# . Futures process.
#===================================================================================================================================
if HasFutures:
    def SGOFProcessFuturesFunctionGradients ( objectiveFunction, x ):
        """Futures function and gradients."""
        g = Array.WithExtent ( len ( x ) ) ; g.Set ( 0.0 )
        f = objectiveFunction.FunctionGradients ( x, g )
        return ( f, g )

#===================================================================================================================================
# . MPI2 process.
#===================================================================================================================================
if HasMPI2:
    class SGOFProcessMPI2:
        """A SGOF process based upon MPI2."""

        def __init__ ( self, parent = None, rank = None, options = {} ):
            """Constructor."""
            self.isActive = False
            self.parent   = parent
            self.rank     = rank
            self.results  = None
            self.__dict__.update ( options )

        def GetFunctionGradientsResults ( self, g ):
            """Get the results of a function/gradients task."""
            results = self.GetResults ( )
            if results is not None:
                if self.useBP:
                    self.parent.processPipe.Recv ( g, source = self.rank )
                    return results[0]
                else:
                    ( fP, gP ) = results
                    gP.CopyTo ( g )
                return fP
            else:
                return None

        def GetResults ( self ):
            """Get results if there are any."""
            results = self.results
            if results is not None:
                self.isActive = False
                self.results  = None
            return results

        def FunctionGradients ( self, x, wait = False ):
            """Function and gradients."""
            if not self.isActive:
                if self.useBP:
                    data = ( ProcessTaskCode_FunctionGradients, [] )
                    if wait:
                        self.parent.processPipe.send ( data, dest = self.rank )
                        self.parent.processPipe.Send ( x   , dest = self.rank )
                    else:
                        self.parent.processPipe.isend ( data, dest = self.rank )
                        self.parent.processPipe.Isend ( x   , dest = self.rank )
                else:
                    data = ( ProcessTaskCode_FunctionGradients, [ x ] )
                    if wait: self.parent.processPipe.send  ( data, dest = self.rank )
                    else:    self.parent.processPipe.isend ( data, dest = self.rank )
                self.isActive = True

        def HasResultsWaiting ( self ):
            """Check to see if there are results waiting."""
            waiting = self.isActive and self.parent.processPipe.Iprobe ( source = self.rank )
            if waiting:
                results = self.parent.processPipe.recv ( source = self.rank )
                if results[0] == ProcessExitCode_Failure:
                    raise ValueError ( "MPI2 process error: {:s}.".format ( results[1] ) )
                else:
                    self.results = results[1:]
            return waiting

        def Terminate ( self ):
            """Termination."""
            self.parent.processPipe.isend ( ( ProcessTaskCode_Exit, [] ), dest = self.rank )
            self.isActive = False
            self.parent   = None
            self.rank     = None

#===================================================================================================================================
# . Multiprocessing process.
#===================================================================================================================================
# . No useBP for the moment.
if HasMultiprocessing:
    import multiprocessing

    class SGOFProcessMultiprocessing:
        """A SGOF process based upon Python multiprocessing."""

        def __init__ ( self, objectiveFunction ):
            """Constructor."""
            ( parentReceiver, childSender ) = multiprocessing.Pipe ( duplex = False )
            ( childReceiver, parentSender ) = multiprocessing.Pipe ( duplex = False )
            process = multiprocessing.Process ( target = SGOFProcessMultiprocessingScript, args = ( objectiveFunction, childReceiver, childSender ) )
            process.start ( )
            self.isActive = False
            self.process  = process
            self.receiver = parentReceiver
            self.results  = None
            self.sender   = parentSender

        def GetFunctionGradientsResults ( self, g ):
            """Get the results of a function/gradients task."""
            results = self.GetResults ( )
            if results is not None:
                ( fP, gP ) = results
                gP.CopyTo ( g )
                return fP
            else:
                return None

        def GetResults ( self ):
            """Get results if there are any."""
            results = self.results
            if results is not None:
                self.isActive = False
                self.results  = None
            return results

        def FunctionGradients ( self, x, wait = False ):
            """Function and gradients."""
            if not self.isActive:
                self.sender.send ( ( ProcessTaskCode_FunctionGradients, [ x ] ) )
                self.isActive = True

        def HasResultsWaiting ( self ):
            """Check to see if there are results waiting."""
            if self.receiver.poll ( ):
                results = self.receiver.recv ( )
                if results[0] == ProcessExitCode_Failure:
                    raise ValueError ( "Multiprocessing process error: {:s}.".format ( results[1] ) )
                else:
                    self.results = results[1:]
            return ( self.results is not None )

        def Terminate ( self ):
            """Termination."""
            self.sender.send ( ( ProcessTaskCode_Exit, [] ) )
            self.receiver.close    ( )
            self.sender.close      ( )
# . This is not needed as it kills the process.
#            self.process.terminate ( )

#===================================================================================================================================
# . Multiprocessing remote script.
#===================================================================================================================================
    def SGOFProcessMultiprocessingScript ( objectiveFunction, receiver, sender ):
        """Multiprocessing script."""
        # . Initialization.
        g = Array.WithExtent ( objectiveFunction.numberOfVariables ) ; g.Set ( 0.0 )
        # . Poll for tasks.
        while True:
            if receiver.poll ( ):
                ( task, arguments ) = receiver.recv ( )
                if task == ProcessTaskCode_Exit:
                    break
                elif task == ProcessTaskCode_FunctionGradients:
                    try:
                        x = arguments[0]
                        f = objectiveFunction.FunctionGradients ( x, g )
                        sender.send ( ( ProcessExitCode_Success, f, g ) )
                    except Exception as error:
                        sender.send ( ( ProcessExitCode_Failure, error[0] ) )
                else:
                    sender.send ( ( ProcessExitCode_Failure, "Unrecognized task." ) )
        # . Finish up.

#===================================================================================================================================
# . Run the script.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

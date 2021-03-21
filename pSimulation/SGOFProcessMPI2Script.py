"""Script for remote SGOF MPI2 processes."""

from  pScientific.Arrays     import Array
from .ParallelizationOptions import HasMPI2
from .SGOFProcess            import ProcessExitCode_Failure, ProcessExitCode_Success, ProcessTaskCode_Exit, ProcessTaskCode_FunctionGradients

#===================================================================================================================================
# . Class.
#===================================================================================================================================
if HasMPI2:
    from mpi4py import MPI

    class RemoteSGOFMPI2Process:

        def __init__ ( self ):
            """Constructor."""
            self.log   = None
            self.pipe  = MPI.Comm.Get_parent ( )
            self.rank  = self.pipe.Get_rank  ( )
            self.size  = self.pipe.Get_size  ( )
            self.useBP = False

        def Error ( self, message ):
            """Signal an error."""
            self.pipe.isend ( ( ProcessExitCode_Failure, message ), dest = 0 )

        def Initialize ( self ):
            """Initialization."""
            ( self.objectiveFunction, options ) = self.pipe.bcast ( None, root = 0 )
            self.__dict__.update ( options )
            self.g = Array.WithExtent ( self.objectiveFunction.numberOfVariables ) ; self.g.Set ( 0.0 )
            if self.useBP:
                self.x = Array.WithExtent ( self.objectiveFunction.numberOfVariables ) ; self.x.Set ( 0.0 )

        def FunctionGradients ( self, *arguments ):
            """Function and gradients."""
            try:
                if self.useBP:
                    x = self.x
                    self.pipe.Recv ( x, source = 0 )
                else:
                    x = arguments[0]
                f = self.objectiveFunction.FunctionGradients ( x, self.g )
                if self.useBP:
                    self.pipe.isend ( ( ProcessExitCode_Success, f ), dest = 0 )
                    self.pipe.Isend ( self.g, dest = 0 )
                else:
                    self.pipe.isend ( ( ProcessExitCode_Success, f, self.g ), dest = 0 )
            except Exception as error:
                self.Error ( error[0] )

        def Receive ( self ):
            """Receive data."""
            return self.pipe.recv ( source = 0 )

        def Terminate ( self ):
            """Termination."""
            self.pipe.Disconnect ( )

else:
    class RemoteSGOFMPI2Process:

        def __init__ ( self ):
            raise ValueError ( "MPI2 SGOF processes unavailable." )

        def Error ( self, message ): pass

        def Initialize ( self ): pass

        def FunctionGradients ( self, x ): pass

        def Receive ( self ): return ( ProcessTaskCode_Exit, tuple ( ) )

        def Terminate ( self ): pass

#===================================================================================================================================
# . Run the script.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Activate a process.
    process = RemoteSGOFMPI2Process ( )
    process.Initialize ( )
    # . Poll for tasks.
    while True:
        ( task, arguments ) = process.Receive ( )
        if   task == ProcessTaskCode_Exit:              break
        elif task == ProcessTaskCode_FunctionGradients: process.FunctionGradients ( *arguments )
        else:                                           process.Error ( "Unrecognized task." )
    # . Finish up.
    process.Terminate ( )

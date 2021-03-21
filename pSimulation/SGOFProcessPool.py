"""System geometry objective function processes."""

from  pCore                  import AttributableObject , \
                                    logFile            , \
                                    LogFileActive      , \
                                    SummarizableObject
from  pMolecule              import SystemGeometryObjectiveFunction
from .ParallelizationOptions import HasFutures         , \
                                    HasMPI2            , \
                                    HasMultiprocessing

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_PoolTypes = {}

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SGOFProcessPool ( AttributableObject ):
    """Base (serial) class for a pool of processes that requires a SystemGeometryObjectiveFunction."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "maximumProcesses"  :    1 ,
                             "objectiveFunction" : None ,
                             "processes"         : None } )

    def Finalize ( self ):
        """Finalization."""
        pass

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction, maximumProcesses = 1 ):
        """Constructor given an objective function."""
        if not isinstance ( objectiveFunction, SystemGeometryObjectiveFunction ): raise ValueError ( "Invalid objective function." )
        self                   = selfClass ( )
        self.maximumProcesses  = max ( 1, maximumProcesses )
        self.objectiveFunction = objectiveFunction
        self.Initialize ( )
        return self

    def FunctionGradients ( self, x, g ):
        """Function and gradients for the pool."""
        of = self.objectiveFunction
        return [ of.FunctionGradients ( xI, gI ) for ( xI, gI ) in zip ( x, g ) ]

    def Initialize ( self ):
        """Initialization."""
        pass

# . Add to the pool types.
_PoolTypes = { "Serial" : SGOFProcessPool }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . This is less efficient as processes are spawned at each method invocation. Nevertheless it is easy and maybe OK for one shot pools.
if HasFutures:
    import concurrent.futures
    from .SGOFProcess import SGOFProcessFuturesFunctionGradients

    class SGOFProcessPoolFutures ( SGOFProcessPool ):
        """SystemGeometryObjectiveFunction process pool based upon Python concurrent futures."""

        _attributable = dict ( SGOFProcessPool._attributable )
        _attributable.update ( { "executor" : None } )

        def Finalize ( self ):
            """Finalization."""
            self.executor.shutdown ( )

        def FunctionGradients ( self, x, g ):
            """Function and gradients for the pool."""
            futures = dict ( ( self.executor.submit ( SGOFProcessFuturesFunctionGradients, self.objectiveFunction, x ), i ) for ( i, x ) in enumerate ( x ) )
            f       = [ 0.0 for i in range ( len ( x ) ) ]
            for future in concurrent.futures.as_completed ( futures ):
                i = futures[future]
                if future.exception ( ) is not None:
                    raise ValueError ( "Futures process error: {:s}.".format ( future.exception ( ) ) )
                else:
                    ( fP, gP ) = future.result ( )
                    f[i] = fP
                    gP.CopyTo ( g[i] )
            return f

        def Initialize ( self ):
            """Initialization."""
            self.executor = concurrent.futures.ProcessPoolExecutor ( max_workers = self.maximumProcesses )

    # . Add to the pool types.
    _PoolTypes["Futures"] = SGOFProcessPoolFutures

#===================================================================================================================================
# . Class.
#===================================================================================================================================
if HasMPI2 or HasMultiprocessing:
    class SGOFProcessPoolSpawning ( SGOFProcessPool ):
        """SystemGeometryObjectiveFunction process pool based upon the spawning of processes."""

        # . "useBP" is a flag to use the buffer protocol (new only).
        # . This seems to work for MPI2 but not Multiprocessing (Python 2.7).
        # . In any case, it doesn't seem to make a huge difference.
        _attributable = dict ( SGOFProcessPool._attributable )
        _processClass = None
        _attributable.update ( { "useBP" : False } )

        def Finalize ( self ):
            """Finalization."""
            for process in self.processes:
                process.Terminate ( )

        def FunctionGradients ( self, x, g ):
            """Function and gradients for the pool."""
            # . Initialization.
            numberOfTasks     = len ( x )
            numberOfProcesses = min ( self.maximumProcesses, numberOfTasks )
            # . Initial loading of processes.
            assigned = []
            for t in range ( numberOfProcesses ):
                self.processes[t].FunctionGradients ( x[t], wait = True )
                assigned.append ( t )
            # . Poll for results.
            f          = [ 0.0 for i in range ( numberOfTasks ) ]
            nextT      = numberOfProcesses
            numberDone = 0
            while numberDone < numberOfTasks:
                for p in range ( numberOfProcesses ):
                    process = self.processes[p]
                    # . Receive data.
                    if process.HasResultsWaiting ( ):
                        t    = assigned[p]
                        f[t] = process.GetFunctionGradientsResults ( g[t] )
                        numberDone += 1
                        # . Send data.
                        if nextT < numberOfTasks:
                            process.FunctionGradients ( x[nextT] )
                            assigned[p] = nextT
                            nextT += 1
            # . Finish up.
            return f

        def Initialize ( self ):
            """Initialization."""
            self.processes = [ self.__class__._processClass ( self.objectiveFunction ) for i in range ( self.maximumProcesses ) ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
if HasMPI2:
    import os, os.path, sys
    from mpi4py       import MPI
    from .SGOFProcess import SGOFProcessMPI2

    _scriptPath = os.path.split ( __file__ )[0] # . In the same directory.

    class SGOFProcessPoolMPI2 ( SGOFProcessPoolSpawning ):
        """SystemGeometryObjectiveFunction process pool based upon MPI2."""

        _attributable  = dict ( SGOFProcessPoolSpawning._attributable )
        _processClass  = SGOFProcessMPI2
        _processScript = os.path.join ( _scriptPath, "SGOFProcessMPI2Script.py" )
        _attributable.update ( { "processPipe" : None } )

        def Finalize ( self ):
            """Finalization."""
            MPI.COMM_WORLD.Barrier ( )
            for process in self.processes: process.Terminate ( )
            self.processPipe.Disconnect ( )
            MPI.COMM_WORLD.Barrier ( )
            self.processPipe = None
            self.processes   = None

        def Initialize ( self ):
            """Initialization."""
            # . Check the MPI environment.
            rank = MPI.COMM_WORLD.Get_rank ( )
            size = MPI.COMM_WORLD.Get_size ( )
            if ( rank != 0 ) or ( size > 1 ): raise ValueError ( "Invalid MPI environment (wrong rank or multiple processes)." )
            # . Spawn processes.
            commands  = self.maximumProcesses * [ sys.executable  ]
            arguments = self.maximumProcesses * [ [ self.__class__._processScript ] ]
            processes = self.maximumProcesses * [ 1 ]
            self.processPipe = MPI.COMM_WORLD.Spawn_multiple ( commands, arguments, processes )
            # . Broadcasting basic data.
            options = { "useBP" : self.useBP }
            self.processPipe.bcast ( ( self.objectiveFunction, options ), root = MPI.ROOT )
            MPI.COMM_WORLD.Barrier ( )
            # . Set up processes.
            self.processes = [ SGOFProcessMPI2 ( parent = self, rank = i, options = options ) for i in range ( self.maximumProcesses ) ]

    # . Add to the pool types.
    _PoolTypes["MPI2"] = SGOFProcessPoolMPI2

#===================================================================================================================================
# . Class.
#===================================================================================================================================
if HasMultiprocessing:
    from .SGOFProcess import SGOFProcessMultiprocessing

    class SGOFProcessPoolMultiprocessing ( SGOFProcessPoolSpawning ):
        """SystemGeometryObjectiveFunction process pool based upon Python multiprocessing."""

        _processClass = SGOFProcessMultiprocessing

    # . Add to the pool types.
    _PoolTypes["Multiprocessing"] = SGOFProcessPoolMultiprocessing

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SGOFProcessPoolFactory ( SummarizableObject ):
    """Factory class for creating process pools."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "SGOF Process Pool Factory"
    _summarizable = dict ( SummarizableObject._summarizable )
    _attributable.update ( { "maximumProcesses" :  1                  ,
                             "poolType"         : "Serial"            } )
    _summarizable.update ( { "maximumProcesses" : "Maximum Processes" ,
                             "poolType"         : "Pool Type"         } )

    def _CheckOptions ( self ):
        """Check options."""
        if  ( self.maximumProcesses <= 0 ) or ( self.poolType not in _PoolTypes.keys ( ) ):
            raise ValueError ( "Illegal SGOF process pool factory option." )

    def PoolFromObjectiveFunction ( self, objectiveFunction ):
        """Create a pool given an objective function."""
        poolClass = _PoolTypes[self.poolType]
        pool      = poolClass.FromObjectiveFunction ( objectiveFunction, maximumProcesses = self.maximumProcesses )
        return pool

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

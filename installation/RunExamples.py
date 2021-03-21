"""Run sets of example scripts."""

from argparse import ArgumentParser
from pCore    import CPUTime, TestSets_Get, TestSets_Run

# . List files option (e.g. -lf/--listFiles).
# . Also run files option (e.g. -s set + targets are script file names).

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
class InstallationOptions:
    """Class for determining installation options."""

    def __init__ ( self ):
        """Constructor."""
        self.maximumThreads = self.MaximumNumberOfThreads  ( )

    @classmethod
    def FromCommandLine ( selfClass ):
        """Constructor by getting arguments from a command line."""
        # . Initialization.
        self = selfClass ( )
        # . Parse the command line.
        parser = ArgumentParser ( epilog = "The arguments are optional but, if present, should be the names of the example sets to run. By default all short scripts in all example sets are run." )
        parser.add_argument     ( "-a" , "--all"     , action = "store_true" , dest = "doAll"    , default = False , help = "run all scripts (long and short)"    )
        parser.add_argument     ( "-l" , "--list"    , action = "store_true" , dest = "listSets" , default = False , help = "list example sets only"              )
        parser.add_argument     (        "--threads" , type   = int          , dest = "threads"  , default = 1     , help = "the number of threads"               )
        parser.add_argument     (        "target"    , nargs  = "*"          ,                                       help = "the name of the example sets to run" )
        arguments = parser.parse_args ( )
        # . Find whether to do long scripts.
        self.doLong   = arguments.doAll
        # . Find out whether to list sets only.
        self.listSets = arguments.listSets
        # . Find the number of threads.
        self.numberOfThreads = min ( max ( arguments.threads, 1 ), self.maximumThreads )
        # . Targets have been specified.
        if len ( arguments.target ) > 0:
            self.targets = sorted ( arguments.target )
        else:
            self.targets = None
        # . Finish up.
        return self

    def MaximumNumberOfThreads ( self ):
        """Find the maximum number of threads."""
        try:
            import multiprocessing
            return multiprocessing.cpu_count ( )
        except:
            return 1

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def RunTests ( ):
    """Main function."""
    # . Get the run options.
    options = InstallationOptions.FromCommandLine ( )
    # . Find the example sets to run.
    ( testSets, unknownSets ) = TestSets_Get ( options.targets )
    # . There are unknown sets.
    if len ( unknownSets ) > 0:
        print ( "\nUnknown example sets were specified:" )
        for name in unknownSets: print ( name )
        print ( "" )
    # . No sets were specified.
    elif len ( testSets ) == 0:
        print ( "\nNo example sets were found.\n" )
    # . There are sets but only for printing.
    elif options.listSets:
        print ( "\nExample sets:\n" )
        for testSet in testSets: print ( "  " + testSet )
        print ( "" )
    # . There are sets.
    else:
        TestSets_Run ( testSets, doLong = options.doLong )

#===================================================================================================================================
# . Run the script.
#===================================================================================================================================
if __name__ == "__main__" :
    RunTests ( )

"""pDynamo installation script."""

import time

from collections      import defaultdict
from InstallUtilities import FindRootDirectory    , \
                             GetRootDirectories   , \
                             InstallationOptions  , \
                             InstallationSummary  , \
                             RemoveBuildDirectory , \
                             WriteShellFiles

# . Should print error if Cython or PyYaml are not installed!
# . Also python development headers.\

# . Also use yaml for all configuration files (dependencies and libraries).

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def InstallpDynamo ( ):
    """Main function."""

    # . Initialization.
    report    = defaultdict ( int )
    startTime = time.time ( )

    # . Get the installation options.
    options = InstallationOptions.FromCommandLine ( )

    # . Find the root directory.
    FindRootDirectory ( )

    # . Get the packages for processing in the correct order.
    if options.doCLibraries or options.doExtensions or options.doPyrex or options.listPackages:
        packages = GetRootDirectories ( options.packageNames, update = options.update )

    # . Package listing only.
    if options.listPackages:
        if len ( packages ) > 0:
            print ( "\nPackages to install in order of dependency:\n" )
            for package in packages: print ( "  " + package.name )
            print ( "" )
        else:
            print ( "\nThere are no packages to install.\n" )

    # . Other processing.
    else:

        try:

            # . Install packages with extensions.
            if options.doCLibraries or options.doExtensions or options.doPyrex:
                for package in packages:
                    if package.HasExtensions ( ):
                        package.ProcessWithOptions ( options, report )

                # . Clear up.
                if options.doClearUp: RemoveBuildDirectory ( )

            # . Shell.
            if options.doShellFiles: WriteShellFiles ( report )

            # . Terminating message.
            message = "Installation terminated successfully."

        except ValueError as e:
            if hasattr ( e, "args" ) and ( len ( e.args ) > 0 ): message = "Exiting: " + e.args[0]
            else:                                                message = "Exiting with error."

        # . Finish up.
        totalTime = time.time ( ) - startTime
        InstallationSummary ( report, totalTime, message )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Run the script.
    InstallpDynamo ( )

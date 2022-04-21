"""Surface crossing test."""

# . This algorithm is very sensitive and needs work (including derivatives)!

import glob, math, os.path

from Definitions                            import dataPath
from pBabel                                 import ImportCoordinates3    , \
                                                   ImportSystem
from pCore                                  import Align                 , \
                                                   CPUTime               , \
                                                   logFile               , \
                                                   LogFileActive         , \
                                                   TestScriptExit_Fail
from pMolecule                              import SEAMObjectiveFunction
from pMolecule.QCModel                      import DIISSCFConverger      , \
                                                   ElectronicState       , \
                                                   QCModelMNDO
from pScientific.ObjectiveFunctionIterators import QuasiNewtonMinimizer

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
checkGradients = False
isOK           = True
molPath        = os.path.join ( dataPath, "mol" )

# . Options - the results are very sensitive to the SCF and optimizer tolerances!
converger = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-8, maximumIterations = 500 )
qcModel   = QCModelMNDO.WithOptions      ( converger = converger, hamiltonian = "am1" )
singlet   = ElectronicState.WithOptions  ( charge = 1, isSpinRestricted = False, multiplicity = 1 ) # . Unrestricted!
triplet   = ElectronicState.WithOptions  ( charge = 1, isSpinRestricted = False, multiplicity = 3 )

# . Optimizer.
optimizer = QuasiNewtonMinimizer.WithOptions ( logFrequency         =     1 ,
                                               maximumIterations    =  2500 ,
                                               rmsGradientTolerance =  0.01 )
optimizer.Summary ( )

# . Set up the system.
system                 = ImportSystem ( os.path.join ( molPath, "phenylCation.mol" ) )
system.electronicState = singlet
system.label           = "Phenyl Cation"
system.DefineQCModel ( qcModel )
system.Summary ( )
system.Energy ( )

# . Check both methods.
numberNotConverged = 0
results            = {}
for method in ( "GP", "PF" ):

    # . Reset coordinates.
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( molPath, "phenylCation.mol" ) )
    system.Energy ( )

    # . Set up the objective function.
    seamOF = SEAMObjectiveFunction.FromSystem ( system, singlet, None, triplet, None, method = method )

    # . Gradients.
    if checkGradients:
        seamOF.RemoveRotationTranslation ( )
        seamOF.TestGradients ( delta = 1.0e-05 ) # . Works with 1.0e-10 density tolerance.

    # . Minimize.
    cpu                = CPUTime ( )
    report             = optimizer.Iterate ( seamOF )
    report["CPU Time"] = cpu.CurrentAsString ( )

    # . Final energies.
    ( f1, f2 )         = seamOF.Energies ( doGradients = True )
    report["Energy 1"] = f1
    report["Energy 2"] = f2
    results[method]    = report
    if not report.get ( "Converged", False ): numberNotConverged += 1

# . Print out a summary of the results.
table = logFile.GetTable ( columns = [ 10, 20, 20, 10, 10, 20 ] )
table.Start   ( )
table.Title   ( "Surface Crossing Optimizations" )
table.Heading ( "Method"    )
table.Heading ( "State Energies", columnSpan = 2 )
table.Heading ( "Converged" )
table.Heading ( "Calls"     )
table.Heading ( "Time"      )
for method in ( "GP", "PF" ):
    report = results[method]
    table.Entry ( method, align = Align.Left )
    table.Entry ( "{:20.1f}".format ( report["Energy 1"      ] ) )
    table.Entry ( "{:20.1f}".format ( report["Energy 2"      ] ) )
    table.Entry ( "{!r}"    .format ( report.get ( "Converged", False ) ) )
    table.Entry ( "{:d}"    .format ( report["Function Calls"] ) )
    table.Entry (                     report["CPU Time"      ]   )
table.Stop ( )

# . Footer.
logFile.Footer ( )
if numberNotConverged != 0: TestScriptExit_Fail ( )

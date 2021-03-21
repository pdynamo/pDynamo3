"""Example 10 - XHTML output."""

from Definitions import *

# . Get a pretty header.
header = __file__[0:-3].replace ( "Example", "Example " ).replace ( "XHTML", " - XHTML Output" )

# . Start output.
logFile = XHTMLLogFileWriter.WithOptions ( title = "pDynamo Book Examples" )
logFile.Header ( "pDynamo " + header )

# . Define the molecule and its QC model.
molecule = ImportSystem ( os.path.join ( xyzPath, "bala_c7eq.xyz" ) )
molecule.DefineQCModel ( QCModelMNDO.WithOptions ( hamiltonian = "am1" ) )
molecule.Summary ( log = logFile )

# . Save a copy of the starting coordinates.
coordinates3 = Clone ( molecule.coordinates3 )

# . Determine the starting energy.
eStart = molecule.Energy ( log = logFile )

# . Optimization.  
ConjugateGradientMinimize_SystemGeometry ( molecule                       ,
                                           log                  = logFile ,
                                           logFrequency         =  100    ,
                                           maximumIterations    = 2000    ,
                                           rmsGradientTolerance =    0.1  )

# . Determine the final energy.
eStop = molecule.Energy ( log = logFile )

# . Determine the RMS coordinate deviation between the optimized and unoptimized structures.
masses = Array.FromIterable ( [ atom.mass for atom in molecule.atoms ] )
coordinates3.Superimpose ( molecule.coordinates3, weights = masses )
rms = coordinates3.RootMeanSquareDeviation ( molecule.coordinates3, weights = masses )

# . Print the results.
table = logFile.GetTable ( columns = [ 30, 30 ] )
table.Start ( )
table.Title ( "Minimization Results" )
table.Entry ( "Energy Change",            align = Align.Left )
table.Entry ( "{:20.4f}".format ( eStop - eStart ) )
table.Entry ( "RMS Coordinate Deviation", align = Align.Left )
table.Entry ( "{:20.4f}".format ( rms ) )
table.Stop ( )

# . Finish output.
logFile.Footer ( )

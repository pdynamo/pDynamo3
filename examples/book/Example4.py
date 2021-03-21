"""Example 4."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the list of structures.
xyzFiles = [ "bala_alpha.xyz", "bala_c5.xyz", "bala_c7ax.xyz", "bala_c7eq.xyz" ]

# . Define a molecule.
xyzFile  = xyzFiles.pop ( )
molecule = ImportSystem ( os.path.join ( xyzPath, xyzFile ) )
molecule.Summary ( )

# . Translate the system to its center of mass.
masses = Array.FromIterable ( [ atom.mass for atom in molecule.atoms ] )
molecule.coordinates3.TranslateToCenter ( weights = masses )

# . Calculate and print the inertia matrix before reorientation.
inertia = molecule.coordinates3.InertiaMatrix ( weights = masses )
labels  = [ "x", "y", "z" ]
ArrayPrint2D ( inertia, columnLabels = labels   ,
                        itemFormat   = "{:.3f}" ,
                        rowLabels    = labels   ,
                        title        = "Inertia Matrix Before Reorientation" )

# . Transform to principal axes.
molecule.coordinates3.ToPrincipalAxes ( weights = masses )

# . Calculate and print the inertia matrix after reorientation.
inertia = molecule.coordinates3.InertiaMatrix ( weights = masses )
ArrayPrint2D ( inertia, columnLabels = labels   ,
                        itemFormat   = "{:.3f}" ,
                        rowLabels    = labels   ,
                        title        = "Inertia Matrix After Reorientation" )

# . Define a table for the results.
table = logFile.GetTable ( columns = [ 20, 10, 10 ] )
table.Start ( )
table.Title ( "RMS Coordinate Deviations" )
table.Heading ( "Structure" )
table.Heading ( "Before"    )
table.Heading ( "After"     )

# . Loop over the remaining structures.
for xyzFile in xyzFiles:
    coordinates3 = ImportCoordinates3 ( os.path.join ( xyzPath, xyzFile ), log = None )
    rms0 = coordinates3.RootMeanSquareDeviation ( molecule.coordinates3, weights = masses )
    coordinates3.Superimpose ( molecule.coordinates3, weights = masses )
    rms1 = coordinates3.RootMeanSquareDeviation ( molecule.coordinates3, weights = masses )
    table.Entry ( xyzFile[0:-4], align = Align.Left )
    table.Entry ( "{:.2f}".format ( rms0 ) )
    table.Entry ( "{:.2f}".format ( rms1 ) )

# . Finish up the table.
table.Stop ( )

# . Footer.
logFile.Footer ( )

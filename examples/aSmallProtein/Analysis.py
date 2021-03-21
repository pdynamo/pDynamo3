"""Analysis of the molecular dynamics trajectories."""

from Definitions import *

# . Header.
logFile.Header ( )

# . PDB files and the number of cations to add.
_PDBPaths = (  "1UAO", "2E4E" )

# . Structures.
_Structures = ( "folded", "unfolded" )

# . Initialization.
results = []

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Loop over folded and unfolded systems.
    for structure in _Structures:

        # . Retrieve the system with the initial coordinates.
        system = Unpickle ( os.path.join ( outPath, pdbPath + "_" + structure + "_solvated.pkl" ) )
        system.Summary ( )
        system.Energy  ( )

        # . Get the atom masses and a selection for protein atoms only.
        masses  = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
        protein = AtomSelection.FromAtomPattern ( system, "A:*:*" )

        # . Calculate the radius of gyration.
        rg0 = system.coordinates3.RadiusOfGyration ( selection = protein, weights = masses )

        # . Save the starting coordinates.
        reference3 = Clone ( system.coordinates3 )

        # . Get the trajectory.
        trajectory = ImportTrajectory ( os.path.join ( outPath, pdbPath + "_" + structure + "_md1.mdcrd" ), system )
        trajectory.ReadHeader ( )

        # . Loop over the frames in the trajectory.
        rg  = []
        rms = []
        while trajectory.RestoreOwnerData ( ):
            system.coordinates3.Superimpose ( reference3, selection = protein, weights = masses )
            rg.append  ( system.coordinates3.RadiusOfGyration        (             selection = protein, weights = masses ) )
            rms.append ( system.coordinates3.RootMeanSquareDeviation ( reference3, selection = protein, weights = masses ) )

        # . Set up the statistics calculations.
        rgStatistics  = Statistics ( rg  )
        rmsStatistics = Statistics ( rms )

        # . Save the results.
        results.append ( ( pdbPath, structure, rg0, rgStatistics.mean , rgStatistics.standardDeviation , rgStatistics.maximum , rgStatistics.minimum  , \
                                                    rmsStatistics.mean, rmsStatistics.standardDeviation, rmsStatistics.maximum, rmsStatistics.minimum ) )

# . Output the results.
table = logFile.GetTable ( columns = 11 * [ 10 ] )
table.Start   ( )
table.Title   ( "Radii of Gyration (r) and RMS Deviations (d)" )
table.Heading ( "Protein" )
table.Heading ( "State"   )
table.Heading ( "r0"      )
table.Heading ( "<r>"     )
table.Heading ( "Std. r"  )
table.Heading ( "Max. r"  )
table.Heading ( "Min. r"  )
table.Heading ( "<d>"     )
table.Heading ( "Std. d"  )
table.Heading ( "Max. d"  )
table.Heading ( "Min. d"  )
for ( pdbPath, structure, r0, rmean, rstd, rmax, rmin, dmean, dstd, dmax, dmin ) in results:
    table.Entry ( pdbPath   )
    table.Entry ( structure )
    table.Entry ( "{:.2f}".format ( r0    ) )
    table.Entry ( "{:.2f}".format ( rmean ) )
    table.Entry ( "{:.2f}".format ( rstd  ) )
    table.Entry ( "{:.2f}".format ( rmax  ) )
    table.Entry ( "{:.2f}".format ( rmin  ) )
    table.Entry ( "{:.2f}".format ( dmean ) )
    table.Entry ( "{:.2f}".format ( dstd  ) )
    table.Entry ( "{:.2f}".format ( dmax  ) )
    table.Entry ( "{:.2f}".format ( dmin  ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )

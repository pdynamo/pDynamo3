"""Refinement of vacuum structures."""

from Definitions import *

# . Header.
logFile.Header ( )

# . PDB files.
_PDBPaths = ( "1UAO", "2E4E" )

# . Parameters for the minimization.
_Iterations   = 1000
_Logfrequency =  100
_Tolerance    =  1.0

# . Loop over structures.
for pdbPath in _PDBPaths:

    # . Retrieve the system.
    system = Unpickle ( os.path.join ( outPath, pdbPath + ".pkl" ) )
    system.Summary ( )
    system.Energy ( )

    # . Get selection for the hydrogens and heavy atoms.
    hydrogens = Selection.FromIterable ( [ index for ( index, atom ) in enumerate ( system.atoms ) if atom.atomicNumber == 1 ] )
    heavies   = hydrogens.Complement ( upperBound = len ( system.atoms ) )

    # . Minimize the coordinates of the hydrogens.
    system.freeAtoms = hydrogens
    ConjugateGradientMinimize_SystemGeometry ( system                               ,
                                               maximumIterations    = _Iterations   ,
                                               logFrequency         = _Logfrequency ,
                                               rmsGradientTolerance = _Tolerance    )
    system.freeAtoms = None

    # . Loop over minimizations.
    for forceConstant in ( 1000.0, 500.0, 250.0, 100.0, 50.0 ):

        # . Harmonically restrain heavy atoms.
        tethers = None
        if ( forceConstant is not None ):
            reference          = Clone ( system.coordinates3 )
            tetherEnergyModel  = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
            tethers            = RestraintModel ( )
            tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                       reference   = reference         , 
                                                                       selection   = heavies           )
        system.DefineRestraintModel ( tethers )

        # . Minimize.
        ConjugateGradientMinimize_SystemGeometry ( system                               ,
                                                   maximumIterations    = _Iterations   ,
                                                   logFrequency         = _Logfrequency ,
                                                   rmsGradientTolerance = _Tolerance    )
        system.Energy ( doGradients = True )

    # . Save the refined coordinates.
    ExportSystem ( os.path.join ( outPath, pdbPath + "_folded.xyz" ), system )

# . Footer.
logFile.Footer ( )

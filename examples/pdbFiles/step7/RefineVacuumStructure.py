"""Refine the vacuum structure for PKA."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions       import outPath
from pCore             import Clone                   , \
                              logFile                 , \
                              Pickle                  , \
                              Selection               , \
                              Unpickle
from pMolecule         import RestraintEnergyModel    , \
                              RestraintModel          , \
                              RestraintMultipleTether
from pMolecule.NBModel import NBModelCutOff
from pSimulation       import ConjugateGradientMinimize_SystemGeometry

# . Header.
logFile.Header ( )

# . Parameters.
_Iterations   = 1000
_Logfrequency =  100
_Tolerance    =  1.0

# . Get the system with an appropriate NB model.
system = Unpickle ( os.path.join ( outPath, "step6.pkl" ) )
system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
system.Summary ( )

# . Get a selection with the indices of the atoms whose coordinates were built.
hydrogens = [ index for ( index, atom ) in enumerate ( system.atoms ) if atom.atomicNumber == 1 ]
built     = Selection.FromIterable ( hydrogens + [ system.sequence.AtomIndex ( "A:LYS.8:CB" ), system.sequence.AtomIndex ( "I:SER.17:OG" ) ] )

# . Minimize the coordinates of these atoms.
system.freeAtoms = built
ConjugateGradientMinimize_SystemGeometry ( system,                               \
                                           maximumIterations    = _Iterations,   \
                                           logFrequency         = _Logfrequency, \
                                           rmsGradientTolerance = _Tolerance     )
system.freeAtoms = None

# . Get a selection corresponding to heavy atoms.
heavies = Selection.FromIterable ( [ index for ( index, atom ) in enumerate ( system.atoms ) if atom.atomicNumber > 1 ] )

# . Loop over minimizations.
for forceConstant in ( 1000.0, 500.0, 250.0, 100.0, 50.0, 10.0, 4.0, 1.0, 0.25, None ):

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
    ConjugateGradientMinimize_SystemGeometry ( system,                               \
                                               maximumIterations    = _Iterations,   \
                                               logFrequency         = _Logfrequency, \
                                               rmsGradientTolerance = _Tolerance     )
    system.Energy ( doGradients = True )

# . Save the system.
Pickle ( os.path.join ( outPath, "step7.pkl" ), system )

# . Footer.
logFile.Footer ( )

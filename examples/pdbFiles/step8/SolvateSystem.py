"""Solvate the system using a water box."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions               import outPath                                                                                                                
from pBabel                    import ExportSystem                                                                                                          
from pCore                     import Clone                                    , \
                                      logFile                                  , \
                                      Pickle                                   , \
                                      Selection                                , \
                                      Unpickle
from pMolecule                 import AtomSelection                            , \
                                      RestraintEnergyModel                     , \
                                      RestraintModel                           , \
                                      RestraintMultipleTether
from pMolecule.NBModel         import NBModelCutOff
from pScientific.RandomNumbers import NormalDeviateGenerator                   , \
                                      RandomNumberGenerator
from pSimulation               import ConjugateGradientMinimize_SystemGeometry , \
                                      LangevinDynamics_SystemGeometry          , \
                                      SolvateSystemBySuperposition                                                            

# . Header.
logFile.Header ( )

# . Parameters.
# . Refine options.
_Optimize = True
_Refine   = True
_Steps    = 2000

# . Define the NB model.
nbModel = NBModelCutOff.WithDefaults ( )

# . Get the solute system.
solute = Unpickle ( os.path.join ( outPath, "step8_a.pkl" ) )
solute.Summary ( )

# . Get the solvent system.
solvent = Unpickle ( os.path.join ( outPath, "waterBox.pkl" ) )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Create the solvated system.
solution       = SolvateSystemBySuperposition ( solute, solvent, reorientSolute = False )
solution.label = "Solvated PKA"
solution.DefineNBModel  ( nbModel )
solution.Summary ( )

# . Refinement.
if _Refine:

    # . Check energies.
    solute.Energy ( )
    for system in ( solvent, solution ):
        system.DefineNBModel ( nbModel )
        system.Energy ( )

    # . Get selections for the protein atoms in the A and I chains.
    protein   = ( AtomSelection.FromAtomPattern ( solution, "A:*:*" ) | AtomSelection.FromAtomPattern ( solution, "I:*:*" ) )
    hydrogens = Selection.FromIterable ( [ index for ( index, atom ) in enumerate ( solution.atoms ) if atom.atomicNumber == 1 ] )
    heavies   = protein - hydrogens

    # . Fix all protein atoms.
    solution.freeAtoms = protein.Complement ( upperBound = len ( solution.atoms ) )
    solution.Summary ( )

    # . Optimization.
    if _Optimize:
        ConjugateGradientMinimize_SystemGeometry ( solution                    ,
                                                   maximumIterations    =  200 ,
                                                   logFrequency         =   10 ,
                                                   rmsGradientTolerance = 10.0 )

    # . Define a normal deviate generator in a given state.
    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 517152 ) )

    # . Dynamics loop - fixed atoms then tether restraints.
    for forceConstant in ( None, 500.0, 100.0, 20.0, 4.0 ):

        # . Harmonically restrain protein heavy atoms.
        tethers = None
        if ( forceConstant is not None ):
            reference         = Clone ( solution.coordinates3 )
            tetherEnergyModel = RestraintEnergyModel.Harmonic ( 0.0, forceConstant )
            tethers           = RestraintModel ( )
            tethers["Tethers"] = RestraintMultipleTether.WithOptions ( energyModel = tetherEnergyModel ,
                                                                       reference   = reference         , 
                                                                       selection   = heavies           )
        solution.DefineRestraintModel ( tethers )

        # . Dynamics.
        LangevinDynamics_SystemGeometry ( solution                        ,
                                          collisionFrequency     =   25.0 ,
                                          logFrequency           =    100 ,
                                          normalDeviateGenerator = normalDeviateGenerator ,
                                          steps                  = _Steps ,
                                          temperature            =  300.0 ,
                                          timeStep               =  0.001 )

        # . Unfix atoms and remove tethers.
        solution.freeAtoms = None
        solution.DefineRestraintModel ( None )

# . Save the system.
Pickle ( os.path.join ( outPath, "step8_b.pkl" ), solution )

# . Print PDB file.
ExportSystem ( os.path.join ( outPath, "step8_b.pdb" ), solution )

# . Footer.
logFile.Footer ( )

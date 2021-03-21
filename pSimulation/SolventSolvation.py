"""Helper functions for constructing solvent boxes and solvating molecules.

At the moment these are purely geometric and need to be combined with MC or
MD equilibration to get the correct density.
"""

# . This needs to use the solvent/counterion database. Also need to include PDB names in the latter so can immediately construct sequence.

import math

from  random                    import sample
from  pCore                     import Align                                , \
                                       Clone                                , \
                                       logFile                              , \
                                       LogFileActive                        , \
                                       Selection
from  pScientific               import Constants                            , \
                                       PeriodicTable                        , \
                                       Units
from  pScientific.Arrays        import Array                                , \
                                       ArrayPrint                           , \
                                       ArrayPrint2D
from  pScientific.Geometry3     import Coordinates3                         , \
                                       CrossPairList_FromSingleCoordinates3 , \
                                       Matrix33                             , \
                                       RegularGrid                          , \
                                       Vector3
from  pScientific.RandomNumbers import RandomNumberGenerator
from  pScientific.Symmetry      import CrystalSystemCubic                   , \
                                       PeriodicBoundaryConditions
from .EditMergePrune            import MergeByAtom                          , \
                                       MergeRepeatByAtom                    , \
                                       PruneByAtom
from .SequenceUtilities         import CreateElementSequence                , \
                                       CreateHomogeneousIsolateSequence     , \
                                       DetermineUniqueEntityLabel           , \
                                       RenumberEntityComponents

#===================================================================================================================================
# . Add counterions to a system.
#===================================================================================================================================
# . Parameters.
# . Merge option.
_DoMerge = True

# . Number of tries for generating grid.
_NTries = 2

# . Print options.
_DoPrint = True

# . The reduction factor for grid point size.
_ReductionFactor = 0.95

# . Local function.
def _PlaceIonsAtGridPoints ( ions, grid, unoccupied ):
    """Place ions at grid points."""
#
# . Placement is currently done randomly. In principle it could also be done using energies using a procedure such as:
#
# - Calculate potentials (1/r,1/r^6,1/r^12) on unoccupied grid points due to protein.
# - Loop over ions (+/- alternate pairs):
#       Determine energy of each grid points as q * pot(1/r) + e * s^3 * pot(1/r^6) + ...
#       Find point with lowest energy and put ion at this point.
#       Update potentials due to the ion at this point.
#
# . Requires a single NB model method that calculates/updates potentials.
#
# . However, a similar result can be achieved by using random placement, doing a Monte Carlo (no waters) with
# . an NBModel with dielectric of 80 and then choosing low energy conformations.
#
    # . Number of ions.
    nions = len ( ions.atoms )

    # . Get a random number of points from the unoccupied sites.
    points = sample ( list ( unoccupied ), nions )

    # . Calculate the coordinates - these are implicitly sorted due to the selection.
    coordinates3 = Coordinates3.FromGrid ( grid, selection = Selection.FromIterable ( points ) )

    # . Shuffle again.
    points = sample ( list ( range ( nions ) ), nions )
    for ( n, p ) in enumerate ( points ):
        for i in range ( 3 ): ions.coordinates3[n,i] = coordinates3[p,i]

# . Function.
def AddCounterIons ( system, nAnions, anion, nCations, cation, boxDimensions, log = logFile ):
    """Add counter ions to a system."""

    # . Create a merged ion system.
    ions       = MergeByAtom ( nAnions * [ anion ] + nCations * [ cation ] )
    ions.label = "Counter Ions"
    if system.sequence is not None:
        entityLabel = DetermineUniqueEntityLabel ( system, label = "X" )
        CreateElementSequence ( ions, entityDescription = "Counter Ions", entityLabel = entityLabel )
    ions.Summary ( )

    # . Determine the size of the grid with a grid spacing in each direction of at least twice the maximum vdW radius (as an initial guess).
    nions               = len ( ions.atoms )
    maximumvdwradius    = max ( [ atom.vdwRadius for atom in ions.atoms ] )
    step                = 2.0 * maximumvdwradius
    ( origin, extents ) = system.coordinates3.EnclosingOrthorhombicBox ( )

    # . Readjust the origin and extents to take into account the box size.
    for ( i, ( o, extent, newsize ) ) in enumerate ( zip ( origin, extents, boxDimensions ) ):
        if newsize < extent: raise ValueError ( "The system is bigger than the requested box size." )
        center     = o + 0.5 * extent
        origin[i]  = center - 0.5 * newsize
        extents[i] = newsize

    # . Find the radii for the system atoms and increment by the maximum ion radius.
    radii = Array.FromIterable ( [ atom.vdwRadius for atom in system.atoms ] )
    radii.Increment ( maximumvdwradius )

    # . Loop over grid construction.
    success = False
    for itry in range ( _NTries ):

        # . Basic grid.
        data    = []
        ntotal  =  1
        for ( i, e ) in enumerate ( extents ):
            npoint = int ( math.floor ( e / step ) )
            ntotal *= npoint
            lower  = origin[i] + 0.5 * ( e - float ( npoint ) * step )
            data.append ( { "bins" : npoint, "binSize" : step, "lower" : lower } )

        # . Unoccupied grid points.
        if ntotal >= nions:
            grid       = RegularGrid.FromDimensionData ( data )
            unoccupied = system.coordinates3.IdentifyUnoccupiedGridPoints ( grid, radii )
            if len ( unoccupied ) >= nions:
                grid.Summary ( )
                _PlaceIonsAtGridPoints ( ions, grid, unoccupied )
                success = True
                break

        # . Try again.
        step *= _ReductionFactor

    # . Check for a grid.
    if not success: raise ValueError ( "Unable to generate a grid with sufficient points to place ions." )

    # . Printing.
    if _DoPrint:
        ArrayPrint   ( origin           , itemFormat = "{:.3f}", title = "Box Origin"      )
        ArrayPrint   ( extents          , itemFormat = "{:.3f}", title = "Box Extents"     )
        ArrayPrint2D ( ions.coordinates3, itemFormat = "{:.3f}", title = "Ion Coordinates" )

    # . Merge.
    newSystem = None
    if _DoMerge:

        # . Merge system and ions.
        newSystem       = MergeByAtom ( [ system, ions ] )
        newSystem.label = system.label
        newSystem.Summary ( )

        # . Check total charge.
        qtotal = sum ( [ atom.formalCharge for atom in newSystem.atoms ] )
        if qtotal != 0: log.Paragraph ( "Warning> The new system has a non-zero total charge - {:.1f}.".format ( float ( qtotal ) ) )

    return newSystem

#===================================================================================================================================
# . Build a solvent box of any crystal class.
#===================================================================================================================================
def BuildSolventBox ( crystalSystem, symmetryParameters, molecule, density, log = logFile, randomNumberGenerator = None ):
    """Build a solvent box.

    No attempt is made to avoid overlapping molecules as it is assumed this will be corrected in a subsequent optimization process.
    """

    # . Initialization.
    solvent = None

    # . Find the number of solvent molecules that fit in the box.
    nmolecules = SolventMoleculeNumber ( molecule, symmetryParameters, density )

    # . Check the number of molecules.
    if nmolecules > 0:

        # . Create the solvent system.
        molecule.coordinates3.ToPrincipalAxes ( )
        solvent = MergeRepeatByAtom ( molecule, nmolecules )
        solvent.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystalSystem )
        solvent.symmetryParameters = symmetryParameters

        # . Check the random number generator.
        if randomNumberGenerator is None: randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )

        # . Do random rotations and translations.
        natoms      = len ( molecule.atoms )
        rotation    = Matrix33.Null ( )
        selection   = Selection.FromIterable ( range ( natoms ) )
        translation = Vector3.Null ( )
        for i in range ( nmolecules ):
            rotation.RandomRotation ( randomNumberGenerator )
            solvent.coordinates3.Rotate  ( rotation,    selection = selection )
            for j in range ( 3 ): translation[j] = randomNumberGenerator.NextReal ( )
            symmetryParameters.H.ApplyTo ( translation )
            solvent.coordinates3.Translate ( translation, selection = selection )
            selection = Selection.FromIterable ( [ ( i + natoms ) for i in selection ] )

        # . Do some printing.
        if LogFileActive ( log ):
            items = [ ( "Molecules"         , "{:d}"  .format ( nmolecules                        ) ) ,
                      ( "Density (kg m^-3)" , "{:.3f}".format ( SystemDensity ( solvent )         ) ) ,
                      ( "Volume"            , "{:.3f}".format ( solvent.symmetryParameters.volume ) ) ]
            log.SummaryOfItems ( items, title = "Solvent Box Summary" )

    # . Return the system.
    return solvent

#===================================================================================================================================
# . Build a cubic solvent box.
# . This function needs a better way of determining the molecule size.
#===================================================================================================================================
def BuildCubicSolventBox ( molecule                        ,
                           nmolecules                      ,
                           doRandomRotation      = True    ,
                           excludeHydrogens      = True    ,
                           log                   = logFile ,
                           moleculesize          = None    ,
                           randomNumberGenerator = None    ,
                           scaleSafety           = 1.1     ):
    """Build a cubic solvent box."""

    # . Get the number of molecules in each direction.
    nlinear = int ( math.ceil ( math.pow ( float ( nmolecules ), 1.0 / 3.0 ) ) )

    # . Get the indices of the occupied sites.
    sites = sample ( list ( range ( nlinear**3 ) ), nmolecules )
    sites.sort ( )

    # . Get the number of atoms in the molecule and an appropriate selection.
    natoms    = len ( molecule.atoms )
    selection = Selection.FromIterable ( range ( natoms ) )

    # . Get a copy of molecule's coordinates (reorientated).
    coordinates3 = Clone ( molecule.coordinates3 )
    coordinates3.ToPrincipalAxes ( )

    # . Get the molecule size depending upon the input options.
    # . A molecule size has been specified.
    if moleculesize is not None:
        size = moleculesize
    # . Determine the size of the molecule as the diagonal distance across its enclosing orthorhombic box.
    else:
        radii = Array.FromIterable ( [ atom.vdwRadius for atom in molecule.atoms ] )
        if excludeHydrogens:
            for ( i, atom ) in enumerate ( molecule.atoms ):
                if atom.atomicNumber == 1: radii[i] = 0.0
        ( origin, extents ) = coordinates3.EnclosingOrthorhombicBox ( radii = radii )
        size                = extents.Norm2 ( ) * scaleSafety

    # . Create the new system - temporarily resetting the coordinates.
    temporary3            = molecule.coordinates3
    molecule.coordinates3 = coordinates3
    solvent               = MergeByAtom ( nmolecules * [ molecule ] )
    molecule.coordinates3 = temporary3

    # . Set the system symmetry.
    solvent.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( CrystalSystemCubic ( ) )
    solvent.symmetryParameters = solvent.symmetry.MakeSymmetryParameters ( a = size * float ( nlinear ) )

    # . Set up for random rotations.
    if doRandomRotation:
        if randomNumberGenerator is None: randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )
        rotation = Matrix33.Null ( )

    # . Loop over the box sites.
    origin      = 0.5 * float ( 1 - nlinear ) * size
    n           = 0
    translation = Vector3.Null ( )
    for i in range ( nlinear ):
        translation[0] = origin + size * float ( i )
        for j in range ( nlinear ):
            translation[1] = origin + size * float ( j )
            for k in range ( nlinear ):
                if len ( sites ) == 0: break
                translation[2] = origin + size * float ( k )
                # . Check for an occupied site.
                if sites[0] == n:
                    sites.pop ( 0 )
                    # . Randomly rotate the coordinates.
                    if doRandomRotation:
                        rotation.RandomRotation ( randomNumberGenerator )
                        solvent.coordinates3.Rotate ( rotation, selection = selection )
                    # . Translate the coordinates.
                    solvent.coordinates3.Translate ( translation, selection = selection )
                    # . Increment the selection for the next molecule.
                    selection = Selection.FromIterable ( [ ( i + natoms ) for i in selection ] )
                n += 1

    # . Do some printing.
    if LogFileActive ( log ):
        items = [ ( "Molecules"         , "{:d}"  .format ( nmolecules                   ) ) ,
                  ( "Density (kg m^-3)" , "{:.3f}".format ( SystemDensity ( solvent )    ) ) ,
                  ( "Box Side"          , "{:.3f}".format ( solvent.symmetryParameters.a ) ) ,
                  ( "Molecule Size"     , "{:.3f}".format ( size                         ) ) ]
        log.SummaryOfItems ( items, title = "Cubic Solvent Box Summary" )

    # . Return the cubic system.
    return solvent

#===================================================================================================================================
# . Calculate solvation parameters for a system.
#===================================================================================================================================
# . Charge equations (only two):
#
# . Charge constraint: Sum_i n_i * q_i + q_T = 0.
# . Ionic strength = Sum_i c_i * q_i^2 / 2.
# . c_i * N_av * 10^-27 * V = n_i.
#
# . For a monovalent anion and a monovalent cation, ionic strength is equivalent to molarity.

def CalculateSolvationParameters ( system, bufferValue = 0.0, geometry = "Orthorhombic", ionicStrength = 0.0, reorientSolute = False ):
    """Check the solvation parameters for a system."""

    # . Get the system charge.
    charge         = sum ( [ atom.formalCharge for atom in system.atoms ] )
    absoluteCharge = abs ( charge )

    # . Get the system extents.
    if reorientSolute:
        masses       = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
        coordinates3 = Clone ( system.coordinates3 )
        coordinates3.ToPrincipalAxes ( weights = masses )
    else:
        coordinates3 = system.coordinates3
    ( origin, extents ) = coordinates3.EnclosingOrthorhombicBox ( )
    extents.Increment ( bufferValue )

    # . Get solvated system dimensions.
    r = x = y = z = 0.0
    if   geometry == "Cubic":
        x = y = z = max ( extents ) ;
    elif geometry == "Orthorhombic":
        x = extents[0] ; y = extents[1] ; z = extents[2]
    elif geometry == "Spherical":
        r = max ( extents ) / 2.0
    elif geometry == "Tetragonal (X=Y)":
        x = y = max ( extents[0], extents[1] ) ; z = extents[2]
    elif geometry == "Tetragonal (Y=Z)":
        y = z = max ( extents[1], extents[2] ) ; x = extents[0]
    elif geometry == "Tetragonal (Z=X)":
        z = x = max ( extents[2], extents[0] ) ; y = extents[1]

    # . Get the minimum values.
    minimumR = max ( 0.0, r - bufferValue )
    minimumX = max ( 0.0, x - bufferValue )
    minimumY = max ( 0.0, y - bufferValue )
    minimumZ = max ( 0.0, z - bufferValue )

    # . Get the volume.
    if geometry == "Spherical": v = 4.0 * math.pi * r**3 / 3.0
    else:                       v = x * y * z

    # . Determine the number of anions and cations given the ionic strength.
    minimumNumberAnions  = 0
    minimumNumberCations = 0
    numberAnions         = 0
    numberCations        = 0
    if ionicStrength > 0.0:
        numberAnions = numberCations = int ( round ( ionicStrength * v * ( Constants.Avogadro_Number * 1.0e-27 ) ) )
    if   charge > 0:
        minimumNumberAnions  = absoluteCharge
        numberAnions        += absoluteCharge
    elif charge < 0:
        minimumNumberCations = absoluteCharge
        numberCations       += absoluteCharge

    # . Finish up.
    solvationParameters = { "minimumNumberAnions"  : minimumNumberAnions  ,
                            "minimumNumberCations" : minimumNumberCations ,
                            "minimumRadius"        : minimumR             ,
                            "minimumX"             : minimumX             ,
                            "minimumY"             : minimumY             ,
                            "minimumZ"             : minimumZ             ,
                            "numberAnions"         : numberAnions         ,
                            "numberCations"        : numberCations        ,
                            "radius"               : r                    ,
                            "x"                    : x                    ,
                            "y"                    : y                    ,
                            "z"                    : z                    }
    return solvationParameters

#===================================================================================================================================
# . Determine solvation parameters for a system.
#===================================================================================================================================
# . Parameters.
_BufferValues    = ( 0.0, 2.5, 5.0, 7.5, 10.0 )    # . Angstroms.
_ChargeTolerance = 1.0e-6
_Molarities      = ( 0.0, 0.01, 0.1, 0.154, 1.0 )  # . M (moles/litre).

# . Solvent data.
_WaterDensity = 1000.0   # . kg/m^3.
_WaterMass    = 18.01528 # . amu.

# . Function.
def DetermineSolvationParameters ( system                         ,
                                   buffervalues   = _BufferValues ,
                                   log            = logFile       ,
                                   molarities     = _Molarities   ,
                                   solventdensity = _WaterDensity ,
                                   solventmass    = _WaterMass    ):
    """Determine some solvation parameters for a system."""

    # . Get the system's extents (before and after reorientation).
    masses                = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
    ( origin0, extents0 ) = system.coordinates3.EnclosingOrthorhombicBox ( )
    coordinates3          = Clone ( system.coordinates3 )
    coordinates3.ToPrincipalAxes ( weights = masses )
    ( origin1, extents1 ) = coordinates3.EnclosingOrthorhombicBox ( )

    # . Get the system's charge and mass.
    charge   = sum ( [ atom.formalCharge for atom in system.atoms ] )
    mmcharge = sum ( system.AtomicCharges ( )                       )
    mass     = sum ( masses )

    # . Convert the mass to milligrams.
    mass *= Units.Mass_AMU_To_Kg * 1000.0

    # . Calculate the solvent volume in A^3.
    solventvolume = solventmass * ( Units.Mass_AMU_To_Kg * 1.0e+30 ) / solventdensity

    # . Basic printing.
    if LogFileActive ( log ):
        if math.fabs ( float ( charge ) - mmcharge ) > _ChargeTolerance: log.Paragraph ( "Total formal and MM charges differ. MM charge = {:.3f}.".format ( mmcharge ) )
        log.Paragraph ( "Total formal charge of system = {:.1f}.".format ( float ( charge ) ) )

    # . Detailed printing.
    # . Initialization.
    nbuffers = len ( buffervalues )
    if nbuffers > 0:

        # . Loop over extents.
        for ( extents, tag ) in ( ( extents0, "Original" ), ( extents1, "Reoriented" ) ):

            # . Get the extents.
            x = extents[0]
            y = extents[1]
            z = extents[2]

            # . Get the buffer data.
            data = [ [] for i in range ( 4 ) ]
            for buffer in buffervalues:
                data[0].append ( x + buffer )
                data[1].append ( y + buffer )
                data[2].append ( z + buffer )
                data[3].append ( ( x + buffer ) * ( y + buffer ) * ( z + buffer ) )

            # . Table header.
            if LogFileActive ( log ):
                table   = log.GetTable ( columns = [ 30 ] + nbuffers * [ 20 ] )
                table.Start ( )
                table.Title ( "Solvation Information for " + tag + " Coordinates in an Orthorhombic Box" )
                table.Heading ( "Property" )
                table.Heading ( "Buffers (Angstroms)", columnSpan = nbuffers )
                table.Heading ( "" )
                for buffer in buffervalues: table.Heading ( "{:.2f}".format ( buffer ) )

                # . Box data.
                for ( values, property ) in zip ( data, ( "X (A)", "Y", "Z", "Volume (A^3)" ) ):
                    table.Entry ( "Box " + property, align = Align.Left )
                    for value in values: table.Entry ( "{:.2f}".format ( value ) )

                # . Number of solvent molecules in solvent-only box.
                table.Entry ( "Molecules in Pure Solvent Box", align = Align.Left )
                for v in data[3]: table.Entry ( "{:d}".format ( int ( round ( v / solventvolume ) ) ) )

                # . System concentration.
                table.Entry ( "System Concentration (mg/l)", align = Align.Left )
                for v in data[3]: table.Entry ( "{:.2f}".format ( mass / ( 1.0e-27 * v ) ) )

                # . Number of ions in box for a given molarity - with the constraint that the total charge (system + ions) is zero.
                for molarity in molarities:
                    table.Entry ( "Ions {:.3f} M".format ( molarity ), align = Align.Left )
                    for v in data[3]:
                        nnegative = npositive = int ( round ( molarity * v * ( Constants.Avogadro_Number * 1.0e-27 ) ) )
                        if   charge > 0: nnegative += abs ( charge )
                        elif charge < 0: npositive += abs ( charge )
                        table.Entry ( "+{:d}/-{:d}".format ( npositive, nnegative ) )

                # . Finish up.
                table.Stop ( )

    # . Spherical systems.
    # . Radii.
    radius0 = max ( extents0 ) / 2.0
    radius1 = max ( extents1 ) / 2.0
    if nbuffers > 0:

        # . Loop over extents.
        for ( radius, tag ) in ( ( radius0, "Original" ), ( radius1, "Reoriented" ) ):

            # . Get the buffer data.
            data = [ [] for i in range ( 2 ) ]
            for buffer in buffervalues:
                r = radius + buffer
                data[0].append ( r )
                data[1].append ( 4.0 * math.pi * r**3 / 3.0 )

            # . Table header.
            if LogFileActive ( log ):
                table   = log.GetTable ( columns = [ 35 ] + nbuffers * [ 20 ] )
                table.Start ( )
                table.Title ( "Solvation Information for " + tag + " Coordinates in a Sphere" )
                table.Heading ( "Property" )
                table.Heading ( "Buffers (Angstroms)", columnSpan = nbuffers )
                table.Heading ( "" )
                for buffer in buffervalues: table.Heading ( "{:.2f}".format ( buffer ) )

                # . Box data.
                for ( values, property ) in zip ( data, ( "Radius (A)", "Volume (A^3)" ) ):
                    table.Entry ( "Sphere " + property, align = Align.Left )
                    for value in values: table.Entry ( "{:.2f}".format ( value ) )

                # . Number of solvent molecules in solvent-only box.
                table.Entry ( "Molecules in Pure Solvent Sphere", align = Align.Left )
                for v in data[1]: table.Entry ( "{:d}".format ( int ( round ( v / solventvolume ) ) ) )

                # . System concentration.
                table.Entry ( "System Concentration (mg/l)", align = Align.Left )
                for v in data[1]: table.Entry ( "{:.2f}".format ( mass / ( 1.0e-27 * v ) ) )

                # . Number of ions in box for a given molarity - with the constraint that the total charge (system + ions) is zero.
                for molarity in molarities:
                    table.Entry ( "Ions {:.3f} M".format ( molarity ), align = Align.Left )
                    for v in data[1]:
                        nnegative = npositive = int ( round ( molarity * v * ( Constants.Avogadro_Number * 1.0e-27 ) ) )
                        if   charge > 0: nnegative += abs ( charge )
                        elif charge < 0: npositive += abs ( charge )
                        table.Entry ( "+{:d}/-{:d}".format ( npositive, nnegative ) )

                # . Finish up.
                table.Stop ( )

#===================================================================================================================================
# . Calculate the dimension of a cubic box given a number of molecules and a density.
#===================================================================================================================================
def SolventCubicBoxDimensions ( molecule, nmolecules, density ):
    """Calculate the dimension of a cubic box given a number of molecules and a density.

    |density| is in kg m^-3.
    """
    #  . Get the volumes of the molecule and the box.
    vmolecule = sum ( [ atom.mass for atom in molecule.atoms ] ) * ( Units.Mass_AMU_To_Kg * 1.0e+30 ) / density
    vbox      = float ( nmolecules ) * vmolecule
    # . Return the box size.
    return vbox**( 1.0 / 3.0 )

#===================================================================================================================================
# . Calculate the number of molecules required to fill a box of a given size and a given density.
#===================================================================================================================================
def SolventMoleculeNumber ( molecule, symmetryParameters, density ):
    """Calculate the number of molecules required to fill a box of a given size.

    |density| is in kg m^-3.
    """
    # . Get the volumes of the molecule and the box.
    vmolecule = sum ( [ atom.mass for atom in molecule.atoms ] ) * ( Units.Mass_AMU_To_Kg * 1.0e+30 ) / density
    vbox      = symmetryParameters.volume
    # . Return the molecule number.
#    return int ( vbox / vmolecule )
    return round ( vbox / vmolecule )

#===================================================================================================================================
# . Solvate a system by superposition.
# . This function needs a better way of determining molecule size.
#===================================================================================================================================
def SolvateSystemBySuperposition ( solute, solvent, log = logFile, excludeHydrogens = True, printPairs = False, reorientSolute = True, scaleRadii = 0.75 ):
    """Solvate a system by superposition."""

    # . Check that the solvent has a full connectivity.
    if len ( solvent.connectivity.bonds ) <= 0: raise ValueError ( "The solvent must have a full connectivity." )

    # . Create a sequence.
    if solute.sequence is not None:
        entityLabel = DetermineUniqueEntityLabel ( solute, label = "W" )
        CreateHomogeneousIsolateSequence ( solvent, entityLabel = entityLabel )

    # . Merge the two systems.
    solution = MergeByAtom ( solute, solvent )

    # . Get the atom masses and selections corresponding to the solute and solvent atoms.
    masses       = Array.FromIterable ( [ atom.mass for atom in solution.atoms ] )
    nsolute      = len ( solute.atoms  )
    nsolvent     = len ( solvent.atoms )
    soluteatoms  = Selection.FromIterable ( range ( nsolute                     ) )
    solventatoms = Selection.FromIterable ( range ( nsolute, nsolute + nsolvent ) )

    # . Reorient the solute coordinates if required.
    if reorientSolute: solution.coordinates3.ToPrincipalAxes ( selection = soluteatoms, weights = masses )

    # . Check the extents of the solute and solvent - to ensure the box is big enough!
    ( origin1, extents1 ) = solution.coordinates3.EnclosingOrthorhombicBox ( selection = soluteatoms  )
    ( origin2, extents2 ) = solution.coordinates3.EnclosingOrthorhombicBox ( selection = solventatoms )
    extents2.Add ( extents1, scale = -1.0 )
    if min ( extents2 ) < 0.0: raise ValueError ( "It is likely that the solvent box is too small for the solute." )

    # . Move the center of the solute atoms to that of the solvent atoms.
    translation = Clone ( origin2 )
    translation.Add ( extents2, scale =  0.5 )
    translation.Add ( origin1 , scale = -1.0 )
    solution.coordinates3.Translate ( translation, selection = soluteatoms )

    # . Reduce the selections if hydrogens are to be excluded.
    if excludeHydrogens:
        indices = []
        for i in soluteatoms:
            if solution.atoms[i].atomicNumber != 1: indices.append ( i )
        soluteatoms = Selection.FromIterable ( indices )
        indices = []
        for i in solventatoms:
            if solution.atoms[i].atomicNumber != 1: indices.append ( i )
        solventatoms = Selection.FromIterable ( indices )

    # . Get the radii.
    radii = Array.FromIterable ( [ atom.vdwRadius for atom in solution.atoms ] )
    radii.Scale ( scaleRadii )

    # . Get the pairlist.
    pairlist = CrossPairList_FromSingleCoordinates3 ( solution.coordinates3, selection1 = solventatoms, selection2 = soluteatoms, radii = radii )
    npairs   = len ( pairlist )

    # . Find the solvent atoms in the pairlist (which are the first indices in the list).
    indices = set ( )
    for ( i, j ) in pairlist: indices.add ( i )

    # . Now find the indices of all the atoms in the corresponding solvent isolates.
    if len ( indices ) > 0:
        toRemove  = solution.connectivity.isolateIndices.UnionOfItems ( Selection.FromIterable ( indices ) )
        nToRemove = len ( toRemove )
        # . Remove unwanted atoms.
        reduced = PruneByAtom ( solution, toRemove.Complement ( upperBound = nsolute + nsolvent ) )
        if solution.sequence is not None: RenumberEntityComponents ( reduced, entityLabels = [ entityLabel ] )
    else:
        nToRemove = 0
        reduced   = solution

    # . Set the symmetry of solution by copying that of solvent.
    if ( solvent.symmetry is not None ) and ( solvent.symmetryParameters is not None ):
        reduced.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( solvent.symmetry.crystalSystem, transformations = solvent.symmetry.transformations )
        reduced.symmetryParameters = Clone ( solvent.symmetryParameters )

    # . Do some printing.
    if LogFileActive ( log ):

        # . Distances.
        if printPairs:
            table   = log.GetTable ( columns = [ 10, 10, 20, 20 ] )
            table.Start ( )
            table.Title ( "Solute/Solvent Pairs" )
            table.Heading ( "Solvent"  )
            table.Heading ( "Solute"   )
            table.Heading ( "Distance" )
            table.Heading ( "CutOff"   )
            for ( i, j ) in pairlist:
                table.Entry ( PeriodicTable.Symbol ( solution.atoms[i].atomicNumber, index = i ) )
                table.Entry ( PeriodicTable.Symbol ( solution.atoms[j].atomicNumber, index = j ) )
                table.Entry ( "{:.3f}".format ( solution.coordinates3.Distance ( i, j ) ) )
                table.Entry ( "{:.3f}".format ( radii[i] + radii[j]                     ) )
            table.Stop ( )

        # . Summary.
        items = [ ( "Solute Atoms"         , "{:d}".format ( nsolute              ) ) ,
                  ( "Solvent Atoms"        , "{:d}".format ( nsolvent - nToRemove ) ) ,
                  ( "Atoms Removed"        , "{:d}".format ( nToRemove            ) ) ,
                  ( "Solute/Solvent Pairs" , "{:d}".format ( npairs               ) ) ]
        log.SummaryOfItems ( items, title = "Solvation By Superposition Summary" )

    # . Finish up.
    return reduced

#===================================================================================================================================
# . Calculate the density of a system.
#===================================================================================================================================
def SystemDensity ( system ):
    """Calculate the density of a system in kg m^-3."""
    # . Check that the system has a volume.
    try:
        volume = system.symmetryParameters.volume
        mass   = sum ( [ atom.mass for atom in system.atoms ] )
        return ( ( mass / volume ) * ( Units.Mass_AMU_To_Kg * 1.0e+30 ) )
    except:
        return 0.0

#===================================================================================================================================
# . Determine a system's extents.
#===================================================================================================================================
def SystemExtents ( system, log = logFile ):
    """Calculate the extents of a system."""
    extents = None
    if system.coordinates3 is not None:
        ( origin, extents ) = system.coordinates3.EnclosingOrthorhombicBox ( )
        if LogFileActive ( log ):
            table = log.GetTable ( columns = [ 5, 20 ] )
            table.Start ( )
            table.Title ( "System Extents" )
            for ( label, extent ) in zip ( [ "x", "y", "z" ], extents ):
                table.Entry ( label, align = Align.Left )
                table.Entry ( "{:.3f}".format ( extent ) )
            table.Stop ( )
    return extents

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

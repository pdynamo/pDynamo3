"""Helper functions for performing normal mode and associated thermodynamics calculations."""

import math

from  enum                      import Enum
from  pCore                     import Align                              , \
                                       AttributableObject                 , \
                                       Clone                              , \
                                       logFile                            , \
                                       LogFileActive
from  pMolecule                 import SystemGeometryObjectiveFunction
from  pScientific               import Constants                          , \
                                       Units
from  pScientific.Arrays        import Array
from  pScientific.Geometry3     import Coordinates3                       , \
                                       Vector3
from  pScientific.LinearAlgebra import EigenPairs
from  pScientific.Symmetry      import Find3DGraphPointGroup              , \
                                       IdentifyIrreducibleRepresentations , \
                                       PrintIrreducibleRepresentations
from .TrajectoryUtilities       import CovarianceMatrix

#
# . RMS displacements for an atom are:
#
#   D(t) = sum_k a_k V_k cos ( w_k t + p_k )
#
#   where a_k = Sqrt ( 2 * k_B * T ) / w_k
#         p_k = arbitrary phase
#         V_k = kth normal mode vector
#         w_k = kth frequency
#
#   < D(t) > = 0 as < cos p_k > = < sin p_k > = 0
#
#   < D(t) D(t') > = k_B * T * sum_k V_k^2 / w_k^2 as < cos^2 p_k > = < sin^2 p_k > = 1/2 and < cos p_k sin p_k > = 0.
#
#   The sum_k runs over all modes and the total RMS displacement is taken by picking out the appropriate components for the atom.
#

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class ModifyOption ( Enum ):
    """Normal mode frequency modify options."""
    Null    = 1
    Project = 2
    Raise   = 3

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Frequency output options.
_FrequencyColumns =  8
_FrequencyFormat  = "{:10.3f}"
_FrequencyWidths  = 12

# . Hessian modify options.
_ModifyOption    = ModifyOption.Null
_RaiseEigenvalue = 50000.0

# . Mode output options.
_ModeColumns        = 6
_ModeFormat         = "{:11.5f}"
_ModeWidthComponent =  3
_ModeWidthElement   = 12
_ModeWidthIndex     =  6
_ModeWidthName      =  4

# . Mode trajectory options.
_LowFrequency       = 0.1

# . The conversion factor from cm^-1 to ps^-1.
_To_HZ = Constants.Speed_Of_Light / 1.0e+10

# . The conversion factor from K to kJ mol^-1.
_To_kJMol = Constants.Avogadro_Number * Constants.Boltzmann * 1.0e-03

# . The conversion factor from internal units to cm^-1.
_To_Wavenumbers = 1.0e+11 / ( 2.0 * math.pi * Constants.Speed_Of_Light )

# . The conversion factor from wavenumbers to Joules.
_Wavenumbers_To_Joules = 1.0e+12 * Constants.Planck * _To_HZ

# . Thermodynamic options.
_ExponentialUnderflow    = 75.0
_RotationLinearTolerance = 1.0e-3
_ThermodynamicProperties = ( "Constant Pressure Heat Capacity" ,
                             "Constant Volume Heat Capacity"   ,
                             "Enthalpy"                        ,
                             "Entropy"                         ,
                             "Gibbs Free Energy"               ,
                             "Helmholtz Free Energy"           ,
                             "Internal Energy"                 ,
                             "Log Partition Function"          )

# . Tolerances.
_QHTolerance = 1.0e-05

#===================================================================================================================================
# . Utility functions.
#===================================================================================================================================
# . Useful for hessians coming from elsewhere (e.g. ORCA).
def MassWeightHessian ( hessian, masses, freeAtoms = None ):
    """Mass-weight a hessian."""
    if freeAtoms is None: freeAtoms = range ( len ( masses ) )
    if ( 3 * len ( freeAtoms ) ) == hessian.rows:
        weights = []
        for i in freeAtoms:
            w = 1.0 / math.sqrt ( masses[i] )
            weights.extend ( 3 * [ w ] )
        weightedHessian = Clone ( hessian )
        for i in range ( hessian.rows ):
            wI = weights[i]
            for j in range ( i + 1 ):
                weightedHessian[i,j] *= ( wI * weights[j] )
        return weightedHessian
    else:
        raise ValueError ( "Incompatible hessian and mass arrays." )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NormalModeState ( AttributableObject ):
    """A class to contain the results of a normal mode calculation."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "covariance"  : None  ,
                             "dimension"   :    0  ,
                             "freeAtoms"   : None  ,
                             "frequencies" : None  ,
                             "hessian"     : None  ,
                             "modes"       : None  ,
                             "nRTModes"    :    0  ,
                             "weights"     : None  } )
        
#===================================================================================================================================
# . Compute infrared intensities.
#===================================================================================================================================
# . Intensity parameters.
_IntensityColumns       =  8
_IntensityFormat        = "{:10.3f}"
_IntensityScalingFactor = 42.2561 # 1 Debye2-angstrom-2-amu-1 (IR intensity unit) = 42.2561 km-mol-1
_IntensityWidths        = 12

def NormalModes_InfraredIntensities ( system, center = None, delta = 1.0e-4, log = logFile, title = "Infrared Intensities" ):
    """Compute IR intensities.

    A previous normal mode calculation for the system must have been done.
    """
    # . Get the state.
    state = system.scratch.Get ( "nmState", None )
    if state is None: raise ValueError ( "Normal modes not found for system." )
    # . Get state-related information.
    if state.freeAtoms == None: freeAtoms = range ( len ( system.atoms ) )
    else:                       freeAtoms = state.freeAtoms
    frequencies  = state.frequencies
    modes        = state.modes
    n            = state.dimension
    # . Initialization.
    coordinates3       = Clone ( system.coordinates3 )
    dipoleDerivatives  = Array.WithExtents ( n, 3 ) ; dipoleDerivatives.Set  ( 0.0 )
    dipoleDerivativesQ = Array.WithExtents ( n, 3 ) ; dipoleDerivativesQ.Set ( 0.0 )
    intensities        = Array.WithExtent  ( n    ) ; intensities.Set        ( 0.0 )
    # . Determine the dipole derivatives by finite differences.
    i = 0
    for f in freeAtoms:
        for c in range ( 3 ):
            r = system.coordinates3[f,c]
            system.coordinates3[f,c] = r - delta
            system.Energy ( log = None )
            dipoleM = system.DipoleMoment ( center = center )
            system.coordinates3[f,c] = r + delta
            system.Energy ( log = None )
            dipoleP = system.DipoleMoment ( center = center )
            system.coordinates3[f,c] = r
            dipoleP.Add ( dipoleM, scale = -1.0 )
            dipoleP.Scale ( 1.0 / ( 2.0 * delta ) )
            for d in range ( 3 ): dipoleDerivatives[i,d] = dipoleP[d]
            i += 1
    # . Determine the transformed derivatives.
    dipoleDerivativesQ.MatrixMultiply ( modes, dipoleDerivatives, xTranspose = True )
    # . Determine the intensities.
    for i in range ( n ):
        d = 0.0
        for c in range ( 3 ): d += ( dipoleDerivativesQ[i,c]**2 )
        intensities[i] = d
    intensities.Scale ( _IntensityScalingFactor )
    # . Do some printing.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = _IntensityColumns * [ _IntensityWidths ] )
        table.Start ( )
        table.Title ( title )
        for f in intensities: table.Entry ( _IntensityFormat.format ( f ) )
        table.Stop ( )
    # . Finish up.
    state.intensities = intensities
    return state

#===================================================================================================================================
# . Irreducible representations.
#===================================================================================================================================
# . Parameters.
_NormalModeDegeneracyTolerance = 5.0 # . cm^-1.

def _CharacterFunction ( modeIndices, modes, weights ):
    """Return a function that can compute characters for a set of modes under a symmetry operation."""
    def GetItemCharacters ( operation ):
        """Get item characters."""
        if operation.order == 1:
            characters = [ 1.0 for i in range ( len ( modeIndices ) ) ]
        else:
            characters = []
            vector1    = Vector3.Uninitialized ( )
            vector2    = Vector3.Uninitialized ( )
            for modeIndex in modeIndices:
                character = 0.0
                for u in range ( len ( weights ) // 3 ):
                    v = operation.mapping[u]
                    # . Unweight modes to obtain normalized vectors.
                    for c in range ( 3 ): vector1[c] = modes[3*u+c,modeIndex] * weights[3*u+c]
                    for c in range ( 3 ): vector2[c] = modes[3*v+c,modeIndex] * weights[3*v+c]
                    operation.ApplyTo ( vector1 )
                    character += vector1.Dot ( vector2 )
                characters.append ( character )
        return characters
    return GetItemCharacters

def NormalModes_IrreducibleRepresentations ( system, degeneracyTolerance = _NormalModeDegeneracyTolerance ,
                                                     log                 = logFile                        ,
                                                     maximumIRs          = 6                              ,
                                                     modeIndices         = None                           ,
                                                     results             = None                           ):
    """Find the irreducible representations of some of the modes."""
    # . Get the state.
    state = system.scratch.Get ( "nmState", None )
    if state is None: raise ValueError ( "Normal modes not found for system." )
    # . Get results if it does not exist.
    if results is None:
        masses  = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
        numbers = [ atom.atomicNumber for atom in system.atoms ]
        results = Find3DGraphPointGroup ( numbers                                ,
                                          system.coordinates3                    ,
                                          doCharacterSymmetryOperations = True   ,
                                          log                           = log    ,
                                          weights                       = masses )
    # . Mode data.
    if modeIndices is None: modeIndices = range ( 3 * len ( system.atoms ) )
    itemValues = [ state.frequencies[i] for i in modeIndices ]
    # . Get the representations and convert to lower case (as is standard for normal modes).
    GetItemCharacters   = _CharacterFunction ( modeIndices, state.modes, state.weights )
    ( iRs, characters ) = IdentifyIrreducibleRepresentations ( results                 ,
                                                               itemValues              ,
                                                               GetItemCharacters       ,
                                                               degeneracyTolerance     ,
                                                               maximumIRs = maximumIRs ) # . For rotational and translation modes.
    iRs = [ iR.lower ( ) for iR in iRs ]
    PrintIrreducibleRepresentations ( results                  ,
                                      itemValues               ,
                                      iRs                      ,
                                      characters               ,
                                      itemName   = "Mode"      ,
                                      itemFormat = "{:.1f}"    ,
                                      log        = logFile     ,
                                      valueName  = "Frequency" )
    return iRs

#===================================================================================================================================
# . Normal modes.
#===================================================================================================================================
# . Note that the input hessian here, if provided, should already be mass-weighted.
def NormalModes_SystemGeometry ( system, hessian = None, log = logFile, modify = _ModifyOption, title = "Harmonic Frequencies (cm^(-1))" ):
    """Determine the normal modes for a system."""
    # . Get the Hessian with mass-weighting.
    of = SystemGeometryObjectiveFunction.FromSystem ( system )
    of.DefineWeights ( )
    n  = of.numberOfVariables
    if hessian is None:
        x = Array.WithExtent ( n )
        x.Set ( 0.0 )
        of.VariablesGet ( x )
        hessian = of.NumericalHessian ( x )
    # . Get the mass-weighted rotation-translation vectors and count their number.
    of.RemoveRotationTranslation ( )
    if of.linearScalars is None: nRTModes = 0
    else:                        nRTModes = len ( of.linearScalars )
    # . Modify the Hessian.
    if   modify == ModifyOption.Project : hessian = hessian.ProjectOutVectors ( of.linearVectors                   )
    elif modify == ModifyOption.Raise   :           hessian.Raise             ( of.linearVectors, _RaiseEigenvalue )
    # . Diagonalization.
    # . Maybe should save hessian here as it is destroyed by the diagonalization.
    eigenValues  = Array.WithExtent  ( n )    ; eigenValues.Set  ( 0.0 )
    eigenVectors = Array.WithExtents ( n, n ) ; eigenVectors.Set ( 0.0 )
    EigenPairs ( hessian, eigenValues, eigenVectors )
    # . Convert eigenvalues to frequencies.
    for ( i, e ) in enumerate ( eigenValues ):
        f = math.sqrt ( math.fabs ( e ) ) * _To_Wavenumbers
        if e < 0.0: f *= -1.0
        eigenValues[i] = f
    # . Un-mass-weight the modes.
    for r in range ( n ):
#        w = 1.0 / of.variableWeights[r]
#        for c in range ( n ): eigenVectors[r,c] *= w
        eigenVectors[r,:].Scale ( 1.0 / of.variableWeights[r] )
    # . Do some printing.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = _FrequencyColumns * [ _FrequencyWidths ] )
        table.Start ( )
        table.Title ( title )
        for f in eigenValues: table.Entry ( _FrequencyFormat.format ( f ) )
        table.Stop ( )
    # . Save all data.
    state = NormalModeState.WithOptions ( dimension = n, freeAtoms = of.freeAtoms, frequencies = eigenValues, modes = eigenVectors, nRTModes = nRTModes, weights = of.variableWeights )
    system.scratch.nmState = state
    # . Finish up.
    return state

#===================================================================================================================================
# . Normal mode printing.
#===================================================================================================================================
def NormalModesPrint_SystemGeometry ( system, log = logFile, modes = None, selection = None, state = None, title = "Normal Mode Eigenvectors" ):
    """Print the normal modes."""
    # . Check for printing.
    if LogFileActive ( log ):
        # . Get the state.
        if state is None: state = system.scratch.nmState
        # . Get state-related information.
        n            = state.dimension
        eigenVectors = state.modes
        frequencies  = state.frequencies
        # . Get the atom selection as the intersection of the input selection and the free atoms.
        # . Input selection.
        if selection == None: atoms = range ( len ( system.atoms ) )
        else:                 atoms = selection
        # . Free atoms.
        if state.freeAtoms == None: freeAtoms = range ( len ( system.atoms ) )
        else:                       freeAtoms = state.freeAtoms
        # . The intersection.
        atomSelection = list ( set ( atoms ).intersection ( set ( freeAtoms ) ) )
        atomSelection.sort ( )
        # . Get the mode selection.
        if modes == None: modeSelection = range ( n )
        else:             modeSelection = modes
        nmodes = len ( modeSelection )
        # . There are atoms and modes.
        if ( len ( atomSelection ) > 0 ) and ( nmodes > 0 ):
            # . Get the free atom positions.
            freeindices = {}
            for ( i, f ) in enumerate ( freeAtoms ): freeindices[f] = i
            # . Get atom columns and names.
            atomcolumns = []
            atomnames   = []
            for atom in atomSelection:
                atomcolumns.append ( 3 * freeindices[atom] )
                atomnames.append ( system.atoms[atom].symbol )
            # . Do the printing.
            for start in range ( 0, nmodes, _ModeColumns ):
                stop    = min ( start + _ModeColumns, nmodes )
                columns = [ _ModeWidthIndex, _ModeWidthName, _ModeWidthComponent ]
                for row in range ( start, stop ): columns.append ( _ModeWidthElement )
                table = log.GetTable ( columns = columns )
                table.Start ( )
                table.Title ( title )
                table.Entry ( "Mode",          align = Align.Left, columnSpan = 3 )
                for i in range ( start, stop ): table.Entry ( "{:d}".format ( modeSelection[i] ) )
                table.Entry ( "Freq. (cm^-1)", align = Align.Left, columnSpan = 3 )
                for i in range ( start, stop ): table.Entry ( _ModeFormat.format ( frequencies[modeSelection[i]]          ) )
                table.Entry ( "Freq. (ps^-1)", align = Align.Left, columnSpan = 3 )
                for i in range ( start, stop ): table.Entry ( _ModeFormat.format ( frequencies[modeSelection[i]] * _To_HZ ) )
                for ( atom, column, name ) in zip ( atomSelection, atomcolumns, atomnames ):
                    table.Entry ( "{:d}".format ( atom ) )
                    table.Entry ( name   )
                    table.Entry ( "x", align = Align.Center )
                    for i in range ( start, stop ): table.Entry ( _ModeFormat.format ( eigenVectors[column  ,modeSelection[i]] ) )
                    table.Entry ( None )
                    table.Entry ( None )
                    table.Entry ( "y", align = Align.Center )
                    for i in range ( start, stop ): table.Entry ( _ModeFormat.format ( eigenVectors[column+1,modeSelection[i]] ) )
                    table.Entry ( None )
                    table.Entry ( None )
                    table.Entry ( "z", align = Align.Center )
                    for i in range ( start, stop ): table.Entry ( _ModeFormat.format ( eigenVectors[column+2,modeSelection[i]] ) )
                table.Stop ( )

#===================================================================================================================================
# . Generate a normal mode trajectory.
#===================================================================================================================================
def NormalModesTrajectory_SystemGeometry ( system, trajectory, cycles = 10, frames = 21, mode = 0, state = None, temperature = 300.0 ):
    """Generate a normal mode trajectory."""
    # . Get the state.
    if state is None: state = system.scratch.nmState
    # . Get state-related information.
    if state.freeAtoms == None: freeAtoms = range ( len ( system.atoms ) )
    else:                       freeAtoms = state.freeAtoms
    frequencies = state.frequencies
    modes       = state.modes
    # . Get the mode frequency.
    omega = math.fabs ( frequencies[mode] )
    # . Calculate the number of frames.
    total = cycles * frames
    # . Check for a calculation.
    if ( mode >= 0 ) and ( mode < state.dimension ) and ( omega > _LowFrequency ) and ( total > 0 ):
        # . Calculate the amplitude (in Angstroms).
        amplitude = math.sqrt ( 2.0 * 1.0e-3 * Constants.Avogadro_Number * Constants.Boltzmann * temperature ) * ( _To_Wavenumbers / omega )
        # . Allocate space for the coordinates and mode.
        coordinates3 = Clone ( system.coordinates3 )
        displacement = Coordinates3.WithExtent ( len ( system.atoms ) )
        displacement.Set ( 0.0 )
        # . Get the displacement.
        for ( i, f ) in enumerate ( freeAtoms ):
            displacement[f,0] = modes[3*i  ,mode]
            displacement[f,1] = modes[3*i+1,mode]
            displacement[f,2] = modes[3*i+2,mode]
        # . Loop over the cycles and frames.
        # . Calculate the displacement prefactor using sine instead of cosine.
        trajectory.WriteHeader ( )
        for c in range ( cycles ):
            for f in range ( frames ):
                factor = amplitude * math.sin ( 2.0 * math.pi * float ( f ) / float ( frames ) )
                coordinates3.CopyTo ( system.coordinates3  )
                system.coordinates3.Add ( displacement, scale = factor )
                trajectory.WriteOwnerData ( )
        # . Finish up.
        trajectory.WriteFooter ( )
        trajectory.Close       ( )

#===================================================================================================================================
# . Quasi-harmonic modes.
#===================================================================================================================================
#
# . For a regular NM analysis, one has H * A = w^2 * M * A or H' * A' = e * A'.
# . Here H' = M^(-1/2) * H * M^(-1/2), A' = M^1/2 * A, and e = w^2. Therefore,
# . w = Sqrt ( 1/e ) and A = M^(-1/2) * A'.
#
# . Here H = T * inverse ( C ) so equations are T * inverse ( C ) * A = w^2 * M * A.
# . Rearranging gives: C * M * A = ( T / w^2 ) * A or C' * A' = e * A' where
# . C' = (M^1/2) * C * (M^1/2), A' = (M^1/2) * A and e = ( T / w^2 ). Therefore,
# . w = Sqrt ( T/e ) and A = M^(-1/2) * A'.
#
def QuasiHarmonic_SystemGeometry ( system, covariance = None, log = logFile, modify = _ModifyOption, temperature = 300.0, title = "Quasi-Harmonic Frequencies (cm^(-1))", trajectories = None ):
    """Determine the quasi-harmonic modes for a system."""
    # . Initialization.
    state = None
    # . Get the covariance matrix.
    if covariance is None:
        covariance = CovarianceMatrix ( trajectories, system, selection = system.freeAtoms )
    # . Proceed with the analysis.
    if covariance is not None:
        # . Mass-weight the covariance matrix.
        # . Weights are square roots of masses.
        of = SystemGeometryObjectiveFunction.FromSystem ( system )
        of.DefineWeights ( )
        n  = of.numberOfVariables
        for i in range ( n ):
            wI = of.variableWeights[i]
            for j in range ( i + 1 ):
                wJ = of.variableWeights[j]
                covariance[i,j] *= ( wI * wJ )
        # . Get the mass-weighted rotation-translation vectors and count their number.
        of.RemoveRotationTranslation ( )
        if of.linearScalars is None: nRTModes = 0
        else:                        nRTModes = len ( of.linearScalars )
        # . Modify the Hessian.
        if   modify == ModifyOption.Project : covariance = covariance.ProjectOutVectors ( of.linearVectors      )
        elif modify == ModifyOption.Raise   :              covariance.Raise             ( of.linearVectors, 0.0 )
        # . Diagonalization.
        eigenValues  = Array.WithExtent  ( n )    ; eigenValues.Set  ( 0.0 )
        eigenVectors = Array.WithExtents ( n, n ) ; eigenVectors.Set ( 0.0 )
        EigenPairs ( covariance, eigenValues, eigenVectors )
        # . Convert eigenValues to frequencies.
        conversionFactor = math.sqrt ( _To_kJMol * temperature ) * _To_Wavenumbers
        numberZero       = 0
        for ( i, e ) in enumerate ( eigenValues ):
            eAbs = math.fabs ( e )
            if eAbs <= _QHTolerance:
                f = 0.0
                numberZero += 1
            else:
                f = math.sqrt ( 1.0 / eAbs ) * conversionFactor
                if e < 0.0: f *= -1.0
            eigenValues[i] = f
        # . Un-mass-weight the modes.
        for r in range ( n ):
#            w = 1.0 / of.variableWeights[r]
#            for c in range ( n ): eigenVectors[r,c] *= w
            eigenVectors[r,:].Scale ( 1.0 / of.variableWeights[r] )
        # . Reverse in place (excluding zero modes).
        temporary = Array.WithExtent ( n )
        for i in range ( ( n - numberZero ) // 2 ):
            # . Indices.
            lower = i + numberZero
            upper = n - i - 1
            # . Eigenvalues.
            e = eigenValues[upper]
            eigenValues[upper] = eigenValues[lower]
            eigenValues[lower] = e
            # . Eigenvectors.
            for j in range ( n ): temporary[j]          = eigenVectors[j,upper]
            for j in range ( n ): eigenVectors[j,upper] = eigenVectors[j,lower]
            for j in range ( n ): eigenVectors[j,lower] = temporary[j]
        # . Do some printing.
        if LogFileActive ( log ):
            table = log.GetTable ( columns = _FrequencyColumns * [ _FrequencyWidths ] )
            table.Start ( )
            table.Title ( title )
            for f in eigenValues: table.Entry ( _FrequencyFormat.format ( f ) )
            table.Stop ( )
        # . Save all data.
        state = NormalModeState.WithOptions ( dimension = n, freeAtoms = of.freeAtoms, frequencies = eigenValues, modes = eigenVectors, nRTModes = nRTModes, weights = of.variableWeights )
        system.scratch.qhState = state
    # . Finish up.
    return state

#===================================================================================================================================
# . Thermodynamical quantities within the RRHO approximation.
#===================================================================================================================================
# . Allowing multiple P and T values would be more efficient.
def ThermodynamicsRRHO_SystemGeometry ( system, pressure = 1.0, state = None, symmetryNumber = 1, temperature = 300.0 ):
    """Determine thermodynamical quantities within the RRHO approximation."""
    # . Get the state.
    if state is None: state = system.scratch.nmState
    # . Calculate pressure, R, RT and volume.
    p  = Units.Pressure_Atmospheres_To_Pascals * pressure
    R  = Constants.Molar_Gas / 1.0e+3
    RT = R * temperature
    v  = Constants.Boltzmann * temperature / p
    # . Get the size of the problem.
    natoms = len ( system.atoms )
    if state.freeAtoms == None: nfree = natoms
    else:                       nfree = len ( state.freeAtoms )
    # . Get the masses of the free atoms.
    if system.freeAtoms is None: indices = range ( natoms )
    else:                        indices = system.freeAtoms
    masses = Array.WithExtent ( natoms )
    masses.Set ( 0.0 )
    for s in indices: masses[s] = system.atoms[s].mass
    totmas = Units.Mass_AMU_To_Kg * sum ( masses )
    # . Electronic contributions.
    # . Assume there is one state with a multiplicity of one.
    electronic = {}
    # . Rotational contributions.
    # . Initialization.
    rotation = {}
    # . Polyatomic system.
    if nfree > 1:
        # . Get the moments of inertia factor (in amu angstroms**2).
        coordinates3  = Clone ( system.coordinates3 )
        coordinates3.ToPrincipalAxes ( weights = masses )
        inertiamatrix = coordinates3.InertiaMatrix ( weights = masses )
        mproduct      = 1.0
        nzeromoments  = 0
        for i in range ( 3 ):
            m = inertiamatrix[i,i]
            if math.fabs ( m ) <= _RotationLinearTolerance: nzeromoments += 1
            else:                                           mproduct     *= m
        # . Calculate some factors.
        factor = ( 8.0 * math.pi * math.pi * Constants.Boltzmann * temperature * Units.Mass_AMU_To_Kg ) / ( Constants.Planck**2 * 1.0e+20 )
        # . Linear molecule.
        if nzeromoments >= 2:
            z = factor * mproduct / float ( symmetryNumber )
            rotation["Log Partition Function"       ] = math.log ( z )
            rotation["Constant Volume Heat Capacity"] = R
            rotation["Entropy"                      ] = R * ( rotation["Log Partition Function"] + 1.0 )
            rotation["Internal Energy"              ] = RT
        # . Non-linear molecule.
        else:
            z = math.sqrt ( math.pi * factor**3 * mproduct ) / float ( symmetryNumber )
            rotation["Log Partition Function"]        = math.log ( z )
            rotation["Constant Volume Heat Capacity"] = 1.5 * R
            rotation["Entropy"                      ] = R * ( rotation["Log Partition Function"] + 1.5 )
            rotation["Internal Energy"              ] = 1.5 * RT
        # . Calculate the common terms.
        rotation["Helmholtz Free Energy"          ] = - RT * rotation["Log Partition Function"]
        rotation["Constant Pressure Heat Capacity"] = rotation["Constant Volume Heat Capacity"]
        rotation["Helmholtz Free Energy"          ] = rotation["Internal Energy"]
        rotation["Gibbs Free Energy"              ] = rotation["Helmholtz Free Energy"] - temperature * rotation["Entropy"]
    # . Translational contributions.
    z = ( math.sqrt ( 2.0 * math.pi * totmas * Constants.Boltzmann * temperature ) / Constants.Planck )**3 * v
    translation                                    = {}
    translation["Log Partition Function"         ] = math.log ( z )
    translation["Helmholtz Free Energy"          ] = - RT * translation["Log Partition Function"]
    translation["Constant Pressure Heat Capacity"] = 2.5 * R
    translation["Constant Volume Heat Capacity"  ] = 1.5 * R
    translation["Helmholtz Free Energy"          ] = 2.5 * RT
    translation["Entropy"                        ] = R * ( translation["Log Partition Function"] + 1.5 )
    translation["Internal Energy"                ] = 1.5 * RT
    translation["Gibbs Free Energy"              ] = translation["Helmholtz Free Energy"] - temperature * translation["Entropy"]
    # . Vibrational contributions.
    # . Initialization.
    vibration = {}
    # . Polyatomic system.
    if natoms > 1:
        # . Get some state-related information.
        frequencies = state.frequencies
        nRTModes    = state.nRTModes
        # . Remove imaginary and rotation-translation modes.
        start = len ( frequencies )
        zero  = []
        for ( i, f ) in enumerate ( frequencies ):
            zero.append ( f )
            if len ( zero ) > nRTModes:
                old = math.fabs ( zero.pop ( 0 ) )
                if f > old:
                    start = i
                    break
        # . Loop over the frequencies.
        cv  = 0.0
        lnz = 0.0
        u   = 0.0
        zpe = 0.0
        for i in range ( start, len ( frequencies ) ):
            omega = _Wavenumbers_To_Joules * frequencies[i]
            hvkt  = omega / ( Constants.Boltzmann * temperature )
            lnz  -= 0.5 * hvkt
            zpe  += 0.5 * 1.0e-3 * Constants.Avogadro_Number * omega
            if hvkt <= _ExponentialUnderflow:
                expm  = math.exp ( - hvkt )
                expp  = 1.0 - expm
                lnz  -= math.log ( expp )
                expp  = expm / expp
                cv   += R * hvkt * hvkt * expp * ( 1.0 + expp )
                u    += RT * hvkt * expp
        # . Assign the terms.
        vibration["Log Partition Function"         ] = lnz
        vibration["Helmholtz Free Energy"          ] = - RT * lnz
        vibration["Constant Pressure Heat Capacity"] = cv
        vibration["Constant Volume Heat Capacity"  ] = cv
        vibration["Internal Energy"                ] = u + zpe
        vibration["Entropy"                        ] = R * lnz + vibration["Internal Energy"] / temperature
        vibration["Gibbs Free Energy"              ] = vibration["Internal Energy"] - temperature * vibration["Entropy"]
        vibration["Helmholtz Free Energy"          ] = vibration["Internal Energy"]
    # . Get the totals.
    tdics = {}
    for key in _ThermodynamicProperties:
        tdics[key] = electronic.get ( key, 0.0 ) + rotation.get ( key, 0.0 ) + translation.get ( key, 0.0 ) + vibration.get ( key, 0.0 )
    # . Add in some extra terms pertaining to the ideal gas.
    factor = math.log ( Constants.Avogadro_Number ) - 1.0
    tdics["Gibbs Free Energy"]     -= RT * factor
    tdics["Entropy"]               -= R  * factor
    tdics["Helmholtz Free Energy"] += RT * factor
    # . Return.
    return tdics

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

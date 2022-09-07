"""Defines the object function class for a system's geometrical variables."""

import math

from pCore                                  import Clone                  , \
                                                   logFile                , \
                                                   LogFileActive
from pScientific                            import Constants              , \
                                                   Units
from pScientific.Arrays                     import Array
from pScientific.Geometry3                  import Coordinates3
from pScientific.RandomNumbers              import NormalDeviateGenerator , \
                                                   RandomNumberGenerator
from pScientific.ObjectiveFunctionIterators import ObjectiveFunction
from pScientific.Symmetry                   import CrystalSystem

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The conversion factor from amu A^2 ps^-2 to kJ mol^-1 (equivalent to AMU_TO_KG * NAVOGADRO * 10^-3 / MS_TO_APS^2 ).
_AMUA2PS2_TO_KJMOL = 1.0e-2

# . The conversion factor from m s^-1 to A ps^-1.
_MS_TO_APS = 1.0e-2

# . The conversion factor from amu A^2 ps^-2 to Kelvin.
_AMUA2PS2_TO_K = Constants.Atomic_Mass / ( Constants.Boltzmann * _MS_TO_APS**2 )

# . The conversion factor from kJ mol^-1 to amu A^2 ps^-2.
_KJMOL_TO_AMUA2PS2 = 1.0e+2

# . The conversion factor from atm. Angstroms^3 to kJ mol^-1.
_PV_TO_KJMOL = Units.Pressure_Atmospheres_To_Pascals * Constants.Avogadro_Number * 1.0e-33

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemGeometryObjectiveFunction ( ObjectiveFunction ):
    """The object function for geometry optimizations and dynamics of a single system."""

    _attributable = dict ( ObjectiveFunction._attributable )
    _attributable.update ( { "atomWeights"          :  None ,
                             "degreesOfFreedom"     :     0 ,
                             "fractional"           :  None ,
                             "freeAtoms"            :  None ,
                             "hasPressure"          : False ,
                             "hasSymmetry"          : False ,
                             "iCoordinates3"        :  None ,
                             "linearScalars"        :  None ,
                             "linearVectors"        :  None ,
                             "log"                  :  None , # . Debugging only.
                             "nCoordinateVariables" :     0 ,
                             "nSymmetryVariables"   :     0 ,
                             "nVariables"           :     0 ,
                             "rtReference"          :  None ,
                             "startingHessian"      :  None ,
                             "system"               :  None ,
                             "variableWeights"      :  None } )

    def __len__ ( self ): return self.numberOfVariables

    def AccelerationConversionFactor ( self ):
        """Return the conversion factor from gradient units to acceleration units."""
        # . K to amu A^2 ps^-2.
        return ( - _KJMOL_TO_AMUA2PS2 )

    def Accelerations ( self, variables, accelerations ):
        """Evaluate the function and accelerations."""
        f = self.FunctionGradients ( variables, accelerations )
        accelerations.Scale ( - _KJMOL_TO_AMUA2PS2 )
        return f

    def AddLinearConstraint ( self, constraint ):
        """Add a linear constraint."""
        if len ( constraint ) != self.nVariables: raise ValueError ( "Invalid number of variables in linear constraint." )
        # . Orthogonalize to existing constraints.
        if self.linearVectors is not None:
            constraint = Clone ( constraint )
            self.linearVectors.ProjectOutOfArray ( constraint )
        # . Check to see if the constraint is valid.
        cNorm2 = constraint.Norm2 ( )
        if cNorm2 > 1.0e-10:
            constraint.Scale ( 1.0 / cNorm2 )
            # . Allocate space for new constraints.
            ncolumns = 1
            if self.linearVectors is not None: ncolumns += self.linearVectors.columns
            newConstraints = Array.WithExtents ( len ( constraint ), ncolumns )
            # . Copy over constraints.
            if self.linearVectors is not None:
                for r in range ( self.linearVectors.rows ):
                    for c in range ( self.linearVectors.columns ):
                        newConstraints[r,c] = self.linearVectors[r,c]
            for r in range ( len ( constraint ) ): newConstraints[r,ncolumns-1] = constraint[r]
            self.linearVectors = newConstraints
            # . Determine the linear scalars.
            self.linearScalars = Array.WithExtent ( self.linearVectors.columns )
            reference          = Array.WithExtent ( self.linearVectors.rows    )
            if self.rtReference is None: self.iCoordinates3.CopyTo        ( reference.iterator )
            else:                        self.rtReference.iterator.CopyTo ( reference.iterator )
            self.linearVectors.VectorMultiply ( reference, self.linearScalars, transpose = True )
            # . Reset the number of degrees of freedom.
            self.degreesOfFreedom = self.nVariables - len ( self.linearScalars )

    def ApplyConstraintsToVector ( self, vector, **keywordArguments ):
        """Apply constraints to a vector."""
        if self.linearVectors is not None: self.linearVectors.ProjectOutOfArray ( vector )

    def DefinePressure ( self ):
        """Set up pressure calculation."""
        try:
            cc = self.system.symmetry.crystalSystem
            if not isinstance ( cc, CrystalSystem ): raise
            self.hasPressure = True
        except:
            pass

    def DefineWeights ( self ):
        """Define atom weights - use masses by default."""
        self.atomWeights = Array.FromIterable ( [ atom.mass for atom in self.system.atoms ] )
        if self.nVariables > 0:
            self.variableWeights = Array.WithExtent ( self.nVariables )
            self.variableWeights.Set ( 0.0 )
            if self.freeAtoms is None: indices = range ( len ( self.system.atoms ) )
            else:                      indices = self.freeAtoms
            n = 0
            for iatom in indices:
                w = math.sqrt ( self.atomWeights[iatom] )
                for j in range ( 3 ):
                    self.variableWeights[n] = w
                    n += 1

    def DistanceSquared ( self, variables1, variables2, gradients ):
        """Calculate the distance squared and its gradients between two points."""
        # . This needs to be updated for crystal cases.
        variables1.CopyTo ( gradients )
        gradients.Add ( variables2, scale = -1.0 )
        d2 = gradients.DotSelf ( )
        gradients.Scale ( 2.0 )
        return d2

    @classmethod
    def FromSystem ( selfClass, system ):
        """Constructor given a system."""
        # . Basic object.
        self = selfClass ( )
        # . Define the system.
        self.system = system
        # . Get the number of atoms.
        nAtoms = len ( system.atoms )
        nFree  = nAtoms
        # . Check for fixed atoms.
        if system.freeAtoms is not None:
            self.freeAtoms = system.freeAtoms
            nFree = len ( self.freeAtoms )
        # . Define the coordinate iterator.
        if self.freeAtoms is None: self.iCoordinates3 = self.system.coordinates3.iterator
        else:                      self.iCoordinates3 = self.system.coordinates3.RowIterator ( selection = self.freeAtoms )
        # . Set the number of variables.
        self.nCoordinateVariables = 3 * nFree
        self.nVariables           = self.nCoordinateVariables
        self.degreesOfFreedom     = self.nVariables
        # . Finish up.
        return self

    def Function ( self, variables ):
        """Evaluate the function."""
        self.VariablesPut ( variables )
        f = self.system.Energy ( log = self.log )
        return f

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        self.VariablesPut ( variables )
        f = self.system.Energy ( doGradients = True, log = self.log )
        self.GradientsGet ( gradients )
        return f

    def GradientsGet ( self, gradients ):
        """Get the gradients."""
        # . Needs to be made more efficient by not redoing each time if there are fixed atoms.
        if self.freeAtoms is None: iGradients3 = self.system.scratch.gradients3.iterator
        else:                      iGradients3 = self.system.scratch.gradients3.RowIterator ( selection = self.freeAtoms )
        if self.hasSymmetry:
            spg = self.system.scratch.symmetryParameterGradients
            spg.MakeFractionalDerivatives ( self.system.symmetryParameters, self.system.coordinates3, self.system.scratch.gradients3 )
            spg.MakeCrystalDerivatives    ( self.system.symmetryParameters )
            iGradients3.CopyTo ( gradients[0:self.nCoordinateVariables].iterator )
            self.system.symmetry.crystalSystem.EmptyGradientsToVector ( gradients, spg, self.nCoordinateVariables )
        else:
            iGradients3.CopyTo ( gradients.iterator )
        if self.variableWeights is not None: gradients.Divide ( self.variableWeights )
        if self.linearVectors is not None: self.linearVectors.ProjectOutOfArray ( gradients )

    def IncludeSymmetryParameters ( self ):
        """Flag the symmetry parameters as variables."""
        if ( self.system.symmetry is not None ) and ( self.system.symmetry.crystalSystem is not None ):
            nVariables = len ( self.system.symmetry.crystalSystem )
            if nVariables > 0:
                self.fractional         = Coordinates3.WithExtent ( self.nCoordinateVariables // 3 )
                self.nSymmetryVariables = nVariables
                self.nVariables        += nVariables
                self.hasSymmetry        = True

    def Pressure ( self, ke ):
        """Calculate the pressure and volume."""
        pressure = 0.0
        volume   = 0.0
        if self.hasPressure:
            cc       = self.system.symmetry.crystalSystem
            crd3     = self.system.scratch.Get ( "coordinates3NB", self.system.coordinates3 ) # . Centered coordinates if there are any.
            sp       = self.system.symmetryParameters
            spg      = self.system.scratch.symmetryParameterGradients
            spg.MakeFractionalDerivatives ( sp, crd3, self.system.scratch.gradients3 )
            spg.MakeCrystalDerivatives    ( sp )
            volume   = sp.volume
            pressure = ( ( 2.0 * ke ) / ( 3.0 * volume ) - cc.GetVolumeDerivative ( sp, spg ) ) / _PV_TO_KJMOL
        return ( pressure, volume )

    def RemoveRotationTranslation ( self, reference = None ):
        """Remove rotation and translational degrees of freedom.

        |reference| is the reference coordinates.

        Nothing is done if there are fixed atoms.

        Weighting is done using |atomWeights| if it is has been defined.
        """
        if self.freeAtoms is None:
            # . Check that RT can be removed.
            doRotation    = ( self.system.symmetry is None ) or ( ( self.system.symmetry is not None ) and ( self.system.symmetry.crystalSystem is None ) )
            doTranslation = True
            # . Remove RT.
            if doRotation or doTranslation:
                # . Assign the reference coordinates and modify the system coordinates if necessary.
                if reference is None:
                    self.rtReference = Clone ( self.system.coordinates3 )
                else:
                    self.rtReference = Clone ( reference )
                    self.system.coordinates3.Superimpose ( self.rtReference, weights = self.atomWeights )
                # . Get the RT vectors.
                ( self.linearVectors, self.linearScalars ) = self.rtReference.RotationTranslationVectors ( doRotation, doTranslation, dimension = self.nVariables, weights = self.atomWeights )
                self.degreesOfFreedom = self.nVariables - len ( self.linearScalars )

    def SetStartingHessian ( self, hessian ):
        """Set the starting hessian."""
        self.startingHessian = hessian

    def StartingHessian ( self, variables ):
        """Get a starting hessian."""
        if self.startingHessian is None: hessian = self.NumericalHessian ( variables )
        else:                            hessian = self.startingHessian
        if self.linearVectors is not None:
            hessian = hessian.ProjectOutVectors ( self.linearVectors          )
            hessian.Raise                       ( self.linearVectors, 3.0e+04 )
        return hessian

    def Temperature ( self, velocities ):
        """Calculate the kinetic energy and temperature."""
        ke           = velocities.DotSelf ( )
        temperature  = _AMUA2PS2_TO_K * ke / float ( self.degreesOfFreedom )
        ke          *= 0.5 * _AMUA2PS2_TO_KJMOL
        return ( ke, temperature )

    def TemperatureConversionFactor ( self ):
        """Return the conversion factor from temperature units to dynamics units."""
        # . K to amu A^2 ps^-2.
        return ( 1.0 / _AMUA2PS2_TO_K )

    def UndefineWeights ( self ):
        """Remove any weighting."""
        if ( self.atomWeights is not None ) or ( self.variableWeights is not None ):
            self.atomWeights     = None
            self.variableWeights = None
            if self.linearVectors is not None: self.RemoveRotationTranslation ( )

    def VariablesAllocate ( self ):
        """Return an object to hold the variables."""
        if self.numberOfVariables > 0:
            variables = Array.WithExtent ( self.numberOfVariables )
            variables.Set ( 0.0 )
            return variables
        else:
            return None

    def VariablesGet ( self, variables ):
        """Fill the variable array."""
        if self.hasSymmetry:
            sp = self.system.symmetryParameters
            self.iCoordinates3.CopyTo       ( self.fractional.iterator )
            self.fractional.Rotate          ( sp.inverseH              )
            self.fractional.iterator.CopyTo ( variables[0:self.nCoordinateVariables].iterator )
            self.system.symmetry.crystalSystem.EmptySymmetryParametersToVector ( variables, sp, self.nCoordinateVariables )
        else:
            self.iCoordinates3.CopyTo       ( variables.iterator       )
        if self.variableWeights is not None: variables.Multiply ( self.variableWeights )

    def VariablesPut ( self, variables ):
        """Empty the variable array."""
        if self.variableWeights is not None: variables.Divide ( self.variableWeights )
        if self.hasSymmetry:
            # . Get new box shape first.
            sp = self.system.symmetryParameters
            self.system.symmetry.crystalSystem.FillSymmetryParametersFromVector ( variables, sp, self.nCoordinateVariables )
            # . Now get coordinates.
            variables[0:self.nCoordinateVariables].iterator.CopyTo ( self.fractional.iterator )
            self.fractional.Rotate          ( sp.H               )
            self.fractional.iterator.CopyTo ( self.iCoordinates3 )
        else:
            variables.iterator.CopyTo       ( self.iCoordinates3 )
        if self.variableWeights is not None: variables.Multiply ( self.variableWeights )

    def VelocitiesAllocate ( self ):
        """Return an object with velocities.

        The velocities for the system must already exist in scratch.
        """
        return self.system.scratch.velocities

    def VelocitiesAssign ( self, temperature, normalDeviateGenerator = None ):
        """Set up the velocities for the system at a particular temperature.

        If |temperature| is None the velocities must already exist.
        """
        # . Check for an existing set of velocities of the correct size.
        velocities = self.system.scratch.Get ( "velocities", None )
        QASSIGN = ( velocities is None ) or ( ( velocities is not None ) and ( len ( velocities ) != self.nVariables ) )
        # . Assign velocities.
        if QASSIGN:
            if temperature is None:
                raise ValueError ( "Velocities need to be assigned for the system but a temperature has not been specified." )
            else:
                # . Get the generator.
                if normalDeviateGenerator is None:
                    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithRandomSeed ( ) )
                # . Get the velocities.
                sigma      = _MS_TO_APS * math.sqrt ( Constants.Boltzmann * temperature / Constants.Atomic_Mass )
                velocities = Array.WithExtent ( self.numberOfVariables )
                normalDeviateGenerator.NextStandardDeviates ( velocities )
                velocities.Scale ( ( sigma ) )
                self.system.scratch.velocities = velocities
        # . Project out linear constraints.
        if self.linearVectors is not None: self.linearVectors.ProjectOutOfArray ( velocities )
        # . Scale velocities if necessary (even for assigned velocities).
        if temperature is not None:
            ( ke, tactual ) = self.Temperature ( velocities )
            velocities.Scale ( math.sqrt ( temperature / tactual ) )

    def VolumeScale ( self, scale ):
        """Scale the volume of the system."""
        if self.hasPressure:
            self.system.symmetry.crystalSystem.ScaleVolume ( self.system.symmetryParameters, scale )

    @property
    def numberOfVariables ( self ):
        return self.nVariables

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

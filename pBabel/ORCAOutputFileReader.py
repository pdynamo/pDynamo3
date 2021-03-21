"""Classes and functions for reading ORCA output files."""

from  pCore                 import logFile            , \
                                   LogFileActive      , \
                                   TextFileReader
from  pMolecule             import System
from  pScientific           import PeriodicTable      , \
                                   Units
from  pScientific.Geometry3 import Coordinates3       , \
                                   Vector3
from  pScientific.Spectra   import ContinuousSpectrum , \
                                   DataRange          , \
                                   DataSet            , \
                                   GaussianLineShape  , \
                                   StickSpectrum
from .ExportImport          import _Importer

import os, os.path, re

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultFWHM        = 1500.0 # . Wavenumbers.
_DefaultSystemLabel = "ORCA Output File"

#===================================================================================================================================
# . Data class.
#===================================================================================================================================
class Data:
    """A class to hold data."""

    def __init__ ( self ):
        """Constructor."""
        self.index = {}

    def GetItem ( self, fullName, default = None ):
        """Get an item from the index."""
        return self.index.get ( fullName, default )

    def SetItem ( self, fullName, item ):
        """Set an item and index it."""
        self.index[fullName] = item
        return item

#===================================================================================================================================
# . ORCA output file reader class.
#===================================================================================================================================
class ORCAOutputFileReader ( TextFileReader ):
    """ORCAOutputFileReader is the class for ORCA output files that are to be read."""

    _attributable = dict ( TextFileReader._attributable )
    _attributable.update ( { "context" : None ,
                             "frames"  : None } )

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            currentFrame = Data ( )
            self.context = None
            self.frames  = [ currentFrame ]
            # . Open the file.
            self.Open ( )
            try:
                while True:
                    # . Get the next line.
                    line = self.GetLine ( signalWarnings = False )
                    # . Cartesian coordinates.
                    if line.startswith ( "CARTESIAN COORDINATES (ANGSTROEM)" ):
                        atomicNumbers = currentFrame.SetItem ( "Atomic Numbers", [] )
                        xyz           = currentFrame.SetItem ( "XYZ"           , [] )
                        self.GetLine ( )
                        while True:
                            tokens = self.GetTokens ( converters =  [ PeriodicTable.AtomicNumber, float, float, float ] )
                            if len ( tokens ) > 0:
                                atomicNumbers.append ( tokens[0] )
                                xyz.append ( [ tokens[1], tokens[2], tokens[3] ] )
                            else:
                                break
                    # . CHELPG charges.
                    elif line == "CHELPG Charges":
                        self.GetLine ( )
                        espCharges = currentFrame.SetItem ( "CHELPG Charges", [] )
                        while True:
                            tokens = self.GetTokens ( )
                            if len ( tokens ) < 4: break
                            else: espCharges.append ( float ( tokens[-1] ) )
                    # . Charge.
                    elif line.startswith ( "Total Charge" ):
                        currentFrame.SetItem ( "Charge" , int ( line.split ( )[-1] ) )
                    # . Convergence OK.
                    elif line.find ( "SCF CONVERGED AFTER" ) >= 0:
                        tokens = line.split ( )
                        currentFrame.SetItem ( "SCF Converged" , True               )
                        currentFrame.SetItem ( "SCF Cycles"    , int ( tokens[-3] ) )
                    # . Convergence not OK.
                    elif line.find ( "SCF NOT CONVERGED AFTER" ) >= 0:
                        tokens          = line.split ( )
                        currentFrame.SetItem ( "SCF Converged" , False              )
                        currentFrame.SetItem ( "SCF Cycles"    , int ( tokens[-3] ) )
                    # . Dipole.
                    elif line.startswith ( "Total Dipole Moment" ):
                        tokens = line.split ( )
                        if not hasattr ( self, "dipole" ): dipole = currentFrame.SetItem ( "Dipole Moment Vector", [] )
                        for ( i, token ) in enumerate ( tokens[-3:] ):
                            dipole.append ( Units.Dipole_Atomic_Units_To_Debyes * float ( token ) )
                    # . Dispersion correction.
                    elif line.startswith ( "Dispersion correction" ):
                        try   : currentFrame.SetItem ( "Dispersion Correction", float ( line.split ( )[-1] ) )
                        except: pass
                    # . Final single point energy (includes everything including dispersion except outlying charge correction if present).
                    elif line.startswith ( "FINAL SINGLE POINT ENERGY" ):
                        currentFrame.SetItem ( "Potential Energy", float ( line.split ( )[-1] ) )
                    # . Energy after SCF with outlying charge correction (excludes dispersion).
                    elif line.startswith ( "Total Energy after outlying charge correction =" ):
                        currentFrame.SetItem ( "Total Energy After Outlying Charge Correction", float ( line.split ( )[-2] ) )
                    # . Energy after SCF (excludes dispersion).
                    elif line.startswith ( "Total Energy       " ):
                        currentFrame.SetItem ( "Total Energy", float ( line.split ( )[3] ) )
                    # . Geometry optimization converged.
                    elif line.find ( "THE OPTIMIZATION HAS CONVERGED" ) >= 0:
                        currentFrame.SetItem ( "Geometry Optimization Converged", True )
                    # . Geometry optimization not converged.
                    elif line.find ( "The optimization has not yet converged" ) >= 0:
                        currentFrame.SetItem ( "Geometry Optimization Converged", False )
                    # . Geometry optimization cycle.
                    elif line.find ( "GEOMETRY OPTIMIZATION CYCLE" ) >= 0:
                        cycle = int ( line.split ( )[-2] )
                        if cycle > 1:
                            currentFrame = Data ( )
                            self.frames.append ( currentFrame )
                        currentFrame.SetItem ( "Geometry Optimization Cycle", cycle )
                    # . Geometry optimization - final energy evaluation.
                    elif line.find ( "FINAL ENERGY EVALUATION AT THE STATIONARY POINT" ) >= 0:
                        currentFrame = Data ( )
                        self.frames.append ( currentFrame )
                        currentFrame.SetItem ( "Geometry Optimization Stationary Point", True )
                        currentFrame.SetItem ( "Geometry Optimization Cycles", int ( self.GetTokens ( )[2] ) )
                    # . Gibbs free enthalpy minus the electronic energy.
                    elif line.startswith ( "G-E(el)" ):
                        currentFrame.SetItem ( "Gibbs Free Energy", float ( line.split ( )[2] ) )
                    # . Loewdin charges.
                    elif line == "LOEWDIN ATOMIC CHARGES":
                        self.GetLine ( )
                        loewdinCharges = currentFrame.SetItem ( "Loewdin Charges", [] )
                        while True:
                            tokens = self.GetTokens ( )
                            if ( len ( tokens ) <= 0 ) or ( tokens[0] == "Sum" ): break
                            else: loewdinCharges.append ( float ( tokens[-1] ) )
                    # . Loewdin charges and spin densities.
                    elif line == "LOEWDIN ATOMIC CHARGES AND SPIN POPULATIONS":
                        self.GetLine ( )
                        loewdinCharges = currentFrame.SetItem ( "Loewdin Charges", [] )
                        loewdinSpins   = currentFrame.SetItem ( "Loewdin Spins"  , [] )
                        while True:
                            tokens = self.GetTokens ( )
                            if ( len ( tokens ) <= 0 ) or ( tokens[0] == "Sum" ): break
                            else:
                                loewdinCharges.append ( float ( tokens[-2] ) )
                                loewdinSpins.append   ( float ( tokens[-1] ) )
                    # .  Mayer bond orders.
                    elif line.startswith ( "Mayer bond orders larger than" ):
                        mayerBondOrders = currentFrame.SetItem ( "Mayer Bond Orders", [] )
                        while True:
                            line = self.GetLine ( )
                            if line.startswith ( "B(" ):
                                items = line[2:].split ( "B(" )
                                for item in items:
                                    tokens = item.split ( ":" )
                                    words  = tokens[0].replace ( ",", " " ).replace ( "-", " " ).split ( )
                                    i      = int ( words[0] )
                                    j      = int ( words[2] )
                                    value  = float ( tokens[1] )
                                    mayerBondOrders.append ( [ i, j, value ] )
                            else:
                                break
                    # . Mulliken charges.
                    elif line == "MULLIKEN ATOMIC CHARGES":
                        self.GetLine ( )
                        mullikenCharges = currentFrame.SetItem ( "Mulliken Charges", [] )
                        while True:
                            tokens = self.GetTokens ( )
                            if ( len ( tokens ) <= 0 ) or ( tokens[0] == "Sum" ): break
                            else: mullikenCharges.append ( float ( tokens[-1] ) )
                    # . Mulliken charges and spin densities.
                    elif line == "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS":
                        self.GetLine ( )
                        mullikenCharges = currentFrame.SetItem ( "Mulliken Charges", [] )
                        mullikenSpins   = currentFrame.SetItem ( "Mulliken Spins"  , [] )
                        while True:
                            tokens = self.GetTokens ( )
                            if ( len ( tokens ) <= 0 ) or ( tokens[0] == "Sum" ): break
                            else:
                                mullikenCharges.append ( float ( tokens[-2] ) )
                                mullikenSpins.append   ( float ( tokens[-1] ) )
                    # . Multiplicity.
                    elif line.startswith ( "Multiplicity" ):
                        currentFrame.SetItem ( "Multiplicity" , int ( line.split ( )[-1] ) )
                    # . Orbital energies.
                    elif line == "ORBITAL ENERGIES":
                        doContinue = True
                        self.GetLine ( )
                        # . Loop over sets.
                        while doContinue:
                            line = self.GetLine ( )
                            if   line.find ( "SPIN UP ORBITALS"   ) >= 0:
                                tag      = "Alpha "
                            elif line.find ( "SPIN DOWN ORBITALS" ) >= 0:
                                doContinue = False
                                tag        = "Beta "
                            else:
                                doContinue = False
                                tag        = ""
                            HOMO = -1
                            LUMO = -1
                            orbitalEnergies    = currentFrame.SetItem ( tag + "Orbital Energies"    , [] )
                            orbitalOccupancies = currentFrame.SetItem ( tag + "Orbital Occupancies" , [] )
                            for i in range ( 3 ): self.GetLine ( )
                            index = 0
                            while True:
                                tokens = self.GetTokens ( )
                                if len ( tokens ) < 3: break
                                else:
                                    occupancy = float ( tokens[1] )
                                    if ( HOMO == -1 ) and ( LUMO == -1 ) and ( occupancy <= 1.0e-6 ):
                                        HOMO = index - 1
                                        LUMO = index
                                    orbitalEnergies.append    ( float ( tokens[2] ) )
                                    orbitalOccupancies.append ( occupancy           )
                                    index += 1
                            currentFrame.SetItem ( tag + "HOMO Index", HOMO )
                            currentFrame.SetItem ( tag + "LUMO Index", LUMO )
                    # . Run type - energy.
                    elif line.find ( "Single Point Calculation" ) >= 0:
                        currentFrame.SetItem ( "Run Type", "Single Point" )
                    # . Run type - geometry optimization.
                    elif line.find ( "Geometry Optimization Run" ) >= 0:
                        currentFrame.SetItem ( "Run Type", "Geometry Optimization" )
                    # . <S**2>.
                    elif line.startswith ( "Expectation value of <S**2>" ):
                        currentFrame.SetItem ( "Spin Squared", float ( line.split ( ":", 1 )[-1] ) )
                    # . Spectrum - absorption with transition electric dipole moments.
                    elif line.startswith ( "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" ):
                        currentFrame.SetItem ( "Absorption Spectrum (Transition Electric Dipole Moments)", self.ParseAbsorptionSpectrum ( ) )
                    # . Spectrum - absorption with transition velocity dipole moments.
                    elif line.startswith ( "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" ):
                        currentFrame.SetItem ( "Absorption Spectrum (Transition Velocity Dipole Moments)", self.ParseAbsorptionSpectrum ( ) )
                    # . Spectrum - absorption.
                    elif line.startswith ( "ABSORPTION SPECTRUM" ):
                        key  = "Absorption Spectrum"
                        if self.context is not None: key += " ({:s})".format ( self.context )
                        item = self.ParseAbsorptionSpectrum ( )
                        if item is not None: currentFrame.SetItem ( key, item )
                    # . Spectrum - CD.
                    elif line.startswith ( "CD SPECTRUM" ):
                        key = "CD Spectrum"
                        if self.context is not None:
                            key += " ({:s})".format ( self.context )
                            self.context = None
                        currentFrame.SetItem ( key, self.ParseCDSpectrum ( ) )
                    # . Wavefunction type.
                    elif line.startswith ( "Kohn-Sham wavefunction type" ):
                        tag = line[-3:]
                        currentFrame.SetItem ( "Is Spin Restricted", ( tag not in ( "UHF", "UKS" ) ) )
                    # . Contexts.
                    elif line.startswith ( "CASSCF UV, CD spectra and dipole moments"                            ): self.context = "CASSCF"
                    elif line.startswith ( "CASSCF (NEVPT2 diagonal energies) UV, CD spectra and dipole moments" ): self.context = "NEVPT2"
                    elif line.startswith ( "CI-EXCITATION SPECTRA"                                               ): self.context = "CI"
                    elif line.startswith ( "TD-DFT/TDA-EXCITATION SPECTRA"                                       ): self.context = "TDDFT"
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.context = None
            self.log     = None
            self.isParsed = True

    def ParseAbsorptionSpectrum ( self ):
        """Parse an absorption spectrum."""
        self.GetLine ( ) # . Header.
        tokens = self.GetTokens ( )
        if ( len ( tokens ) >= 1 ) and ( tokens[0] in ( "State", "States" ) ):
            for i in range ( 2 ): self.GetLine ( ) # . Remaining lines in header.
            fields = [ "Energy (cm^-1)"            ,
                       "Wavelength (nm)"           ,
                       "Oscillator Strength"       ,
                       "Transition Moment Squared" ,
                       "Transition Moment X"       ,
                       "Transition Moment Y"       ,
                       "Transition Moment Z"       ]
            values = []
            while True:
                tokens = self.GetTokens ( ) #converters =  [ None, float, float, float, float, float, float, float ] )
                if len ( tokens ) >= 8:
                    try   : values.append ( [ float ( token ) for token in tokens[-7:] ] )
                    except: self.Warning ( "Unable to convert absorption spectrum tokens.", True )
                else: break
            return { "Fields" : fields, "Values" : values }
        else:
            return None

    def ParseCDSpectrum ( self ):
        """Parse a CD spectrum."""
        for i in range ( 4 ): self.GetLine ( ) # . Top of table.
        fields = [ "Energy (cm^-1)"      ,
                   "Wavelength (nm)"     ,
                   "Rotatory Strength"   ,
                   "Transition Moment X" ,
                   "Transition Moment Y" ,
                   "Transition Moment Z" ]
        values = []
        while True:
            tokens = self.GetTokens ( ) #converters =  [ None, float, float, float, float, float, float ] )
            if len ( tokens ) >= 7:
                try   : values.append ( [ float ( token ) for token in tokens[-6:] ] )
                except: self.Warning ( "Unable to convert CD spectrum tokens.", True )
            else: break
        return { "Fields" : fields, "Values" : values }

    @classmethod
    def PathToAbsorptionSpectrum ( selfClass                                         ,
                                   path                                              ,
                                   context                            = None         ,
                                   frameIndex                         = -1           ,
                                   fwhm                               = _DefaultFWHM ,
                                   log                                = logFile      ,
                                   useTransitionVelocityDipoleMoments = False        ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToAbsorptionSpectrum ( context                            = context                            ,
                                             frameIndex                         = frameIndex                         ,
                                             fwhm                               = fwhm                               ,
                                             useTransitionVelocityDipoleMoments = useTransitionVelocityDipoleMoments )

    @classmethod
    def PathToCoordinates3 ( selfClass, path, frameIndex = -1, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( frameIndex = frameIndex )

    @classmethod
    def PathToSystem ( selfClass, path, frameIndex = -1, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToSystem ( frameIndex = frameIndex )

    def ToAbsorptionSpectrum ( self, context = None, frameIndex = -1, fwhm = _DefaultFWHM, useTransitionVelocityDipoleMoments = False ):
        """Return an absorption spectrum using the same defaults as ORCA."""
        continuousSpectrum = None
        if self.isParsed:
            if context is None:
                if useTransitionVelocityDipoleMoments: key = "Absorption Spectrum (Transition Velocity Dipole Moments)"
                else:                                  key = "Absorption Spectrum (Transition Electric Dipole Moments)"
            else: key = "Absorption Spectrum ({:s})".format ( context )
            rawData = self.frames[frameIndex].GetItem ( key )
            if rawData is not None:
                fields  = rawData["Fields"]
                peaks   = rawData["Values"]
                xField  = fields.index ( "Energy (cm^-1)"      )
                yField  = fields.index ( "Oscillator Strength" )
                # . Generate the stick spectrum.
                stickSpectrum = StickSpectrum.WithSize ( len ( peaks ) )
                xData         = stickSpectrum.abscissae.data
                yData         = stickSpectrum.ordinates.data
                for ( p, peak ) in enumerate ( peaks ):
                    xData[p] = peak[xField]
                    yData[p] = peak[yField]
                # . Generate the continuous spectrum in the appropriate units (wavenumbers and extinction coefficient).
                abscissae          = DataSet.FromDataRange ( DataRange.FromSizeStartStop ( 1024, 5000.0, 40000.0 ) )
                lineShape          = GaussianLineShape     ( center = 0.0, fwhm = fwhm )
                continuousSpectrum = ContinuousSpectrum.FromStickSpectrum ( stickSpectrum, abscissae, lineShape )
                continuousSpectrum.ScaleOrdinates ( 1.0 / 4.33e-09 )
                continuousSpectrum.abscissae.units = "cm^-1"
                continuousSpectrum.ordinates.units = None
        return continuousSpectrum

    def ToCoordinates3 ( self, frameIndex = -1 ):
        """Return a coordinates3 object."""
        coordinates3 = None
        if self.isParsed:
            frame = self.frames[frameIndex]
            xyz   = frame.GetItem ( "XYZ" )
            if xyz is not None:
                coordinates3 = Coordinates3.WithExtent ( len ( xyz ) )
                for ( i, ( x, y, z ) ) in enumerate ( xyz ):
                    coordinates3[i,0] = x
                    coordinates3[i,1] = y
                    coordinates3[i,2] = z
                frame.SetItem ( "Coordinates3", coordinates3 )
            else:
                coordinates3 = frame.GetItem ( "Coordinates3" )
        return coordinates3

    def ToElectronicState ( self, frameIndex = -1 ):
        """Return an electronic state object."""
        if self.isParsed:
            frame = self.frames[frameIndex]
            return { "Charge"       : frame.GetItem ( "Charge"       , default = 0 ) ,
                     "Multiplicity" : frame.GetItem ( "Multiplicity" , default = 1 ) }
        else:
            return None

    def ToSystem ( self, frameIndex = -1 ):
        """Return a system."""
        system = None
        if self.isParsed:
            frame         = self.frames[frameIndex]
            atomicNumbers = frame.GetItem ( "Atomic Numbers" )
            if atomicNumbers is not None:
                ( head, tail )      = os.path.split ( self.path )
                system              = System.FromAtoms ( atomicNumbers )
                system.label        = _DefaultSystemLabel + " " + tail
                system.coordinates3 = self.ToCoordinates3 ( )
                system.Set ( "Electronic State", self.ToElectronicState ( ), attribute = "electronicState" )
        return system

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : ORCAOutputFileReader.PathToCoordinates3 ,
                         System       : ORCAOutputFileReader.PathToSystem       } ,
                       [ "olog", "oout", "OLOG", "OOUT" ], "ORCA Output", defaultFunction = ORCAOutputFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

"""The DFTB+ QC model."""

import glob, math, os, os.path, subprocess, re

from   pCore                     import logFile           , \
                                        LogFileActive     , \
                                        NotInstalledError
from   pScientific               import PeriodicTable     , \
                                        Units
from   pScientific.Arrays        import Array
from   pScientific.Geometry3     import Coordinates3      , \
                                        Matrix33          , \
                                        Vector3
from   pScientific.RandomNumbers import RandomString
from  .QCModel                   import QCModel           , \
                                        QCModelState
from  .QCModelError              import QCModelError

# . DFTB+ has many options with a structured input. Only the basic options are processed here.
# . Modifications can be made to the basic options using DFTB+'s extended format.
# . To see what the input file will be like, call "WriteInputFile" with the path keyword.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Orbital quantum number labels.
_OrbitalQuantumNumberLabels = { 0 : "s" ,
                                1 : "p" ,
                                2 : "d" ,
                                3 : "f" ,
                                4 : "g" }

# . Default error suffix.
_DefaultErrorPrefix = "error_"

# . Command environment variable.
_DFTBCommand = "PDYNAMO3_DFTBCOMMAND"

# . Scratch directory.
_DFTBScratch = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), "dftbPlusScratch" )

# . File names - these are invariant.
_AutoTestFile    = "autotest.tag"
_ChargesFile     = "charges.bin"
_DetailedFile    = "detailed.out"
_InputFile       = "dftb_in.hsd"
_OutputFile      = "dftb_out.log"
_PinFile         = "dftb_pin.hsd"
_PointChargeFile = "pointCharges.dat"

# . Input blocks.
_ElectricFieldInput = \
"""  ElectricField = {{
    PointCharges = {{
      GaussianBlurWidth = {:.1f}
      CoordsAndCharges [Angstrom] = DirectRead {{
        Records = {:d}
	File = \"{:s}\"
      }}
    }}
  }}
"""

# . Key/path dictionary.
_ErrorKeys = ( "Detailed" , "Input" , "Output" , "PointCharge" )
_KeyNames  = ( ( "AutoTest"    , _AutoTestFile    ) ,
               ( "Charges"     , _ChargesFile     ) ,
               ( "Detailed"    , _DetailedFile    ) ,
               ( "Input"       , _InputFile       ) ,
               ( "Output"      , _OutputFile      ) ,
               ( "PinFile"     , _PinFile         ) ,
               ( "PointCharge" , _PointChargeFile ) )

# . Occupancy tolerance for HOMO/LUMO determination.
_OccupancyTolerance = 1.0e-6

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelDFTBState ( QCModelState ):
    """A QC model state."""

    _attributable = dict ( QCModelState._attributable )
    _attributable.update ( { "deleteJobFiles" : True , 
                             "paths"          : None } )

    def __del__ ( self ):
        """Deallocation."""
        self.DeleteJobFiles ( )

    def DeleteJobFiles ( self ):
        """Delete job files."""
        if self.deleteJobFiles:
            try:
                scratch  = self.paths["Scratch"]
                jobFiles = glob.glob ( os.path.join ( scratch, "*" ) )
                for jobFile in jobFiles: os.remove ( jobFile )
                os.rmdir ( scratch ) # . Always remove as job files have constant names.
            except:
                pass

    def DeterminePaths ( self, scratch, deleteJobFiles = True, randomScratch = False ):
        """Determine the paths needed by a DFTB job."""
        paths = {}
        if randomScratch:
            scratch = os.path.join ( scratch, RandomString ( ) )
        paths["Scratch"] = scratch
        if not os.path.exists ( scratch ): os.mkdir ( scratch )
        for ( key, name ) in _KeyNames:
            paths[key] = os.path.join ( scratch, name )
        # . Finish up.
        self.deleteJobFiles = deleteJobFiles
        self.paths          = paths

    def SaveErrorFiles ( self, message ):
        """Save the input and output files for inspection if there is an error."""
        for key in _ErrorKeys:
            path = self.paths[key]
            if os.path.exists ( path ):
                ( head, tail ) = os.path.split ( path )
                os.rename ( path, os.path.join ( head, _DefaultErrorPrefix + tail ) )
        raise QCModelError ( message + "\nCheck the files \"{:s}\".".format ( os.path.join ( self.paths["Scratch"], "{:s}*".format ( _DefaultErrorPrefix ) ) ) )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelDFTB ( QCModel ):
    """The DFTB+ QC model class."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "DFTB QC Model"
    _stateObject  = QCModelDFTBState
    _summarizable = dict ( QCModel._summarizable )
    _attributable.update ( { "deleteJobFiles"       : True         ,
                             "extendedInput"        : None         ,
                             "fermiTemperature"     : 10.0         ,
                             "gaussianBlurWidth"    :  0.0         ,
                             "maximumSCCIterations" : 300          ,
                             "randomScratch"        : False        ,
                             "sccTolerance"         : 1.0e-8       ,
                             "scratch"              : _DFTBScratch ,
                             "skfPath"              : "."          ,
                             "useSCC"               : False        } )
    _summarizable.update ( { "deleteJobFiles"       :   "Delete Job Files"                    ,
                             "fermiTemperature"     : ( "Fermi Temperature"      , "{:.3f}" ) ,
                             "gaussianBlurWidth"    : ( "Gaussian Blur Width"    , "{:.3f}" ) ,
                             "maximumSCCIterations" :   "Maximum SCC Iterations"              ,
                             "randomScratch"        :   "Random Scratch"                      ,
                             "sccTolerance"         : ( "SCC Tolerance"          , "{:g}"   ) ,
                             "useSCC"               :   "Use SCC"                           } )

    def AtomicCharges ( self, target, chargeModel = None ):
        """Atomic charges."""
        return target.scratch.dftbOutputData.get ( "Charges", None )

    def AtomicSpins ( self, target, chargeModel = None ):
        """Atomic spins."""
        return target.scratch.dftbOutputData.get ( "Spins", None )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelDFTB, self ).BuildModel ( target, qcSelection = qcSelection )
        state.DeterminePaths ( self.scratch                         ,
                               deleteJobFiles = self.deleteJobFiles ,
                               randomScratch  = self.randomScratch  )
        return state

    def DipoleMoment ( self, target, center = None ):
        """Dipole Moment."""
        return target.scratch.dftbOutputData.get ( "Dipole", None )

    def Energy ( self, target ):
        """Calculate the quantum chemical energy."""
        doGradients    = target.scratch.doGradients
        dftbOutputData = {}
        state          = getattr ( target, self.__class__._stateName )
        target.scratch.dftbOutputData = dftbOutputData
        self.WriteInputFile ( target, doGradients, ( target.nbModel is not None ), target.scratch.qcCoordinates3 )
        isOK = self.Execute ( state )
        if not isOK: state.SaveErrorFiles ( "Error executing program." )
        isOK = self.ReadAutoTestFile ( target, dftbOutputData )
        if not isOK: state.SaveErrorFiles ( "Error reading autotest file." )
        isOK = self.ReadDetailedFile ( target, dftbOutputData )
        if not isOK: state.SaveErrorFiles ( "Error reading detailed file." )
        target.scratch.energyTerms["DFTB QC"] = ( dftbOutputData["Energy"] * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def Execute ( self, state ):
        """Execute the job."""
        # . DFTB+ jobs must be run in the directory where all the DFTB+ files reside.
        try:
            outFile = open ( state.paths["Output"], "w" )
            subprocess.check_call ( [ self.command, _InputFile ], cwd = state.paths["Scratch"], stderr = outFile, stdout = outFile )
            outFile.close ( )
            return True
        except:
            return False

    def ReadAutoTestFile ( self, target, dftbOutputData ):
        """Read the autotest output file."""
        # . Energy, gradients, stress and volume.
        # . Input data is in atomic units throughout.
        state = getattr ( target, self.__class__._stateName )
        try:
            atFile = open ( state.paths["AutoTest"], "r" )
            while True:
                line = next ( atFile )
                if   line.startswith ( "cell_volume" ):
                    dftbOutputData["Volume"] = float ( next ( atFile ) )
                elif line.startswith ( "forces " ): # . Note trailing space to distinguish from next item.
                    gradients3 = target.scratch.qcGradients3AU
                    for r in range ( gradients3.rows ):
                        for ( c, v ) in enumerate ( [ float ( t ) for t in next ( atFile ).split ( ) ] ):
                            gradients3[r,c] = v   
                    gradients3.Scale ( - 1.0 )
                elif line.startswith ( "forces_ext_charges" ):
                    gradientsMM = Coordinates3.WithExtent ( int ( line.split ( "," )[-1] ) )
                    for r in range ( gradientsMM.rows ):
                        for ( c, v ) in enumerate ( [ float ( t ) for t in next ( atFile ).split ( ) ] ):
                            gradientsMM[r,c] = v
                    gradientsMM.Scale ( - Units.Length_Angstroms_To_Bohrs * Units.Energy_Hartrees_To_Kilojoules_Per_Mole ) # . Convert to pDynamo units.
                    dftbOutputData["MM Gradients"] = gradientsMM
                elif line.startswith ( "mermin_energy" ):
                    dftbOutputData["Energy"] = float ( next ( atFile ) )
                elif line.startswith ( "stress" ):
                    data = Matrix33.Null ( )
                    for r in range ( 3 ):
                        for ( c, v ) in enumerate ( [ float ( t ) for t in next ( atFile ).split ( ) ] ):
                            data[r,c] = v
                    dftbOutputData["Stress"] = data
            atFile.close ( )
            return True
        except StopIteration:
            atFile.close ( )
            return True
        except:
            return False

    def ReadDetailedFile ( self, target, dftbOutputData ):
        """Read the detailed output file."""
        state = getattr ( target, self.__class__._stateName )
        try:
            isSpinPolarized = False
            n               = len ( state.atomicNumbers )
            orbitalTag      = ""
            scratch         = { "Is Successful" : False }
            dFile           = open ( state.paths["Detailed"], "r" )
            while True:
                try:
                    line = next ( dFile ).strip ( )
                    # . Atomic charges.
                    if line == "Atomic gross charges (e)":
                        data = Array.WithExtent ( n )
                        line = next ( dFile )
                        for i in range ( n ):
                            data[i] = float ( next ( dFile ).split ( )[-1] )
                        scratch["Charges"] = data
                    # . Component/spin down and up flags.
                    elif line.startswith ( "COMPONENT = down" ) or line.startswith ( "Spin  down" ):
                        isSpinPolarized = True
                        orbitalTag      = "Beta "
                    elif line.startswith ( "COMPONENT = up" ) or line.startswith ( "Spin  up" ):
                        isSpinPolarized = True
                        orbitalTag      = "Alpha "
                    # . Dipole - Debyes using pDynamo's own conversion factor.
                    elif line.startswith ( "Dipole moment" ) and line.endswith ( "au" ):
                        data = Vector3.Null ( )
                        for ( i, word ) in enumerate ( line.split ( )[2:5] ):
                            data[i] = Units.Dipole_Atomic_Units_To_Debyes * float ( word )
                        scratch["Dipole"] = data
                    # . Fermi energy - atomic units.
                    elif line.startswith ( "Fermi level:" ):
                        scratch[orbitalTag + "Fermi Energy"] = float ( line.split ( )[2] )
                    # . Lattice derivatives - atomic units.
                    elif line.startswith ( "Total lattice derivs" ):
                        data = Matrix33.Null ( )
                        for r in range ( 3 ):
                            for ( c, v ) in enumerate ( [ float ( t ) for t in next ( dFile ).split ( ) ] ):
                                data[r,c] = v
                        scratch["Lattice Derivatives"] = data
                    # . Non-SCC calculation.
                    elif line.find ( "Non-SCC calculation" ) >= 0:
                        scratch["Is Successful"] = True
                    # . Orbital energies (atomic units) and occupancies.
                    elif line.startswith ( "Eigenvalues (H) and fillings (e)" ):
                        energies    = []
                        occupancies = []
                        while True:
                            tokens = next ( dFile ).split ( )
                            if len ( tokens ) == 3:
                                energies.append    ( float ( tokens[1] ) )
                                occupancies.append ( float ( tokens[2] ) )
                            else: break
                        scratch[orbitalTag + "Orbital Energies"   ] = energies
                        scratch[orbitalTag + "Orbital Occupancies"] = occupancies
                        orbitalTag = ""
                    # . Pressure - atomic units.
                    elif line.startswith ( "Pressure" ):
                        scratch["Pressure"] = float ( line.split ( )[1] )
                    # . SCC calculation converged.
                    elif line.find ( "SCC converged" ) >= 0:
                        scratch["Is Successful"] = True
                    # . SCC calculation not converged.
                    elif line.find ( "SCC is NOT converged" ) >= 0:
                        scratch["Is Successful"] = False
                    # . Spins.
                    elif line.startswith ( "Atom populations (down)" ) and isSpinPolarized:
                        data = scratch["Spins"]
                        line = next ( dFile )
                        for i in range ( n ):
                            data[i] -= float ( next ( dFile ).split ( )[-1] )
                    elif line.startswith ( "Atom populations (up)"   ) and isSpinPolarized:
                        data = Array.WithExtent ( n )
                        line = next ( dFile )
                        for i in range ( n ):
                            data[i] = float ( next ( dFile ).split ( )[-1] )
                        scratch["Spins"] = data
                except StopIteration:
                    break
            dFile.close ( )
            dftbOutputData.update ( scratch )
            return scratch["Is Successful"]
        except Exception as e:
            return False

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModel, self ).SummaryItems ( )
        items.append ( ( "Extended Input", repr ( self.extendedInput is not None ) ) )
        return items

    def WriteInputFile ( self, target, doGradients, doQCMM, coordinates3, path = None ):
        """Write an input file."""
        # . Initialization.
        state = getattr ( target, self.__class__._stateName )
        if doGradients: gTag = "Yes"
        else:           gTag = "No"
        if path is None: path = state.paths["Input"] # . This option is for testing.
        inFile = open ( path, "w" )
        # . Element data.
        uniqueNumbers = sorted ( set ( state.atomicNumbers ) )
        uniqueSymbols = [ PeriodicTable[n].symbol                      for n in uniqueNumbers ]
        uniqueQNs     = [ PeriodicTable[n].maximumOrbitalQuantumNumber for n in uniqueNumbers ]
        # . Symmetry.
        sp           = target.symmetryParameters
        if sp is None: sTag = "C"
        else:          sTag = "S"
        # . Geometry block (Angstroms).
        inFile.write ( "Geometry = GenFormat {{\n{:d} {:s}\n".format ( len ( state.atomicNumbers ), sTag ) )
        inFile.write ( "  {:s}\n".format ( " ".join ( uniqueSymbols ) ) )
        for ( i, p ) in enumerate ( state.atomicNumbers ):
            inFile.write ( "  {:6d} {:6d} {:25.15f}{:25.15f}{:25.15f}\n".format ( i+1, uniqueNumbers.index ( p )+1 ,
                                                                                               coordinates3[i,0]   ,
                                                                                               coordinates3[i,1]   ,
                                                                                               coordinates3[i,2] ) )
        if sp is not None:
            inFile.write ( "  0.0 0.0 0.0\n" )
            for i in range ( 3 ):
                inFile.write ( "  {:25.15f}{:25.15f}{:25.15f}\n".format ( sp.H[0,i], sp.H[1,i], sp.H[2,i] ) )
        inFile.write ( "}\n" )
        # . Start Hamiltonian block - basic options only.
        inFile.write ( "Hamiltonian = DFTB {\n" )
        inFile.write ( "  Charge             = {:d}\n".format ( target.electronicState.charge ) )
        inFile.write ( "  Eigensolver        = QR {}\n" )
        inFile.write ( "  Filling            = Fermi {{ Temperature [Kelvin] = {:.1f} }}\n".format ( self.fermiTemperature ) )
        # . SCC.
        if self.useSCC:
            if os.path.exists ( state.paths["Charges"] ): inFile.write ( "  ReadInitialCharges = Yes\n" )
            inFile.write ( "  SCC                = Yes\n" )
            inFile.write ( "  MaxSCCIterations   = {:d}\n".format ( self.maximumSCCIterations ) )
            inFile.write ( "  SCCTolerance       = {:g}\n".format ( self.sccTolerance         ) )
        else:
            inFile.write ( "  SCC = No\n" )
        # . Maximum angular momenta.
        inFile.write ( "  MaxAngularMomentum = {\n" )
        for ( s, q ) in zip ( uniqueSymbols, uniqueQNs ): inFile.write ( "    {:2s} = \"{:s}\"\n".format ( s, _OrbitalQuantumNumberLabels[q] ) )
        inFile.write ( "  }\n" )
        # . Slater-Koster files (type 2 only).
        inFile.write ( "  SlaterKosterFiles = Type2FileNames {\n" )
        inFile.write ( "    Prefix    = \"{:s}\"\n    Separator = \"-\"\n    Suffix    = \".skf\"\n  }}\n".format ( os.path.join ( self.skfPath, "" ) ) )
        # . User specified options using extended format.
        if self.extendedInput is not None: inFile.write ( "{:s}\n".format ( self.extendedInput ) )
        # . QC/MM - there are assumptions here about the QC/MM model.
        if doQCMM: 
            target        = target
            nRecords      = len ( target.mmState.pureMMAtoms )
            bpCoordinates = target.scratch.Get ( "bpCoordinates3", None )
            if bpCoordinates is not None: nRecords += bpCoordinates.rows
            inFile.write ( _ElectricFieldInput.format ( self.gaussianBlurWidth, nRecords, _PointChargeFile ) )#state.paths["PointCharge"] ) )
        # . Symmetry.
        if sp is not None: inFile.write ( "  KPointsAndWeights = { 0.0 0.0 0.0 1.0 }\n" )
        # . Unrestricted calculation - this is incomplete as spin constants need to be specified!
        if target.electronicState.multiplicity != 1:
            inFile.write ( "  SpinPolarisation = Colinear {{\n    UnpairedElectrons = {:d}\n  }}\n".format ( target.electronicState.multiplicity - 1 ) )
        # . Finish Hamiltonian block.
        inFile.write ( "}\n" )
        # . Other blocks.
        inFile.write ( "Driver        = {}\n" )
        inFile.write ( "Analysis      = {{\n  CalculateForces = {:s}\n}}\n".format ( gTag ) )
        inFile.write ( "Options       = {\n  WriteAutotestTag   = Yes\n  }\n" )
        inFile.write ( "ParserOptions = {\n  ParserVersion = 5\n}\n" )
        # . Finish up.
        inFile.close ( )

    @property
    def command ( self ):
        """Get the command to execute the program."""
        command = self.__dict__.get ( "_command", None )
        if command is None:
            command = os.getenv ( _DFTBCommand )
            # . Command must point to an executable file.
            if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
                raise NotInstalledError ( "DFTB executable not found." )
            else:
                self.__dict__["_command"] = command
        return command

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass


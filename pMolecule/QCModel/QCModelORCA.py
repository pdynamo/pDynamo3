"""The ORCA QC model."""

import glob, math, os, os.path, subprocess, re

from  pCore                     import logFile           , \
                                       LogFileActive     , \
                                       NotInstalledError
from  pScientific               import PeriodicTable     , \
                                       Units
from  pScientific.Arrays        import Array
from  pScientific.Geometry3     import Coordinates3      , \
                                       Vector3
from  pScientific.RandomNumbers import RandomString
from .QCDefinitions             import ChargeModel
from .QCModel                   import QCModel           , \
                                       QCModelState
from .QCModelError              import QCModelError

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Default error suffix.
_DefaultErrorPrefix = "error_"

# . Default job name.
_DefaultJobName = "orcaJob"

# . Command environment variable.
_ORCACommand = "PDYNAMO3_ORCACOMMAND"

# . Scratch directory.
_ORCAScratch = os.path.join ( os.getenv ( "PDYNAMO3_SCRATCH" ), "orcaScratch" )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelORCAState ( QCModelState ):
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
                jobFiles = glob.glob ( os.path.join ( self.paths["Glob"] + "*" ) )
                for jobFile in jobFiles: os.remove ( jobFile )
                scratch  = self.paths.get ( "Scratch", None )
                if scratch is not None: os.rmdir ( scratch ) # . Only deleted if random.
            except:
                pass

    def DeterminePaths ( self, scratch, deleteJobFiles = True, randomJob = False, randomScratch = False ):
        """Determine the paths needed by an ORCA job."""
        paths = {}
        if randomJob: job = RandomString ( )
        else:            job = _DefaultJobName
        if randomScratch:
            scratch          = os.path.join ( scratch, RandomString ( ) )
            paths["Scratch"] = scratch # . Only set if random.
        if not os.path.exists ( scratch ): os.mkdir ( scratch )
        jobRoot       = os.path.join ( scratch, job )
        paths["Glob"] = jobRoot
        for ( key, ext ) in ( ( "EnGrad" , "engrad" ) ,
                              ( "Input"  , "inp"    ) ,
                              ( "Output" , "log"    ) ,
                              ( "PCGrad" , "pcgrad" ) ,
                              ( "PC"     , "pc"     ) ):
            paths[key] = "{:s}.{:s}".format ( jobRoot, ext )
        # . Finish up.
        self.deleteJobFiles = deleteJobFiles
        self.paths          = paths

    def SaveErrorFiles ( self, message ):
        """Save the input and output files for inspection if there is an error."""
        for key in ( "EnGrad" , "Input" , "Output" , "PCGrad" , "PC" ):
            path = self.paths[key]
            if os.path.exists ( path ):
                ( head, tail ) = os.path.split ( path )
                os.rename ( path, os.path.join ( head, _DefaultErrorPrefix + tail ) )
        ( head, tail ) = os.path.split ( self.paths["Glob"] )
        raise QCModelError ( message + "\nCheck the files \"{:s}*\".".format ( os.path.join ( head, _DefaultErrorPrefix + tail ) ) )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelORCA ( QCModel ):
    """The ORCA QC model class."""

    _attributable = dict ( QCModel._attributable )
    _classLabel   = "ORCA QC Model"
    _stateObject  = QCModelORCAState
    _summarizable = dict ( QCModel._summarizable )
    _attributable.update ( { "deleteJobFiles" : True                           ,
                             "keywords"       : [ "HF", "6-31G*", "TIGHTSCF" ] ,
                             "randomJob"      : False                          ,
                             "randomScratch"  : False                          ,
                             "scratch"        : _ORCAScratch                   } )
    _summarizable.update ( { "deleteJobFiles" : "Delete Job Files"             ,
                             "randomJob"      : "Random Job"                   ,
                             "randomScratch"  : "Random Scratch"               } )

    def AtomicCharges ( self, target, chargeModel = ChargeModel.Mulliken ):
        """Atomic charges."""
        source = target.scratch.orcaOutputData
        if   chargeModel is ChargeModel.CHelpG : return source.get ( "Chelpg Charges"   , None )
        elif chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Charges"  , None )
        else:                                    return source.get ( "Mulliken Charges" , None )

    def AtomicSpins ( self, target, chargeModel = "Mulliken" ):
        """Atomic spins."""
        source = target.scratch.orcaOutputData
        if chargeModel is ChargeModel.Loewdin: return source.get ( "Loewdin Spins"  , None )
        else:                                  return source.get ( "Mulliken Spins" , None )

    def BondOrders ( self, target, chargeModel = None ):
        """Bond Orders - Mayer only."""
        return target.scratch.orcaOutputData.get ( "Mayer Bond Orders", None )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelORCA, self ).BuildModel ( target, qcSelection = qcSelection )
        state.DeterminePaths ( self.scratch                         ,
                               deleteJobFiles = self.deleteJobFiles ,
                               randomJob      = self.randomJob      ,
                               randomScratch  = self.randomScratch  )
        return state

    def DipoleMoment ( self, target, center = None ):
        """Dipole Moment."""
        return target.scratch.orcaOutputData.get ( "Dipole", None )

    def Energy ( self, target ):
        """Calculate the quantum chemical energy."""
        doGradients    = target.scratch.doGradients
        orcaOutputData = {}
        state          = getattr ( target, self.__class__._stateName )
        target.scratch.orcaOutputData = orcaOutputData
        self.WriteInputFile ( target, doGradients, ( target.nbModel is not None ), target.scratch.qcCoordinates3AU )
        isOK = self.Execute ( state )
        if not isOK: state.SaveErrorFiles ( "Error executing program." )
        if doGradients:
            isOK = self.ReadEngradFile ( target, orcaOutputData, target.scratch.qcGradients3AU )
            if not isOK: state.SaveErrorFiles ( "Error reading engrad file." )
        isOK = self.ReadOutputFile ( target, orcaOutputData ) # . Returns whether converged or an error.
        if not isOK: state.SaveErrorFiles ( "Error reading output file." )
        target.scratch.energyTerms["ORCA QC"] = ( orcaOutputData["Energy"] * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def Execute ( self, state ):
        """Execute the job."""
        try:
            outFile = open ( state.paths["Output"], "w" )
            subprocess.check_call ( [ self.command, state.paths["Input"] ], stderr = outFile, stdout = outFile )
            outFile.close ( )
            return True
        except:
            return False

    # . Alpha/beta?
    def OrbitalEnergies ( self, target ):
        """Orbital energies and HOMO and LUMO indices."""
        return ( target.scratch.orcaOutputData.get ( "Orbital Energies", None ) ,
                 target.scratch.orcaOutputData.get ( "HOMO"            , -1   ) ,
                 target.scratch.orcaOutputData.get ( "LUMO"            , -1   ) )

    def ReadEngradFile ( self, target, orcaOutputData, gradients3 ):
        """Read an engrad file."""
        # . The energy and gradients are in atomic units.
        state = getattr ( target, self.__class__._stateName )
        try:
            egFile = open ( state.paths["EnGrad"], "r" )
            for i in range ( 7 ): next ( egFile )     # . The number of atoms section and the energy header.
            orcaOutputData["Energy"] = float ( next ( egFile ) )
            for i in range ( 3 ): ( next ( egFile ) ) # . The gradients header.
            for i in range ( len ( state.atomicNumbers ) ):
                for j in range ( 3 ):
                    gradients3[i,j] = float ( ( next ( egFile ) ) )
            egFile.close ( )
            return True
        except:
            return False

    def ReadOutputFile ( self, target, orcaOutputData ):
        """Read an output file."""
        state  = getattr ( target, self.__class__._stateName )
        try:
            n       = len ( state.atomicNumbers )
            scratch = { "Is Converged" : False }
            outFile = open ( state.paths["Output"], "r" )
            while True:
                try:
                    line = next ( outFile ).strip ( )
                    # . CHELPG charges.
                    if line == "Chelpg Charges":
                        data = Array.WithExtent ( n )
                        line = next ( outFile )
                        for i in range ( n ):
                            words   = next ( outFile ).split ( ":", 1 )
                            data[i] = float ( words[-1] )
                        scratch["CHELPG Charges"] = data
                    # . Convergence OK.
                    elif line.find ( "SCF CONVERGED AFTER" ) >= 0:
                        words                   = line.split ( )
                        scratch["Cycles"]       = int ( words[-3] )
                        scratch["Is Converged"] = True
                    # . Convergence not OK.
                    elif line.find ( "SCF NOT CONVERGED AFTER" ) >= 0:
                        words             = line.split ( )
                        scratch["Cycles"] = int ( words[-3] )
                    # . Dipole.
                    elif line.startswith ( "Total Dipole Moment" ):
                        data  = Vector3.Null ( )
                        words = line.split ( )
                        for ( i, word ) in enumerate ( words[-3:] ):
                            data[i] = Units.Dipole_Atomic_Units_To_Debyes * float ( word )
                        scratch["Dipole"] = data
                    # . Energy.
                    elif line.startswith ( "FINAL SINGLE POINT ENERGY" ):
                        scratch["Energy"] = float ( line.split ( )[-1] )
                    # . Loewdin charges.
                    elif line == "LOEWDIN ATOMIC CHARGES":
                        data = Array.WithExtent ( n )
                        next ( outFile )
                        for i in range ( n ):
                            words   = next ( outFile ).split ( ":", 1 )
                            data[i] = float ( words[-1] )
                        scratch["Loewdin Charges"] = data 
                    # . Loewdin charges and spin densities.
                    elif line.startswith ( "LOEWDIN ATOMIC CHARGES AND SPIN " ):
                        data = Array.WithExtent ( n )
                        datb = Array.WithExtent ( n )
                        next ( outFile )
                        for i in range ( n ):
                            words = next ( outFile ).split ( )
                            data[i] = float ( words[-2] )
                            datb[i] = float ( words[-1] )
                        scratch["Loewdin Charges"] = data 
                        scratch["Loewdin Spins"  ] = datb
                    # . Mayer bond orders.
                    elif line.startswith ( "Mayer bond orders larger than" ):
                        data = []
                        while True:
                            line = next ( outFile )
                            atompairs=re.findall(r'(\d+)-\w+\s*,\s*(\d+)-\w+', line)
                            bondOrders=p=re.findall(r'(\d\.\d+)', line)
                            if atompairs:
                                for (pair, bond) in zip(atompairs, bondOrders):
                                    i=int(pair[0])
                                    j=int(pair[1])
                                    data.append ( ( i, j, float ( bond ) ) )
                            else:
                                break
                        scratch["Mayer Bond Orders"] = data
                    # . Mulliken charges.
                    elif line == "MULLIKEN ATOMIC CHARGES":
                        data = Array.WithExtent ( n )
                        line = next ( outFile )
                        for i in range ( n ):
                            words = next ( outFile ).split ( ":", 1 )
                            data[i] = float ( words[-1] )
                        scratch["Mulliken Charges"] = data
                    # . Mulliken charges and spin densities.
                    elif line.startswith ( "MULLIKEN ATOMIC CHARGES AND SPIN " ):
                        data = Array.WithExtent ( n )
                        datb = Array.WithExtent ( n )
                        next ( outFile )
                        for i in range ( n ):
                            words = next ( outFile ).split ( )
                            data[i] = float ( words[-2] )
                            datb[i] = float ( words[-1] )
                        scratch["Mulliken Charges"] = data 
                        scratch["Mulliken Spins"  ] = datb
                    # . Orbital energies.
                    elif line == "ORBITAL ENERGIES":
                        HOMO = -1
                        LUMO = -1
                        orbitalEnergies = []
                        for i in range ( 3 ): next ( outFile )
                        index = 0
                        while True:
                            tokens = next ( outFile ).split ( )
                            if len ( tokens ) < 3: break
                            else:
                                occupancy = float ( tokens[1] )
                                if ( HOMO == -1 ) and ( LUMO == -1 ) and ( occupancy <= 1.0e-6 ):
                                    HOMO = index - 1
                                    LUMO = index
                                orbitalEnergies.append ( float ( tokens[2] ) )
                                index += 1
                        energies = Array.WithExtent ( len ( orbitalEnergies ) )
                        energies.Set ( 0.0 )
                        for ( i, e ) in enumerate ( orbitalEnergies ): energies[i] = e
                        scratch["HOMO"] = HOMO
                        scratch["LUMO"] = LUMO
                        scratch["Orbital Energies"] = energies
                    # . <S**2>.
                    elif line.startswith ( "Expectation value of <S**2>" ):
                        scratch["Spin Squared"] = float ( line.split ( ":", 1 )[-1] )
                except StopIteration:
                    break
            outFile.close ( )
            orcaOutputData.update ( scratch )
            return scratch["Is Converged"]
        except Exception as e:
            return False

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelORCA, self ).SummaryItems ( )
        if self.keywords is not None:
            n = len ( self.keywords )
            items.append ( ( "Keywords", "{:s}".format ( "/".join ( self.keywords[0:min(2,n)] ) ) ) )
        return items

    # . Coordinates are written in atomic units (keywords must contain "Bohrs").
    def WriteInputFile ( self, target, doGradients, doQCMM, coordinates3 ):
        """Write an input file."""
        state  = getattr ( target, self.__class__._stateName )
        inFile = open ( state.paths["Input"], "w" )
        inFile.write ( "#\n" )
        inFile.write ( "# ORCA Job.\n" )
        inFile.write ( "#\n" )
        if doQCMM     : inFile.write ( '%pointcharges "' + state.paths["PC"] + '"\n' )
        if doGradients: mode = "ENGRAD"
        else          : mode = "ENERGY"
        inFile.write ( "! " + mode + " BOHRS " + " ".join ( self.keywords ) + "\n" )
        inFile.write ( "* xyz {:d} {:d}\n".format ( target.electronicState.charge, target.electronicState.multiplicity ) )
        for ( i, n ) in enumerate ( state.atomicNumbers ):
            inFile.write ( "{:<12s}{:20.10f}{:20.10f}{:20.10f}\n".format ( PeriodicTable.Symbol ( n ) ,
                                                                           coordinates3[i,0]          ,
                                                                           coordinates3[i,1]          ,
                                                                           coordinates3[i,2]        ) )
        inFile.write ( "*\n" )
        inFile.close ( )

    @property
    def command ( self ):
        """Get the command to execute the program."""
        command = self.__dict__.get ( "_command", None )
        if command is None:
            command = os.getenv ( _ORCACommand )
            # . Command must point to an executable file.
            if  ( command is None ) or not ( os.path.isfile ( command ) and os.access ( command, os.X_OK ) ):
                raise NotInstalledError ( "ORCA executable not found." )
            else:
                self.__dict__["_command"] = command
        return command

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass


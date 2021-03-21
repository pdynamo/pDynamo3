"""MCModelGMCT is a class for the calculation of protonation state probabilities in GMCT."""

import os, subprocess

from  pCore                 import logFile                      , \
                                   LogFileActive
from  pScientific           import Constants                    , \
                                   Units
from .CEModelError          import CEModelError
from .GMCTOutputFileReader  import GMCTOutputFileReader
from .InputFileWriter       import WriteInputFile

_DefaultPathGMCT           = "/usr/local/bin"
_DefaultDoubleFlip         = 2.
_DefaultTripleFlip         = 3.
_DefaultProductionScans    = 20000
_DefaultEquilibrationScans = 500

_DefaultSetupGMCT = """
blab        1
nconfflip   10
tlimit      {:f}
itraj       0
nmcfull     {:d}
temp        {:f}
icorr       0
limit       {:f}
nmcequi     {:d}
nmu         1
mu          {:f}  {:f}  0.0  0  0
"""


class MCModelGMCT:
    """A parent class."""

    _attributable = {
        "owningModel"         :  None                       ,
        "doubleFlip"          :  _DefaultDoubleFlip         ,
        "tripleFlip"          :  _DefaultTripleFlip         ,
        "productionScans"     :  _DefaultProductionScans    ,
        "equilibrationScans"  :  _DefaultEquilibrationScans ,
        "pathGMCT"            :  _DefaultPathGMCT           ,
            }

    def __init__ (self, *arguments, **options):
        """Constructor."""
        for (key, value) in self.__class__._attributable.items (): setattr (self, key, value)
        for (key, value) in                 options.items (): setattr (self, key, value)


    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ): log.SummaryOfItems ( self.SummaryItems ( ), title = "GMCT Monte Carlo Sampling Model" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Equilibration scans"    , "{:d}".format   ( self.equilibrationScans ) ) ,
                 ( "Production scans"       , "{:d}".format   ( self.productionScans    ) ) ,
                 ( "Limit for double moves" , "{:.1f}".format ( self.doubleFlip         ) ) ,
                 ( "Limit for triple moves" , "{:.1f}".format ( self.tripleFlip         ) ) ]

    def PrintPairs (self, log=logFile):
        """Compatibility method for MCModelDefault."""
        pass


    def Initialize (self, meadModel):
        """Link the Monte Carlo model to its owning energy model."""
        self.owningModel = meadModel


    def CalculateOwnerProbabilities (self, pH=7.0, dryRun=False, logFrequency=None, log=logFile):
        """Calculate probabilities of the owning model."""
        owner = self.owningModel
        if not owner.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")

        sites       = None
        project     = "job"
        potential   = - Constants.Molar_Gas * Units.Energy_Kilojoules_Per_Mole_To_Kilocalories_Per_Mole * owner.temperature * Constants.Ln10 * pH
        fileContent = _DefaultSetupGMCT.format ( self.tripleFlip, self.productionScans, owner.temperature, self.doubleFlip, self.equilibrationScans, potential, potential )

        # Prepare input files and directories for GMCT
        dirConf   = os.path.join ( owner.pathScratch , "gmct" , "conf"                 )
        dirCalc   = os.path.join ( owner.pathScratch , "gmct" , "{:s}".format ( pH   ) )
        fileGint  = os.path.join ( dirConf ,           "{:s}.gint".format  ( project ) )
        fileInter = os.path.join ( dirConf ,           "{:s}.inter".format ( project ) )
        fileConf  = os.path.join ( dirCalc ,           "{:s}.conf".format  ( project ) )
        fileSetup = os.path.join ( dirCalc ,           "{:s}.setup".format ( project ) )
        linkname  = os.path.join ( dirCalc ,           "conf"                          )
        if not os.path.exists ( dirConf   ): os.makedirs      ( dirConf                  )
        if not os.path.exists ( dirCalc   ): os.makedirs      ( dirCalc                  )
        if not os.path.exists ( fileGint  ): owner.WriteGintr ( fileGint  , precision=8  )
        if not os.path.exists ( fileInter ): owner.WriteW     ( fileInter , precision=8  )
        if not os.path.exists ( fileConf  ): WriteInputFile   ( fileConf  , ["conf  0.0  0.0  0.0\n"] )
        if not os.path.exists ( fileSetup ): WriteInputFile   ( fileSetup , fileContent  )
        if not os.path.exists ( linkname  ): os.symlink       ( "../conf" , linkname     )

        if not dryRun:
            output = os.path.join (dirCalc, "{:s}.gmct-out".format ( project ) )
            error  = os.path.join (dirCalc, "{:s}.gmct-err".format ( project ) )

            if os.path.exists (os.path.join (dirCalc, output)):
                pass
            else:
                command = [os.path.join (self.pathGMCT, "gmct"), project]
                try:
                    out = open (output, "w")
                    err = open (error,  "w")
                    subprocess.check_call (command, stderr=err, stdout=out, cwd=dirCalc)
                    out.close ()
                    err.close ()
                except:
                    raise CEModelError ("Failed running command: {:s}".format ( " ".join (command)) )

            # Read probabilities from the output file
            reader = GMCTOutputFileReader.FromPath (output)
            reader.Parse (temperature=owner.temperature)

            # Copy probabilities from the reader to the owner
            for site in owner.sites:
                for instance in site.instances:
                    key                  = "conf_{:s}_{:s}{:d}_{:s}".format ( site.segName, site.resName, site.resSerial, instance.label )
                    probability          = reader.probabilities[key][0]
                    instance.probability = probability


#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass

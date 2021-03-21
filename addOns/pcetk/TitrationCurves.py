"""TitrationCurves is a class for calculating titration curves.

CurveThread is a class for running parallel calculations of titration curves."""

import os, threading

from  collections     import defaultdict
from  pCore           import Align         , \
                             logFile       , \
                             LogFileActive
from .CEModelError    import CEModelError
from .InputFileWriter import WriteInputFile
from .MCModelGMCT     import MCModelGMCT

_DefaultDirectory  = "curves"
_DefaultSampling   =   .5
_DefaultStart      =  0.
_DefaultStop       = 14.


class CurveThread (threading.Thread):
    """Calculate each pH-step in a separate thread."""

    def __init__ (self, curves, pH):
        threading.Thread.__init__ (self)
        self.curves = curves
        self.sites  = None
        self.pH     = pH

    def run (self):
        """The method that runs the calculations."""
        curves     = self.curves
        model      = curves.owningModel
        self.sites = model.CalculateProbabilities (pH=self.pH, log=None, isCalculateCurves=True, unfolded=curves.unfolded)


#-------------------------------------------------------------------------------
class TitrationCurves:
    """Titration curves."""

    _attributable = {
        "curveSampling"  :  _DefaultSampling  ,
        "curveStart"     :  _DefaultStart     ,
        "curveStop"      :  _DefaultStop      ,
        "unfolded"       :  False             ,
            }

    def __init__ (self, meadModel, log=logFile, *arguments, **options):
        """Constructor."""
        for (key, value) in self.__class__._attributable.items (): setattr (self, key, value)
        for (key, value) in                 options.items (): setattr (self, key, value)

        if not meadModel.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")
        self.owningModel   =  meadModel
        self.nsteps        =  int ((self.curveStop - self.curveStart) / self.curveSampling + 1)
        self.steps         =  None
        self.halves        =  None
        self.isHalves      =  False
        self.isCalculated  =  False


    #===============================================================================
    def _ProbabilitiesSave (self):
        """Backup the probabilities existing in the model."""
        meadModel     = self.owningModel
        energyModel   = meadModel.energyModel
        probabilities = [0.] * meadModel.ninstances

        for site in meadModel.sites:
            for instance in site.instances:
                probabilities[instance._instIndexGlobal] = energyModel.GetProbability (instance._instIndexGlobal)
        return probabilities


    def _ProbabilitiesRestore (self, probabilities):
        """Restore the original probabilities to the model."""
        meadModel   = self.owningModel
        energyModel = meadModel.energyModel

        for site in meadModel.sites:
            for instance in site.instances:
                energyModel.SetProbability (instance._instIndexGlobal, probabilities[instance._instIndexGlobal])


    #===============================================================================
    def CalculateCurves (self, forceSerial=True, printTable=False, log=logFile):
        """Calculate titration curves."""
        if not self.isCalculated:
            owner   = self.owningModel
            restore = False
            if owner.isProbability:
                probabilities  = self._ProbabilitiesSave ()
                restore        = True
            nthreads  =  owner.nthreads
            steps     =  []
            tab       =  None

            if LogFileActive (log):
                if nthreads < 2 or forceSerial:
                    log.Paragraph ( "Starting serial run." )
                else:
                    log.Paragraph ( "Starting parallel run on {:d} CPUs.".format ( nthreads ) )

                if printTable:
                    tab = log.GetTable (columns=[10, 10])
                    tab.Start ()
                    tab.Heading ("Step")
                    tab.Heading ("pH")

            # Serial run?
            if nthreads < 2 or forceSerial:
                for step in range (self.nsteps):
                    sites = owner.CalculateProbabilities (pH=(self.curveStart + step * self.curveSampling), log=None, isCalculateCurves=True, unfolded=self.unfolded)
                    steps.append (sites)
                    if tab:
                        tab.Entry ( "{:d}".format   ( step ) )
                        tab.Entry ( "{:.2f}".format ( pH   ) )
            # Parallel run?
            else:
                limit   = nthreads - 1
                batches = []
                batch   = []
                for step in range (self.nsteps):
                    batch.append (CurveThread (self, self.curveStart + step * self.curveSampling))
                    if len (batch) > limit:
                        batches.append (batch)
                        batch = []
                if batch:
                    batches.append (batch)

                # If GMCT is to be used, first perform a dry run in serial mode to create directories and files
                owner   = self.owningModel
                sampler = owner.sampler
                if isinstance (sampler, MCModelGMCT):
                    for step in range (self.nsteps):
                        sampler.CalculateOwnerProbabilities (pH=(self.curveStart + step * self.curveSampling), dryRun=True, log=None)

                step = 0
                for batch in batches:
                    for thread in batch: thread.start ()
                    for thread in batch: thread.join ()
                    for thread in batch:
                        steps.append (thread.sites)
                        if tab:
                            tab.Entry ( "{:d}".format   ( step ) )
                            tab.Entry ( "{:.2f}".format ( self.curveStart + step * self.curveSampling ) )
                        step = step + 1

            if LogFileActive (log):
                if tab:
                    tab.Stop ()
                log.Paragraph ( "Calculating titration curves complete." )

            # Finalize
            if restore : self._ProbabilitiesRestore (probabilities)
            else       : owner.isProbability = False
            self.isCalculated = True
            self.steps        = steps


    #===============================================================================
    def WriteCurves (self, directory=_DefaultDirectory, log=logFile):
        """Write calculated curves to a directory."""
        if self.isCalculated:
            # Initialize output directory
            if not os.path.exists (directory):
                try:
                    os.makedirs (directory)
                except:
                    raise CEModelError ( "Cannot create directory {:s}".format ( directory ) )

            # For each instance of each site, write a curve file
            meadModel = self.owningModel
            phdata    = self.steps

            for site in meadModel.sites:
                for instance in site.instances:
                    # Collect instance data
                    lines = []
                    for step in range (self.nsteps):
                        lines.append ( "{:f} {:f}\n".format ( self.curveStart + step * self.curveSampling, phdata[step][site.siteIndex][instance.instIndex] ) )
                    # Write instance data
                    filename = os.path.join ( directory, "{:s}_{:s}.dat".format ( site.label, instance.label ) )
                    WriteInputFile (filename, lines)

            if LogFileActive (log):
                log.Paragraph ( "Writing curve files complete." )


    #===============================================================================
    def CalculateHalfpKs (self):
        """Scan the calculated curves to find pK1/2 values."""
        if self.isCalculated:
            owner   =  self.owningModel
            data    =  self.steps
            nsteps  =  self.nsteps
            halves  =  []
            for site in owner.sites:
                instances = []
                for instance in site.instances:
                    pa  = data[0][site.siteIndex][instance.instIndex]
                    pKs = []
                    for step in range (1, nsteps):
                        pb = data[step][site.siteIndex][instance.instIndex]
                        if ((pa - .5) * (pb - .5)) < 0.:
                            a  = self.curveStart + (step - 1) * self.curveSampling
                            b  = self.curveStart + (step    ) * self.curveSampling
                            pK = a + (.5 - pa) * (b - a) / (pb - pa)
                            pKs.append (pK)
                        pa = pb
                    instances.append (pKs)
                halves.append (instances)

            self.halves   = halves
            self.isHalves = True


    #===============================================================================
    def _GetEntry ( self, site, decimalPlaces = 2 ):
        entry   = []
        fFormat = "{:" + ".{:d}".format ( decimalPlaces ) + "f}"
        if self.isHalves:
            for instance in site.instances:
                pKs     = self.halves[site.siteIndex][instance.instIndex]
                nvalues = len ( pKs )
                entry.append ( instance.label )
                if nvalues < 1:
                    entry.append ( "n/a" )
                else:
                    for pK in pKs:
                        entry.append ( fFormat.format ( pK ) )
        return entry

    #===============================================================================
    def PrintHalfpKs (self, decimalPlaces=2, sortSites=False, log=logFile):
        """Print pK1/2 values."""
        if LogFileActive (log):
            if self.isCalculated:
                if not self.isHalves: self.CalculateHalfpKs ( )
                table   = []
                lengths = defaultdict ( list )
                for site in self.owningModel.sites:
                    items = self._GetEntry ( site, decimalPlaces )
                    for ( i, item ) in enumerate ( items ):
                        lengths[i].append ( len ( item ) )
                    table.append ( [ site.segName, site.resName, site.resSerial, items ] )
                if sortSites: table.sort ( key = lambda k: ( k[0], k[1], k[2] ) )
                columns = [ 6, 6, 6 ]
                for i in sorted ( lengths.keys ( ) ): columns.append ( max ( max ( lengths[i] ), 7 ) )
                nItems  = len ( columns ) - 3
                tab = log.GetTable ( columns = columns )
                tab.Start ( )
                tab.Heading ( "Site", columnSpan = 3 )
                tab.Heading ( "pK1/2 values of instances", columnSpan = nItems )
                for ( segName, resName, resSerial, items ) in table:
                    tab.Entry ( segName )
                    tab.Entry ( resName )
                    tab.Entry ( "{:d}".format ( resSerial ) )
                    for item in items: tab.Entry ( item, align = Align.Right )
                    if len ( items ) < nItems: tab.EndRow ( )
                tab.Stop ()

#===============================================================================
# Testing
#===============================================================================
if __name__ == "__main__": pass

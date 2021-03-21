import os, os.path, time

from  pBabel          import ExportSystem
from  pCore           import logFile           , \
                             LogFileActive     , \
                             NotInstalledError , \
                             Selection         , \
                             YAMLUnpickle
from .CEModel         import CEModel
from .CEModelError    import CEModelError
from .Constants       import YAMLPATHIN
from .InputFileWriter import WriteInputFile
from .InstanceThread  import InstanceThread
from .SiteMEAD        import SiteMEAD

_DEFAULT_THREADS      = 1
_DEFAULT_PATH_SCRATCH = os.getenv ( "PDYNAMO3_SCRATCH"  )
_MEADPath             = "PDYNAMO3_MEADPATH"

class CEModelMEAD (CEModel):
    """A class to represent a continuum electrostatic model based on MEAD."""
    _attributable = {
        "nthreads"             :   _DEFAULT_THREADS       ,
        "pathMEAD"             :   None                   ,
        "pathScratch"          :   _DEFAULT_PATH_SCRATCH  ,
        "deleteJobFiles"       :   False                  ,
        "splitToDirectories"   :   True                   ,
        }
    _attributable.update (CEModel._attributable)


    defaultAttributeNames = {
        "Threads"              :  "nthreads"              ,
        "Split Directories"    :  "splitToDirectories"    ,
        "Delete Job Files"     :  "deleteJobFiles"        ,
        }
    defaultAttributeNames.update (CEModel.defaultAttributeNames)

    @property
    def label (self):
        return "MEAD"

    #-------------------------------------------------------------------------------
    def __init__ (self, system, customFiles=None, log=logFile, **options):
        """Constructor."""
        super (CEModelMEAD, self).__init__ (system, customFiles=customFiles, log=log, **options)

        # . Prepare filenames (do not write actual files)
        generate = (
                ("pathPqrProtein" ,  "protein.pqr"),
                ("pathPqrBack"    ,  "back.pqr"   ),
                ("pathFptSites"   ,  "site.fpt"   ), )
        for attribute, filename in generate:
            setattr (self, attribute, os.path.join (self.pathScratch, filename))

        # . Mead path.
        if self.pathMEAD is None:
            path = os.getenv ( _MEADPath )
            if ( path is None ) or not ( os.path.isdir ( path ) ): raise NotInstalledError ( "MEAD path not found." )
            else: self.pathMEAD = path

    #-------------------------------------------------------------------------------
    def _CreateSite (self, **options):
        """Create a site and its instances specific to the MEAD-based CE model."""
        newSite = SiteMEAD (
            parent            =  self                                       ,
            siteIndex         =  options  [ "siteIndex"        ]   ,
            segName           =  options  [ "segName"          ]   ,
            resName           =  options  [ "resName"          ]   ,
            resSerial         =  options  [ "resSerial"        ]   ,
            siteAtomIndices   =  options  [ "siteAtomIndices"  ]   ,
            modelAtomIndices  =  options  [ "modelAtomIndices" ]   ,)
        # . Initialize instances
        libSite            = options [ "libSite"         ]
        instIndexGlobal    = options [ "instIndexGlobal" ]
        updatedIndexGlobal = newSite._CreateInstances (libSite.instances, instIndexGlobal)

        # . Calculate center of geometry
        newSite._CalculateCenter (centralAtom=libSite.center)

        # . Add the site to the list of sites
        self.sites.append (newSite)

        # . Finalize
        return updatedIndexGlobal


    #-------------------------------------------------------------------------------
    def CalculateElectrostaticEnergies (self, calculateETA=False, asymmetricTolerance=0.05, asymmetricSummary=False, log=logFile):
        if self.isFilesWritten:
            ninstances = self.ninstances
            times      = []
            tab        = None

            if LogFileActive (log):
                if self.nthreads < 2:
                    log.Paragraph ( "Starting serial run." )
                else:
                    log.Paragraph ( "Starting parallel run on {:d} CPUs.".format ( self.nthreads ) )

                heads = [ ("Instance of a site" , 4),
                          ("Gborn_model"        , 0),
                          ("Gback_model"        , 0),
                          ("Gborn_protein"      , 0),
                          ("Gback_protein"      , 0),
                          ("Gmodel"             , 0),
                          ("Gintr"              , 0),]
                columns = [6, 6, 6, 6, 16, 16, 16, 16, 16, 16]
                if calculateETA:
                    heads.append (("ETA", 0))
                    columns.append (16)
                tab = log.GetTable (columns = columns)
                tab.Start ()
                for head, span in heads:
                    if span > 0:
                        tab.Heading (head, columnSpan = span)
                    else:
                        tab.Heading (head)


            if self.nthreads < 2:
                for site in self.sites:
                    for instance in site.instances:
                        time0 = time.time ()
                        instance.CalculateModelCompound (log)
                        instance.CalculateProtein (log)
                        instance.CalculateGintr (log)

                        if calculateETA:
                            times.append (time.time () - time0)
                            averageTimePerInstance = sum (times) / len (times)
                            ninstances = ninstances - 1
                            instance._TableEntry (tab, secondsToCompletion = averageTimePerInstance * ninstances)
                        else:
                            instance._TableEntry (tab)
            else:
                batches = []
                threads = []
                limit   = self.nthreads - 1

                for site in self.sites:
                    for instance in site.instances:
                        if len (threads) > limit:
                            batches.append (threads)
                            threads = []
                        thread = InstanceThread (instance, log)
                        threads.append (thread)

                if threads:
                    batches.append (threads)

                for batch in batches:
                    for thread in batch: thread.start ()
                    for thread in batch: thread.join ()

                    secondsToCompletion = None
                    if calculateETA:
                        # . Collect times of execution
                        nthreads = len (batch)
                        for thread in batch:
                            times.append (thread.time)

                        averageTimePerInstance = sum (times) / len (times) / nthreads
                        ninstances = ninstances - nthreads
                        secondsToCompletion = averageTimePerInstance * ninstances

                    # . Print the results at the end of each batch, otherwise they come in random order
                    for thread in batch:
                        instance = thread.instance
                        instance._TableEntry (tab, secondsToCompletion = secondsToCompletion)
            if tab:
                tab.Stop ( )
                log.Paragraph ( "Calculating electrostatic energies complete." )

            # . Check for symmetricity of the matrix of interactions
            self._CheckIfSymmetric (tolerance=asymmetricTolerance, printSummary=asymmetricSummary, log=log)

            # . Symmetrize interaction energies inside the matrix of interactions
            self.energyModel.SymmetrizeInteractions (log=log)

            # . Finalize
            self.isCalculated = True


    #-------------------------------------------------------------------------------
    def _CheckIfSymmetric (self, tolerance=0.05, printSummary=False, log=logFile):
        """This method is a wrapper for the EnergyModel's CheckIfSymmetric method.

        The wrapper is able to print summaries."""
        isSymmetric, maxDeviation = self.energyModel.CheckIfSymmetric (tolerance=tolerance)

        if LogFileActive (log):
            if isSymmetric:
                log.Paragraph ("Interactions are symmetric within the given tolerance ({:.4f} kcal/mol).".format ( tolerance ) )
            else:
                if not printSummary:
                    log.Paragraph ( "WARNING: Maximum deviation of interactions is {:.4f} kcal/mol.".format ( maxDeviation ) )
                else:
                    heads = [("Instance of a site A" , 4),
                             ("Instance of a site B" , 4),
                             ("Deviation"            , 0),]
                    columns = (7, 7, 7, 7, 7, 7, 7, 7, 12)
                    gaps = ("{:7s}", "{:7s}", "{:7d}", "{:7s}")

                    tab = log.GetTable (columns=columns)
                    tab.Start ()
                    tab.Title ("Deviations of interactions")
                    for head, span in heads:
                        if span > 0:
                            tab.Heading (head, columnSpan=span)
                        else:
                            tab.Heading (head)

                    # . This fragment should be rewritten to work faster
                    report = []
                    for rowSite in self.sites:
                        for rowInstance in rowSite.instances:
                            for columnSite in self.sites:
                                for columnInstance in columnSite.instances:
                                    deviation = self.energyModel.GetDeviation (rowInstance._instIndexGlobal, columnInstance._instIndexGlobal)
                                    if abs (deviation) > tolerance:
                                        report.append ([rowInstance, columnInstance, deviation])

                    for ainstance, binstance, deviation in report:
                        asite = ainstance.parent
                        for gap, content in zip (gaps, (asite.segName, asite.resName, asite.resSerial, ainstance.label)):
                            tab.Entry ( gap.format ( content ) )
                        bsite = binstance.parent
                        for gap, content in zip (gaps, (bsite.segName, bsite.resName, bsite.resSerial, binstance.label)):
                            tab.Entry ( gap.format ( content ) )

                        tab.Entry ( "{:.4f}".format ( deviation ) )
                    tab.Stop ()
        return isSymmetric


    #-------------------------------------------------------------------------------
    def WriteJobFiles (self, log=logFile):
        """Write files: PQR, FPT, OGM and MGM."""
        if self.isInitialized:
            # . Get atomic charges and radii for the system
            system        = self.system
            systemCharges = system.AtomicCharges ()
            systemRadii   = []
            systemTypes   = [ system.mmState.atomTypes[i] for i in system.mmState.atomTypeIndices ]
            radii         = YAMLUnpickle ( os.path.join ( YAMLPATHIN, "radii.yaml" ) )

            for atomType in systemTypes:
                if atomType in radii:
                    radius = radii[atomType]
                else:
                    generalAtomType = "{:s}*".format ( atomType[0] )
                    if generalAtomType in radii:
                        radius = radii[generalAtomType]
                    else:
                        raise CEModelError ( "Cannot find atomic radius for atom type {:s}".format ( atomType  ) )
                systemRadii.append (radius)

            # . Prepare scratch space
            if not os.path.exists (self.pathScratch):
                try:
                    os.makedirs (self.pathScratch)
                except:
                    raise CEModelError ( "Cannot create scratch directory {:s}".format ( self.pathScratch ) )

            # . Create subdirectories, if necessary
            if self.splitToDirectories:
                for site in self.sites:
                    sitePqr   = site.instances[0].sitePqr
                    directory = os.path.dirname (sitePqr)
                    if not os.path.exists (directory):
                        try:
                            os.makedirs (directory)
                        except:
                            raise CEModelError ( "Cannot create directory {:s}".format ( directory ) )

            # . Write PQR, OGM and MGM files of all instances of all sites
            for site in self.sites:
                site._WriteMEADFiles (system, systemCharges, systemRadii)

            # . Write background PQR file
            ExportSystem (self.pathPqrBack, system, selection=Selection (self.backAtomIndices), charges=systemCharges, radii=systemRadii)

            # . Write full-protein PQR file (to be used as eps2set_region)
            ExportSystem (self.pathPqrProtein, system, selection=Selection (self.proteinAtomIndices), charges=systemCharges, radii=systemRadii)

            # . Write FPT-file
            lines = []
            for siteIndex, site in enumerate (self.sites):
                for instanceIndex, instance in enumerate (site.instances):
                    for atomIndex, charge in zip (site.siteAtomIndices, instance.charges):
                        x, y, z = system.coordinates3[atomIndex]
                        line    = "{:d} {:d} {:f} {:f} {:f} {:f}\n".format ( siteIndex, instanceIndex, x, y, z, charge )
                        lines.append (line)
            WriteInputFile (self.pathFptSites, lines)

            self.isFilesWritten = True


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

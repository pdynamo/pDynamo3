import os

from  pCore             import Align            , \
                               logFile          , \
                               LogFileActive    , \
                               Selection        , \
                               YAMLUnpickle
from  pMolecule         import System
from .CEModelError      import CEModelError
from .Constants         import TERM_REMOVE      , \
                               PROTEIN_RESIDUES , \
                               NEXT_RESIDUE_GLY , \
                               NEXT_RESIDUE_PRO , \
                               NEXT_RESIDUE     , \
                               PREV_RESIDUE     , \
                               REMOVE_RESIDUES
from .EnergyModel       import EnergyModel
from .InputFileWriter   import WriteInputFile
from .MCModelDefault    import MCModelDefault
from .MCModelGMCT       import MCModelGMCT
from .TemplatesLibrary  import TemplatesLibrary

_DEFAULT_FOCUSING_STEPS    =  ((121, 2.),  (101, 1.),  (101, .5),  (101, .25))
_DEFAULT_EXCLUDE_SEGMENTS  =  ("WATA", "WATB", "WATC", "WATD", )

_DEFAULT_TEMPERATURE       =  300.      # 300 K
_DEFAULT_IONIC_STRENGTH    =     .1     # 100 mM = 0.1 M
_DEFAULT_EPSILON_WATER     =   80.
_DEFAULT_EPSILON_PROTEIN   =    4.

# CEAtom = collections.namedtuple ("CEAtom", "label  x  y  z  charge  radii")


class CEModel:
    """Base class for continuum electrostatic models.

    This class should not be used directly."""
    _attributable = {
        "temperature"      :   _DEFAULT_TEMPERATURE      ,
        "ionicStrength"    :   _DEFAULT_IONIC_STRENGTH   ,
        "focusingSteps"    :   _DEFAULT_FOCUSING_STEPS   ,
        "epsilonWater"     :   _DEFAULT_EPSILON_WATER    ,
        "epsilonProtein"   :   _DEFAULT_EPSILON_PROTEIN  ,
        }

    defaultAttributeNames = {
        "Temperature"           :   "temperature"     ,
        "Ionic Strength"        :   "ionicStrength"   ,
        "Initialized"           :   "isInitialized"   ,
        "Files Written"         :   "isFilesWritten"  ,
        "Calculated"            :   "isCalculated"    ,
        "Calculated Prob."      :   "isProbability"   ,
        "Focusing Steps"        :   "focusingSteps"   ,
        "Water   Diel. Const."  :   "epsilonWater"    ,
        "Protein Diel. Const."  :   "epsilonProtein"  ,
        }

    @property
    def nsites (self):
        if hasattr (self, "sites"):
            return len (self.sites)
        return 0

    @property
    def ninstances (self):
        if hasattr (self, "sites"):
            counter = 0
            for site in self.sites:
                counter += site.ninstances
            return counter
        return 0

    @property
    def label (self):
        return "Base model"


    #-------------------------------------------------------------------------------
    def __init__ (self, system, customFiles=None, log=logFile, **options):
        """Constructor."""
        # . Perform initial checks
        if not isinstance (system, System):
            raise CEModelError ("Cannot assign system.")

        if "CHARMM" not in system.mmModel.__class__._classLabel:
            raise CEModelError ("Cannot use a non-CHARMM energy model.")


        # . Set system and library
        self.system  = system
        self.library = TemplatesLibrary ( customFiles=customFiles, log=log )

        # . Set attributes
        attributes = self.__class__._attributable
        for (key, val) in attributes.items ():
            setattr (self, key, val)

        for (key, val) in options.items ():
            if key not in attributes:
                raise CEModelError ( "Unknown attribute: {:s}".format ( key ) )
            setattr (self, key, val)

        # . Set status attributes
        for attribute in ("isInitialized", "isFilesWritten", "isCalculated", "isProbability"):
            setattr (self, attribute, False)


    #-------------------------------------------------------------------------------
    def __del__ (self):
        """Deallocation."""
        pass


    #-------------------------------------------------------------------------------
    def Initialize (self, excludeSegments=_DEFAULT_EXCLUDE_SEGMENTS, excludeResidues=None, includeTermini=False, log=logFile):
        """Decompose the system into model compounds, sites and a background charge set.

        |excludeSegments| is a sequence of segment names to exclude from the model, usually segments of water molecules.
        |excludeResidues| is a sequence of three-element sequences (segmentName, residueName, residueSerial).

        It is possible to leave some of the elements blank, for example ("PRTA", "CYS", "") means exclude all cysteines in segment PRTA.
        """
        if not self.isInitialized:
            # . Perform a dry run to calculate the numbers of sites and instances
            totalSites, totalInstances = self._SplitModel (excludeSegments=excludeSegments, excludeResidues=excludeResidues, includeTermini=includeTermini, dryRun=True, log=log)

            # . Allocate arrays of Gmodels, protons, intrinsic energies, interaction energies and probabilities
            self.energyModel = EnergyModel (self, totalSites, totalInstances)

            # . Perform the actual initialization
            self._SplitModel (excludeSegments=excludeSegments, excludeResidues=excludeResidues, includeTermini=includeTermini, dryRun=False, log=None)

            # . Complete the initialization of the energy model
            self.energyModel.Initialize ()

            # . Construct the background set of charges and the protein (to be used as eps2set_region)
            self._SetupBackground ()

            # . Finish up
            self.isInitialized = True

    #-------------------------------------------------------------------------------
    def CalculateElectrostaticEnergies (self, **options):
        """
        Calculate for each instance of each site:
        - self (Born) energy in the model compound
        - interaction energy   between the site and the background charge set of the model compound

        - self (Born) energy in the protein
        - interaction energy   between the site and the background charge set of the protein
        - interaction energies between the site and the other sites in their different protonation forms

        Finally, use the calculated heterotransfer energies to calculate Gintr from Gmodel.
        """
        pass


    #-------------------------------------------------------------------------------
    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            if hasattr ( self, "label" ): label = self.label
            else:                         label = "Unknown"
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Continuum Electrostatic Model ({:s})".format ( label ) )

    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( "Number Of Sites"    , "{:d}".format ( self.nsites     ) ) ,
                   ( "Number Of Instances", "{:d}".format ( self.ninstances ) ) ]
        for ( name, key ) in self.__class__.defaultAttributeNames.items ( ):
            value = getattr ( self, key )
            if   isinstance ( value, bool ):
                convert = "True" if value else "False"
            elif isinstance ( value, float ):
                convert = "{:g}".format ( value )
            elif isinstance ( value, str ):
                convert = value
            elif isinstance ( value, tuple ):
                convert = "{:d}".format ( len ( value ) )
            else:
                convert = str ( value )
            items.append ( ( name, convert ) )
        return items

    #-------------------------------------------------------------------------------
    def SummarySites (self, log=logFile):
        """List titratable residues."""
        if LogFileActive (log):
            if self.isInitialized:
                tab = log.GetTable (columns=[8, 8, 8, 8, 10, 10, 10, 10])
                # . Write header
                tab.Start ()
                tab.Heading ("SiteID"),
                tab.Heading ("Site", columnSpan=3),
                tab.Heading ("Instances"),
                tab.Heading ("Center", columnSpan=3),
                # . Write sites
                for site in self.sites:
                    items = (
                        ( "{:d}"  .format ( site.siteIndex + 1 ) ),
                        ( "{:s}"  .format ( site.segName       ) ),
                        ( "{:s}"  .format ( site.resName       ) ),
                        ( "{:d}"  .format ( site.resSerial     ) ),
                        ( "{:d}"  .format ( site.ninstances    ) ),
                        ( "{:.3f}".format ( site.center[0]     ) ),
                        ( "{:.3f}".format ( site.center[1]     ) ),
                        ( "{:.3f}".format ( site.center[2]     ) ), )
                    for item in items:
                        tab.Entry (item)
                tab.Stop ()


    #-------------------------------------------------------------------------------
    def SummaryProbabilities (self, reportOnlyUnusual=False, maxProbThreshold=0.75, log=logFile):
        """List probabilities of occurance of instances."""
        unusualProtonations = {
            "HIS" : ("HSE", "HSD", "fd"),
            "ARG" : ("d", ),
            "ASP" : ("p", ),
            "CYS" : ("d", ),
            "GLU" : ("p", ),
            "LYS" : ("d", ),
            "TYR" : ("d", ), }

        if LogFileActive (log):
            if self.isProbability:
                maxinstances = 0
                for site in self.sites:
                    ninstances = len (site.instances)
                    if ninstances > maxinstances: maxinstances = ninstances

                tab = log.GetTable (columns=[6, 6, 6] + [8, 8] * maxinstances)
                tab.Start ()
                tab.Heading ("Site", columnSpan=3)
                tab.Heading ("Probabilities of instances", columnSpan=maxinstances * 2)
                for site in self.sites:
                    maxProb = 0.
                    for instance in site.instances:
                        if instance.probability > maxProb:
                            maxLabel = instance.label
                            maxIndex = instance.instIndex
                            maxProb  = instance.probability

                    skipSite = False
                    if reportOnlyUnusual:
                        if site.resName in unusualProtonations:
                            labels = unusualProtonations[site.resName]
                            if maxLabel not in labels:
                                skipSite = True

                    # . If the maximum probability is lower than a certain threshold, do not skip the site in the report
                    if maxProb < maxProbThreshold:
                        skipSite = False
                    if not skipSite:
                        tab.Entry ( site.segName )
                        tab.Entry ( site.resName )
                        tab.Entry ( "{:d}".format ( site.resSerial ) )
                        for instance in site.instances:
                            if instance.instIndex == maxIndex:
                                label = "*{:s}".format ( instance.label )
                            else:
                                label = instance.label
                            tab.Entry ( label )
                            tab.Entry ( "{:.4f}".format ( instance.probability ) )
                        for filler in range (maxinstances - site.ninstances):
                            tab.Entry ( "" )
                            tab.Entry ( "" )
                tab.Stop ()


    #-------------------------------------------------------------------------------
    def SedScript_FromProbabilities (self, filename="his_repl.sed", overwrite=False, putPath=False, log=logFile):
        """Generate a sed script for substituting histidines in the source PDB file based on the calculated probabilities."""
        if not self.isProbability:
            raise CEModelError ("First calculate probabilities.")

        if not overwrite:
            if os.path.exists (filename):
                if LogFileActive (log):
                    log.Paragraph ( "File {:s} already exists, skipping.".format ( filename ) )
                    return
        lines = []
        if putPath:
            lines.append ("# {:s}\n".format ( os.path.abspath (filename) ) )

        # . First write histidines
        for site in self.sites:
            if site.resName in ("HIS", "HSP"):
                mostProbValue, mostProbIndex, mostProbLabel = site.GetMostProbableInstance ()
                lines.append ("/H.. .{:4d}/  s/H../{:3s}/  # {:.4f}\n".format ( site.resSerial, mostProbLabel, mostProbValue ) )
        # . Then everything else
        unusualProtonations = {
            "ARG" : ("d", ),
            "ASP" : ("p", ),
            "CYS" : ("d", ),
            "GLU" : ("p", ),
            "LYS" : ("d", ),
            "TYR" : ("d", ), }
        translateLabels = {
            "p" :   "protonated",
            "d" : "deprotonated", }
        warnings = []

        for site in self.sites:
            if site.resName not in ("HIS", "HSP"):
                mostProbValue, mostProbIndex, mostProbLabel = site.GetMostProbableInstance ()

                if site.resName in unusualProtonations:
                    unusualLabels = unusualProtonations[site.resName]
                    if mostProbLabel in unusualLabels:
                        if   site.resName == "ASP":
                            lines.append ("# patch ASPP {:4s} {:4d} setup  ! {:.4f}\n".format ( site.segName, site.resSerial, mostProbValue ) )
                        elif site.resName == "GLU":
                            lines.append ("# patch GLUP {:4s} {:4d} setup  ! {:.4f}\n".format ( site.segName, site.resSerial, mostProbValue ) )
                        else:
                            warnings.append ("# Warning: {:4s} {:3s} {:4d} is {:s} with probability of {:.4f}\n".format ( site.segName, site.resName, site.resSerial, translateLabels[mostProbLabel], mostProbValue ) )
                else:
                    warnings.append ("# Unknown residue: {:4s} {:3s} {:4d} (most probable instance is \"{:s}\" with probability of {:.4f})\n".format ( site.segName, site.resName, site.resSerial, mostProbLabel, mostProbValue ) )
        lines.extend (warnings)

        # . Write a sed script
        WriteInputFile (filename, lines)
        # . Summarize
        if LogFileActive (log):
            log.Paragraph ( "Wrote file: {:s}".format ( filename ) )


    #-------------------------------------------------------------------------------
    def PrintInteractions ( self, log = logFile ):
        """Print the matrix of interactions in a readable form.

        It only remains readable as long as there are a few titratable sites with two instances each."""
        if self.isCalculated:
            if LogFileActive ( log ):
                # . Find labels.
                instanceLabels = []
                siteLabels     = []
                for asite in self.sites:
                    siteLabels.append ( ( asite.label, len ( asite.instances ) ) )
                    for ainstance in asite.instances: instanceLabels.append ( ainstance.label )
                # . Start table.
                table = log.GetTable ( columns = [ 12, 4 ] + len ( instanceLabels ) * [ 10 ] )
                table.Start ( )
                table.Title ( "Site Interaction Matrix" )
                # . Headings.
                table.Heading ( "" )
                table.Heading ( "" )
                for ( label, span ) in siteLabels: table.Heading ( label, columnSpan = span )
                table.Heading ( "" )
                table.Heading ( "" )
                for label in instanceLabels: table.Heading ( label )
                # . Double loop over sites.
                for asite in self.sites:
                    for ainstance in asite.instances:
                        table.Entry ( asite.label    , align = Align.Left )
                        table.Entry ( ainstance.label )
                        for bsite in self.sites:
                            for binstance in bsite.instances:
                                Wij  = self.energyModel.GetInteractionSymmetric ( ainstance._instIndexGlobal, binstance._instIndexGlobal )
                                table.Entry ( "{:.2f}".format ( Wij ) )
                table.Stop ( )

    #-------------------------------------------------------------------------------
    def WriteW (self, filename="W.dat", precision=3, log=logFile):
        """Write a GMCT-compatible matrix of interactions."""
        if self.isCalculated:
            items = (
                ( "idSiteA"  ,   8  ,   0 ),
                ( "idInstA"  ,   8  ,   0 ),
                ( "labSiteA" ,  14  ,  -1 ),
                ( "labInstA" ,  10  ,  -1 ),
                ( "idSiteB"  ,   8  ,   0 ),
                ( "idInstB"  ,   8  ,   0 ),
                ( "labSiteB" ,  14  ,  -1 ),
                ( "labInstB" ,  10  ,  -1 ),
                ( "Wij_symm" ,  10  ,   4 ),
                ( "Wij"      ,  10  ,   4 ),
                ( "Wij_err"  ,  10  ,   4 ), )
            # . Prepare header
            header = "# "
            for label, width, digits in items:
                header = "{:s}{:s}".format ( header, label.center (width) )
            lines  = ["{:s}\n".format ( header ), ]
            # . Prepare formating string
            form = ""
            for label, width, digits in items:
                if   digits < 0:
                    form += "{:" + "{:d}s".format ( width ) + "}"
                elif digits < 1:
                    form += "{:" + "{:d}d".format ( width ) + "}"
                else:
                    form += "{:" + "{:d}.{:d}f".format ( width, digits ) + "}"
            form += "\n"
            # . Outer loop
            for asite in self.sites:
                for ainstance in asite.instances:
                    # . Inner loop
                    for bsite in self.sites:
                        for binstance in bsite.instances:
                            interSymmetric = self.energyModel.GetInteractionSymmetric (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            interaction    = self.energyModel.GetInteraction          (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            deviation      = self.energyModel.GetDeviation            (ainstance._instIndexGlobal, binstance._instIndexGlobal)
                            lines.append (form.format (asite.siteIndex + 1, ainstance.instIndex + 1, asite.label, ainstance.label, bsite.siteIndex + 1, binstance.instIndex + 1, bsite.label, binstance.label, interSymmetric, interaction, deviation))
            # . Write to a file
            WriteInputFile (filename, lines)


    #-------------------------------------------------------------------------------
    def WriteGintr (self, filename="gintr.dat", precision=3, log=logFile):
        """Write a GMCT-compatible file containing intrinsic energies of each instance of each site."""
        if self.isCalculated:
            items = (
                ( "siteID"    ,   8  ,   0 ),
                ( "instID"    ,   8  ,   0 ),
                ( "siteLabel" ,  14  ,  -1 ),
                ( "instLabel" ,  10  ,  -1 ),
                ( "Gintr"     ,  12  ,   4 ),
                ( "protons"   ,   8  ,   0 ), )
            # . Prepare header
            header = "# "
            for label, width, digits in items:
                header = "{:s}{:s}".format (header, label.center (width))
            lines  = ["{:s}\n".format ( header ), ]
            # . Prepare formating string
            form = ""
            for label, width, digits in items:
                if   digits < 0:
                    form += "{:" + "{:d}s".format ( width ) + "}"
                elif digits < 1:
                    form += "{:" + "{:d}d".format ( width ) + "}"
                else:
                    form += "{:" + "{:d}.{:d}f".format ( width, digits ) + "}"
            form += "\n"
            # . Prepare instances
            for site in self.sites:
                for instance in site.instances:
                    lines.append ( form.format (site.siteIndex + 1, instance.instIndex + 1, site.label, instance.label, instance.Gintr, instance.protons))
            # . Write to a file
            WriteInputFile (filename, lines)


    #-------------------------------------------------------------------------------
    def DefineMCModel (self, sampler, log=logFile):
        """Assign a Monte Carlo model to the continuum electrostatic model."""
        if sampler is not None:
            # . Perform initial checks
            checks = (
                isinstance (sampler, MCModelDefault) ,
                isinstance (sampler, MCModelGMCT)    ,)
            if not any (checks):
                raise CEModelError ("Cannot define MC model.")

            # . Assign and initialize MC model
            self.sampler = sampler
            self.sampler.Initialize (self)
            self.sampler.PrintPairs (log=log)
            self.sampler.Summary    (log=log)


    #-------------------------------------------------------------------------------
    def CalculateMicrostateEnergy (self, stateVector, pH=7.0):
        """Calculate microstate energy."""
        return self.energyModel.CalculateMicrostateEnergy (stateVector, pH=pH)


    #-------------------------------------------------------------------------------
    def CalculateProbabilities (self, pH=7.0, unfolded=False, isCalculateCurves=False, logFrequency=-1, trajectoryFilename="", log=logFile):
        """Calculate probabilities.

        Setting |trajectoryFilename| will cause writing energies of sampled states to a file."""
        nstates = -1
        sites   = None

        if       hasattr (self, "sampler")  and     unfolded:
            raise CEModelError ("Monte Carlo sampling of unfolded proteins unsupported.")
        elif     hasattr (self, "sampler")  and not unfolded:
            self.sampler.CalculateOwnerProbabilities (pH=pH, logFrequency=logFrequency, trajectoryFilename=trajectoryFilename, log=log)
        elif not hasattr (self, "sampler")  and     unfolded:
            if (trajectoryFilename != ""):
                raise CEModelError ("Writing trajectories of unfolded proteins unsupported.")
            nstates = self.energyModel.CalculateProbabilitiesAnalyticallyUnfolded (pH=pH)
        elif not hasattr (self, "sampler")  and not unfolded:
            # TODO !!!
            if (trajectoryFilename != ""):
                raise CEModelError ("Writing trajectories unsupported.")
            nstates = self.energyModel.CalculateProbabilitiesAnalytically (pH=pH)

        if isCalculateCurves:
            sites = []
            for site in self.sites:
                instances = []
                for instance in site.instances:
                    instances.append (instance.probability)
                sites.append (instances)
        if nstates > 0:
            if LogFileActive (log):
                log.Paragraph ( "Calculated {:d} protonation states.".format ( nstates ) )

        self.isProbability = True
        return sites


    #-------------------------------------------------------------------------------
    def _GetResidueInfo (self, residue):
        system        = self.system
        ParseLabel    = system.sequence.ParseLabel
        # . Get segment label
        segment       = residue.parent
        segmentName   = segment.label
        # . Get residue label and serial
        residueName, residueSerial = ParseLabel (residue.label, fields=2)
        residueSerial = int (residueSerial)

        return (segmentName, residueName, residueSerial)


    #-------------------------------------------------------------------------------
    def _GetIndices (self, residue, atomLabels, check=True):
        missingLabels = []
        atomIndices   = []
        atoms         = residue.children
        for label in atomLabels:
            index = -1
            for atom in atoms:
                if label == atom.label:
                    index = atom.index
                    break
            if index >= 0:
                atomIndices.append (index)
            else:
                missingLabels.append (label)
            if check and missingLabels:
                segmentName, residueName, residueSerial = self._GetResidueInfo (residue)
                raise CEModelError ( "Cannot include residue {:s} {:s} {:d} because of missing atoms: {:s}".format ( segmentName, residueName, residueSerial, " ".join (missingLabels) ) )
        return atomIndices


    #-------------------------------------------------------------------------------
    def _SetupBackground (self):
        allSiteAtomIndices  = []
        for site in self.sites:
            allSiteAtomIndices.extend (site.siteAtomIndices)
        system              = self.system
        segments            = system.sequence.children
        backAtomIndices     = []
        proteinAtomIndices  = []

        #============ Iterate segments ============
        for segment in segments:
            residues = segment.children
        
            #============ Iterate residues ============
            for residue in residues:
                foo, residueName, residueSerial = self._GetResidueInfo (residue)
        
                # . Remove residues not defined in PROTEIN_RESIDUES, usually waters and ions
                if residueName not in REMOVE_RESIDUES:
                    atoms = residue.children
        
                    #============ Iterate atoms ============
                    for atom in atoms:
                        proteinAtomIndices.append (atom.index)
                        if atom.index not in allSiteAtomIndices:
                            backAtomIndices.append (atom.index)
        self.proteinAtomIndices = proteinAtomIndices
        self.backAtomIndices    = backAtomIndices


    #-------------------------------------------------------------------------------
    def _SetupSites (self, residue, prevResidue=None, nextResidue=None, terminal=None, log=logFile):
        segmentName, residueName, residueSerial = self._GetResidueInfo (residue)
        setupSites    = []
        if terminal:
            if   terminal == "N":
                if   residueName == "GLY":
                    libTerm = self.library["NGL"]
                elif residueName == "PRO":
                    libTerm = self.library["NPR"]
                else:
                    libTerm = self.library["NTR"]
            elif terminal == "C":
                libTerm = self.library["CTR"]
            termIndices = self._GetIndices (residue, libTerm.atomLabels)
            setupSites.append ([terminal, libTerm, termIndices, termIndices])

        # . Check for a titratable residue
        if residueName in self.library:
            libSite      = self.library[residueName]
            siteIndices  = self._GetIndices (residue, libSite.atomLabels)
            modelIndices = []
            for atom in residue.children:
                if terminal:
                    if (atom.label in TERM_REMOVE) or (atom.index in termIndices):
                        continue
                modelIndices.append (atom.index)

            if not terminal:
                prevIndices  = []
                if prevResidue:
                    foo, prevResidueName, prevResidueSerial = self._GetResidueInfo (prevResidue)
                    if prevResidueName in PROTEIN_RESIDUES:
                        prevIndices = self._GetIndices (prevResidue, PREV_RESIDUE, check=False)
                    modelIndices = prevIndices + modelIndices
    
                nextIndices  = []
                if nextResidue:
                    foo, nextResidueName, nextResidueSerial = self._GetResidueInfo (nextResidue)
                    if nextResidueName in PROTEIN_RESIDUES:
                        if   nextResidueName == "GLY":
                            nextLabels = NEXT_RESIDUE_GLY
                        elif nextResidueName == "PRO":
                            nextLabels = NEXT_RESIDUE_PRO
                        else:
                            nextLabels = NEXT_RESIDUE
                        nextIndices = self._GetIndices (nextResidue, nextLabels, check=False)
                    modelIndices = modelIndices + nextIndices
            setupSites.append (["SITE", libSite, siteIndices, modelIndices])

            if terminal == "C":
                setupSites.reverse ()
        return setupSites


    #-------------------------------------------------------------------------------
    def _CreateSite (self, **options):
        """Create a site and its instances specific to the CE model."""
        pass


    #-------------------------------------------------------------------------------
    def _SplitModel (self, excludeSegments=_DEFAULT_EXCLUDE_SEGMENTS, excludeResidues=None, includeTermini=False, dryRun=True, log=logFile):
        siteIndex        =  0
        instIndexGlobal  =  0
        totalInstances   =  0

        excluded = []
        if not dryRun:
            self.sites = []
        system    = self.system
        segments  = system.sequence.children

        #============ Iterate segments ============
        for segment in segments:
            segmentName = segment.label

            # . Include segment?
            if segmentName not in excludeSegments:
                residues  = segment.children
                nresidues = len (residues)

                #============ Iterate residues ============
                for residueIndex, residue in enumerate (residues):
                    foo, residueName, residueSerial = self._GetResidueInfo (residue)
                    includeResidue = self._CheckResidue ( excludeResidues, segmentName, residueName, residueSerial )
                    if not includeResidue:
                        excluded.append ( ( segmentName, residueName, residueSerial ) )
                        continue
                    if residueIndex > 0:
                        prevResidue = residues[residueIndex - 1]
                    else:
                        prevResidue = None
                    if residueIndex < (nresidues - 1):
                        nextResidue = residues[residueIndex + 1]
                    else:
                        nextResidue = None
                    if   residueIndex < 1:
                        terminal = "N"
                    elif residueIndex > (nresidues - 2):
                        terminal = "C"
                    else:
                        terminal = None

                    setupSites = self._SetupSites (residue, prevResidue, nextResidue, terminal=(terminal if (includeTermini and (residueName in PROTEIN_RESIDUES)) else None), log=log)
                    for siteType, libSite, siteAtomIndices, modelAtomIndices in setupSites:
                        if not dryRun:
                            if   siteType == "N":
                                updatedSerial = 998
                                updatedName   = libSite.label
                            elif siteType == "C":
                                updatedSerial = 999
                                updatedName   = libSite.label
                            else:
                                updatedSerial = residueSerial
                                updatedName   = residueName

                            # . Create a new site in the CE model
                            instIndexGlobal = self._CreateSite (
                                siteIndex        = siteIndex         ,
                                segName          = segmentName       ,
                                resName          = updatedName       ,
                                resSerial        = updatedSerial     ,
                                siteAtomIndices  = siteAtomIndices   ,
                                modelAtomIndices = modelAtomIndices  ,

                                libSite          = libSite           ,
                                instIndexGlobal  = instIndexGlobal   ,
                                )
                        # . Update the numbers of sites and instances
                        siteIndex      += 1
                        totalInstances += len (libSite.instances)

        # . Printing.
        if LogFileActive ( log ):
            if len ( excluded ) > 0:
                table = log.GetTable ( columns = [ 10, 10, 10 ] )
                table.Start ( )
                table.Title ( "CE Model Excluded Residues" )
                for ( segmentName, residueName, residueSerial ) in excluded:
                    table.Entry ( segmentName )
                    table.Entry ( residueName )
                    table.Entry ( "{:d}".format ( residueSerial ) )
                table.Stop ( )
            else:
                log.Paragraph ( "No residues were excluded from CE model." )

        # . Return the numbers of sites and instances
        totalSites = siteIndex
        return (totalSites, totalInstances)


    #-------------------------------------------------------------------------------
    def _CheckResidue ( self, excludeResidues, segmentName, residueName, residueSerial ):
        """Check if the residue should be included."""
        includeResidue = True
        if excludeResidues:
            for exclSegmentName, exclResidueName, exclResidueSerial in excludeResidues:
                if   (    exclSegmentName) and (    exclResidueName) and (    exclResidueSerial):
                    if exclSegmentName == segmentName and exclResidueName == residueName and exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (    exclSegmentName) and (    exclResidueName) and (not exclResidueSerial):
                    if exclSegmentName == segmentName and exclResidueName == residueName:
                        includeResidue = False
                        break

                elif (    exclSegmentName) and (not exclResidueName) and (    exclResidueSerial):
                    if exclSegmentName == segmentName and exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (    exclSegmentName) and (not exclResidueName) and (not exclResidueSerial):
                    if exclSegmentName == segmentName:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (    exclResidueName) and (    exclResidueSerial):
                    if exclResidueName == residueName and exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (    exclResidueName) and (not exclResidueSerial):
                    if exclResidueName == residueName:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (not exclResidueName) and (    exclResidueSerial):
                    if exclResidueSerial == residueSerial:
                        includeResidue = False
                        break

                elif (not exclSegmentName) and (not exclResidueName) and (not exclResidueSerial):
                    includeResidue = False
                    break

        return includeResidue

#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

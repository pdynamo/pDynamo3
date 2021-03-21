"""Substate is a class for calculating energies of substates of protonation states."""

from  pCore           import logFile, LogFileActive
from .CEModelError    import CEModelError
from .InputFileWriter import WriteInputFile
from .StateVector     import StateVector

class Substate:
    """Substate of a protonation state."""

    def __init__ (self, meadModel, selectedSites, pH=7.0, log=logFile):
        """Construct a substate defined by |selectedSites|.

        |selectedSites| is a sequence of three-element sequences (segmentName, residueName, residueSerial)"""
        pairs = []
        indicesOfSites = []

        for ( selectedSegment, selectedResidueName, selectedResidueSerial ) in selectedSites:
            foundSite = False
            for siteIndex, site in enumerate ( meadModel.sites ):
#                print ( "Sites>", selectedSegment, selectedResidueName, selectedResidueSerial, site.segName, site.resSerial, type ( site.resSerial ) )
                if ( site.segName == selectedSegment ) and ( site.resSerial == selectedResidueSerial ):
                    indicesOfSites.append (siteIndex)
                    foundSite = True
                    break

            if not foundSite:
                raise CEModelError ( "Site {:s} {:s} {:d} not found.".format ( selectedSegment, selectedResidueName, selectedResidueSerial ) )
            pairs.append ([selectedSegment, selectedResidueSerial])

        self.indicesOfSites = indicesOfSites
        self.isCalculated   = False
        self.substates      = None
        self.owningModel    = meadModel
        self.pH             = pH
        self.vector         = self._DetermineLowestEnergyVector ()
        self.vector.DefineSubstate (pairs)

        if LogFileActive (log):
            nsites = len (indicesOfSites)
            log.Paragraph ( "Substate is initialized with {:d} site{:s}.".format ( nsites, "s" if nsites > 1 else "" ) )


    def _DetermineLowestEnergyVector (self):
        """Find a vector representing the lowest energy protonation state."""
        owner   = self.owningModel
        restore = False
        if owner.isProbability:
            backup  = self._ProbabilitiesSave ()
            restore = True
        owner.CalculateProbabilities (pH=self.pH, log=None)
        vector = StateVector_FromProbabilities (owner)
        if restore:
            self._ProbabilitiesRestore (backup)
        else:
            owner.isProbabilty = False
        return vector


    def _ProbabilitiesSave (self):
        """Backup the probabilities of the owning model."""
        meadModel     = self.owningModel
        probabilities = None

        if meadModel.isProbability:
            energyModel   = meadModel.energyModel
            probabilities = [0.] * meadModel.ninstances
            for site in meadModel.sites:
                for instance in site.instances:
                    probabilities[instance._instIndexGlobal] = energyModel.GetProbability (instance._instIndexGlobal)
        return probabilities


    def _ProbabilitiesRestore (self, probabilities):
        """Restore the original probabilities to the owning model."""
        if probabilities:
            meadModel   = self.owningModel
            energyModel = meadModel.energyModel
            for site in meadModel.sites:
                for instance in site.instances:
                    energyModel.SetProbability (instance._instIndexGlobal, probabilities[instance._instIndexGlobal])


    def CalculateSubstateEnergies (self, log=logFile):
        """Calculate microstate energies for a substate."""
        if not self.isCalculated:
            indicesOfSites = self.indicesOfSites
            increment      = True
            substates      = []
            vector         = self.vector
            owner          = self.owningModel
            vector.ResetSubstate ()

            while increment:
                Gmicro = owner.CalculateMicrostateEnergy (vector, pH=self.pH)
                indicesOfInstances = []
                for siteIndex in indicesOfSites:
                    indicesOfInstances.append (vector[siteIndex])
                substates.append ([Gmicro, indicesOfInstances])

                increment = vector.IncrementSubstate ()

            substates.sort ()
            lowestSubstate = substates[0]
            lowestEnergy   = lowestSubstate[0]

            self.substates    = substates
            self.zeroEnergy   = lowestEnergy
            self.isCalculated = True

            if LogFileActive (log):
                log.Paragraph ( "Calculating substate energies at pH={:.1f} complete.".format ( self.pH ) )


    def Summary (self, relativeEnergy=True, roundCharge=True, title="", log=logFile):
        """Summarize calculated substate energies in a table."""
        if self.isCalculated:
            indicesOfSites = self.indicesOfSites
            zeroEnergy     = self.zeroEnergy
            substates      = self.substates
            owner          = self.owningModel
            nsites         = len (indicesOfSites)

            if LogFileActive (log):
                tab = log.GetTable (columns = [6, 9, 8, 8] + [14] * nsites)
                tab.Start ()
                if title: tab.Title (title)
                tab.Heading ("State")
                tab.Heading ("Gmicro")
                tab.Heading ("Charge")
                tab.Heading ("Protons")

                for siteIndex in indicesOfSites:
                    site = owner.sites[siteIndex]
                    tab.Heading ("{:s} {:s} {:d}".format ( site.segName, site.resName, site.resSerial ) )

                for stateIndex, (energy, indicesOfInstances) in enumerate (substates):
                    tab.Entry ("{:6d}".format ( stateIndex + 1 ) )
                    if relativeEnergy:
                        energy = energy - zeroEnergy
                    tab.Entry ( "{:.2f}".format ( energy ) )

                    nprotons = 0
                    charge   = 0.
                    labels   = []
                    for siteIndex, instanceIndex in zip (indicesOfSites, indicesOfInstances):
                        site     = owner.sites     [siteIndex]
                        instance = site.instances  [instanceIndex]
                        nprotons = nprotons + instance.protons
                        charge   = charge   + sum (instance.charges)
                        labels.append (instance.label)

                    # Charges should ALWAYS sum up to integer values
                    if roundCharge:
                        tab.Entry ("{:d}".format ( round ( charge ) ) )
                    else:
                        tab.Entry ("{:.1f}".format ( charge ) )
                    tab.Entry ("{:d}".format ( nprotons ) )
                    for label in labels:
                        tab.Entry ( label )
                tab.Stop ()


    def Summary_ToLatex (self, filename="table.tex", relativeEnergy=True, includeSegment=False, translateToLatex=None):
        """Summarize calculated substate energies in a Latex table."""
        transl = {
            "ASP" : {"p" : "0", "d" : "($-$)"},
            "GLU" : {"p" : "0", "d" : "($-$)"},
            "HIS" : {"HSP" : "$\\epsilon$, $\\delta$(+)", "HSE" : "$\\epsilon$", "HSD" : "$\\delta$", "fd" : "($-$)"},
                 }
        if translateToLatex:
            transl.update (translateToLatex)

        if self.isCalculated:
            indicesOfSites = self.indicesOfSites
            zeroEnergy     = self.zeroEnergy
            substates      = self.substates
            owner          = self.owningModel
            nsites         = len (indicesOfSites)

            lines  = [ "\\begin{tabular}{@{\\extracolsep{2mm}}cc" + "{:s}".format ( "l" * nsites ) + "c}" ]
            header = "State & $\\Delta E$ (kcal/mol) & "
            for siteIndex in indicesOfSites:
                site = owner.sites[siteIndex]
                if includeSegment:
                    header = "{:s} {:s} {:s}{:d} &".format (header, site.segName, site.resName.capitalize (), site.resSerial)
                else:
                    header = "{:s} {:s}{:d} &".format  ( header , site.resName.capitalize (), site.resSerial)
            header = "{:s} No. of protons \\\\".format ( header )
            lines.append (header)
            lines.append ("\\hline\\noalign{\\smallskip}")

            for substateCount, (energy, indicesOfInstances) in enumerate (substates, 1):
                if relativeEnergy:
                    energy = energy - zeroEnergy
                line = "{:3d} & {:5.1f} & ".format ( substateCount, energy )
                nprotons = 0

                for siteIndex, instanceIndex in zip (indicesOfSites, indicesOfInstances):
                    site     = owner.sites     [siteIndex]
                    instance = site.instances  [instanceIndex]
                    if site.resName in transl:
                        dic  = transl[site.resName]
                        line = "{:s}  {:26s} &".format ( line, dic[instance.label] )
                    else:
                        line = "{:s}  {:26s} &".format ( line, instance.label )
                    nprotons = nprotons + instance.protons
                line = "{:s}  {:1d} \\\\".format ( line, nprotons )
                lines.append (line)
            lines.append ("\\hline\\noalign{\\smallskip}")
            lines.append ("\\end{tabular}")

            WriteInputFile (filename, lines, addLineBreaks=True)


#-------------------------------------------------------------------------------
class MEADSubstate (Substate):
    """A class to maintain backward compatibility."""
    pass


#===============================================================================
# . Helper functions
#===============================================================================
def StateVector_FromProbabilities (meadModel):
    """Create a state vector from the previously calculated probabilities."""
    if meadModel.isProbability:
        vector = StateVector (meadModel)
        for siteIndex, site in enumerate (meadModel.sites):
            pairs = []
            for instanceIndex, instance in enumerate (site.instances):
                pair = (instance.probability, instanceIndex)
                pairs.append (pair)
            maxProbPair = max (pairs)
            probability, instanceIndex = maxProbPair
            vector[siteIndex] = instanceIndex
        return vector
    else:
        raise CEModelError ("First calculate probabilities.")


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

import os, subprocess

from  pCore                import logFile, LogFileActive
from .CEModelError         import CEModelError
from .Instance             import Instance
from .MEADOutputFileReader import MEADOutputFileReader

class InstanceMEAD (Instance):
    """A class to represent a MEAD type of instance."""

    _attributable = {
        }
    _attributable.update (Instance._attributable)
    # modelPqr  modelLog  modelGrid  sitePqr  siteLog  siteGrid

    def __init__ (self, **options):
        """Constructor."""
        super (InstanceMEAD, self).__init__ (**options)


    #-------------------------------------------------------------------------------
    def CalculateModelCompound (self, log=logFile):
        """Calculate Gborn and Gback of a site in a model compound."""
        site  = self.parent
        model = site.parent

        if model.isFilesWritten:
            if os.path.exists (self.modelLog):
                pass
            else:
                instancePqr        , ext = os.path.splitext (self.sitePqr)
                modelBackgroundPqr , ext = os.path.splitext (self.modelPqr)
                command = [
                    os.path.join ( model.pathMEAD, "my_2diel_solver"   ) , 
                    "-T"        , "{:f}".format ( model.temperature    ) , 
                    "-ionicstr" , "{:f}".format ( model.ionicStrength  ) , 
                    "-epsin"    , "{:f}".format ( model.epsilonProtein ) , 
                    "-epsext"   , "{:f}".format ( model.epsilonWater   ) , 
                    instancePqr, 
                    modelBackgroundPqr
                    ]
                try:
                    outFile = open (self.modelLog, "w")
                    subprocess.check_call (command, stderr=outFile, stdout=outFile)
                    outFile.close ()
                except:
                    raise CEModelError ( "Failed running command: {:s}".format ( " ".join (command) ) )
            reader = MEADOutputFileReader.FromPath (self.modelLog)
            reader.Parse ()

            checks = (hasattr (reader, "born"), hasattr (reader, "back"), )
            if not all (checks):
                raise CEModelError ("Output file {:s} empty or corrupted. Empty the scratch directory and start anew.".format ( self.modelLog ) )
            self.Gborn_model = reader.born
            self.Gback_model = reader.back


    #-------------------------------------------------------------------------------
    def CalculateProtein (self, log=logFile):
        """Calculate Gborn, Gback and Wij of a site in protein environment."""
        site  = self.parent
        model = site.parent

        if model.isFilesWritten:
            if os.path.exists (self.siteLog):
                pass
            else:
                # . Assign removing extensions, otherwise MEAD does not work
                sitesFpt             , ext = os.path.splitext (model.pathFptSites)
                proteinPqr           , ext = os.path.splitext (model.pathPqrProtein)
                proteinBackgroundPqr , ext = os.path.splitext (model.pathPqrBack)
                instancePqr          , ext = os.path.splitext (self.sitePqr)

                # . epsin1 is never used but must be given

                # . eps2set defines the whole protein
                command = [
                    os.path.join ( model.pathMEAD, "my_3diel_solver" ), 
                    "-T"        , "{:f}".format ( model.temperature    ) ,
                    "-ionicstr" , "{:f}".format ( model.ionicStrength  ) ,
                    "-epsin1"   , "{:f}".format ( 1.0                  ) ,
                    "-epsin2"   , "{:f}".format ( model.epsilonProtein ) ,
                    "-epsext"   , "{:f}".format ( model.epsilonWater   ) ,
                    "-eps2set"  , "{:s}".format ( proteinPqr           ) ,
                    "-fpt"      , "{:s}".format ( sitesFpt             ) ,
                    instancePqr , 
                    proteinBackgroundPqr
                    ]
                try:
                    outFile = open (self.siteLog, "w")
                    subprocess.check_call (command, stderr=outFile, stdout=outFile)
                    outFile.close ()
                except:
                    raise CEModelError ("Failed running command: {:s}".format ( " ".join (command) ) )
            reader = MEADOutputFileReader.FromPath (self.siteLog)
            reader.Parse ()

            checks = (hasattr (reader, "born"), hasattr (reader, "back"), hasattr (reader, "interactions"), )
            if not all (checks):
                raise CEModelError ("Output file {:s} empty or corrupted. Empty the scratch directory and start anew.".format ( self.modelLog ) )
            self.Gborn_protein = reader.born
            self.Gback_protein = reader.back

            # . Create a list of interactions
            interactions    = []
            instances       = []
            siteIndexOld    = 99999
            parentSiteIndex = self.parent.siteIndex

            for siteIndex, instanceIndex, energy in reader.interactions:
                if siteIndex > siteIndexOld:
                    interactions.append (instances)
                    instances = []
                # . Set the interaction energy to zero if the site is interacting with itself
                if siteIndex == parentSiteIndex:
                    energy = 0.
                instances.append (energy)
                siteIndexOld = siteIndex

            if instances:
                interactions.append (instances)

            # . Copy the interactions to the centralized array
            indexGlobal = 0
            for site in interactions:
                for instance in site:
                    energy = instance
                    model.energyModel.SetInteraction (self._instIndexGlobal, indexGlobal, energy)
                    indexGlobal = indexGlobal + 1


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

from  pCore        import logFile, LogFileActive
from .CEModelError import CEModelError

class Instance:
    """Base class for an instance of a titratable site.

    This class should not be used directly."""
    _attributable = {
        }
    # parent  instIndex  _instIndexGlobal  label  charges  Gborn_model  Gback_model  Gborn_protein  Gback_protein

    def __init__ (self, **options):
        """Constructor."""
        for (key, val) in options.items ():
            setattr (self, key, val)


    def _GetEnergyModel (self):
        if hasattr (self, "parent"):
            site        = self.parent
            ceModel     = site.parent
            energyModel = ceModel.energyModel
        else:
            raise CEModelError ("Energy model is undefined.")
        return energyModel


    @property
    def Gmodel (self):
        em = self._GetEnergyModel ()
        return em.GetGmodel (self._instIndexGlobal)

    @Gmodel.setter
    def Gmodel (self, value):
        em = self._GetEnergyModel ()
        em.SetGmodel (self._instIndexGlobal, value)

    @property
    def Gintr (self):
        em = self._GetEnergyModel ()
        return em.GetGintr (self._instIndexGlobal)

    @Gintr.setter
    def Gintr (self, value):
        em = self._GetEnergyModel ()
        em.SetGintr (self._instIndexGlobal, value)

    @property
    def protons (self):
        em = self._GetEnergyModel ()
        return em.GetProtons (self._instIndexGlobal)

    @protons.setter
    def protons (self, value):
        em = self._GetEnergyModel ()
        em.SetProtons (self._instIndexGlobal, value)

    @property
    def probability (self):
        em = self._GetEnergyModel ()
        return em.GetProbability (self._instIndexGlobal)

    @probability.setter
    def probability (self, value):
        em = self._GetEnergyModel ()
        em.SetProbability (self._instIndexGlobal, value)

    @property
    def interactions (self):
        # . Return a list of interactions of the current instance with other instances
        if not hasattr (self, "parent"):
            raise CEModelError ("Energy model is undefined.")
        site        = self.parent
        ceModel     = site.parent
        energyModel = ceModel.energyModel
        energies    = []
        for indexOther in range (ceModel.ninstances):
            energies.append (energyModel.GetInteractionSymmetric (self._instIndexGlobal, indexOther))
        return energies


    def CalculateModelCompound (self, log=logFile):
        """Calculate Gborn and Gback of a site in a model compound."""
        pass


    def CalculateProtein (self, log=logFile):
        """Calculate Gborn, Gback and Wij of a site in protein environment."""
        pass


    def CalculateGintr (self, log=logFile):
        """Calculate Gintr of an instance of a site in a protein."""
        checks = (
            hasattr (self , "Gborn_protein") ,
            hasattr (self , "Gback_protein") ,
            hasattr (self , "Gborn_model"  ) ,
            hasattr (self , "Gback_model"  ) ,)
        if all (checks):
            self.Gintr = self.Gmodel + (self.Gborn_protein - self.Gborn_model) + (self.Gback_protein - self.Gback_model)


    def PrintInteractions (self, sort=False, log=logFile):
        """Print interactions of an instance of a site with other instances of other sites."""
        if LogFileActive (log):
            site         = self.parent
            model        = site.parent
            interactions = self.interactions

            if model.isCalculated:
                instances = []
                for site in model.sites:
                    for instance in site.instances:
                        wij = interactions[instance._instIndexGlobal]
                        instances.append ([wij, site.segName, site.resName, site.resSerial, instance.label])
                if sort:
                    instances.sort ()

                tab = log.GetTable (columns=[6, 6, 6, 6, 16])
                tab.Start ()
                tab.Heading ("Instance of a site", columnSpan=4)
                tab.Heading ("Wij")
                for wij, segName, resName, resSerial, label in instances:
                    entries = ( ( "{:s}"  .format ( segName   ) ) ,
                                ( "{:s}"  .format ( resName   ) ) ,
                                ( "{:d}"  .format ( resSerial ) ) ,
                                ( "{:s}"  .format ( label     ) ) ,
                                ( "{:.4f}".format ( wij       ) ) , )
                    for entry in entries:
                        tab.Entry (entry)
                tab.Stop ()


    def _TableEntry (self, tab=None, secondsToCompletion=None):
        """Report calculated energies in a table.

        Optionally, include Estimated Time for Accomplishment (ETA).

        ETA has to be calculated outside of this method."""
        if tab:
            site = self.parent
            entries = (( "{:s}"  .format (  site.segName        ) ),
                       ( "{:s}"  .format (  site.resName        ) ),
                       ( "{:d}"  .format (  site.resSerial      ) ),
                       ( "{:s}"  .format (  self.label          ) ),
                       ( "{:.4f}".format (  self.Gborn_model    ) ),
                       ( "{:.4f}".format (  self.Gback_model    ) ),
                       ( "{:.4f}".format (  self.Gborn_protein  ) ),
                       ( "{:.4f}".format (  self.Gback_protein  ) ),
                       ( "{:.4f}".format (  self.Gmodel         ) ),
                       ( "{:.4f}".format (  self.Gintr          ) ), )
            for entry in entries:
                tab.Entry (entry)

            if isinstance (secondsToCompletion, float):
                minutes, seconds = divmod (secondsToCompletion, 60)
                hours, minutes   = divmod (minutes, 60)
                tab.Entry ( "{:d}:{:0>2}:{:0>2}".format ( hours, minutes, seconds ) )


#===============================================================================
# . Main program
#===============================================================================
if __name__ == "__main__": pass

from  pCore        import logFile       , \
                          LogFileActive
from .CEModelError import CEModelError

DEF ANALYTIC_STATES = 67108864

cdef class EnergyModel:
    """A class defining the energy model."""

    def __dealloc__ (self):
        """Deallocate."""
        EnergyModel_Deallocate (&self.cObject)

    # Setters
    def SetGmodel (self, CInteger instIndexGlobal, CReal value):
        EnergyModel_SetGmodel (self.cObject, instIndexGlobal, value)

    def SetGintr (self, CInteger instIndexGlobal, CReal value):
        EnergyModel_SetGintr (self.cObject, instIndexGlobal, value)

    def SetProtons (self, CInteger instIndexGlobal, CInteger value):
        EnergyModel_SetProtons (self.cObject, instIndexGlobal, value)

    def SetProbability (self, CInteger instIndexGlobal, CReal value):
        EnergyModel_SetProbability (self.cObject, instIndexGlobal, value)

    def SetInteraction (self, CInteger instIndexGlobalA, CInteger instIndexGlobalB, CReal value):
        EnergyModel_SetInteraction (self.cObject, instIndexGlobalA, instIndexGlobalB, value)

    # Getters
    def GetGmodel (self, CInteger instIndexGlobal):
        cdef CReal Gmodel = EnergyModel_GetGmodel (self.cObject, instIndexGlobal)
        return Gmodel

    def GetGintr (self, CInteger instIndexGlobal):
        cdef CReal Gintr = EnergyModel_GetGintr (self.cObject, instIndexGlobal)
        return Gintr

    def GetProtons (self, CInteger instIndexGlobal):
        cdef CInteger protons = EnergyModel_GetProtons (self.cObject, instIndexGlobal)
        return protons

    def GetProbability (self, CInteger instIndexGlobal):
        cdef CReal prob = EnergyModel_GetProbability (self.cObject, instIndexGlobal)
        return prob

    def GetInteraction (self, CInteger instIndexGlobalA, CInteger instIndexGlobalB):
        cdef CReal interac = EnergyModel_GetInteraction (self.cObject, instIndexGlobalA, instIndexGlobalB)
        return interac

    def GetInteractionSymmetric (self, CInteger instIndexGlobalA, CInteger instIndexGlobalB):
        cdef CReal interac = EnergyModel_GetInterSymmetric (self.cObject, instIndexGlobalA, instIndexGlobalB)
        return interac

    def GetDeviation (self, CInteger instIndexGlobalA, CInteger instIndexGlobalB):
        cdef CReal deviate = EnergyModel_GetDeviation (self.cObject, instIndexGlobalA, instIndexGlobalB)
        return deviate


    def __init__ (self, ceModel, CInteger totalSites, CInteger totalInstances):
        """Constructor."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = EnergyModel_Allocate (totalSites, totalInstances, &cStatus)
        if cStatus != CStatus_OK: raise CEModelError ( "Cannot allocate energy model." )
        self.cObject.temperature = ceModel.temperature
        self.cObject.ninstances  = totalInstances
        self.isOwner             = False
        self.owningModel         = ceModel


    def Initialize (self):
        """Set up sites and calculate the number of possible states."""
        cdef CInteger nstates, ninstances, totalSites
        cdef CInteger indexSite, indexDown, indexUp, index
        cdef CStatus  cStatus = CStatus_OK
        totalSites = self.cObject.vector.nsites
        ceModel    = self.owningModel
        nstates    = 1
        for indexSite from 0 <= indexSite < totalSites:
            ninstances = 0
            indexUp    = 0
            indexDown  = 9999
            site       = ceModel.sites[indexSite]
            for instance in site.instances:
                ninstances += 1
                index = instance._instIndexGlobal
                if index < indexDown : indexDown = index
                if index > indexUp   : indexUp   = index
            StateVector_SetSite ( self.cObject.vector, indexSite, indexDown, indexUp, &cStatus )
            if nstates <= ANALYTIC_STATES:
                nstates = nstates * ninstances
        self.cObject.nstates = nstates

    def CheckIfSymmetric (self, CReal tolerance=0.05):
        """Check if the matrix of interactions is symmetric within the given threshold."""
        cdef CReal    maxDeviation
        cdef CBoolean isSymmetric
        isSymmetric = EnergyModel_CheckInteractionsSymmetric (self.cObject, tolerance, &maxDeviation)
        return (isSymmetric, maxDeviation)


    def SymmetrizeInteractions (self, log=logFile):
        """Symmetrize the matrix of interactions."""
        cdef CStatus cStatus = CStatus_OK
        EnergyModel_SymmetrizeInteractions (self.cObject, &cStatus)
        if LogFileActive (log): log.Paragraph ( "Symmetrizing interactions complete." )


    def ResetInteractions (self):
        """Set all interactions to zero."""
        EnergyModel_ResetInteractions (self.cObject)


    def ScaleInteractions (self, CReal scale):
        """Scale interactions."""
        EnergyModel_ScaleInteractions (self.cObject, scale)


    def CalculateMicrostateEnergy (self, StateVector vector, CReal pH=7.0):
        """Calculate energy of a protonation state (=microstate)."""
        cdef CReal Gmicro
        ceModel = self.owningModel

        if not ceModel.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")
        Gmicro = EnergyModel_CalculateMicrostateEnergy (self.cObject, vector.cObject, pH)
        return Gmicro


    def CalculateProbabilitiesAnalytically (self, CReal pH=7.0):
        """Calculate probabilities of protonation states analytically."""
        cdef CStatus cStatus = CStatus_OK
        ceModel = self.owningModel
        if self.cObject.nstates > ANALYTIC_STATES:
            raise CEModelError ( "Maximum number of states for analytic treatment ({:d}) exceeded.".format ( ANALYTIC_STATES ) )
        if not ceModel.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")
        EnergyModel_CalculateProbabilitiesAnalytically (self.cObject, pH, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ("Cannot allocate Boltzmann factors.")
        return self.cObject.nstates


    def CalculateProbabilitiesAnalyticallyUnfolded (self, CReal pH=7.0):
        """Calculate probabilities of protonation states analytically (unfolded protein)."""
        cdef CStatus cStatus = CStatus_OK
        ceModel = self.owningModel

        if self.cObject.nstates > ANALYTIC_STATES:
            raise CEModelError ( "Maximum number of states for analytic treatment ({:d}) exceeded.".format ( ANALYTIC_STATES ) )
        if not ceModel.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")

        EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded (self.cObject, pH, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ("Cannot allocate Boltzmann factors.")
        return self.cObject.nstates


    def CalculateMicrostateEnergyUnfolded (self, StateVector vector, CReal pH=7.0):
        """Calculate energy of a protonation state (=microstate) in an unfolded protein."""
        cdef CReal Gmicro
        ceModel = self.owningModel
        if not ceModel.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")
        Gmicro = EnergyModel_CalculateMicrostateEnergyUnfolded (self.cObject, vector.cObject, pH)
        return Gmicro


    def CalculateZfolded (self, CReal Gzero=0.0, CReal pH=7.0):
        """Calculate partition function of a folded protein."""
        cdef CStatus  cStatus = CStatus_OK
        cdef CReal    Zfolded
        if self.cObject.nstates > ANALYTIC_STATES:
            raise CEModelError ( "Maximum number of states for analytic treatment ({:d}) exceeded.".format ( ANALYTIC_STATES ) )

        Zfolded = EnergyModel_CalculateZfolded (self.cObject, pH, Gzero, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ("Cannot allocate Boltzmann factors.")
        return Zfolded


    def CalculateZunfolded (self, CReal Gzero=0.0, CReal pH=7.0):
        """Calculate partition function of an unfolded protein."""
        cdef CStatus  cStatus = CStatus_OK
        cdef CReal    Zunfolded
        if self.cObject.nstates > ANALYTIC_STATES:
            raise CEModelError ( "Maximum number of states for analytic treatment ({:d}) exceeded.".format ( ANALYTIC_STATES ) )

        Zunfolded = EnergyModel_CalculateZunfolded (self.cObject, pH, Gzero, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ("Cannot allocate Boltzmann factors.")
        return Zunfolded

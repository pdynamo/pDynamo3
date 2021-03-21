from  pCore        import logFile       , \
                          LogFileActive
from .CEModelError import CEModelError

DEF ANALYTIC_STATES = 67108864

cdef class StateVector:
    """A class defining the state vector."""

    def __len__ (self):
        """Return the size of the vector."""
        return self.cObject.nsites

    def __dealloc__ (self):
        """Deallocate."""
        StateVector_Deallocate (&self.cObject)

    def __copy__ (self):
        """Copying."""
        pass

    def __deepcopy__ (self):
        """Copying."""
        pass

    def __getitem__ (self, CInteger index):
        """Get an item."""
        cdef CInteger item
        cdef CStatus  cStatus = CStatus_OK
        item = StateVector_GetItem (self.cObject, index, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ( "Index ({:d}) out of range.".format ( index ) )
        return item


    def __setitem__ (self, CInteger index, CInteger value):
        """Set an item."""
        cdef CStatus cStatus = CStatus_OK
        StateVector_SetItem (self.cObject, index, value, &cStatus)
        if cStatus != CStatus_OK:
            if cStatus == CStatus_IndexOutOfRange:
                raise CEModelError ( "Index ({:d}) out of range.".format ( index ) )
            else: # Status_ValueError
                raise CEModelError ( "Invalid value ({:d}).".format ( value ) )


    def GetActualItem (self, CInteger index):
        """Get an actual item."""
        cdef CInteger item
        cdef CStatus  cStatus = CStatus_OK
        item = StateVector_GetActualItem (self.cObject, index, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ( "Index ({:d}) out of range.".format ( index ) )
        return item


    def SetActualItem (self, CInteger index, CInteger value):
        """Set an actual item."""
        cdef CStatus cStatus = CStatus_OK
        StateVector_SetActualItem (self.cObject, index, value, &cStatus)
        if cStatus != CStatus_OK:
            if cStatus == CStatus_IndexOutOfRange:
                raise CEModelError ( "Index ({:d}) out of range.".format ( index ) )
            else: # Status_ValueError
                raise CEModelError ( "Invalid value ({:d}).".format ( value ) )


    def __init__ (self, ceModel):
        """Constructor."""
        cdef CInteger numberOfSites
        cdef CInteger indexSite, indexDown, indexUp, index
        cdef CStatus  cStatus = CStatus_OK
        numberOfSites = ceModel.nsites
        self.cObject  = StateVector_Allocate (numberOfSites, &cStatus)

        if cStatus != CStatus_OK:
            raise CEModelError ("Cannot allocate state vector.")

        for indexSite from 0 <= indexSite < numberOfSites:
            indexUp   = 0
            indexDown = 9999
            site      = ceModel.sites[indexSite]
            for instance in site.instances:
                index = instance._instIndexGlobal
                if index < indexDown : indexDown = index
                if index > indexUp   : indexUp   = index
            StateVector_SetSite (self.cObject, indexSite, indexDown, indexUp, &cStatus)

        self.isOwner     = False
        self.owningModel = ceModel


    def CopyTo (StateVector self, StateVector other):
        """In-place copy of one state vector to another."""
        pass

    def Increment (self):
        """Generate a new state incrementing the state vector.

        Return False after the incrementation has finished."""
        cdef CBoolean moreStates = StateVector_Increment (self.cObject)
        if moreStates != CTrue:
            return False
        else:
            return True


    def Reset (self):
        """Set all components of the vector to their minimum values (formerly zeros)."""
        StateVector_Reset (self.cObject)


    def ResetToMaximum (self):
        """Set all components of the vector to their maximum values."""
        StateVector_ResetToMaximum (self.cObject)


    def Print (self, verbose=True, title=None, log=logFile):
        """Print the state vector."""
        cdef CInteger indexSite, indexSiteSubstate, indexActive
        cdef CBoolean isSubstate
        cdef CStatus  cStatus = CStatus_OK

        if LogFileActive (log):
            tab = log.GetTable (columns=[7, 7, 7, 8, 8, 2])
            tab.Start ()
            if title : tab.Title (title)
            else     : tab.Title ("State vector")
            tab.Heading ("Site"     , columnSpan=3)
            tab.Heading ("Instance" , columnSpan=3)

            for indexSite from 0 <= indexSite < self.cObject.nsites:
                indexActive = StateVector_GetItem (self.cObject, indexSite, &cStatus)
                isSubstate  = StateVector_IsSubstate (self.cObject, indexSite, &cStatus)
                ceModel     = self.owningModel
                site        = ceModel.sites[indexSite]
                instance    = site.instances[indexActive]
                substate    = "@" if isSubstate == CTrue else " "

                tab.Entry (site.segName)
                tab.Entry (site.resName)
                tab.Entry ( "{:d}".format ( site.resSerial ) )
                tab.Entry (instance.label)
                tab.Entry ( "{:d}".format ( instance.instIndex ) )
                tab.Entry (substate)
            tab.Stop ()


# Try to further improve this method with data types from C
    def DefineSubstate (self, selectedSites):
        """Define a substate.

        |selectedSites| is a sequence of two-element sequences (segmentName, residueSerial)"""
        cdef CBoolean foundSite
        cdef CInteger siteIndex, substateSiteIndex, nselected, nsites
        cdef CStatus  cStatus = CStatus_OK
        ceModel    = self.owningModel
        nsites     = len (ceModel.sites)
        nselected  = len (selectedSites)

        StateVector_AllocateSubstate (self.cObject, nselected, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ("Memory allocation failure.")

        for substateSiteIndex from 0 <= substateSiteIndex < nselected:
            selectedSegment, selectedSerial = selectedSites[substateSiteIndex]
            foundSite = CFalse

            for siteIndex from 0 <= siteIndex < nsites:
                site = ceModel.sites[siteIndex]
                if site.segName == selectedSegment and site.resSerial == selectedSerial:
                    StateVector_SetSubstateItem (self.cObject, siteIndex, substateSiteIndex, &cStatus)
                    foundSite = CTrue
                    break

            if foundSite == CFalse:
                raise CEModelError ( "Site {:s} {:d} not found.".format ( selectedSegment, selectedSerial ) )


    def IncrementSubstate (self):
        """Generate a new substate.

        Return False after the incrementation has finished."""
        cdef CBoolean moreStates = StateVector_IncrementSubstate (self.cObject)
        if moreStates != CTrue:
            return False
        else:
            return True


    def ResetSubstate (self):
        """Set all components of the substate to their minimum values (formerly zeros)."""
        StateVector_ResetSubstate (self.cObject)

from  pCore        import Align         , \
                          logFile       , \
                          LogFileActive
from .CEModelError import CEModelError

_DefaultDoubleFlip         = 2.
_DefaultProductionScans    = 20000
_DefaultEquilibrationScans = 500


cdef class MCModelDefault:
    """A class defining the default Monte Carlo model."""

    def __dealloc__ (self):
        """Deallocate."""
        MCModelDefault_Deallocate (&self.cObject)


    def __init__ (self, CReal doubleFlip=_DefaultDoubleFlip, CInteger nequil=_DefaultEquilibrationScans, CInteger nprod=_DefaultProductionScans, CInteger randomSeed=-1):
        """Constructor."""
        cdef CStatus cStatus = CStatus_OK
        self.isOwner  = True
        self.cObject  = MCModelDefault_Allocate (doubleFlip, nequil, nprod, randomSeed, &cStatus)
        if cStatus != CStatus_OK: raise CEModelError ("Cannot allocate Monte Carlo model.")


    def Initialize (self, ceModel):
        """Link Monte Carlo and energy models."""
        cdef EnergyModel energyModel = ceModel.energyModel
        cdef CStatus      cStatus      = CStatus_OK
        cdef CInteger     npairs

        MCModelDefault_LinkToEnergyModel (self.cObject, energyModel.cObject, &cStatus)
        if cStatus != CStatus_OK:
            raise CEModelError ("Cannot initialize Monte Carlo model.")

        npairs = MCModelDefault_FindPairs (self.cObject, -1, &cStatus)
        if npairs > 0:
            MCModelDefault_FindPairs (self.cObject, npairs, &cStatus)
            if cStatus != CStatus_OK:
                raise CEModelError ("Cannot allocate pairs.")
        self.isOwner     = False
        self.owningModel = ceModel


    def CalculateOwnerProbabilities (self, CReal pH=7.0, CInteger logFrequency=0, trajectoryFilename="", log=logFile):
        """Calculate probabilities of the owning model."""
        cdef CMCModelDefault   *mcModel
        cdef CInteger  moves, movesAcc, flips, flipsAcc
        cdef CInteger  nmoves, i, j
        cdef CReal     scale, Gmicro
        cdef CBoolean  active, writing

        if not self.owningModel.isCalculated:
            raise CEModelError ("First calculate electrostatic energies.")

        if (logFrequency <= 0) and (trajectoryFilename == ""):
            # . Do a quiet run
            active = CTrue if (LogFileActive (log)) else CFalse

            MCModelDefault_Equilibration (self.cObject, pH)
            if active:
                log.Paragraph ( "Completed {:d} equilibration scans.".format ( self.cObject.nequil ) )

            MCModelDefault_Production (self.cObject, pH)
            if active:
                log.Paragraph ( "Completed {:d} production scans.".format ( self.cObject.nprod ) )
        else:
            # . Do logging or trajectory writing
            mcModel = self.cObject
            nmoves  = mcModel.vector.nsites + mcModel.vector.npairs
            active  = CTrue if (LogFileActive (log) and (logFrequency > 0)) else CFalse
    
            if active:
                j         =  0
                moves     =  0
                flips     =  0
                movesAcc  =  0
                flipsAcc  =  0
                table  = log.GetTable (columns=[8, 10, 10, 10, 10])
                table.Start ()
                table.Heading ( "Scan"  )
                table.Heading ( "Moves" )
                table.Heading ( "M-Acc" )
                table.Heading ( "Flips" )
                table.Heading ( "F-Acc" )
    
            # . Do equilibration phase
            StateVector_Randomize (mcModel.vector, mcModel.generator)

            for i from 0 <= i < mcModel.nequil:
                MCModelDefault_MCScan (mcModel, pH, nmoves, &moves, &movesAcc, &flips, &flipsAcc)
                if active:
                    j += 1
                    if j >= logFrequency:
                        table.Entry ( "{:d}".format ( i + 1    ) )
                        table.Entry ( "{:d}".format ( moves    ) )
                        table.Entry ( "{:d}".format ( movesAcc ) )
                        table.Entry ( "{:d}".format ( flips    ) )
                        table.Entry ( "{:d}".format ( flipsAcc ) )
                        j         =  0
                        moves     =  0
                        flips     =  0
                        movesAcc  =  0
                        flipsAcc  =  0
            if active:
                j        =  0
                moves    =  0
                flips    =  0
                movesAcc =  0
                flipsAcc =  0
                table.Entry ("** Equilibration phase done **".center (8 + 10 * 4), columnSpan=5)
   
            # . Do production phase
            writing = CFalse
            if (trajectoryFilename != ""):
                writing = CTrue
                output = open (trajectoryFilename, "w")
            RealArray1D_Set (mcModel.energyModel.probabilities, 0.)
 
            for i from 0 <= i < mcModel.nprod:
                Gmicro = MCModelDefault_MCScan (mcModel, pH, nmoves, &moves, &movesAcc, &flips, &flipsAcc)
                if writing:
                    output.write ( "{:f}\n".format ( Gmicro ) )

                MCModelDefault_UpdateProbabilities (mcModel)
                if active:
                    j += 1
                    if j >= logFrequency:
                        table.Entry ( "{:d}".format ( i + 1    ) )
                        table.Entry ( "{:d}".format ( moves    ) )
                        table.Entry ( "{:d}".format ( movesAcc ) )
                        table.Entry ( "{:d}".format ( flips    ) )
                        table.Entry ( "{:d}".format ( flipsAcc ) )
                        j         =  0
                        moves     =  0
                        flips     =  0
                        movesAcc  =  0
                        flipsAcc  =  0
            if active:
                table.Entry ("** Production phase done **".center (8 + 10 * 4), columnSpan=5)
                table.Stop ()

            # . Finalize
            if writing:
                output.close ()
            scale = 1. / mcModel.nprod
            RealArray1D_Scale (mcModel.energyModel.probabilities, scale)


    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ): log.SummaryOfItems ( self.SummaryItems ( ), title = "Default Monte Carlo Sampling Model" )

    def SummaryItems ( self ):
        """Summary items."""
        cdef CInteger nequil     = 0
        cdef CInteger nprod      = 0
        cdef CReal    doubleFlip = 0.0
        if self.cObject != NULL:
            nequil     = self.cObject.nequil
            nprod      = self.cObject.nprod
            doubleFlip = self.cObject.limit
        return [ ( "Equilibration scans"    , "{:d}".format   ( nequil     ) ) ,
                 ( "Production scans"       , "{:d}".format   ( nprod      ) ) ,
                 ( "Limit for double moves" , "{:.1f}".format ( doubleFlip ) ) ]

    def PrintPairs ( self, log = logFile ):
        """Summary of strongly interacting pairs."""
        cdef CInteger indexPair, indexSiteA, indexSiteB
        cdef CInteger npairs
        cdef CReal    Wmax
        cdef CStatus  cStatus = CStatus_OK
        npairs = self.cObject.vector.npairs
        if LogFileActive ( log ):
            if npairs == 0:
                log.Paragraph ( "There were no pairs of strongly interacting sites." )
            else:
                ceModel = self.owningModel
                sites   = ceModel.sites
                table   = log.GetTable ( columns = [ 8, 8, 6, 4, 8, 8, 6, 12 ] )
                table.Start   ( )
                table.Title   ( "Pairs of Strongly Interacting Sites" )
                table.Heading ( "Site 1", columnSpan = 3 )
                table.Heading ( "" )
                table.Heading ( "Site 2", columnSpan = 3 )
                table.Heading ( "Energy" )
                for indexPair from 0 <= indexPair < npairs:
                    StateVector_GetPair ( self.cObject.vector, indexPair, &indexSiteA, &indexSiteB, &Wmax, &cStatus )
                    siteFirst  = sites[indexSiteA]
                    siteSecond = sites[indexSiteB]
                    table.Entry ( siteFirst.segName, align = Align.Left  )
                    table.Entry ( siteFirst.resName, align = Align.Left  )
                    table.Entry ( "{:d}".format ( siteFirst.resSerial  ) )
                    table.Entry ( "" )
                    table.Entry ( siteSecond.segName, align = Align.Left )
                    table.Entry ( siteSecond.resName, align = Align.Left )
                    table.Entry ( "{:d}".format ( siteSecond.resSerial ) )
                    table.Entry ( "{:.2f}".format ( Wmax ) )
                table.Stop ( )

"""The base class for QC/MM electrostatic models."""

from   pCore                 import logFile               , \
                                    LogFileActive
from   pScientific.Arrays    import Array
from   pScientific.Geometry3 import Coordinates3
from  .NBModelError          import NBModelError
from ..EnergyModel           import EnergyClosurePriority , \
                                    EnergyModel           , \
                                    EnergyModelState

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelState ( EnergyModelState ):
    """A QC/MM electrostatic model state."""

    _attributable = dict ( EnergyModelState._attributable )
    _attributable.update ( { "baMMPartners" : None  ,
                             "bpCharges"    : None  } )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.bpCharges is None: n = 0
        else:                      n = len ( self.bpCharges )
        return [ ( "Number of BP Atoms" , "{:d}".format ( n ) ) ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModel ( EnergyModel ):
    """The base class for QC/MM electrostatic models."""

    _attributable             = dict ( EnergyModel._attributable )
    _classLabel               = "QC/MM Electrostatic Model"
    _pairwiseInteractionClass = None
    _stateName                = "qcmmState"
    _stateObject              = QCMMElectrostaticModelState
    _summarizable             = dict ( EnergyModel._summarizable )
    _attributable.update ( { "dielectric"          : 1.0   ,
                             "pairwiseInteraction" : None  ,
                             "rdLACharges"         : False } )
    _summarizable.update ( { "dielectric"          : ( "Dielectric", "{:.3f}" ) ,
                             "pairwiseInteraction" : None                       ,
                             "rdLACharges"         : "RD Approximation"       } )

    def BuildModel ( self, target ):
        """Build the model."""
        state = super ( QCMMElectrostaticModel, self ).BuildModel ( target )
        if ( target.mmModel is None ) or ( target.qcModel is None ):
            raise NBModelError ( "Missing QC or MM model in QC/MM model definition." )
        baMMPartners = getattr ( target.qcState, "baMMPartners", None )
        if ( baMMPartners is not None ) and ( len ( baMMPartners ) > 0 ):
            state.baMMPartners = baMMPartners
            state.bpCharges    = self.GetBPCharges ( baMMPartners, target.mmState.charges )
        return state

    def _CheckOptions ( self ):
        """Check options."""
        pwiClass = self.__class__._pairwiseInteractionClass
        if pwiClass is not None:
            if self.pairwiseInteraction is None:
                self.pairwiseInteraction = pwiClass ( )
            elif not isinstance ( self.pairwiseInteraction, pwiClass ):
                raise TypeError ( "Invalid pairwise interaction attribute." ) 

    def ClearScratch ( self, scratch ):
        """Clear scratch."""
        scratch.Clear ( )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.EnergyInitialize ( target )
        def b ( ): self.QCMMPotentials   ( target )
        def c ( ): self.QCMMGradients    ( target )
        def d ( ): self.EnergyFinalize   ( target )
        return [ ( EnergyClosurePriority.QCInitialization, a, "QC/MM Electrostatic Initialization" ) ,
                 ( EnergyClosurePriority.QCIntegrals     , b, "QC/MM Electrostatic Potentials"     ) ,
                 ( EnergyClosurePriority.QCGradients     , c, "QC/MM Electrostatic Gradients"      ) ,
                 ( EnergyClosurePriority.QCFinalization  , d, "QC/MM Electrostatic Finalization"   ) ]

    def EnergyFinalize ( self, target ):
        """Energy finalization."""
        scratch = target.scratch
        if scratch.doGradients:
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            target.qcModel.ScatterQCGradients3 ( target.qcState, coordinates3, scratch.qcGradients3QCMM, scratch.gradients3 )
            state     = getattr ( target, self.__class__._stateName )
            bpCharges = getattr ( state, "bpCharges", None )
            if bpCharges is not None:
                self.ScatterBPGradients3 ( state.baMMPartners, scratch.bpGradients3, scratch.gradients3 )

    def EnergyInitialize ( self, target ):
        """Energy initialization."""
        scratch      = target.scratch
        coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
        n            = len ( target.qcState.atomicNumbers )
        crd3         = scratch.Get ( "qcCoordinates3QCMM", None )
        if crd3 is None:
            crd3 = Coordinates3.WithExtent ( n )
            scratch.qcCoordinates3QCMM = crd3
        target.qcModel.GatherQCCoordinates3 ( target.qcState, coordinates3, crd3 )
        if scratch.doGradients:
           grd3 = scratch.Get ( "qcGradients3QCMM", None )
           if grd3 is None:
               grd3 = Coordinates3.WithExtent ( n )
               scratch.qcGradients3QCMM = grd3
           grd3.Set ( 0.0 )
        state = getattr ( target, self.__class__._stateName )
        if state.bpCharges is not None:
            n    = len ( state.bpCharges )
            crd3 = scratch.Get ( "bpCoordinates3", None )
            if crd3 is None:
                crd3 = Coordinates3.WithExtent ( n )
                scratch.bpCoordinates3 = crd3
            self.GatherBPCoordinates3 ( state.baMMPartners, coordinates3, crd3 )
            if scratch.doGradients:
                grd3 = scratch.Get ( "bpGradients3", None )
                if grd3 is None:
                    grd3 = Coordinates3.WithExtent ( n )
                    scratch.bpGradients3 = grd3
                grd3.Set ( 0.0 )

    def GatherBPCoordinates3 ( self, baMMPartners, coordinates3, bpCoordinates3 ):
        """Gather the BP coordinates without unit conversion."""
        n = 0
        for ( b, mm ) in baMMPartners:
            xB = coordinates3[b,0] ; yB = coordinates3[b,1] ; zB = coordinates3[b,2]
            for p in mm:
                xP = coordinates3[p,0] ; yP = coordinates3[p,1] ; zP = coordinates3[p,2]
                bpCoordinates3[n,0] = 0.5 * ( xB + xP )
                bpCoordinates3[n,1] = 0.5 * ( yB + yP )
                bpCoordinates3[n,2] = 0.5 * ( zB + zP )
                n += 1
            if self.rdLACharges:
                for p in mm:
                    xP = coordinates3[p,0]   ; yP = coordinates3[p,1]   ; zP = coordinates3[p,2]
                    bpCoordinates3[n,0] = xP ; bpCoordinates3[n,1] = yP ; bpCoordinates3[n,2] = zP
                    n += 1

    def GetBPCharges ( self, baMMPartners, charges ):
        """Return the charges of the BP atoms."""
        bpCharges = None
        n         = 0
        for ( b, mm ) in baMMPartners: n += len ( mm )
        if self.rdLACharges: n *= 2
        if n > 0:
            bpCharges = Array.WithExtent ( n ) ; bpCharges.Set ( 0.0 )
            n = 0
            for ( b, mm ) in baMMPartners:
                nP = len ( mm )
                dQ = charges[b] / float ( nP )
                if self.rdLACharges:
                    for a in range ( nP ):
                        bpCharges[n] = 2.0 * dQ ; n += 1
                    for a in range ( nP ):
                        bpCharges[n] =     - dQ ; n += 1
                else:
                    for a in range ( nP ):
                        bpCharges[n] = dQ ; n += 1
        return bpCharges

    def QCMMGradients ( self, target ):
        """Calculate the QC/MM electrostatic gradients."""
        pass

    def QCMMPotentials ( self, target ):
        """Calculate the QC/MM electrostatic potentials."""
        pass

    def ScatterBPGradients3 ( self, baMMPartners, bpGradients3, gradients3 ):
        """Scatter the BP gradients."""
        n = 0
        for ( b, mm ) in baMMPartners:
            for p in mm:
                gX = bpGradients3[n,0] ; gY = bpGradients3[n,1] ; gZ = bpGradients3[n,2]
                gX *= 0.5e+00 ; gY *= 0.5e+00 ; gZ *= 0.5e+00
                gradients3[b,0] += gX ; gradients3[b,1] += gY ; gradients3[b,2] += gZ
                gradients3[p,0] += gX ; gradients3[p,1] += gY ; gradients3[p,2] += gZ
                n += 1
            if self.rdLACharges:
                for p in mm:
                    gX = bpGradients3[n,0] ; gY = bpGradients3[n,1] ; gZ = bpGradients3[n,2]
                    gradients3[p,0] += gX ; gradients3[p,1] += gY ; gradients3[p,2] += gZ
                    n += 1

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

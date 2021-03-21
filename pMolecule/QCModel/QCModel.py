"""The base class for QC models."""

import math

from   pCore                 import DataType              , \
                                    logFile               , \
                                    LogFileActive         , \
                                    Selection
from   pScientific           import PeriodicTable         , \
                                    Units
from   pScientific.Arrays    import Array
from   pScientific.Geometry3 import Coordinates3
from  .ElectronicState       import ElectronicState
from  .QCModelError          import QCModelError
from ..EnergyModel           import EnergyClosurePriority , \
                                    EnergyModel           , \
                                    EnergyModelState

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelState ( EnergyModelState ):
    """A QC model state."""

    _attributable = dict ( EnergyModelState._attributable )
    _attributable.update ( { "atomicNumbers" : None ,
                             "baMMPartners"  : None ,
                             "baQCPartners"  : None ,
                             "boundaryAtoms" : None ,
                             "pureQCAtoms"   : None ,
                             "qcAtoms"       : None } )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.target is not None:
            if self.atomicNumbers is None: n = 0
            else:                          n = len ( self.atomicNumbers )
            items.append ( ( "Number of QC Atoms" , "{:d}".format ( n ) ) )
            if self.boundaryAtoms is None: n = 0
            else:                          n = len ( self.boundaryAtoms )
            items.append ( ( "Number of Boundary Atoms", "{:d}".format ( n ) ) )
        return items

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModel ( EnergyModel ):
    """The base class for quantum chemical models.

    This class should not be used directly.
    """

    _attributable = dict ( EnergyModel._attributable )
    _classLabel   = "QC Model"
    _stateName    = "qcState"
    _stateObject  = QCModelState
    _summarizable = dict ( EnergyModel._summarizable )
    _attributable.update ( { "relativeLADistance" : True } )
    _summarizable.update ( { "relativeLADistance" : "Relative Link Atom Distance" } )

    def AtomicCharges ( self, target, chargeModel = None ):
        """Atomic charges for QC atoms only."""
        return None

    def AtomicSpins ( self, target, chargeModel = None ):
        """Atomic spins."""
        return None

    def BondOrders ( self, target, chargeModel = None ):
        """Bond orders."""
        return None

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModel, self ).BuildModel ( target )
        # . Default electronic state.
        if target.electronicState is None: target.electronicState = ElectronicState ( )
        # . Atom indices.
        if qcSelection is None:
            atomIndices = Selection.FromIterable ( range ( len ( target.atoms ) ) )
        else:
            atomIndices = qcSelection
            if   len ( atomIndices ) <= 0:
                raise QCModelError ( "There are no QC atoms in the QC model." )
            elif max ( atomIndices ) >= len ( target.atoms ):
                raise QCModelError ( "Atom index out of range in QC atom selection." )
        # . Boundary atom dictionary of BA index: BA partners, MM partners, QC partners.
        # . Boundary atoms with multiple QC partners are not permitted.
        if target.mmModel is None: boundaryAtoms = {}
        else:                      boundaryAtoms = target.mmState.IdentifyBoundaryAtoms ( atomIndices )
        for ( b, partners ) in boundaryAtoms.items ( ):
            if len ( partners[2] ) > 1: raise QCModelError ( "A QC/MM boundary atom - {:d} - has multiple QC partners - {:s}.".format ( b, repr ( partners[2] ) ) )
        # . Get selections.
        state.boundaryAtoms = Selection.FromIterable ( boundaryAtoms.keys ( ) )
        state.pureQCAtoms   = atomIndices
        if len ( state.boundaryAtoms ) == 0: state.qcAtoms = atomIndices
        else:                                state.qcAtoms = Selection.FromIterable ( set ( atomIndices ).union ( set ( state.boundaryAtoms ) ) )
        # . Get boundary atom information.
        state.baMMPartners = []
        state.baQCPartners = []
        for b in state.boundaryAtoms:
            ( ba, mm, qc ) = boundaryAtoms[b]
            m   = list ( ba.union ( mm ) ) ; m.sort ( )
            q   = qc.pop ( )
            nQ  = target.atoms[q].atomicNumber
            eQ  = PeriodicTable[nQ]
            rQL = eQ.GetSingleBondDistance (  1 )
            if rQL is None: raise QCModelError ( "Unknown link atom distance for elements with atomic numbers {:d} and 1.".format ( nQ ) )
            if self.relativeLADistance:
                nM  = target.atoms[b].atomicNumber
                rQM = eQ.GetSingleBondDistance ( nM )
                if rQM is None: raise QCModelError ( "Unknown link atom distance for elements with atomic numbers {:d} and {:d}.".format ( nQ, nM ) )
                d = rQL / rQM
            else:
                d = rQL
            state.baMMPartners.append ( ( b, m    ) )
            state.baQCPartners.append ( ( b, q, d ) )
        # . Set atomic numbers - boundary link atoms are always hydrogens for the moment.
        n = len ( state.pureQCAtoms )
        state.atomicNumbers = Array.WithExtent ( n + len ( state.boundaryAtoms ), dataType = DataType.Integer )
        for ( i, s ) in enumerate ( state.pureQCAtoms ): state.atomicNumbers[i]   = target.atoms[s].atomicNumber
        for i in range ( len ( state.boundaryAtoms )  ): state.atomicNumbers[i+n] = 1
        # . Adjust the MM model.
        if target.mmModel is not None:
            target.mmState.DeactivateQCAtoms ( state.pureQCAtoms, state.boundaryAtoms )
            target.mmState.CheckActiveAtomTotalCharge ( )
        # . Finish up.
        return state

    def ClearScratch ( self, scratch ):
        """Clear scratch."""
        scratch.Clear ( )

    def DipoleMoment ( self, target, center = None, dipole = None ):
        """The dipole moment in Debyes."""
        return None

    def Energy ( self, target ):
        """Calculate the QC energy."""
        pass

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.EnergyInitialize ( target )
        def b ( ): self.Energy           ( target )
        def c ( ): self.EnergyFinalize   ( target )
        return [ ( EnergyClosurePriority.QCInitialization, a, "QC Initialization" ) ,
                 ( EnergyClosurePriority.QCEnergy        , b, "QC Energy"         ) ,
                 ( EnergyClosurePriority.QCFinalization  , c, "QC Finalization"   ) ]

    def EnergyFinalize ( self, target ):
        """Energy finalization."""
        scratch = target.scratch
        if scratch.doGradients:
            state        = getattr ( target, self.__class__._stateName )
            qcGradients3 = scratch.qcGradients3AU
            qcGradients3.Scale ( Units.Length_Angstroms_To_Bohrs * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
            self.ScatterQCGradients3 ( state, target.coordinates3, qcGradients3, scratch.gradients3 )

    def EnergyInitialize ( self,target ):
        """Energy initialization."""
        scratch = target.scratch
        state   = getattr ( target, self.__class__._stateName )
        n       = len ( state.atomicNumbers )
        crd3    = scratch.Get ( "qcCoordinates3"  , None )
        crd3AU  = scratch.Get ( "qcCoordinates3AU", None )
        if crd3 is None:
            crd3 = Coordinates3.WithExtent ( n )
            scratch.qcCoordinates3 = crd3
        if crd3AU is None:
            crd3AU = Coordinates3.WithExtent ( n )
            scratch.qcCoordinates3AU = crd3AU
        self.GatherQCCoordinates3 ( state, target.coordinates3, crd3 )
        crd3.CopyTo  ( crd3AU )
        crd3AU.Scale ( Units.Length_Angstroms_To_Bohrs )
        if scratch.doGradients:
           grd3 = scratch.Get ( "qcGradients3AU", None )
           if grd3 is None:
               grd3 = Coordinates3.WithExtent ( n )
               scratch.qcGradients3AU = grd3
           grd3.Set ( 0.0 )

    def GatherQCCoordinates3 ( self, state, coordinates3, qcCoordinates3 ):
        """Gather the QC coordinates without unit conversion."""
        qcCoordinates3.Gather ( coordinates3, selection = state.pureQCAtoms )
        n = len ( state.pureQCAtoms )
        for ( b, q, d ) in state.baQCPartners:
            x  = coordinates3[q,0]     ; y  = coordinates3[q,1]     ; z  = coordinates3[q,2]
            dX = coordinates3[b,0] - x ; dY = coordinates3[b,1] - y ; dZ = coordinates3[b,2] - z
            if not self.relativeLADistance: d /= math.sqrt ( dX*dX + dY*dY + dZ*dZ )
            qcCoordinates3[n,0] = ( x + d * dX )
            qcCoordinates3[n,1] = ( y + d * dY )
            qcCoordinates3[n,2] = ( z + d * dZ )
            n += 1

    def ModifyElectronicState ( self, target ):
        """Modify the electronic state."""
        target.scratch.Clear ( )

    def ScatterQCGradients3 ( self, state, coordinates3, qcGradients3, gradients3 ):
        """Scatter the QC gradients."""
        qcGradients3.ScatterAdd ( 1.0, gradients3, selection = state.pureQCAtoms )
        n = len ( state.pureQCAtoms )
        for ( b, q, d ) in state.baQCPartners:
            gX = qcGradients3[n,0] ; gY = qcGradients3[n,1] ; gZ = qcGradients3[n,2]
            gradients3[q,0] += gX  ; gradients3[q,1] += gY  ; gradients3[q,2] += gZ
            # . Ratio of distances.
            if self.relativeLADistance:
                gX *= d ; gY *= d ; gZ *= d
            # . Absolute distance.
            else:
                dX  = coordinates3[b,0] - coordinates3[q,0]
                dY  = coordinates3[b,1] - coordinates3[q,1]
                dZ  = coordinates3[b,2] - coordinates3[q,2]
                r   = math.sqrt ( dX*dX + dY*dY + dZ*dZ ) ;
                dX /= r ;
                dY /= r ;
                dZ /= r ;
                gX *= ( d / r )
                gY *= ( d / r )
                gZ *= ( d / r )
                dG  = dX * gX + dY * gY + dZ * gZ
                gX -= dG * dX
                gY -= dG * dY
                gZ -= dG * dZ
            gradients3[b,0] += gX ; gradients3[b,1] += gY ; gradients3[b,2] += gZ
            gradients3[q,0] -= gX ; gradients3[q,1] -= gY ; gradients3[q,2] -= gZ
            n += 1

    def UnbuildModel ( self, target ):
        """Unbuild the model."""
        super ( QCModel, self ).UnbuildModel ( target )
        if target.mmModel is not None:
            target.mmState.ActivateAllAtoms ( )
            if target.freeAtoms is not None:
                target.mmState.FixAtoms ( target.freeAtoms )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

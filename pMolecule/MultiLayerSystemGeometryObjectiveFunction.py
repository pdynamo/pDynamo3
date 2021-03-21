"""Multi-layer system geometry objective function."""

from  pCore                           import Clone, logFile, LogFileActive
from  pScientific                     import PeriodicTable
from  pScientific.Geometry3           import Coordinates3
from .System                          import System
from .SystemGeometryObjectiveFunction import SystemGeometryObjectiveFunction

#===================================================================================================================================
#
# . Layers comprise non-interacting subsets of the original system arranged hierarchically.
#
# . Steps:
#
# - Define full system as normal.
# - Define each layer with atom selection and energy model.
#
# . At the moment, subSystems must have QC energy models only and the first-layer atoms of the original system must also be all
# . in the QC region. This avoids problems of having to construct an MM or QC/MM subSystem with link atoms (when an MM model
# . which allows building is unavailable).
#
# . Probably the easiest way around this (at least initially) is to construct MM or QC/MM subSystems by hand and then add them to
# . the list of subSystems explicitly.
#
#===================================================================================================================================

#===================================================================================================================================
# . Boundary atom class.
#===================================================================================================================================
class BoundaryAtom:
    """A boundary atom between layers."""

    # . Link-ratio method used (bond length is ratio of partner-H/partner-boundary atom lengths).

    def __init__ ( self, index, linkRatio, rIndex, rPartner ):
        """Constructor."""
        self.index     = index
        self.linkRatio = linkRatio
        self.rIndex    = rIndex
        self.rPartner  = rPartner

    def GradientsGet ( self, gradients3, reference3 ):
        """Put the gradients into those of the reference."""
        b = self.index
        d = self.linkRatio
        e = 1.0 - d
        l = self.rIndex
        u = self.rPartner
        gx = gradients3[b,0]
        gy = gradients3[b,1]
        gz = gradients3[b,2]
        reference3[l,0] += d * gx
        reference3[l,1] += d * gy
        reference3[l,2] += d * gz
        reference3[u,0] += e * gx
        reference3[u,1] += e * gy
        reference3[u,2] += e * gz

    def VariablesPut ( self, coordinates3, reference3 ):
        """Get the coordinates from those of the reference."""
        b = self.index
        d = self.linkRatio
        l = self.rIndex
        u = self.rPartner
        xu  = reference3[u,0]
        yu  = reference3[u,1]
        zu  = reference3[u,2]
        xlu = reference3[l,0] - xu
        ylu = reference3[l,1] - yu
        zlu = reference3[l,2] - zu
        coordinates3[b,0] = xu + xlu * d
        coordinates3[b,1] = yu + ylu * d
        coordinates3[b,2] = zu + zlu * d

#===================================================================================================================================
# . Subsystem class.
#===================================================================================================================================
# . This should work no matter what the subSystem's energy model.
class MultiLayerSubsystem:
    """A multilayer subSystem."""

    def __init__ ( self, system, weight, selection, boundaryAtoms = [] ):
        """Constructor."""
        self.boundaryAtoms = boundaryAtoms
        self.selection     = selection
        self.system        = system
        self.weight        = weight

    def Energy ( self, doGradients = False, log = logFile ):
        """Energy and gradient calculation."""
        return ( self.weight * self.system.Energy ( doGradients = doGradients, log = log ) )

    def GradientsGet ( self, gradients3 ):
        """Put the gradients into those of system."""
        grd3 = self.system.scratch.gradients3
        grd3.Scale ( self.weight ) # . Scale - OK in place as unlikely to need subSystem gradients again.
        # . Non-boundary atoms.
        for ( i, s ) in enumerate ( self.selection ):
            for c in range ( 3 ): gradients3[s,c] += grd3[i,c] #. Increment.
        # . Boundary atoms.
        for atom in self.boundaryAtoms: atom.GradientsGet ( grd3, gradients3 )

    def VariablesPut ( self, coordinates3 ):
        """Get the variables of the subSystem from those of the complete system."""
        crd3 = self.system.coordinates3
        # . Non-boundary atoms.
        crd3.Gather ( coordinates3, selection = self.selection )
        # . Boundary atoms.
        for atom in self.boundaryAtoms: atom.VariablesPut ( crd3, coordinates3 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MultiLayerSystemGeometryObjectiveFunction ( SystemGeometryObjectiveFunction ):
    """A multilayer system geometry objective function."""

    _attributable = dict ( SystemGeometryObjectiveFunction._attributable )
    _attributable.update ( { "subSystems" : list } )

    def DefineQCLayer ( self, selection, qcModel, electronicState = None ):
        """Define a QC layer.

        A layer consists of two subSystems with the same atomic composition.
        One uses the new energy model and a weight of +1 whereas the second
        has the old energy model and a weight of -1.
        """
        # . Get the latest defined system.
        if len ( self.subSystems ) > 0:
            number    = len ( self.subSystems )
            oldSystem = self.subSystems[-1].system
        else:
            number    = 0
            oldSystem = self.system
        # . Get the QC atom indices of the full system (excluding boundary atoms).
        try:    oldQCSet = set ( oldSystem.qcState.pureQCAtoms )
        except: raise ValueError ( "Unable to retrieve the previous system's QC atom indices." )
        # . Map a subSystem's indices to those of the full system.
        if len ( self.subSystems ) > 0:
            mapping  = self.subSystems[-1].selection
            unmapped = oldQCSet
            oldQCSet = set ( )
            for i in unmapped: oldQCSet.add ( mapping[i] )
        # . Check that the selection is a subset of the previous one.
        newQCSet = set ( selection )
        if not newQCSet.issubset ( oldQCSet ): raise ValueError ( "The atoms in the new QC layer must be a subset of the QC atoms in the underlying system." )
        # . Get the old QC model (which at this stage is known to exist).
        oldModel = oldSystem.qcModel
        # . Find the atomic numbers.
        atomicNumbers = [ self.system.atoms[s].atomicNumber for s in selection ]
        # . Find boundary atoms.
        boundaryAtoms = self.FindBoundaryAtoms ( selection, len ( atomicNumbers ) )
        # . Extend atomicNumbers by the appropriate number of boundary atom hydrogens.
        atomicNumbers.extend ( len ( boundaryAtoms ) * [ 1 ] )
        # . Define a basic QC system.
        system0 = System.FromAtoms ( atomicNumbers )
        if electronicState is not None: system0.electronicState = electronicState
        system0.coordinates3 = Coordinates3.WithExtent ( len ( system0.atoms ) )
        system0.coordinates3.Set ( 0.0 )
        # . Define a basic QC system - it is cleaner to do this for a QC system from scratch without pruning, etc.
        system0 = System.FromAtoms ( atomicNumbers )
        system0.coordinates3 = Coordinates3.WithExtent ( len ( system0.atoms ) )
        system0.coordinates3.Set ( 0.0 )
        # . Define the subsystems of the layer.
        system1 = Clone ( system0 )
        for ( i, ( system, model, weight ) ) in enumerate ( ( ( system0, oldModel, -1.0 ), ( system1, qcModel, 1.0 ) ) ):
            system.label = "MultiLayer Objective Function Subsystem {:d} with Weight {:.1f}".format ( number+i, weight )
            system.DefineQCModel ( model )
            subSystem = MultiLayerSubsystem ( system, weight, selection, boundaryAtoms = boundaryAtoms )
            subSystem.VariablesPut ( self.system.coordinates3 )
            self.subSystems.append ( subSystem )

    def FindBoundaryAtoms ( self, selection, increment ):
        """Find the boundary atoms for a selection."""
        # . Initialization.
        boundaryAtoms = []
        # . Get basic boundary atom data.
        # . Use connectivity data.
        if len ( self.system.connectivity.bonds ) > 0:
            baData = {}
            self.system.connectivity.IdentifyBoundaryAtoms ( selection, baData )
        # . Use MM model data.
        else:
            if self.system.mmModel is None: raise ValueError ( "Unable to identify system's boundary atoms due to absence of both connectivity and MM model data." )
            baData = self.system.mmModel.IdentifyBoundaryAtoms ( selection )
        # . Process data.
        indices = list ( baData.keys ( ) )
        indices.sort ( )
        for ( i, b ) in enumerate ( indices ):
            # . Get the partner index.
            upartners = baData[b][2]
            if len ( upartners ) > 1: raise ValueError ( "A lower-layer atom - {:d} - has multiple partners with those of the upper-layer - {!r}.".format ( b, upartners ) )
            p = upartners.pop ( )
            # . Get the link ratio.
            nB  = self.system.atoms[b].atomicNumber
            nP  = self.system.atoms[p].atomicNumber
            eP  = PeriodicTable[nP]
            rPH = eP.GetSingleBondDistance (  1 )
            rPB = eP.GetSingleBondDistance ( nB )
            if ( rPB is None ) or ( rPH is None ):
                raise ValueError ( "Unable to determine boundary-atom link-ratio for elements " + PeriodicTable.Symbol ( nB ) + ", " + PeriodicTable.Symbol ( nP ) + " and H." )
            linkRatio = rPH / rPB # . Ratio of standard bond lengths (p-H/p-b).
            # . Store the data.
            boundaryAtoms.append ( BoundaryAtom ( increment + i, linkRatio, b, p ) )
        # . Finish up.
        return boundaryAtoms

    def Function ( self, variables ):
        """Evaluate the function."""
        # . Full system.
        self.VariablesPut ( variables )
        f = self.system.Energy ( log = None )
        # . Subsystems.
        for subSystem in self.subSystems:
            subSystem.VariablesPut ( self.system.coordinates3 )
            f += subSystem.Energy ( log = None )
        return f

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        # . Full system.
        self.VariablesPut ( variables )
        f = self.system.Energy ( doGradients = True, log = None )
        # . Subsystems.
        for subSystem in self.subSystems:
            subSystem.VariablesPut ( self.system.coordinates3 )
            f += subSystem.Energy ( doGradients = True, log = None )
            subSystem.GradientsGet ( self.system.scratch.gradients3 )
        # . Gradients.
        self.GradientsGet ( gradients )
        return f

    def SubsystemSummary ( self, log = logFile ):
        """Write out a summary of the subSystems."""
        if LogFileActive ( log ):
           table = log.GetTable ( columns = [ 8, 16, 12, 12, 12, 40 ] )
           table.Start ( )
           table.Title ( "Multi-Layer Objective Function Subsystems" )
           table.Heading ( "Index"        )
           table.Heading ( "Total Atoms"  )
           table.Heading ( "Link Atoms"   )
           table.Heading ( "Weight"       )
           table.Heading ( "QC Charge"    )
           table.Heading ( "Energy Model" )
           # . System.
           table.Entry ( "0" )
           table.Entry ( "{:d}".format ( len ( self.system.atoms ) ) )
           table.Entry ( "0" )
           table.Entry ( "{:.3f}".format ( 1.0 ) )
           table.Entry ( "{:d}".format ( self.system.electronicState.charge ) )
           table.Entry ( self.system.energyModelLabel )
           # . Subsystems.
           for ( i, subSystem ) in enumerate ( self.subSystems ):
               table.Entry ( "{:d}".format ( i+1 ) )
               table.Entry ( "{:d}".format ( len ( subSystem.system.atoms  ) ) )
               table.Entry ( "{:d}".format ( len ( subSystem.boundaryAtoms ) ) )
               table.Entry ( "{:.3f}".format ( subSystem.weight ) )
               table.Entry ( "{:d}".format ( subSystem.system.electronicState.charge ) )
               table.Entry ( subSystem.system.energyModelLabel )
           # . Finish up the table.
           table.Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

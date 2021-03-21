"""Full QC/MM electrostatic density model for Gaussian basis QC models."""

from   pCore                             import Clone
from   pMolecule.QCModel                 import GaussianBasisQCMMEvaluator
from   pScientific                       import Units
from   pScientific.Geometry3             import Coordinates3
from  .QCMMElectrostaticModelDensityBase import QCMMElectrostaticModelDensityBase
from ..EnergyModel                       import EnergyClosurePriority

# . No pairlist used here as the electron-MM terms are three-body interactions.
# . Instead the pure MM atom selection is used directly to omit BP atom terms.

# . Also eventually best to harmonize different QC/MM electrostatic modules so
# . that same units used in each (either all standard or all atomic units).

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMElectrostaticModelDensityFullGaussianBasis ( QCMMElectrostaticModelDensityBase ):
    """A full QC/MM electrostatic density model for Gaussian basis QC models."""

    _attributable = dict ( QCMMElectrostaticModelDensityBase._attributable )
    _classLabel   = "Gaussian Basis Full Density QC/MM Electrostatic Model"
    _attributable.update ( { "evaluator" : GaussianBasisQCMMEvaluator } )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.QCBPCore ( target )
        def b ( ): self.QCMMCore ( target )
        closures = super ( QCMMElectrostaticModelDensityFullGaussianBasis, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCIntegrals, a, "QC/BP Electrostatic Core" ) ,
                            ( EnergyClosurePriority.QCIntegrals, b, "QC/MM Electrostatic Core" ) ] )
        return closures

    def EnergyFinalize ( self, target ):
        """Energy finalization."""
        scratch = target.scratch
        if scratch.doGradients:
            gScale = Units.Length_Angstroms_To_Bohrs * Units.Energy_Hartrees_To_Kilojoules_Per_Mole
            scratch.qcGradients3QCMM.Add ( scratch.qcGradients3QCMMAU, scale = gScale )
            scratch.gradients3.Add       ( scratch.mmGradients3AU    , scale = gScale )
            state     = getattr ( target, self.__class__._stateName )
            bpCharges = getattr ( state, "bpCharges", None )
            if bpCharges is not None:
                scratch.bpGradients3.Add ( scratch.bpGradients3AU, scale = gScale )
        super ( QCMMElectrostaticModelDensityFullGaussianBasis, self ).EnergyFinalize ( target )

    def EnergyInitialize ( self, target ):
        """Energy initialization."""
        super ( QCMMElectrostaticModelDensityFullGaussianBasis, self ).EnergyInitialize ( target )
        scratch = target.scratch
        state   = getattr ( target, self.__class__._stateName )
        if scratch.Get ( "qcCoordinates3QCMMAU", None ) is None:
            scratch.qcCoordinates3QCMMAU = Clone ( scratch.qcCoordinates3QCMM )
        else:
            scratch.qcCoordinates3QCMM.CopyTo ( scratch.qcCoordinates3QCMMAU )
        scratch.qcCoordinates3QCMMAU.Scale ( Units.Length_Angstroms_To_Bohrs )
        mmCoordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
        if scratch.Get ( "mmCoordinates3AU", None ) is None:
            scratch.mmCoordinates3AU = Clone ( mmCoordinates3 )
        else:
            mmCoordinates3.CopyTo ( scratch.mmCoordinates3AU )
        scratch.mmCoordinates3AU.Scale ( Units.Length_Angstroms_To_Bohrs )
        if state.bpCharges is not None:
            if scratch.Get ( "bpCoordinates3AU", None ) is None:
                scratch.bpCoordinates3AU = Clone ( scratch.bpCoordinates3 )
            else:
                scratch.bpCoordinates3.CopyTo ( scratch.bpCoordinates3AU )
            scratch.bpCoordinates3AU.Scale ( Units.Length_Angstroms_To_Bohrs )
        if scratch.doGradients:
            grd3 = scratch.Get ( "qcGradients3QCMMAU", None )
            if grd3 is None:
                grd3 = Coordinates3.WithExtent ( scratch.qcGradients3QCMM.rows )
                scratch.qcGradients3QCMMAU = grd3
            grd3.Set ( 0.0 )
            grd3 = scratch.Get ( "mmGradients3AU", None )
            if grd3 is None:
                grd3 = Coordinates3.WithExtent ( scratch.gradients3.rows )
                scratch.mmGradients3AU = grd3
            grd3.Set ( 0.0 )
            if state.bpCharges is not None:
                grd3 = scratch.Get ( "bpGradients3AU", None )
                if grd3 is None:
                    grd3 = Coordinates3.WithExtent ( scratch.bpGradients3.rows )
                    scratch.bpGradients3AU = grd3
                grd3.Set ( 0.0 )

    def QCBPCore ( self, target ):
        """Calculate the QC/BP core energy and gradients."""
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if bpCharges is not None:
            scratch = target.scratch
            if scratch.doGradients:
                qcGradients3 = scratch.qcGradients3QCMMAU
                bpGradients3 = scratch.bpGradients3AU
            else:
                qcGradients3 = None
                bpGradients3 = None
            energy = self.evaluator.Core ( target.qcState.nuclearCharges ,
                                           bpCharges                     ,
                                           scratch.qcCoordinates3QCMMAU  ,
                                           scratch.bpCoordinates3AU      ,
                                           None                          ,
                                           qcGradients3                  ,
                                           bpGradients3                  )
            scratch.energyTerms["QC/BP Core"] = energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole

    def QCBPGradients ( self, target ):
        """Calculate the QC/BP electrostatic gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            state     = getattr ( target, self.__class__._stateName )
            bpCharges = getattr ( state, "bpCharges", None )
            if bpCharges is not None:
                self.evaluator.Gradients ( target.qcState.orbitalBases  ,
                                           bpCharges                    ,
                                           scratch.qcCoordinates3QCMMAU ,
                                           scratch.bpCoordinates3AU     ,
                                           None                         ,
                                           scratch.onePDMP.density      ,
                                           scratch.qcGradients3QCMMAU   ,
                                           scratch.bpGradients3AU       )

    def QCBPPotentials ( self, target ):
        """Calculate the QC/BP electrostatic potentials."""
        state     = getattr ( target, self.__class__._stateName )
        bpCharges = getattr ( state, "bpCharges", None )
        if bpCharges is not None:
            scratch = target.scratch
            self.evaluator.Integrals ( target.qcState.orbitalBases  ,
                                       bpCharges                    ,
                                       scratch.qcCoordinates3QCMMAU ,
                                       scratch.bpCoordinates3AU     ,
                                       None                         ,
                                       scratch.qcmmPotentials       )

    def QCMMCore ( self, target ):
        """Calculate the QC/MM core energy and gradients."""
        scratch   = target.scratch
        if scratch.doGradients:
            qcGradients3 = scratch.qcGradients3QCMMAU
            mmGradients3 = scratch.mmGradients3AU
        else:
            qcGradients3 = None
            mmGradients3 = None
        energy = self.evaluator.Core ( target.qcState.nuclearCharges ,
                                       target.mmState.charges        ,
                                       scratch.qcCoordinates3QCMMAU  ,
                                       scratch.mmCoordinates3AU      ,
                                       target.mmState.pureMMAtoms    ,
                                       qcGradients3                  ,
                                       mmGradients3                  )
        scratch.energyTerms["QC/MM Core"] = energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole

    def QCMMGradients ( self, target ):
        """Calculate the QC/MM electrostatic gradients."""
        scratch = target.scratch
        if scratch.doGradients:
            self.evaluator.Gradients ( target.qcState.orbitalBases  ,
                                       target.mmState.charges       ,
                                       scratch.qcCoordinates3QCMMAU ,
                                       scratch.mmCoordinates3AU     ,
                                       target.mmState.pureMMAtoms   ,
                                       scratch.onePDMP.density      ,
                                       scratch.qcGradients3QCMMAU   ,
                                       scratch.mmGradients3AU       )

    def QCMMPotentials ( self, target ):
        """Calculate the QC/MM electrostatic potentials."""
        scratch = target.scratch
        self.evaluator.Integrals ( target.qcState.orbitalBases  ,
                                   target.mmState.charges       ,
                                   scratch.qcCoordinates3QCMMAU ,
                                   scratch.mmCoordinates3AU     ,
                                   target.mmState.pureMMAtoms   ,
                                   scratch.qcmmPotentials       )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

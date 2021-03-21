"""Defines the DFT-D2 QC dispersion model."""

import math, os, os.path

from   pCore              import Clone                 , \
                                 logFile               , \
                                 LogFileActive         , \
                                 YAMLUnpickle
from   pScientific        import Units
from   pScientific.Arrays import Array
from  .QCDispersionDFTD2  import QCDispersionDFTD2_Energy
from  .QCModelError       import QCModelError
from ..EnergyModel        import EnergyClosurePriority , \
                                 EnergyModel

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultYAMLPath = "dftD2.yaml"

#===================================================================================================================================
# . State class.
#===================================================================================================================================
class QCDispersionModelState ( EnergyModelState ):
    """D2 DFT dispersion model state."""

    _attributable = dict ( EnergyModelState._attributable )
    _attributable.update ( { "sqrtC6" : None ,
                             "r0"     : None } )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCDispersionModelDFTD2 ( EnergyModel ):
    """D2 DFT dispersion model."""

    # . Defaults.
    _attributable = dict ( EnergyModel._attributable )
    _classLabel   = "DFT-D2 QC Dispersion Model Summary"
    _stateName    = "qcDispersionState"
    _stateObject      = QCDispersionModelState
    _summarizable = dict ( EnergyModel._summarizable )
    _attributable.update ( { "sqrtC6"   : None ,
                             "dR"       : 0.0  ,
                             "elements" : dict ,
                             "r0"       : None ,
                             "s6"       : 0.0  ,
                             "sR"       : 0.0  } )
    _summarizable.update ( { "dR" : ( "dR", "{:10.5f}" ) ,
                             "sR" : ( "sR", "{:10.5f}" ) ,
                             "s6" : ( "s6", "{:10.5f}" ) }

    def BuildModel ( self, target ):
        """Build the model."""
        qcState = target.qcState
        if qcState is None:
            raise QCModelError ( "This QC dispersion model requires a pre-defined QC model." )
        else:
            state        = super ( QCDispersionModelDFTD2, self ).BuildModel ( target )
            n            = len ( qcState.atomicNumbers )
            state.sqrtC6 = Array.WithExtent ( n )
            state.r0     = Array.WithExtent ( n )
            factor       = 1000.0 * math.pow ( Units.Length_Angstroms_To_Bohrs, 6 ) / Units.Energy_Hartrees_To_Kilojoules_Per_Mole
            for ( i, n ) in enumerate ( qcState.atomicNumbers ):
                if n in self.elements:
                    ( c6, r0 )      = self.elements[n]
                    state.sqrtC6[i] = math.sqrt ( c6 * factor )
                    state.r0    [i] = Units.Length_Angstroms_To_Bohrs * r0
                else:
                    raise QCModelError ( "DFT-D2 dispersion parameters not found for element {:d}.".format ( n ) )

    def Energy ( self, target ):
        """Energy."""
        coordinates3 = target.scratch.qcCoordinates3AU
        gradients3   = target.scratch.Get ( "qcGradients3AU", None )
        state        = getattr ( target, self.__class__._stateName )
        energy       = QCDispersionDFTD2_Energy ( self.s6      ,
                                                  self.sR      ,
                                                  self.dR      ,
                                                  state.sqrtC6 ,
                                                  state.r0     ,
                                                  coordinates3 ,
                                                  gradients3   )
        return { "DFT-D2 Dispersion" : energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole }

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ):
            results = self.Energy ( target )
            target.scratch.energyTerms.update ( results )
        return [ ( EnergyClosurePriority.QCEnergy, a, "QC Dispersion Energy" ) ]

    @classmethod
    def FromYAMLFile ( selfClass, path = None ):
        """Constructor from YAML file."""
        self = selfClass ( )
        self.LoadYAMLFile ( path = path )
        return self

    def LoadYAMLFile ( self, path = None ):
        ""'Load parameters from a YAML file."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "qcDispersion", _DefaultYAMLPath )
        parameters = YAMLUnpickle ( path )
        # . Element parameters.
        self.elements = { n : ( c6, r0 ) for ( n, c6, r0 ) in parameters["Element Values"] }
        # . General parameters.
        generalParameters = parameters["Parameters"]
        self.dR = generalParameters["dR"]
        self.sR = generalParameters["sR"]
        self.s6 = generalParameters["s6"]

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCDispersionModelDFTD2, self ).SummaryItems ( )
        items.append ( ( "Number Of Elements" , "{:d}".format ( len ( self.elements ) ) )
        return items

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

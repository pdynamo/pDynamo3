"""MM parameter assigner."""

import glob, os, os.path

from  pCore                                import Align                    , \
                                                  AttributableObject       , \
                                                  Clone                    , \
                                                  DataType                 , \
                                                  logFile                  , \
                                                  LogFileActive            , \
                                                  YAMLMappingFile_ToObject , \
                                                  YAMLPickleFileExtension
from  pScientific.Arrays                   import Array
from .CMAPDihedralParameterContainer       import CMAPDihedralParameterContainer
from .CosineParameterContainer             import CosineAngleParameterContainer      , \
                                                  CosineDihedralParameterContainer   , \
                                                  CosineOutOfPlaneParameterContainer
from .CosineTermContainer                  import CosineAngleContainer               , \
                                                  CosineDihedralContainer            , \
                                                  CosineOutOfPlaneContainer
from .FourierDihedralContainer             import FourierDihedralContainer
from .FourierDihedralParameterContainer    import FourierDihedralParameterContainer
from .FourierOutOfPlaneParameterContainer  import FourierOutOfPlaneParameterContainer
from .HarmonicAngleContainer               import HarmonicAngleContainer
from .HarmonicAngleParameterContainer      import HarmonicAngleParameterContainer
from .HarmonicBondContainer                import HarmonicBondContainer
from .HarmonicBondParameterContainer       import HarmonicBondParameterContainer
from .HarmonicOutOfPlaneParameterContainer import HarmonicOutOfPlaneParameterContainer
from .LennardJonesParameterContainer       import LennardJonesParameterContainer
from .LJParameterContainer                 import LJParameterContainer
from .MMModelError                         import MMModelError
from .UreyBradleyParameterContainer        import UreyBradleyParameterContainer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
# . Mapping from files names to classes.
_PathClassMapping = { "cmapDihedralParameters"       : CMAPDihedralParameterContainer       ,
                      "cosineAngleParameters"        : CosineAngleParameterContainer        ,
                      "cosineDihedralParameters"     : CosineDihedralParameterContainer     ,
                      "cosineOutOfPlaneParameters"   : CosineOutOfPlaneParameterContainer   ,
                      "fourierDihedralParameters"    : FourierDihedralParameterContainer    ,
                      "fourierOutOfPlaneParameters"  : FourierOutOfPlaneParameterContainer  ,
                      "harmonicAngleParameters"      : HarmonicAngleParameterContainer      ,
                      "harmonicBondParameters"       : HarmonicBondParameterContainer       ,
                      "harmonicOutOfPlaneParameters" : HarmonicOutOfPlaneParameterContainer ,
                      "lennardJones14Parameters"     : LennardJonesParameterContainer       ,
                      "lennardJonesParameters"       : LennardJonesParameterContainer       ,
                      "ureyBradleyParameters"        : UreyBradleyParameterContainer        }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMParameterAssigner ( AttributableObject ):
    """Assign MM parameters."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "lennardJonesParameters"   : None ,
                             "lennardJones14Parameters" : None ,
                             "lennardJonesScale14"      :  1.0 ,
                             "parameterContainers"      : None ,
                             "parameterFactories"       : None ,
                             "path"                     : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( MMParameterAssigner, self )._CheckOptions ( )
        self.ReadData ( )

    def AssignParameters ( self, connectivity, atomTypes, atomCharges, mmState, log ):
        """Assign parameters to the model."""
        # . Initialization.
        missingParameters = set ( )
        uniqueAtomTypes   = self.FindUniqueAtomTypes ( atomTypes )
        # . Atom data.
        self.MakeAtomData ( atomTypes, atomCharges, uniqueAtomTypes, mmState )
        # . LJ parameters.
        for ( localAttribute, modelAttribute ) in ( ( "lennardJonesParameters"  , "ljParameters"   ) ,
                                                    ( "lennardJones14Parameters", "ljParameters14" ) ):
            container = getattr ( self, localAttribute, None )
            if container is not None:
                ( parameters, localMissingParameters ) = container.MakeParameterContainer ( uniqueAtomTypes )
                missingParameters.update ( localMissingParameters )
                if parameters is not None: setattr ( mmState, modelAttribute, parameters )
        # . Exclusions.
        ( mmState.exclusions, mmState.interactions14 ) = connectivity.Make1234And14PairLists ( )
        mmState.exclusions.label     = "Exclusions"
        mmState.interactions14.label = "1-4 Interactions"
        # . MM terms.
        for container in self.parameterContainers:
            ( localMMTerms, localMissingParameters ) = container.MakeMMTermsFromConnectivity ( atomTypes, connectivity )
            missingParameters.update ( localMissingParameters )
            mmState.mmTerms.extend   ( localMMTerms           )
        # . Finish up.
        self.CheckMissingParameters ( missingParameters, log )

    def CheckMissingParameters ( self, missingParameters, log ):
        """Print any missing parameters and raise an error."""
        # . This has been enlarged so as to also print ambiguous and duplicate parameters.
        numberMissing = len ( missingParameters )
        if numberMissing > 0:
            if LogFileActive ( log ):
                # . Sort.
                missingParameters = list ( missingParameters )
                missingParameters.sort ( )
                # . Find label data.
                labelCount  = 0
                labelLength = 0
                tagLength   = 0
                for ( tag, labels ) in missingParameters:
                    labelCount = max ( len ( labels ), labelCount )
                    tagLength  = max ( len ( tag    ), tagLength  )
                    for label in labels:
                        labelLength = max ( len ( label ), labelLength )
                # . Output.
                table = log.GetTable ( columns = [ tagLength + 2 ] + labelCount * [ max ( 10, labelLength + 2 ) ] )
                table.Start  ( )
                table.Title  ( "Missing Force Field Parameters" )
                for ( tag, labels ) in missingParameters:
                    table.Entry ( tag, align = Align.Left )
                    for label in labels: table.Entry ( label )
                    if len ( labels ) < labelCount: table.EndRow ( )
                table.Stop ( )
            raise MMModelError ( "There are {:d} ambiguous, duplicate or missing force field parameters.".format ( numberMissing ) )

    def FindUniqueAtomTypes ( self, atomTypes ):
        """Find a sorted list of unique atom type labels."""
        temporary       = set  ( atomTypes )
        uniqueAtomTypes = list ( temporary )
        uniqueAtomTypes.sort ( )
        return uniqueAtomTypes

    @classmethod
    def FromPath ( selfClass, path, **options ):
        """Constructor from path and other options."""
        options         = dict ( options )
        options["path"] = path
        return selfClass.WithOptions ( **options )

    def MakeAtomData ( self, atomTypes, atomCharges, uniqueAtomTypes, mmState ):
        """Make the atom-specific data."""
        # . Create a type label index.
        labelIndex = {}
        for ( i, label ) in enumerate ( uniqueAtomTypes ):
            labelIndex[label] = i
        # . Create the data.
        n = len ( atomCharges )
        atomTypeIndices = Array.WithExtent ( n, dataType = DataType.Integer )
        charges         = Array.WithExtent ( n )
        for ( i, ( charge, label ) ) in enumerate ( zip ( atomCharges, atomTypes ) ):
            atomTypeIndices[i] = labelIndex[label]
            charges        [i] = charge
        # . Assign the data. 
        mmState.atomTypeIndices = atomTypeIndices
        mmState.atomTypes       = uniqueAtomTypes
        mmState.charges         = charges
        mmState.ljTypeIndices   = Clone ( atomTypeIndices )

    def MakeLennardJones14Parameters ( self ):
        """Make an appropriate container of LJ Parameters for the 1-4 interactions."""
        lj    = self.lennardJonesParameters
        scale = self.lennardJonesScale14
        # . Remove 1-4 parameters if the scaling is zero.
        if scale == 0.0:
            lj14 = None
        # . Processing only if LJs are present too.
        elif lj is not None:
            # . Clone and scale.
            lj14           = Clone ( lj )
            lj14.termLabel = "1-4 Lennard-Jones"
            lj14.ScaleEnergies ( scale )
            # . Update by any predefined 1-4 values.
            oldLJ14 = self.lennardJones14Parameters
            if oldLJ14 is not None:
                lj14.UpdateParameters ( oldLJ14 )
        # . Reset the 14 parameters.
        self.lennardJones14Parameters = lj14

    def ReadData ( self ):
        """Read the parameter data."""
        if self.path is not None:
            self.parameterContainers = []
            paths                    = glob.glob ( os.path.join ( self.path, "*Parameters" + YAMLPickleFileExtension ) )
            for path in paths:
                ( head, tail )             = os.path.split ( path )
                cKey                       = os.path.splitext ( tail )[0]
                container                  = YAMLMappingFile_ToObject ( path, _PathClassMapping[cKey] )
                container.parameterFactory = self.parameterFactories.get ( cKey, None )
                if isinstance ( container, LennardJonesParameterContainer ):
                    if path.find ( "14Parameters" ) >= 0: self.lennardJones14Parameters = container
                    else:                                 self.lennardJonesParameters   = container
                else: self.parameterContainers.append ( container )
            self.MakeLennardJones14Parameters ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

"""Utilities for manipulating a system's Lennard-Jones parameters."""

# . In principle can do the same for all parameter/term containers. E.g. HBPC/HBTC.
# . Raise error if not analytic form. Have hasAnalyticForm property in LJPC.

from  pCore                import Align                , \
                                  logFile              , \
                                  LogFileActive
from  pScientific.Arrays   import ArrayPrint
from .LJParameterContainer import LJForm               , \
                                  LJParameterContainer
from .MMModelError         import MMModelError

#===================================================================================================================================
# . Private functions.
#===================================================================================================================================
def _ChangeEpsilonSigma ( ljParameters, label, epsilon, sigma ):
    """Change an LJ parameter."""
    # . Get data.
    state    = ljParameters.__getstate__ ( )
    labels   = state["parameterKeys"]
    epsilons = state["epsilons"     ]
    sigmas   = state["sigmas"       ]
    # . Change data.
    index           = labels.index ( label )
    epsilons[index] = epsilon
    sigmas  [index] = sigma
    # . Reset data.
    ljParameters               = LJParameterContainer.FromEpsilonSigma ( len ( labels ), epsilons, sigmas, analyticForm = state["analyticForm"] )
    ljParameters.parameterKeys = labels
    return ljParameters

def _GetSystemLJParameters ( system, mode ):
    """Get the LJ parameters for a system corresponding to a given mode."""
    # . System data structures.
    mmModel = getattr ( system, "mmModel", None )
    mmState = getattr ( system, "mmState", None )
    isOK    = ( mmModel is not None ) and ( mmState is not None )
    # . Options.
    do14  = ( mode != "1-5+" ) and isOK
    do15P = ( mode != "1-4"  ) and isOK
    # . Get data.
    data = []
    if do14:
        scale        = getattr ( mmModel, "lennardJonesScale14", 1.0  )
        ljParameters = getattr ( mmState, "ljParameters14"     , None )
        data.append ( ( "lennardJonesScale14", ljParameters, scale ) )
    if do15P:
        ljParameters = getattr ( mmState, "ljParameters", None )
        data.append ( ( "ljParameters", ljParameters, 1.0 ) )
    # . Finish up.
    if len ( data ) > 0: return data
    else:                return None

#===================================================================================================================================
# . Public functions.
#===================================================================================================================================
# . A fudge and needs checking!
def AddNonZeroTIP3PHydrogenLJParameter ( system, label = "HW", mode = None ):
    """Change the current TIP3P hydrogen parameters to their non-zero CHARMM values."""
    epsilon = 0.192464                               # = 0.046 * 4.184.
    if   analyticForm == LJForm.Amber: sigma = 0.449 # =   2.0 * 0.2245
    elif analyticForm == LJForm.OPLS : sigma = 0.400 # = ( 2.0 * 0.2245 ) / ( 2.0 ** ( 1.0/6.0 ) ).
    else: raise MMModelError ( "Invalid analytic form: {:s}.".format ( str ( analyticForm ) ) )
    ChangeLJParameter ( system, label, epsilon, sigma, mode = mode )

def LJUtilities_ChangeParameter ( system, label, epsilon, sigma, mode = None ):
    """Change the LJ parameters for a given LJ type."""
    data = _GetSystemLJParameters ( system, mode )
    if data is not None:
        for ( k, p, s ) in data:
            q = _ChangeEpsilonSigma ( p, label, s * epsilon, sigma )
            setattr ( system.mmState, k, q )

def LJUtilities_ChangeSingleAtomParameter ( system, atomIndex, label, epsilon, sigma, mode = None ):
    """Change an LJ parameter for a single atom."""
    data = _GetSystemLJParameters ( system, mode )
    if data is not None:
        atomTypeIndices = system.mmModel["Atom Type Indices"]
        atomTypes       = system.mmModel["Atom Types"       ]
        ljTypeIndices   = system.mmModel["LJ Type Indices"  ]
        for ( k, p, s ) in data:
            # . Get data.
            ljState    = ljParameters.__getstate__ ( )
            ljLabels   = ljState["parameterKeys"]
            ljEpsilons = ljState["epsilons"     ]
            ljSigmas   = ljState["sigmas"       ]
            # . Change data.
            ljLabels.append   ( label       )
            ljEpsilons.append ( s * epsilon )
            ljSigmas.append   ( sigma       )
            # . Set data.
            ljParameters               = LJParameterContainer.FromEpsilonSigma ( len ( ljLabels ), ljEpsilons, ljSigmas, analyticForm = ljState["analyticForm"] )
            ljParameters.parameterKeys = ljLabels
            system.mmModel[k] = ljParameters
        atomTypes.append ( label )
        atomTypeIndices[atomIndex] = ( len ( atomTypes ) - 1 )
        ljTypeIndices  [atomIndex] = ( len ( ljLabels  ) - 1 ) # . Many assumptions here.

def LJUtilities_PrintParameters ( system, log = logFile, mode = None ):
    """Printing."""
    if LogFileActive ( log ):
        data = _GetSystemLJParameters ( system, mode )
        if data is not None:
            log.Paragraph ( "Lennard-Jones parameters for System \"{:s}\":".format ( system.label ) )
            for ( k, p, s ) in data:
                state    = p.__getstate__ ( )
                form     = state["analyticForm" ]
                label    = p.label
                labels   = state["parameterKeys"]
                epsilons = state["epsilons"     ]
                sigmas   = state["sigmas"       ]
                table    = logFile.GetTable ( columns = [ 10, 20, 20, 20 ] )
                table.Start  ( )
                table.Title  ( "{:s} with Form {:s} and Scaling {:.2f}".format ( label, form.name, s ) )
                table.Heading ( "Index"   )
                table.Heading ( "Label"   )
                table.Heading ( "Epsilon" )
                table.Heading ( "Sigma"   )
                for ( i, ( l, e, s ) ) in enumerate ( zip ( labels, epsilons, sigmas ) ):
                    table.Entry ( "{:d}".format   ( i ) , align = Align.Left )
                    table.Entry (                   l   , align = Align.Left )
                    table.Entry ( "{:.5f}".format ( e ) )
                    table.Entry ( "{:.5f}".format ( s ) )
                table.Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

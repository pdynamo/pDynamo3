"""Functions for the Monte Carlo simulation of molecular systems."""

from pCore                     import logFile               , \
                                      LogFileActive
from pMolecule                 import System
from pMolecule.NBModel         import NBModelMonteCarlo
from pScientific               import Constants             , \
                                      Units
from pScientific.RandomNumbers import RandomNumberGenerator
# . It might be better to switch the use of connectivity isolates to "MM isolates".
# . This could easily be done, along with the use of "free MM isolates" to allow fixed atoms.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The number of random numbers required per Monte Carlo move.
_NRandomIsolate = 7
_NRandomVolume  = 2

# . The simulation options.
_SimulationOptions = { "acceptanceRatio"       :     0.4 ,
                       "adjustFrequency"       :    1000 ,
                       "blocks"                :      10 ,
                       "log"                   : logFile ,
                       "logFrequency"          :       1 ,
                       "moves"                 :   10000 ,
                       "pressure"              :     1.0 ,
                       "randomNumberGenerator" :    None ,
                       "rotation"              :    15.0 ,
                       "temperature"           :   300.0 ,
                       "trajectories"          :    None ,
                       "translation"           :    0.15 ,
                       "volumeChange"          :   400.0 ,
                       "volumeFrequency"       :     500 }

# . The simulation option names.
_SimulationOptionNames = { "Acceptance Ratio"   : "acceptanceRatio" ,
                           "Adjust Frequency"   : "adjustFrequency" ,
                           "Blocks"             : "blocks"          ,
                           "Input RNG"          : "QRNG"            ,
                           "Log Frequency"      : "logFrequency"    ,
                           "Max. Rotation"      : "rotation"        ,
                           "Max. Translation"   : "translation"     ,
                           "Max. Volume Change" : "volumeChange"    ,
                           "Moves"              : "moves"           ,
                           "Pressure"           : "pressure"        ,
                           "Temperature"        : "temperature"     ,
                           "Trajectories"       : "QTRAJECTORIES"   ,
                           "Volume Frequency"   : "volumeFrequency" }

#===================================================================================================================================
# . Private functions.
#===================================================================================================================================
def _CheckSystem ( system ):
    """Check to see if a system is set up for a Monte Carlo calculation."""
    if not ( isinstance ( system        , System            ) and \
             hasattr    ( system        , "nbModel"         ) and \
             isinstance ( system.nbModel, NBModelMonteCarlo ) ):
        raise MonteCarloError ( "System not set up for Monte Carlo calculation." )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MonteCarloError ( Exception ):
    pass

#===================================================================================================================================
# . Calculate the interaction energy between one isolate and the remainder.
#===================================================================================================================================
def MonteCarlo_IsolateInteractionEnergy ( system, isolate, log = None ):
    """Calculate the interaction energy between one isolate and all the others."""
    _CheckSystem ( system )
    system.EnergyInitialize ( False, log )
    system.scratch.energyTerms.update ( system.nbModel.EnergyOfIsolate ( system, isolate ) )
    system.EnergyFinalize ( )
    return system.scratch.energyTerms["Potential Energy"]

#===================================================================================================================================
# . Scale the interaction parameters for an isolate.
#===================================================================================================================================
def MonteCarlo_ScaleIsolateInteractionParameters ( system, isolate, chargeScale = 1.0, epsilonScale = 1.0, sigmaScale = 1.0, log = None ):
    """Scale the interaction parameters for an isolate."""
    _CheckSystem ( system )
    system.nbModel.pairwiseInteraction.ScaleIsolateInteractionParameters ( isolate, chargeScale, epsilonScale, sigmaScale, log = log )

#===================================================================================================================================
# . The simulation function.
#===================================================================================================================================
def MonteCarlo_SystemGeometry ( system, **options ):
    """Metropolis Monte Carlo simulation of a system."""
    cdef IntegerArray1D                ljTypes
    cdef Coordinates3                  coordinates3
    cdef LJParameterContainer          ljParameters
    cdef PairwiseInteractionMonteCarlo pairwiseInteraction
    cdef RealArray1D                   charges
    cdef SelectionContainer            isolates
    cdef SymmetryParameters            symmetryParameters
    cdef CInteger                      blocks, i, iBlock, iMove, moves, tMove
    cdef CMonteCarloSystemGeometry    *cObject    = NULL
    cdef CReal                         dielectric
    cdef CStatus                       status     = CStatus_OK
    # . Get the options.
    unknowns = [ key for key in options if key not in _SimulationOptions ]
    if len ( unknowns ) > 0: raise MonteCarloError ( "Invalid keyword arguments: " + ", ".join ( sorted ( unknowns ) ) + "." )
    for ( key, value ) in _SimulationOptions.items ( ):
        if key not in options: options[key] = value
    # . Get various counters.
    adjustFrequency = options["adjustFrequency"]
    blocks          = options["blocks"]
    moves           = options["moves" ]
    volumeFrequency = options["volumeFrequency"]
    # . Check for printing.
    log             = options["log"]
    logFrequency    = options["logFrequency"]
    doPrinting      = LogFileActive ( log ) and ( logFrequency > 0 ) and ( logFrequency <= blocks )
    if not doPrinting: options["logFrequency"] = 0
    # . Check for a random number generator.
    randomNumberGenerator = options["randomNumberGenerator"]
    options["QRNG"]       = ( randomNumberGenerator is not None )
    if randomNumberGenerator is None: randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )
    # . Check for trajectories.
    logging      = []
    trajectories = options["trajectories"]
    if trajectories is not None:
        for ( trajectory, saveFrequency ) in trajectories:
            if ( trajectory is not None ) and ( saveFrequency > 0 ) and ( saveFrequency <= blocks * moves ):
                logging.append ( ( saveFrequency, trajectory ) )
                trajectory.WriteHeader ( )
    doLogging = ( len ( logging ) > 0 )
    options["QTRAJECTORIES"] = doLogging
    # . Print a summary.
    if doPrinting:
        items = []
        for key in _SimulationOptionNames.keys ( ):
            value = options[_SimulationOptionNames[key]]
            if   isinstance ( value, str ):
                valueString = value
            elif isinstance ( value, bool ):
                valueString = repr ( value )
            elif isinstance ( value, float ):
                valueString = "{:.6g}".format ( value )
            else:
                valueString = "{:s}".format ( repr ( value ) )
            items.append ( ( key, valueString ) )
        log.SummaryOfItems ( items, title = "Monte Carlo Simulation Options" )

    # . Make sure the system is set up correctly.
    # . Basic checks.
    _CheckSystem ( system )
    if system.freeAtoms is not None: raise MonteCarloError ( "Monte Carlo simulations do not currently work with fixed atoms." ) # . This can easily be changed.

    # . Check that there are moves.
    if ( blocks * moves ) > 0:

        # . Set up the Monte Carlo state.
        # . Allocation.
        cObject = MonteCarloSystemGeometry_Allocate ( len ( system.atoms ), max ( _NRandomIsolate, _NRandomVolume ) )
        # . Options.
        cObject.acceptanceRatio = options["acceptanceRatio"]
        cObject.blocks          = blocks
        cObject.moves           = moves
        cObject.rMax            = options["rotation"]
        cObject.tMax            = options["translation"]
        cObject.vMax            = options["volumeChange"]
        # . Aliases.
        charges             = system.mmState.charges             ; cObject.charges             = charges.cObject
        coordinates3        = system.coordinates3                ; cObject.coordinates3        = coordinates3.cObject
        dielectric          = system.nbModel.dielectric          ; cObject.dielectric          = dielectric
        isolates            = system.connectivity.isolateIndices ; cObject.isolates            = isolates.cObject
        ljParameters        = system.mmState.ljParameters        ; cObject.ljParameters        = ljParameters.cObject
        ljTypes             = system.mmState.ljTypeIndices       ; cObject.ljTypes             = ljTypes.cObject   
        pairwiseInteraction = system.nbModel.pairwiseInteraction ; cObject.pairwiseInteraction = pairwiseInteraction.cObject
        symmetryParameters  = system.symmetryParameters          ; cObject.symmetryParameters  = symmetryParameters.cObject
        # . Constants.
        cObject.beta     = 1.0e+03 / ( Constants.Molar_Gas * options["temperature"] )
        cObject.pressure = Units.Pressure_Atmospheres_To_Kilojoules_Per_Mole * options["pressure"]
        cObject.tFactor  = float ( cObject.isolates.capacity ) / cObject.beta
        # . The initial energy and volume.
        cObject.eCurrent = PairwiseInteractionMonteCarlo_MMMMEnergy ( cObject.pairwiseInteraction  ,
                                                                      cObject.charges              ,
                                                                      cObject.ljTypes              ,
                                                                      cObject.ljParameters         ,
                                                                      1.0e+00 / cObject.dielectric ,
                                                                      1.0e+00                      ,
                                                                      cObject.isolates             ,
                                                                      NULL                         ,
                                                                      cObject.coordinates3         ,
                                                                      cObject.symmetryParameters   ,
                                                                      NULL                         ,
                                                                      NULL                         ,
                                                                      &status                      )
        if status != CStatus_OK: raise MonteCarloError ( "Monte Carlo energy error." )
        cObject.volume   = SymmetryParameters_Volume ( cObject.symmetryParameters )
        # . Initialize statistics.
        MonteCarloSystemGeometry_StatisticsStart ( cObject )
        if doPrinting:
            table = log.GetTable ( columns = [ 8, 16, 16, 16, 16, 16, 16, 16 ] )
            table.Start   ( )
            table.Title   ( "Monte Carlo Run Statistics" )
            table.Heading ( "Block"   )
            table.Heading ( "Accept." )
            table.Heading ( "<E>"     )
            table.Heading ( "<dE^2>"  )
            table.Heading ( "<H>"     )
            table.Heading ( "<dH^2>"  )
            table.Heading ( "<V>"     )
            table.Heading ( "<dV^2>"  )

        # . Initialize the total move counter.
        tMove = 0
        # . Loop over the blocks.
        for iBlock from 0 <= iBlock < blocks:
            # . Initialize the block statistics.
            MonteCarloSystemGeometry_StatisticsBlockStart ( cObject )

            # . Loop over the moves per block.
            for iMove from 0 <= iMove < moves:
                # . Increment the total move counter.
                tMove = tMove + 1
                # . Check for a volume move.
                doVolumeMove = ( volumeFrequency > 0 ) and ( ( tMove % volumeFrequency ) == 0 )
                # . Perform a move.
                if doVolumeMove:
                    for i from 0 <= i < _NRandomVolume:
                        cObject.random[i] = randomNumberGenerator.NextReal ( )
                    status = MonteCarloSystemGeometry_MoveVolume  ( cObject )
                else:
                    for i from 0 <= i < _NRandomIsolate:
                        cObject.random[i] = randomNumberGenerator.NextReal ( )
                    status = MonteCarloSystemGeometry_MoveIsolate ( cObject )
                if status != CStatus_OK: raise MonteCarloError ( "Monte Carlo move error." )
                # . Accumulate statistics for the configuration.
                MonteCarloSystemGeometry_StatisticsBlockAccumulate ( cObject )
                # . Adjust move sizes.
                if ( adjustFrequency > 0 ) and ( ( tMove % adjustFrequency ) == 0 ):
                    MonteCarloSystemGeometry_AdjustMoveSizes ( cObject )
                # . Save the configuration.
                if doLogging:
                    for ( saveFrequency, trajectory ) in logging:
                        if ( tMove % saveFrequency ) == 0: trajectory.WriteOwnerData ( )

            # . Do the statistics for the block.
            MonteCarloSystemGeometry_StatisticsBlockStop ( cObject )
            if doPrinting and ( ( iBlock % logFrequency ) == 0 ):
                table.Entry ( "{:d}".format ( iBlock ) )
                table.Entry ( "{:16.4f}".format ( float ( moves - cObject.nReject ) / float ( moves ) ) )
                table.Entry ( "{:16.4f}".format ( cObject.eAv  ) )
                table.Entry ( "{:16.4f}".format ( cObject.eAv2 ) )
                table.Entry ( "{:16.4f}".format ( cObject.hAv  ) )
                table.Entry ( "{:16.4f}".format ( cObject.hAv2 ) )
                table.Entry ( "{:16.4f}".format ( cObject.vAv  ) )
                table.Entry ( "{:16.4f}".format ( cObject.vAv2 ) )

        # . Finish up.
        # . Do the statistics for the run.
        MonteCarloSystemGeometry_StatisticsStop ( cObject )
        if doPrinting:
            if ( blocks > 1 ):
                # . Complete run.
                table.Entry ( "Run" )
                table.Entry ( "{:16.4f}".format ( float ( blocks * moves - cObject.nRejectT ) / float ( blocks * moves ) ) )
                table.Entry ( "{:16.4f}".format ( cObject.eTot  ) )
                table.Entry ( "{:16.4f}".format ( cObject.eTot2 ) )
                table.Entry ( "{:16.4f}".format ( cObject.hTot  ) )
                table.Entry ( "{:16.4f}".format ( cObject.hTot2 ) )
                table.Entry ( "{:16.4f}".format ( cObject.vTot  ) )
                table.Entry ( "{:16.4f}".format ( cObject.vTot2 ) )
                # . Block statistics.
                table.Entry ( "Block" )
                table.Entry ( "-" )
                table.Entry ( "{:16.4f}".format ( cObject.eTotB  ) )
                table.Entry ( "{:16.4f}".format ( cObject.eTotB2 ) )
                table.Entry ( "{:16.4f}".format ( cObject.hTotB  ) )
                table.Entry ( "{:16.4f}".format ( cObject.hTotB2 ) )
                table.Entry ( "{:16.4f}".format ( cObject.vTotB  ) )
                table.Entry ( "{:16.4f}".format ( cObject.vTotB2 ) )
            table.Stop ( )
        # . Write out the final move sizes.
        if doPrinting and ( adjustFrequency > 0 ):
            items = [ ( "Max. Rotation"   , "{:16.6g}".format ( cObject.rMax ) ) ,
                      ( "Max. Translation", "{:16.6g}".format ( cObject.tMax ) ) ,
                      ( "Max. Vol. Move"  , "{:16.6g}".format ( cObject.vMax ) ) ]
            log.SummaryOfItems ( items, title = "Adjusted Move Sizes" )
    # . Deactivate logging.
    if doLogging:
        for ( saveFrequency, trajectory ) in logging:
            trajectory.WriteFooter ( )
            trajectory.Close ( )
        doLogging = False

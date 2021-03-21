"""Titration curve methods for defensin and lysozyme."""

# . It would be good to verify some of the results here ...

import os, os.path

from addOns.pcetk import CEModelMEAD               , \
                         MCModelDefault            , \
                         MEADSubstate              , \
                         StateVector               , \
                         TitrationCurves
from Definitions  import dataPath                  , \
                         outPath
from pBabel       import CHARMMParameterFileReader , \
                         ImportCoordinates3        , \
                         ImportSystem
from pCore        import logFile                   , \
                         TestScriptExit_Fail

# . Header.
logFile.Header ( )

# . Options.
_DefensinSites      = [ ( "PRTA", "ASP" ,  2 ), ( "PRTA", "GLU",  14 ), ( "PRTA", "ARG",  15 ) ]
_LysozymeExclusions = [ ( "PRTA", "CYS",   6 ), ( "PRTA", "CYS", 127 ), ( "PRTA", "CYS",  30 ) ,
                        ( "PRTA", "CYS", 115 ), ( "PRTA", "CYS",  64 ), ( "PRTA", "CYS",  80 ) ,
                        ( "PRTA", "CYS",  76 ), ( "PRTA", "CYS",  94 ), ( "PRTA", "ARG",   0 ) ]
_LysozymeSites      = [ ( "PRTA", "GLU" , 35 ), ( "PRTA", "ASP" , 52 ) ]

# . Paths.
charmmIn = os.path.join ( dataPath, "charmm" )

# . Charmm parameters.
parameters = CHARMMParameterFileReader.PathsToParameters ( [ os.path.join ( charmmIn, "par_all27_prot_na.inp" ), ] )

# . Loop over the systems.
for ( name, includeTermini, exclusions, sites ) in ( ( "defensin" , False, []                 , _DefensinSites ) ,
                                                     ( "lysozyme" , True , _LysozymeExclusions, _LysozymeSites ) ):

    # . Create the system.
    system              = ImportSystem       ( os.path.join ( charmmIn, "{:s}.psfx".format ( name ) ) , isXPLOR = True, parameters = parameters )
    system.coordinates3 = ImportCoordinates3 ( os.path.join ( charmmIn, "{:s}.chm" .format ( name ) ) )
    system.label        = name.capitalize ( )
    system.Summary ( )

    # . Paths.
    resultsPath = os.path.join ( outPath, name )
    if not os.path.exists ( resultsPath ): os.mkdir ( resultsPath )

    # . Create the continuum electrostatic model.
    cem = CEModelMEAD ( system = system, pathScratch = os.path.join ( resultsPath, "mead" ), nthreads = 1 )
    cem.Initialize ( excludeResidues = exclusions, includeTermini = includeTermini )
    cem.Summary       ( )
    cem.SummarySites  ( )
    cem.WriteJobFiles ( )
    cem.CalculateElectrostaticEnergies ( )
    #cem.PrintInteractions ( )

    # . Calculate probabilities and curves analytically.
    logFile.Paragraph ( "*** Calculating protonation states at pH=7 analytically ***" )
    cem.CalculateProbabilities ( pH = 7.0 )
    cem.SummaryProbabilities ( )
    cem.SedScript_FromProbabilities ( filename = os.path.join ( resultsPath, "probabilities_ph7_analytic.sed" ), overwrite = True )
    logFile.Paragraph ( "*** Calculating titration curves analytically ***" )
    curves = TitrationCurves ( cem )
    curves.CalculateCurves   ( )
    curves.WriteCurves  ( directory     = os.path.join ( resultsPath, "curves_analytic" ) )
    curves.PrintHalfpKs ( decimalPlaces = 1 )

    # . Calculate probabilities and curves using Monte Carlo.
    logFile.Paragraph ( "*** Calculating protonation states at pH=7 using Monte Carlo sampling ***" )
    sampling = MCModelDefault  ( )
    cem.DefineMCModel          ( sampling )
    cem.CalculateProbabilities ( pH = 7.0 )
    cem.SummaryProbabilities   ( )
    cem.SedScript_FromProbabilities ( filename = os.path.join ( resultsPath, "probabilities_ph7_MonteCarlo.sed" ), overwrite = True )
    logFile.Paragraph ( "*** Calculating titration curves using Monte Carlo sampling ***" )
    mcc = TitrationCurves      ( cem )
    mcc.CalculateCurves        ( )
    mcc.WriteCurves  ( directory     = os.path.join ( resultsPath, "curves_mc" ) )
    mcc.PrintHalfpKs ( decimalPlaces = 1 )

    # . Calculate the substate energies.
    logFile.Paragraph ( "*** Calculating substate energies at pH=7 using Monte Carlo sampling ***" )
    substate = MEADSubstate ( cem, sites )
    substate.CalculateSubstateEnergies ( )
    substate.Summary ( )
    substate.Summary_ToLatex ( filename = os.path.join ( resultsPath, "substateEnergies_ph7_MonteCarlo.text" ), includeSegment = True )

    # . Calculate the energies of the first few states.
    logFile.Paragraph ( "*** Calculating energies of the first few state vectors (maximum 10) ***" )
    table = logFile.GetTable ( columns = [ 6, 16 ] )
    table.Start ( )
    table.Heading ( "State"  )
    table.Heading ( "Energy" )
    vector = StateVector ( cem )
    vector.Reset ()
    for i in range  ( 10 ):
        e    = cem.CalculateMicrostateEnergy ( vector, pH = 7.0 )
        isOK = vector.Increment ( )
        if not isOK: break
        table.Entry ( "{:d}"  .format ( i + 1 ) )
        table.Entry ( "{:.1f}".format ( e     ) )
    table.Stop ( )

# . Footer.
isOK = True
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )

"""Helper classes and functions for calculating free energies using the WHAM technique."""

import math

from pBabel                                 import SystemRestraintTrajectoryDataHandler
from pCore                                  import Clone                                , \
                                                   logFile                              , \
                                                   LogFileActive
from pScientific                            import Constants                            , \
                                                   Units
from pScientific.Arrays                     import Array
from pScientific.ObjectiveFunctionIterators import ConjugateGradientMinimizer           , \
                                                   LBFGSMinimizer                       , \
                                                   ObjectiveFunction 
from pScientific.Statistics                 import Statistics

# . Minimization methods from F. Zhu and G. Hummer, JCC 33, 453-465 (2012).
#   Other methods can be found there but are not yet implemented.
#   To be done includes:
#   - Error estimation.
#   - Independence of input data points.
#   - Mean force calculations (integration or llsq fit to get free energies).

# . Current bootstrapping method put in by Ramon Crehuet following
# . R. W. Johnson, Teaching Statistics 23, 49-54 (2001).

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Default values.
_Default_FreeEnergyTolerance  = 1.0e-6
_Default_LogFrequency         = 10
_Default_MaximumIterations    = 10000
_Default_Pressure             = 1.0
_Default_ProbabilityTolerance = 1.0e-10
_Default_Resamples            = 100
_Default_Temperature          = 300.0

#===================================================================================================================================
# . WHAM objective function.
#===================================================================================================================================
class WHAMObjectiveFunction ( ObjectiveFunction ):
    """The WHAM objective function."""

    # . The values in |energies| are exp ( F_j / k*T0 ) where F_j is the free energy for window j (note: + not - F_j).


    _attributable = dict ( ObjectiveFunction._attributable )
    _clearable    = { "energies"             : None ,
                      "histogram"            : None ,
                      "pmf"                  : None ,
                      "pointsPerBin"         : None ,
                      "pointsPerWindow"      : None ,
                      "pressure"             : None ,
                      "probabilities"        : None ,
                      "reducedHistogram"     : None ,
                      "reducedProbabilities" : None ,
                      "temperature"          : None ,
                      "temperatureFactor"    : None ,
                      "weightMatrix"         : None ,
                      "workE"                : None ,
                      "workP"                : None }
    _attributable.update ( { "handler"        : None ,
                             "variableOption" : 1    } )
    _attributable.update ( _clearable )

    def Allocate ( self ):
        """Allocate arrays."""
        if ( self.handler is not None ) and ( self.histogram is not None ):
            self.energies      = Array.WithExtent ( self.handler.numberOfWindows ) ; self.energies.Set      ( 1.0 )
            self.probabilities = Array.WithExtent ( self.histogram.bins          ) ; self.probabilities.Set ( 1.0 )
            self.workE         = Array.WithExtent ( self.handler.numberOfWindows ) ; self.workE.Set         ( 0.0 )
            self.workP         = Array.WithExtent ( self.histogram.bins          ) ; self.workP.Set         ( 0.0 )
            self.histogram.Normalize ( self.probabilities )

    def ClearAll ( self ):
        """Clear the function."""
        for ( key, value ) in self.__class__._clearable.items ( ): setattr ( self, key, value )

    def ConvertEnergies ( self ):
        """Convert self.energies to energies."""
        if self.energies is not None:
            self.energies.NaturalLogarithm ( )
            self.energies.Scale ( self.temperatureFactor )

    def Finalize ( self ):
        """Finalization after optimization."""
        self.ConvertEnergies ( )
        self.ReduceHistogram ( )
        self.GetPMF          ( )

    @classmethod
    def FromHandler ( selfClass, handler ):
        """Constructor given a handler."""
        self                   = selfClass ( )
        self.handler           = handler
        return self

    def Function ( self, variables ):
        """Evaluate the function."""
        self.VariablesPut ( variables )
        # . Calculate p.
        self.ProbabilityDistribution ( )
        # . Calculate the negative log-likelihood.
        if self.variableOption != 0: variables.CopyTo     ( self.workE )
        else:                        self.energies.CopyTo ( self.workE ) ; self.workE.NaturalLogarithm ( )
        self.probabilities.CopyTo ( self.workP ) ; self.workP.NaturalLogarithm ( )
        return ( - self.workE.Dot ( self.pointsPerWindow ) - self.workP.Dot ( self.pointsPerBin ) )

    def FunctionGradients ( self, variables, gradients ):
        """Evaluate the function and gradients."""
        likelihood = self.Function ( variables )
        if self.variableOption == 0:
            gradients.Set ( -1.0 )
            gradients.Divide ( self.energies )
            self.weightMatrix.VectorMultiply ( self.probabilities, gradients, beta = 1.0 )
        else:
            self.weightMatrix.VectorMultiply ( self.probabilities, gradients )
            gradients.Multiply ( self.energies )
            gradients.Increment ( -1.0 )
        gradients.Multiply ( self.pointsPerWindow )
        return likelihood

    def GetPMF ( self ):
        """Calculate the PMF."""
        if ( self.reducedHistogram is not None ) and ( self.reducedProbabilities is not None ):
            pmf = Array.WithExtent ( self.reducedHistogram.bins ) ; pmf.Set ( 0.0 )
            self.reducedProbabilities.CopyTo ( pmf )
            pmf.NaturalLogarithm ( )
            pmf.Scale ( -self.temperatureFactor )
            v = min ( pmf )
            pmf.Increment ( -v )
            self.pmf = pmf

    def HistogramData ( self, bins ):
        """Histogram the data."""
        self.ClearAll ( )
        if self.handler is not None:
            n = self.handler.rank
            if   isinstance ( bins, int ): binsArray = n * [ bins ]
            elif isinstance ( bins, ( list, tuple ) ) and ( len ( bins ) == n ): binsArray = bins
            else: raise ValueError ( "Invalid input bins specification." )
            self.histogram = self.handler.HistogramData ( binsArray )

    def Initialize ( self, pressure = _Default_Pressure, temperature = _Default_Temperature ):
        """Create the arrays needed for the WHAM equations."""
        # . Initialization.
        handler   = self.handler
        histogram = self.histogram
        if ( handler is not None ) and ( histogram is not None ):
            p0        = pressure * Units.Pressure_Atmospheres_To_Kilojoules_Per_Mole
            t0        = temperature * Constants.Molar_Gas * 1.0e-3
            # . Gather pressures and temperatures.
            pressures    = []
            temperatures = []
            for ( p, t ) in zip ( handler.pressures, handler.temperatures ):
                if p is None: p = pressure
                if t is None: t = temperature
                pressures.append    ( p * Units.Pressure_Atmospheres_To_Kilojoules_Per_Mole )
                temperatures.append ( t * Constants.Molar_Gas * 1.0e-3                       )
            # . Create c which in the most general case is Exp ( - ( (Bj - B0) * U_i + (Bj*Pj - B0*P0) * V_i + Bj * SC_ji ) ).
            c = Array.WithExtents ( handler.numberOfWindows, histogram.bins ) ; c.Set ( 0.0 )
            for ( i, midPoint ) in enumerate ( histogram.BinMidPointIterator ( ) ):
                # . Energy.
                if handler.hasEnergy:
                    e = midPoint.pop ( 0 )
                    for ( j, tj ) in enumerate ( temperatures ):
                        c[j,i] += e * ( 1.0 / tj - 1.0 / t0 )
                # . Volume.
                if handler.hasVolume:
                    v = midPoint.pop ( 0 )
                    for ( j, pj ) in enumerate ( pressures ):
                        c[j,i] += v * ( pj / tj - p0 / t0 )
                # . Restraints.
                for ( j, tj ) in enumerate ( temperatures ):
                    for ( s, model ) in zip ( midPoint, handler.energyModels[j] ):
                        c[j,i] += model.Energy ( s )[0] / tj
            c.Scale ( -1.0 )
            c.Exponential ( )
            # . Create the points arrays.
            m = Array.WithExtent ( histogram.bins          ) ; m.Set ( 0.0 )
            n = Array.WithExtent ( handler.numberOfWindows ) ; n.Set ( 0.0 )
            for ( i, v ) in enumerate ( histogram.counts     ): m[i] = float ( v )
            for ( i, v ) in enumerate ( handler.windowPoints ): n[i] = float ( v )
            # . Save all.
            self.pointsPerBin      = m
            self.pointsPerWindow   = n
            self.pressure          = pressure
            self.temperature       = temperature
            self.temperatureFactor = t0
            self.weightMatrix      = c

    def ProbabilityDistribution ( self ):
        """Calculate the probability distribution."""
        self.energies.CopyTo ( self.workE )
        self.workE.Multiply ( self.pointsPerWindow )
        self.weightMatrix.VectorMultiply ( self.workE, self.probabilities, transpose = True )
        self.probabilities.Reciprocate ( )
        self.probabilities.Multiply ( self.pointsPerBin )
#        self.histogram.Normalize ( self.probabilities )

    def ReduceHistogram ( self ):
        """Calculate the reduced histogram and probabilities."""
        if ( self.handler is not None ) and ( self.histogram is not None ) and ( self.probabilities is not None ):
            self.histogram.Normalize ( self.probabilities ) # . Here?
            d        = 0
            toReduce = []
            if self.handler.hasEnergy:
                toReduce.append ( d )
                d += 1
            if self.handler.hasVolume:
                toReduce.append ( d )
            if len ( toReduce ) == 0:
                self.reducedHistogram     = self.histogram
                self.reducedProbabilities = self.probabilities
            else:
                self.reducedHistogram     = self.histogram.MakeReducedHistogram ( toReduce )
                self.reducedProbabilities = self.histogram.ReduceData           ( toReduce, self.probabilities )
                self.reducedHistogram.Normalize ( self.reducedProbabilities )

    def RehistogramData ( self, histogram ):
        """Histogram the data using an existing histogram."""
        self.ClearAll ( )
        if self.handler is not None:
            self.histogram = histogram
            self.handler.RehistogramData ( histogram )

    def SummaryFull ( self, log = logFile ):
        """A full summary."""
        if LogFileActive ( log ):
            # . General information.
            text = "PMF generated at a temperature of {:.1f} K".format ( self.temperature )
            if self.handler.hasVolume: text += " and a pressure of {:.1f} atm.".format ( self.pressure )
            else:                      text += "."
            log.Paragraph ( text )
            # . Window information.
            table = log.GetTable ( columns = [ 10, 20, 20 ] )
            table.Start ( )
            table.Title ( "Window Information" )
            table.Heading ( "Window"      )
            table.Heading ( "Points"      )
            table.Heading ( "Free Energy" )
            for ( i, ( p, v ) ) in enumerate ( zip ( self.handler.windowPoints, self.energies ) ):
                table.Entry ( "{:d}"  .format ( i ) )
                table.Entry ( "{:d}"  .format ( p ) )
                table.Entry ( "{:.6g}".format ( v ) )
            table.Stop ( )
            # . PMF.
            table = log.GetTable ( columns = ( len ( self.reducedHistogram.dimensions ) + 2 ) * [ 20 ] )
            table.Start ( )
            table.Title ( "Potential of Mean Force" )
            for dimension in self.reducedHistogram.dimensions: table.Heading ( dimension.label )
            table.Heading ( "PMF" )
            table.Heading ( "PDF" )
            for ( i, midPoint ) in enumerate ( self.reducedHistogram.BinMidPointIterator ( ) ):
                for v in midPoint: table.Entry ( "{:.6g}".format ( v ) )
                table.Entry ( "{:.6g}".format ( self.pmf[i] ) )
                table.Entry ( "{:.6g}".format ( self.reducedProbabilities[i] ) )
            table.Stop ( )

    def VariablesAllocate ( self ):
        """Return an object to hold the variables."""
        variables = Array.WithExtent ( self.numberOfVariables )
        self.VariablesGet ( variables )
        return variables

    def VariablesGet ( self, variables ):
        """Fill the variable array."""
        if self.energies is not variables: self.energies.CopyTo ( variables )
        if self.variableOption != 0: variables.NaturalLogarithm ( )

    def VariablesPut ( self, variables ):
        """Empty the variable array."""
        if self.energies is not variables: variables.CopyTo ( self.energies )
        if self.variableOption != 0: self.energies.Exponential ( )

    @property
    def numberOfVariables ( self ):
        """Return the number of variables."""
        return self.handler.numberOfWindows

#===================================================================================================================================
# . Helper functions for solving the WHAM equations by direct iteration or by minimization.
#===================================================================================================================================
def _BootstrappingSummary ( histogram, log, resamples, pmf, standardError, lowerLimit, upperLimit ):
    """A bootstrapping summary."""
    if LogFileActive ( log ):
        # . General information.
        log.Paragraph ( "Bootstrapping performed with {:d} resamples.".format ( resamples ) )
        # . PMF.
        table = log.GetTable ( columns = ( len ( histogram.dimensions ) + 4 ) * [ 20 ] )
        table.Start ( )
        table.Title ( "Potential of Mean Force" )
        for dimension in histogram.dimensions: table.Heading ( dimension.label )
        table.Heading ( "PMF"          )
        table.Heading ( "Sample Error" )
        table.Heading ( "95% Confidence Interval", columnSpan = 2 )
        for ( i, midPoint ) in enumerate ( histogram.BinMidPointIterator ( ) ):
            for v in midPoint: table.Entry ( "{:.6g}".format ( v ) )
            table.Entry ( "{:.6g}".format ( pmf          [i] ) )
            table.Entry ( "{:.6g}".format ( standardError[i] ) )
            table.Entry ( "{:.6g}".format ( lowerLimit   [i] ) )
            table.Entry ( "{:.6g}".format ( upperLimit   [i] ) )
        table.Stop ( )

_Alpha = 0.05 # . Constant for 95% confidence limit.
def WHAM_Bootstrapping ( paths, **options ):
    """Solve the WHAM equations and estimate errors using bootstrapping."""
    # . Get specific options.
    log       = options.get ( "log", logFile )
    minimizer = options.pop ( "minimizerFunction", WHAM_ConjugateGradientMinimize )
    resamples = options.pop ( "resamples"        , _Default_Resamples             )
    # . Initial minimization - it must succeed to proceed.
    results = minimizer ( paths, **dict ( options ) )
    if results["Converged"]:
        # . Reset options for resampling.
        options["handler"  ] =         results.pop ( "Handler"             )
        options["histogram"] = Clone ( results.pop ( "Unreduced Histogram" ) )
        options["log"      ] = None
        # . Resampling minimizations.
        pmfs = []
        for r in range ( resamples ):
            localResults = minimizer ( paths, **dict ( options ) )
            if localResults["Converged"]: pmfs.append ( localResults.pop ( "PMF" ) )
        # . Process results.
        numberOfResamples = len ( pmfs )
        results.update ( { "Number Of Resamples" : numberOfResamples ,
                           "PMFs"                : pmfs              } )
        if numberOfResamples > 1:
            # . Gather some counters.
            lowerIndex        = int ( math.ceil  (         0.5 * _Alpha   * float ( numberOfResamples ) ) ) - 1
            upperIndex        = int ( math.floor ( ( 1.0 - 0.5 * _Alpha ) * float ( numberOfResamples ) ) ) - 1
            # . Do statistics.
            originalPMF   = results["PMF"]
            lowerLimit    = Array.WithExtent ( len ( originalPMF ) )
            standardError = Array.WithExtent ( len ( originalPMF ) )
            upperLimit    = Array.WithExtent ( len ( originalPMF ) )
            work          = Array.WithExtent ( numberOfResamples   )
            for b in range ( len ( originalPMF ) ):
                for ( r, pmf ) in enumerate ( pmfs ): work[r] = pmf[b]
                work.Sort ( )
                m                = originalPMF[b]
                lowerLimit[b]    = 2.0 * m - work[upperIndex]
                upperLimit[b]    = 2.0 * m - work[lowerIndex]
                standardError[b] = Statistics ( work ).standardDeviation
            # . Printing.
            data = ( ( "Standard Error"                , standardError ) ,
                     ( "Lower 95% Confidence Interval" , lowerLimit    ) ,
                     ( "Upper 95% Confidence Interval" , upperLimit    ) )
            _BootstrappingSummary ( results["Histogram"], log, numberOfResamples, originalPMF, standardError, lowerLimit, upperLimit )
            # . Finish up.
            for ( key, value ) in data: results[key] = value
    return results

#-----------------------------------------------------------------------------------------------------------------------------------
def WHAM_ConjugateGradientMinimize ( paths, **options ):
    """Solve the WHAM equations by conjugate gradient minimization."""
    ( of, options, log ) = _InitializeSolution ( paths, ConjugateGradientMinimizer._attributable, {}, options )
    optimizer = ConjugateGradientMinimizer.WithOptions ( **options )
    optimizer.Summary ( log = log )
    report = optimizer.Iterate ( of, log = log )
    return _FinalizeSolution ( of, log, report )

#-----------------------------------------------------------------------------------------------------------------------------------
def WHAM_DirectIteration ( paths, **options ):
    """Solve the WHAM equations by direct iteration."""
    defaultOptions       = { "freeEnergyTolerance"  : _Default_FreeEnergyTolerance  ,
                             "logFrequency"         : _Default_LogFrequency         ,
                             "maximumIterations"    : _Default_MaximumIterations    ,
                             "probabilityTolerance" : _Default_ProbabilityTolerance }
    ( of, options, log ) = _InitializeSolution ( paths, defaultOptions, {}, options )
    report = _WHAMDirectIterator ( of.histogram       ,
                                   of.weightMatrix    ,
                                   of.pointsPerWindow ,
                                   of.pointsPerBin    ,
                                   of.energies        ,
                                   of.probabilities   ,
                                   of.workE           ,
                                   of.workP           ,
                                   **options          )
    return _FinalizeSolution ( of, log, report )

#-----------------------------------------------------------------------------------------------------------------------------------
def WHAM_LBFGSMinimize ( paths, **options ):
    """Solve the WHAM equations by LBFGS minimization."""
    ( of, options, log ) = _InitializeSolution ( paths, LBFGSMinimizer._attributable, {}, options )
    optimizer = LBFGSMinimizer ( **options )
    optimizer.Summary ( log = log )
    report = optimizer.Iterate ( of, log = log )
    return _FinalizeSolution ( of, log, report )

#-----------------------------------------------------------------------------------------------------------------------------------
def WHAM_TestGradients ( paths, **options ):
    """Test the derivatives of the WHAM equations."""
    ( of, options, log ) = _InitializeSolution ( paths, {}, {}, options )
    delta     = options.get ( "delta"     , 1.0e-4 )
    tolerance = options.get ( "tolerance" , 1.0e-4 )
    of.TestGradients ( delta = delta, log = log, tolerance = tolerance )

#===================================================================================================================================
# . Helper function for cleaning up after solution.
#===================================================================================================================================
def _FinalizeSolution ( of, log, report ):
    """Finalize after solution."""
    of.Finalize    ( )
    of.SummaryFull ( log = log )
    report.update ( { "Free Energies"       : of.energies             ,
                      "Handler"             : of.handler              ,
                      "Histogram"           : of.reducedHistogram     ,
                      "PMF"                 : of.pmf                  ,
                      "Probabilities"       : of.reducedProbabilities ,
                      "Unreduced Histogram" : of.histogram            } )
    return report

#===================================================================================================================================
# . Helper function for setting up solvers.
#===================================================================================================================================
def _InitializeSolution ( paths, defaultOptions, defaultsToChange, inputOptions ):
    """Initialize before solution."""
    # . Get the options.
    options = dict ( defaultOptions   )
    options.update ( defaultsToChange )
    options.update ( inputOptions     )
    # . Get some non-optimizer options.
    bins        = options.pop ( "bins"        , None                 )
    handler     = options.pop ( "handler"     , None                 ) # . Resample an existing handler.
    histogram   = options.pop ( "histogram"   , None                 ) # . Force the use of an existing histogram.
    labels      = options.pop ( "labels"      , None                 )
    log         = options.pop ( "log"         , logFile              )
    pressure    = options.pop ( "pressure"    , _Default_Pressure    )
    temperature = options.pop ( "temperature" , _Default_Temperature )
    # . Make sure label is not None.
    if options.get ( "label", None ) is None: options.pop ( "label", None )
    # . Gather all data.
    if handler is None:
        handler = SystemRestraintTrajectoryDataHandler.FromTrajectoryPaths ( paths )
        if labels is not None: handler.Prune ( labels )
    else:
        handler.ResampleData ( )
    handler.TranslateEnergyVolumeData ( )
    # . Create the objective function.
    of = WHAMObjectiveFunction.FromHandler ( handler )
    # . Histogram the data.
    if histogram is None: of.HistogramData   ( bins      )
    else:                 of.RehistogramData ( histogram )
    # . Create the arrays needed for the WHAM equations.
    of.Allocate   ( )
    of.Initialize ( pressure = pressure, temperature = temperature )
    # . Finish up.
    return ( of, options, log )

#===================================================================================================================================
# . WHAM equation solver using direct iteration.
#===================================================================================================================================
def _WHAMDirectIterator ( histogram, c, m, n, f, p, fOld, pOld , **options ):
    """Solve the WHAM equations via direct iteration."""
    # . Get options.
    freeEnergyTolerance  = options.pop ( "freeEnergyTolerance"  , _Default_FreeEnergyTolerance  )
    log                  = options.pop ( "log"                  , logFile                       )
    logFrequency         = options.pop ( "logFrequency"         , _Default_LogFrequency         )
    maximumIterations    = options.pop ( "maximumIterations"    , _Default_MaximumIterations    )
    probabilityTolerance = options.pop ( "probabilityTolerance" , _Default_ProbabilityTolerance )
    if len ( options ) > 0: raise ValueError ( "Unrecognized keyword arguments: {:s}.".format ( ", ".join ( options.keys ( ) ) ) )
    # . Check for printing.
    doPrint = ( logFrequency > 0 ) and ( logFrequency < maximumIterations ) and LogFileActive ( log )
    if doPrint:
        table = log.GetTable ( columns = [ 10, 20, 20, 20 ] )
        table.Start ( )
        table.Title ( "WHAM Equation Solver" )
        table.Heading ( "Iteration"   )
        table.Heading ( "Likelihood"  )
        table.Heading ( "Change in F" )
        table.Heading ( "Change in P" )
    # . Initialization.
    isConverged = False
    iterations  = 0
    # . Perform the iterations.
    while ( iterations < maximumIterations ) and ( not isConverged ):
        # . Save old data.
        f.CopyTo ( fOld )
        p.CopyTo ( pOld )
        # . Calculate the distribution function.
        # . p_i = n_i / sum_j f_j m_j c_ji
        f.Multiply ( m )
        c.VectorMultiply ( f, p, transpose = True )
        p.Reciprocate ( )
        p.Multiply ( n )
        histogram.Normalize ( p )
        # . Calculate the window free energies.
        # . 1/f_j = sum_i c_ji p_i
        c.VectorMultiply ( p, f )
        f.Reciprocate ( )
        # . Check for convergence.
        fOld.Add ( f, scale = -1.0 )
        pOld.Add ( p, scale = -1.0 )
        fDifference = fOld.AbsoluteMaximum ( )
        pDifference = pOld.AbsoluteMaximum ( )
        isConverged = ( fDifference < freeEnergyTolerance ) and ( pDifference < probabilityTolerance )
        # . Printing.
        if doPrint and ( iterations % logFrequency == 0 ):
            # . Likelihood.
            f.CopyTo ( fOld ) ; fOld.NaturalLogarithm ( )
            p.CopyTo ( pOld ) ; pOld.NaturalLogarithm ( )
            likelihood = - fOld.Dot ( m ) - pOld.Dot ( n )
            # . Output.
            table.Entry ( "{:d}"  .format ( iterations  ) )
            table.Entry ( "{:.6g}".format ( likelihood  ) )
            table.Entry ( "{:.6g}".format ( fDifference ) )
            table.Entry ( "{:.6g}".format ( pDifference ) )
        # . End of loop.
        iterations += 1
    # . Stop printing.
    if doPrint:
        table.Stop ( )
        if isConverged: log.Paragraph ( "WHAM procedure converged after {:d} iterations.".format              ( iterations ) )
        else:           log.Paragraph ( "Warning: WHAM procedure not converged after {:d} iterations.".format ( iterations ) )
    # . Finish up.
    return { "Converged" : isConverged, "Iterations" : iterations }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

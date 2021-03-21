"""Saddle point refinement used in Conjugate Peak Refinment"""

import math

from  pCore                                  import logFile
from  pScientific.Arrays                     import Array
from  pScientific.ObjectiveFunctionIterators import MultiDimensionalMinimizer      , \
                                                    MultiDimensionalMinimizerState , \
                                                    UniDimensionalObjectiveFunction
from .MoreThuenteLineSearchWithMax           import MoreThuenteLineSearcherWithMax

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default parameters for the line searcher.
_LineSearcherParameters = { "bisectionFactor"              :   0.66  ,
                            "functionTolerance"            : 1.0e-04 ,
                            "gradientTolerance"            :   0.1   ,  # . For better maxima.
                            "initialStep"                  :   0.1   ,
                            "logFrequency"                 :   -1    ,
                            "lowerStepExtrapolationFactor" :   1.1   ,
                            "maximumIterations"            :   20    ,
                            "safeguardFactor"              :  0.66   ,
                            "upperStepExtrapolationFactor" :  4.0    ,
                            "maximize"                     :  False  ,
                            "variableTolerance"            : 1.0e-15 ,
                            "rmsdMaximumTolerance"         :  1.5    }

# . A large number for the upper bound of the line search.
_LargeNumber = 1.0e+20

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LineSearchObjectiveFunction ( UniDimensionalObjectiveFunction ):
    """The line search objective function."""

    _attributable = dict ( UniDimensionalObjectiveFunction._attributable )
    _attributable.update ( { "d"                 : None ,
                             "f"                 : None ,
                             "f0"                : None ,
                             "fB"                : None ,
                             "g"                 : None ,
                             "g0"                : None ,
                             "gB"                : None ,
                             "objectiveFunction" : None ,
                             "vB"                : None ,
                             "x"                 : None ,
                             "x0"                : None ,
                             "xB"                : None } )

    @classmethod
    def FromOptions (selfClass, **options):
        """Constructor from options."""
        self = selfClass ()
        for ( key, value ) in options.items ( ): setattr ( self, key, value )
        return self

    def FunctionGradient (self, alpha):
        """Evaluate the function and gradient."""
        # . Save the explicit variable.
        self.variable = alpha
        # . Set up the implicit variables.
        self.x0.CopyTo (self.x)
        self.x.Add ( self.d, scale = alpha )
        # . Function and gradients.
        self.f = self.objectiveFunction.FunctionGradients (self.x, self.g)
        return ( self.f, self.d.Dot ( self.g ) )

    def GetBounds (self): return (-0.5 , _LargeNumber)

    def GetValuesAtOrigin (self):
        """Get the values at the origin."""
        return (0.0, self.f0, self.g0)

    def SetValuesAtOrigin (self, f0, g0):
        """Set the values at the origin."""
        self.f0 = f0
        self.g0 = g0

    def RestoreBestToCurrent (self):
        """Make the best point the current one."""
        self.f = self.fB
        self.variable = self.vB
        self.gB.CopyTo (self.g)
        self.xB.CopyTo (self.x)

    def StoreCurrentAsBest (self):
        """Save the current point as the best point."""
        self.fB = self.f
        self.vB = self.variable
        self.g.CopyTo (self.gB)
        self.x.CopyTo (self.xB)

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CPRSaddlePointRefinementState ( MultiDimensionalMinimizerState ):
    """CPR saddle point refinement state"""

    _attributable = dict ( MultiDimensionalMinimizerState._attributable )
    _attributable.update ( { "g0"                          : None ,
                             "gB"                          : None ,
                             "gDotG"                       : 0.0  ,
                             "lineSearchObjectiveFunction" : None ,
                             "s"                           : None ,
                             "theta"                       : 1.0  ,
                             "f00"                         : 0.0  ,
                             "x0"                          : None ,
                             "y"                           : None ,
                             "h"                           : None ,
                             "saddle"                      : False,
                             "stationary"                  : False,
                             "successiveLowGradients"      :  0   ,
                             "s0"                          : None } )

    def SetUp ( self ):
       """Set up the state."""
       # . Arrays.
       super ( CPRSaddlePointRefinementState, self ).SetUp ( )
       self.g0  = Array.WithExtent ( self.numberOfVariables ) ; self.g0.Set  ( 0.0 )
       self.gB  = Array.WithExtent ( self.numberOfVariables ) ; self.gB.Set  ( 0.0 )
       self.s   = Array.WithExtent ( self.numberOfVariables ) ; self.s.Set   ( 0.0 )
       self.x0  = Array.WithExtent ( self.numberOfVariables ) ; self.x0.Set  ( 0.0 )
       self.y   = Array.WithExtent ( self.numberOfVariables ) ; self.y.Set   ( 0.0 )
       self.s0  = Array.WithExtent ( self.numberOfVariables ) ; self.s0.Set  ( 0.0 )
       self.h   = Array.WithExtent ( self.numberOfVariables ) ; self.h.Set   ( 0.0 )
       self.tmp = Array.WithExtent ( self.numberOfVariables ) ; self.tmp.Set ( 0.0 )
       # . Line-search objective function.
       # . s and y are used here as they are only used for workspace and not for persistent storage.
       self.lineSearchObjectiveFunction = LineSearchObjectiveFunction.FromOptions ( d  = self.d                                ,
                                                                                    g  = self.g                                ,
                                                                                    gB = self.gB                               ,
                                                                                    objectiveFunction = self.objectiveFunction ,
                                                                                    x  = self.x                                ,
                                                                                    x0 = self.x0                               ,
                                                                                    xB = self.y                                )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CPRSaddlePointRefinement ( MultiDimensionalMinimizer ):
    """CPR saddle point refinement class"""

    _attributable = dict ( MultiDimensionalMinimizer._attributable )
    _classLabel   = "CPR Saddle Point Refinement"
    _stateObject  = CPRSaddlePointRefinementState
    _summarizable = dict ( MultiDimensionalMinimizer._summarizable )
    _attributable.update ( { "finiteDifferenceStep"           :-1.0e-3         , # . This step should be > than 1st derivative or numerical accuracy of QC method.
                             "finiteDifference"               : True           ,
                             "orthogonal"                     : False          ,
                             "initialStep"                    : 0.1            ,
                             "lineSearcher"                   : None           ,
                             "numberOfSuccessiveLowGradients" : 2              ,
                             "maximalMinimizationSteps"       : 20             ,
                             "onlyMax"                        : True           ,
                             "tauTol"                         : 0.15           ,
                             "tauReached"                     : False          ,
                             "breakIfTauReached"              : False          ,
                             "stationaryPoint"                : False          ,
                             "betaType"                       : "PRP"          ,
                             "outputLevel"                    :  1             ,
                             "rmsGradientTolerance"           :  0.05          , # . kJ/mol
                             "rmsdMaximumTolerance"           :  1.5           } )
    _summarizable.update ( { "betaType"                       :"Beta Type"     ,
                             "initialStep"                    :"Initial Step"  ,
                             "lineSearcher"                   :"Line Searcher" } )

    def LogHeader (self, log=logFile):
        """Logging header information of CPRSaddle"""
        self.table = log.GetTable (columns=[ 10, 22, 24, 14, 14, 32])
        self.table.Start   ()
        self.table.Title   ( "RMS gradient tolerance: {:.3f}   Successive low gradients for saddle: {:d}   Tau tolerance: {:.2f} \n Options: Finite difference: {:s}   Break if tau reached: {:s}".format (self.rmsGradientTolerance, self.numberOfSuccessiveLowGradients, self.tauTol, str(self.finiteDifference), str(self.breakIfTauReached)))
        self.table.Heading ( "Iteration"   )
        self.table.Heading ( "Energy"      )
        self.table.Heading ( "Rel. Energy" )
        self.table.Heading ( "Tau"         )
        self.table.Heading ( "RMSGradient" )
        self.table.Heading ( "Successive low grad." )

    def LogIteration (self, state, tau, log=logFile):
        """Logging information of each CPRSaddle step"""
        self.table.Entry ("{:d}   ".format (state.numberOfIterations,))
        self.table.Entry ("{: .5f}    ".format (state.f,))
        self.table.Entry ("{: .2f}         ".format ( state.f - state.f00 ) )
        if type (tau) == float or type (tau) == int:
           self.table.Entry ("{:.6f}   ".format (tau,))
        else:
           self.table.Entry ("{:.8s}   ".format (tau,))
        self.table.Entry ("{: .6f}  ".format (state.rmsGradient,))
        self.table.Entry ("  {:d}            ".format (state.successiveLowGradients,))

    def LogTableStop (self, log=logFile):
        self.table.Stop ()

    def Continue (self, state):
        """Check if the refinement should continue"""
        if state.rmsGradient < self.rmsGradientTolerance:
             state.successiveLowGradients += 1
        else:
             state.successiveLowGradients = 0  # . Low gradients are no longer successive, start from 0.
        if self.onlyMax and state.numberOfIterations > 0:
              if (not state.isConverged):
                  self.LogTableStop ()
              return False
        elif (not state.isConverged) and state.numberOfIterations > 1 and  state.rmsGradient >= self.rmsGradientTolerance:
              state.x0.CopyTo (state.x)  # . Refuse the last line search.
              state.g0.CopyTo (state.g)
              state.f = state.f0
              self.LogTableStop ()
              return False  # . Do not go on after a crashed line search.
        elif ((state.numberOfIterations > self.maximalMinimizationSteps) or (state.numberOfIterations > (state.numberOfVariables - 1))):
              self.LogTableStop ()
              return False
        elif state.numberOfIterations > 1:
              tau = ((state.g.Dot (state.s0)) ** 2) / (state.g.DotSelf ( ) * state.s0.DotSelf ( ) )
              #tau = abs( math.sqrt (state.numberOfVariables/3) * (state.g.Dot (state.s0)) / (math.sqrt(state.g.Dot (state.g)) * math.sqrt(state.s0.Dot (state.s0))))
              if (self.outputLevel > 0) and (state.numberOfIterations == 0):
                 self.LogHeader ()
              if (self.outputLevel > 0) and (state.numberOfIterations % 50 == 0):
                 self.table.Title   ("\nRMS gradient tolerance: {:.3f}   Successive low gradients for saddle: {:d}   Tau tolerance: {:.2f} \n Options: Finite difference: {:s}   Break if tau reached: {:s}".format (self.rmsGradientTolerance, self.numberOfSuccessiveLowGradients, self.tauTol, str(self.finiteDifference), str(self.breakIfTauReached)))
                 self.table.Heading ("Iteration")
                 self.table.Heading ("Energy")
                 self.table.Heading ("Rel. Energy")
                 self.table.Heading ("Tau")
                 self.table.Heading ("RMSGradient")
                 self.table.Heading ("Successive low grad.")
              if (self.outputLevel > 0): self.LogIteration (state, tau)
              if tau > self.tauTol:
                  self.tauReached = True
                  #if (self.outputLevel > 0): logFile.Paragraph ("\nSearch direction is no longer conjugate to the path.") # Go on and minimize to a stationary point instead.
                  #self.LogTableStop ()
                  #return True
              if self.tauReached and ( state.successiveLowGradients == 0 ):
                  if self.breakIfTauReached:
                      logFile.Paragraph ( "\nTau was reached." )
                      self.LogTableStop ()
                      return False
                  else:
                      self.stationaryPoint = True
              if (self.tauReached and state.successiveLowGradients > 0):
                  self.stationaryPoint = True
              if (not state.isConverged) and (state.rmsGradient < self.rmsGradientTolerance) and state.successiveLowGradients >= (self.numberOfSuccessiveLowGradients / 10):
                  if (self.stationaryPoint):
                      state.saddle = False  # . It is not worth to minimize further because some plateau is reached. With the gradient below the threshold this image state is marked as a saddle point.
                      state.stationary = True
                      if (self.outputLevel > 0): logFile.Paragraph ("\nPoint seems to be on a plateau, but tau was reached, cannot minimize further, marking as stationary point.")
                  elif (not self.stationaryPoint):
                      state.saddle = True
                      if (self.outputLevel > 0): logFile.Paragraph ("\nPoint seems to be on a plateau, cannot minimize further, marking as saddle point.")
                  self.LogTableStop ()
                  return False
              elif state.successiveLowGradients >= (self.numberOfSuccessiveLowGradients):
                  if (self.stationaryPoint):
                      state.saddle = False
                      state.stationary = True                        
                      if (self.outputLevel > 0): logFile.Paragraph ("\nStationary point found instead of a saddle point, because tau was reached.")
                  elif (not self.stationaryPoint):
                      state.saddle = True
                      if (self.outputLevel > 0): logFile.Paragraph ("\nSaddle-point found.")
                  self.LogTableStop ()
                  return False 
              else:
                  return True              
        else:
              # . Do at least one iteration.
              if (self.outputLevel > 0):
                if state.numberOfIterations == 0:  
                  self.LogHeader ()
                self.LogIteration (state, "n.a.")  # . "n.a." for logging maximization, because tau is not known.
              return True

    def Initialize (self, state):
        """Initialization before iteration."""
        super (CPRSaddlePointRefinement, self).Initialize (state) 
        if (state.numberOfIterations== 0): 
          state.f0 = state.f
          state.f00 = state.f  # . Store the very initial energy.
          state.g.CopyTo (state.g0)
          state.x.CopyTo (state.x0)
        else:
          if (self.outputLevel > 1): logFile.Paragraph ("\nCPRSaddle: Initialize: state.f={:f}, state.f0={:f}".format (state.f, state.f0))
        state.successiveLowGradients = 0
        state.saddle = False
        state.stationary = False
        self.tauReached = False
        self.stationaryPoint = False
        parameters = _LineSearcherParameters
        parameters["initialStep"] = self.initialStep
        parameters["rmsdMaximumTolerance"] = self.rmsdMaximumTolerance
        if self.lineSearcher is None: self.lineSearcher = MoreThuenteLineSearcherWithMax.WithOptions ( **parameters )

    def InitialSearchDirection (self, state):
        """Determine an initial search direction."""
        state.s0.CopyTo (state.d)
        state.d.Normalize ()  # . Initial direction is passed in state.d.
        state.alpha = self.initialStep
        state.gDotG = state.g.DotSelf ( )

    def Iteration (self, state):
        """Perform an iteration."""
        try:
            #logFile.Paragraph ( "\nTau Fischer: %f" % (abs( math.sqrt (state.numberOfVariables/3) * (state.g.Dot (state.s0)) / ((math.sqrt(state.g.Dot (state.g)) * math.sqrt(state.s0.Dot (state.s0)))))))
            #logFile.Paragraph ( "Tau SinclairFletcher: %f" % (((state.g.Dot (state.s0)) ** 2) / (state.g.Dot (state.g) * state.s0.Dot (state.s0))))
            #logFile.Paragraph ( "Skalar projection (g onto s_0): %f" % (abs(((state.g.Dot (state.s0)) ** 1) / (math.sqrt(state.s0.Dot (state.s0))))))
            #logFile.Paragraph ( "Cosine: %f" % abs((state.g.Dot (state.s0)) / ((math.sqrt(state.g.Dot (state.g)) * math.sqrt(state.s0.Dot (state.s0))))) )
            self.NewSearchDirection (state)
            if state.numberOfIterations == 0:
                 self.lineSearcher.maximize = True
            else:
                 self.lineSearcher.maximize = False
            self.LineSearch (state)            
        except Exception as error:
            state.error = error.args[0]
            import traceback, sys
            traceback.print_exc (file=sys.stdout)
        state.numberOfIterations += 1

    def Iterate (self, state):
        """Apply the algorithm to a function."""
        self.Initialize   (state)
        while (self.Continue (state)):
            self.Iteration    (state)
        return state.Finalize ()
   
    def LineSearch (self, state):
        """Do a line search."""
        # . Do a line search for the maximum and return the best point.
        if(self.outputLevel > 1): logFile.Paragraph ("\nCPRSaddle: state.d.Norm2() = {:f}, |g|= {:f}, |g0| = {:f}".format (state.d.Norm2(), state.g.Norm2(), state.g0.Norm2()))
        if(self.outputLevel > 1): logFile.Paragraph ("\nCPRSaddle: Entered the LineSearch function in CPRSaddle module, state.f0={:.5f}, state.d.Dot ( state.g0 ) = {:f}".format (state.f0, state.d.Dot (state.g0)))
        self.lineSearcher.initialStep = state.alpha
        state.lineSearchObjectiveFunction.SetValuesAtOrigin (state.f0, state.d.Dot (state.g0))
        report = self.lineSearcher.Iterate (state.lineSearchObjectiveFunction, log=None)
        # . Get the results.
        state.alpha = report["Variable"      ]
        state.f = report["Function Value"]
        state.numberOfFunctionCalls += report.get ("Function Calls", 0)
        state.rmsGradient = state.g.RootMeanSquare ()
        # . Determine the line search step type (success, partial, error).
        state.isConverged = report["Converged"]  # . Convergence status of last minimization is kept in state.isConverge.
        if report["Converged"]:
            statusFlag = "s"
        elif state.f < state.f0:
            statusFlag = "p"
        else:
            statusFlag = "e"
            state.error = "Line search error: " + report["Status Message"] 
        state.stepType = "L" + repr ( report.get ("Iterations", 0) ) + statusFlag

    def NewSearchDirection (self, state):
        """Determine a new search direction."""
        if (state.numberOfIterations == 0) :
             self.InitialSearchDirection (state)
        else:
           if (state.numberOfIterations == 1):  # . First round is different from all others.
               # . Calculate the Hessian part for all rounds.
               if self.orthogonal:  # . Hessian for orthogonal directions.
                   state.s0.CopyTo (state.h)
                   state.h.Scale ( 1.0 / (state.s0.Norm2() ** 2))
               else:
                  if self.finiteDifference or (not state.isConverged):  # . Use finite difference to calculate Hessian.
                     if not state.isConverged:
                        state.x0.CopyTo (state.x)
                        state.g0.CopyTo (state.g)
                        state.f = state.f0
                     else:
                        state.x.CopyTo (state.x0)
                     state.x0.Add ( state.d, scale = self.finiteDifferenceStep )
                     state.f0 = state.objectiveFunction.FunctionGradients (state.x0, state.g0)
                  # . If detlaS > 0, initial state given before maximization will be used.
                  state.g.CopyTo (state.h)  # . Hessian approximation.
                  state.h.Add( state.g0, scale = -1.0 )
                  deltaGs0 = state.h.Dot (state.s0)
                  if deltaGs0 == 0 :
                     logFile.Paragraph ("H.s0 is 0 after unique line maximization!")
                     return
                  deltaGs0 = 1.0 / deltaGs0
                  state.h.Scale (deltaGs0)
               # . Calculate the first minimization direction.
               hg1 = state.h.Dot (state.g)
               state.s0.CopyTo (state.s)            
               state.s.Scale (hg1)
               state.s.Add ( state.g, scale = -1.0 )
           # . Directions for all other rounds.
           else:
               hgi = state.h.Dot (state.g)
               if self.betaType == "FR":
                  gDotg = state.g.DotSelf ( )
                  g0Dotg0 = state.g0.DotSelf ( ) 
               elif self.betaType == "PRP":
                  state.g.CopyTo (state.tmp)
                  state.tmp.Add ( state.g0, scale = -1.0 )
                  gDotg = state.g.Dot (state.tmp)
                  g0Dotg0 = state.g0.DotSelf ( ) 
               elif self.betaType == "HS":
                  state.g.CopyTo (state.tmp) 
                  state.tmp.Add ( state.g0, scale = -1.0 )
                  gDotg = state.g.Dot (state.tmp)
                  g0Dotg0 = state.s.Dot (state.tmp) 
               state.s.Scale (gDotg / g0Dotg0)
               state.s.Add ( state.s0, scale = hgi  )
               state.s.Add ( state.g , scale = -1.0 )
           # . Save data.
           state.f0 = state.f
           state.g.CopyTo (state.g0)
           state.x.CopyTo (state.x0)
           # . Set new search direction to si.
           state.s.CopyTo (state.d)
           state.d.Normalize ()

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

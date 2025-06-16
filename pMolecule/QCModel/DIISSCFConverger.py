"""Defines classes for the DIIS SCF converger."""

import collections, sys, traceback

from math                                   import fabs, pow, sqrt
from pCore                                  import AttributableObject        , \
                                                   logFile                   , \
                                                   LogFileActive
from pScientific.Arrays                     import Array                     , \
                                                   StorageType
from pScientific.LinearAlgebra              import LinearEquations
from pScientific.ObjectiveFunctionIterators import ObjectiveFunctionIterator

# . Eventually keep state so can reuse if DIIS, orbitals, etc. remain the same sizes ...

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DIISSCFConvergerState ( AttributableObject ):
    """Class for the DIIS SCF converger state."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "currentFrame"          : None  ,
                             "damp"                  : 0.0   ,
                             "dEAverage"             : 0.0   ,
                             "densitiesAreValid"     : False ,
                             "dEOld"                 : 0.0   ,
                             "diisError"             : 0.0   ,
                             "diisOn"                : False ,
                             "energy"                : None  ,
                             "eOld"                  : 0.0   ,
                             "error"                 : None  ,
                             "history"               : 0     ,
                             "inverseOrthogonalizer" : None  ,
                             "isConverged"           : False ,
                             "log"                   : None  ,
                             "numberOfFunctionCalls" : 0     ,
                             "numberOfIterations"    : 0     ,
                             "numberOfSpins"         : 0     ,
                             "overlapMatrix"         : None  ,
                             "orthogonalizer"        : None  ,
                             "printTraceback"        : False ,
                             "rcaMu"                 : 0.0   ,
                             "rmsDifference"         : 0.0   ,
                             "statusMessage"         : None  ,
                             "storeDamping"          : None  ,
                             "storeDIIS"             : None  ,
                             "storeODA"              : None  ,
                             "table"                 : None  ,
                             "target"                : None  ,
                             "workMa"                : None  ,
                             "workMb"                : None  ,
                             "workMc"                : None  ,
                             "workS"                 : None  } )

    def _Allocate ( self, m, n, maximumHistory, useODA ):
        """Allocate space for the state."""
        # . m is extent of error matrix, n is extent of density and Fock matrices, with m <= n.
        # . Stores.
        if useODA:
            self.storeODA = tuple ( [ ( Array.WithExtent ( n, storageType = StorageType.Symmetric ) ,
                                        Array.WithExtent ( n, storageType = StorageType.Symmetric ) ) for s in range ( self.numberOfSpins ) ] )
        else:
            self.storeDamping = collections.deque ( maxlen = 2 )
            for i in range ( 2 ): self.storeDamping.append ( tuple ( [ Array.WithExtent ( n, storageType = StorageType.Symmetric ) for s in range ( self.numberOfSpins ) ] ) )
        self.storeDIIS = collections.deque ( maxlen = maximumHistory )
        for i in range ( maximumHistory ):
            self.storeDIIS.append ( tuple ( [ ( Array.WithExtent ( n, storageType = StorageType.Symmetric     ) ,
                                                Array.WithExtent ( n, storageType = StorageType.Symmetric     ) ,
                                                Array.WithExtent ( m, storageType = StorageType.Antisymmetric ) , [] ) for s in range ( self.numberOfSpins ) ] ) )
        # . Work.
        self.workMa = Array.WithExtents ( n, n )
        self.workMb = Array.WithExtents ( n, n )
        self.workMc = Array.WithExtents ( n, n )
        self.workS  = Array.WithExtent  ( n, storageType = StorageType.Symmetric )

    def DensitiesToTarget ( self ):
        """Copy the densities to the target."""
        scratch = self.target.scratch
        scratch.onePDMP.isValid = self.densitiesAreValid
        if not self.target.electronicState.isSpinRestricted:
            densityA = self.currentFrame[0][0]
            densityB = self.currentFrame[1][0]
            densityP = scratch.onePDMP.density
            densityQ = scratch.onePDMQ.density
            densityA.CopyTo ( densityP ) ; densityP.Add ( densityB               )
            densityA.CopyTo ( densityQ ) ; densityQ.Add ( densityB, scale = -1.0 )
            scratch.onePDMQ.isValid = self.densitiesAreValid

    def Finalize ( self ):
        """Finalization."""
        report = { "Converged"      : self.isConverged           ,
                   "Function Calls" : self.numberOfFunctionCalls ,
                   "Function Value" : self.energy                ,
                   "Iterations"     : self.numberOfIterations    }
        if self.error is not None: report["Error"] = self.error
        return report

    def FockFromTarget ( self ):
        """Copy the Fock matrices from the target."""
        if not self.target.electronicState.isSpinRestricted:
            scratch = self.target.scratch
            fockA   = self.currentFrame[0][1]
            fockB   = self.currentFrame[1][1]
            fockP   = scratch.onePDMP.fock
            fockQ   = scratch.onePDMQ.fock
            fockP.CopyTo ( fockA ) ; fockA.Add ( fockQ               )
            fockP.CopyTo ( fockB ) ; fockB.Add ( fockQ, scale = -1.0 )

    @classmethod
    def FromTarget ( selfClass, target, maximumHistory, useODA ):
        """Constructor given an objective function."""
        self        = selfClass ( )
        self.target = target
        qcModel     = target.qcModel
        scratch     = target.scratch
        extent      = scratch.onePDMP.numberOrbitals
        # . Gather density and orbital data.
        qcModel.SetUpOrbitals ( target )
        # . Spin-restricted cases.
        if target.electronicState.isSpinRestricted:
            self.currentFrame = ( ( scratch.onePDMP.density ,
                                    scratch.onePDMP.fock    ,
                                    scratch.orbitalsP       ) , )
            self.densitiesAreValid = scratch.onePDMP.isValid
            self.numberOfSpins     = 1
        # . Spin-unrestricted cases.
        else:
            densityP = scratch.onePDMP.density
            densityQ = scratch.onePDMQ.density
            self.currentFrame = ( ( Array.WithExtent ( extent, storageType = StorageType.Symmetric ) ,
                                    Array.WithExtent ( extent, storageType = StorageType.Symmetric ) ,
                                    scratch.orbitalsP ) , 
                                  ( Array.WithExtent ( extent, storageType = StorageType.Symmetric ) ,
                                    Array.WithExtent ( extent, storageType = StorageType.Symmetric ) ,
                                    scratch.orbitalsQ ) )
            densityA = self.currentFrame[0][0]
            densityB = self.currentFrame[1][0]
            densityP.CopyTo ( densityA ) ; densityA.Add ( densityQ               ) ; densityA.Scale ( 0.5 )
            densityP.CopyTo ( densityB ) ; densityB.Add ( densityQ, scale = -1.0 ) ; densityB.Scale ( 0.5 )
            self.densitiesAreValid = ( scratch.onePDMP.isValid and scratch.onePDMQ.isValid )
            self.numberOfSpins     = 2
        # . Set the overlap matrix and orthogonalizers.
        self.inverseOrthogonalizer = scratch.Get ( "inverseOrthogonalizer", None )
        self.orthogonalizer        = scratch.Get ( "orthogonalizer"       , None )
        self.overlapMatrix         = scratch.Get ( "overlapMatrix"        , None )
        if ( self.inverseOrthogonalizer is None ) or \
           ( self.orthogonalizer        is None ) or \
           ( self.overlapMatrix         is None ):
            self.inverseOrthogonalizer = None
            self.orthogonalizer        = None
            self.overlapMatrix         = None
        # . Allocate all necessary space.
        if self.orthogonalizer is None: m = extent
        else:                           m = self.orthogonalizer.shape[1]
        self._Allocate ( m, extent, maximumHistory, useODA )
        #self._Allocate ( extent, extent, maximumHistory, useODA )
        return self

    def HandleError ( self, error ):
        """Handle an error."""
        self.error = error.args[0]
        if self.printTraceback: traceback.print_exc ( file = sys.stdout )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class DIISSCFConverger ( ObjectiveFunctionIterator ):
    """Class for the DIIS SCF converger."""

    _attributable = dict ( ObjectiveFunctionIterator._attributable )
    _classLabel   = "DIIS SCF Converger"
    _stateObject  = DIISSCFConvergerState
    _summarizable = dict ( ObjectiveFunctionIterator._summarizable )
    _attributable.update ( { "dampEnergyTolerance"        : 2.0e-04                        , 
                             "dampExtrapolationFrequency" : 15                             , 
                             "dampMaximum"                : 256.0                          , 
                             "dampScaleFactor"            : 16.0                           , 
                             "dampTolerance"              : 1.0                            , 
                             "densityTolerance"           : 1.0e-12                        , 
                             "diisDeviation"              : 1.0e-06                        , 
                             "diisOnset"                  : 0.2                            , 
                             "energyTolerance"            : 2.0e-4                         , 
                             "logFrequency"               : 1                              , 
                             "maximumHistory"             : 10                             , 
                             "maximumIterations"          : 100                            , 
                             "minimumMu"                  : 1.0e-02                        , 
                             "rcaOnset"                   : 0.8                            , 
                             "useODA"                     : True                           } )
    _summarizable.update ( { "dampEnergyTolerance"        : "Damp Energy Tolerance"        ,
                             "dampExtrapolationFrequency" : "Damp Extrapolation Frequency" ,
                             "dampMaximum"                : "Damp Maximum"                 ,
                             "dampScaleFactor"            : "Damp Scale Factor"            ,
                             "dampTolerance"              : "Damp Tolerance"               ,
                             "densityTolerance"           : "Density Tolerance"            ,
                             "diisDeviation"              : "Diis Deviation"               ,
                             "diisOnset"                  : "Diis Onset"                   ,
                             "energyTolerance"            : "Energy Tolerance"             ,
                             "logFrequency"               : "Log Frequency"                ,
                             "maximumHistory"             : "Maximum History"              ,
                             "maximumIterations"          : "Maximum Iterations"           ,
                             "minimumMu"                  : "Minimum Mu"                   ,
                             "rcaOnset"                   : "ODA Onset"                    ,
                             "useODA"                     : "Use ODA"                      } )

    def Continue ( self, state ):
        """Check to see if the calculation should continue."""
        state.isConverged = ( state.numberOfIterations >  0                     ) and \
                            ( state.rmsDifference      <= self.densityTolerance ) and \
                            ( state.dEOld              <= self.energyTolerance  )
        if state.isConverged:
            state.statusMessage = "SCF converged."
        else:
            if state.numberOfIterations >= self.maximumIterations: state.error = "Too many iterations."
            if state.error is not None: state.statusMessage = "SCF error: " + state.error
        return ( state.statusMessage is None )

    def DavidsonDampingFactor ( self, iteration, dE, dEOld, dEAverage, damp ):
        """Get a value for the Davidson damping factor and the average energy difference."""
        if iteration > 0:
            if iteration == 1: dEAverage =   fabs ( dE )
            else:              dEAverage = ( fabs ( dE ) + fabs ( dEOld ) + ( 0.2 * dEAverage ) ) / 2.2
            biggest = max ( damp, dEAverage )
            if dE > 0.0:
                damp = 4.0 * biggest
                if dEOld > 0.0:
                    if   dE >= ( 4.0  * dEOld ): damp  = self.dampScaleFactor * biggest
                    elif dE <= ( 0.25 * dEOld ): damp /= self.dampScaleFactor
                    else:                        damp  = pow ( dE / dEOld, 2 ) * biggest
                else:
                    if   dE           > ( 0.5 * dEAverage ): damp *= self.dampScaleFactor
                    if ( dE - dEOld ) < ( 0.2 * dEAverage ): damp /= self.dampScaleFactor
            else:
                if dEOld > 0.0:
                    damp = 4.0 * biggest
                    if   -dE           > dEAverage: damp *= self.dampScaleFactor
                    if ( -dE + dEOld ) < dEAverage: damp /= self.dampScaleFactor
                else:
                    if dE > dEOld:
                        if dE <= ( 0.25 * dEOld ): damp  = pow ( dE / dEOld, 2 ) * biggest
                        else:                      damp /= self.dampScaleFactor
                    else:
                        if   fabs ( dE ) >= ( 2.0 * dEAverage ): damp  = self.dampScaleFactor * biggest
                        elif fabs ( dE ) <= ( 0.5 * dEAverage ): damp /= self.dampScaleFactor
            damp = min ( damp, self.dampMaximum )
        return ( damp, dEAverage )

    def DavidsonDampingIterate ( self, state ):
        """Perform a Davidson damping iteration."""
        damp            = state.damp
        doExtrapolation = ( self.dampExtrapolationFrequency > 0 ) and \
                          ( state.numberOfIterations        > 1 ) and \
                          ( state.numberOfIterations % self.dampExtrapolationFrequency == 0 )
        for ( s, ( _, f0, _ ) ) in enumerate ( state.currentFrame ):
            f1 = state.storeDamping[0][s]
            f0.Add ( f1, scale = damp ) ; f0.Scale ( 1.0 / ( 1.0 + damp ) )
            if doExtrapolation:
                f2 = state.storeDamping[1][s]
                f3 = state.workS
                f1.CopyTo ( f3 )
                f2.Add ( f1, scale = -1.0 ) # f2 - f1
                f3.Add ( f0, scale = -1.0 ) # f1 - f0
                f11 = f3.Dot ( f3 ) # . Dot is fine here (as not using iterators).
                f12 = f3.Dot ( f2 )
                f22 = f2.Dot ( f2 )
                eps = ( f11 - f21 ) / ( f21 - f22 )
                if ( f21 * f21 ) < ( 0.5 * f11 * f22 ): eps = 0.0
                eps = min ( 0.5, eps )
                f0.Add ( f1, scale = -eps ) ; f0.Scale ( 1.0 / ( 1.0 - eps ) )

    def DavidsonDampingSave ( self, state ):
        """Save the Fock matrix."""
        state.storeDamping.rotate ( 1 ) # 1 -> 2
        for ( s, ( _, f, _ ) ) in enumerate ( state.currentFrame ):
            f.CopyTo ( state.storeDamping[0][s] ) # 0 -> 1

    def DIISCurrentFrame ( self, state ):
        """Process the current frame and calculate the error vectors."""
        # . Newest data is at the front of the store.
        # . Set history counters.
        state.history = min ( state.history + 1, self.maximumHistory )
        state.storeDIIS.rotate ( 1 )
        # . Array aliases.
        s = state.overlapMatrix
        u = state.workMa
        v = state.workMb
        w = state.workMc
        x = state.orthogonalizer
        y = state.inverseOrthogonalizer
        # . Store the matrices and create the error vectors.
        # . When S == I:
        #     e = F * D - D * F
        # . When S != I there are various possibilities (only the last is used):
        #     e = F * D * S - S * D * F                 - AO basis
        #     e = X^T * ( F * D * S - S * D * F ) * X   - MO basis
        #     e = X^T * F * D * Y - Y^T * D * F * X     - MO basis but saving one matrix multiplication
        diisError = 0.0
        for ( c, ( d, f, _ ) ) in enumerate ( state.currentFrame ):
            ( d0, f0, e0, c0 ) = state.storeDIIS[0][c]
            d.CopyTo ( d0 )
            f.CopyTo ( f0 )
            if s is None: e0.MakeCommutatorSS   ( f, d, mA = u, mB = v, mC = w )
            else:         e0.MakeCommutatorXSSY ( f, d, x, y, u, v, w, xTranspose = True )
#            else:         e0.MakeCommutatorSSS   ( f, d, s )
#            else:         e0.MakeCommutatorXSSSX ( f, d, s, x, u, v, w )
            c0.clear ( )
            c0.extend ( [ - e0.TraceOfProduct ( state.storeDIIS[h][c][2] ) for h in range ( state.history ) ] )
            diisError = max ( diisError, e0.AbsoluteMaximum ( ) )
        return diisError

    def DIISIterate ( self, state ):
        """Apply the DIIS procedure."""
        # . Set up and solve the linear equations.
        if state.history > 1:
            # . Loop until there is success or there is no more history.
            while True:
                # . Create A and B.
                n = state.history + 1
                a = Array.WithExtent ( n, storageType = StorageType.Symmetric ) ; a.Set ( 0.0 )
                b = Array.WithExtent ( n ) ; b.Set ( 0.0 )
                for i in range ( state.history ):
                    a[state.history,i] = -1.0
                    for s in range ( state.numberOfSpins ):
                        c = state.storeDIIS[i][s][3]
                        for ( j, v ) in zip ( range ( i, state.history ), c ):
                            a[i,j] += v
                b[state.history] = -1.0
                # . Solve the linear equations.
                isOK = False
                try:
                    LinearEquations ( a, b, preserveInput = False )
                    b[state.history] = -1.0
                    isOK = ( fabs ( sum ( b ) ) < self.diisDeviation )
                except:
                    pass
                if isOK: break
                # . An ill-conditioned solution.
                state.history -= 1
                if state.history <= 1: return
            # . Calculate the new density and Fock matrices.
            for ( s, ( d, f, _ ) ) in enumerate ( state.currentFrame ):
                d.Set ( 0.0 )
                f.Set ( 0.0 )
                for h in range ( state.history ):
                    c = b[h]
                    d.Add ( state.storeDIIS[h][s][0], scale = c )
                    f.Add ( state.storeDIIS[h][s][1], scale = c )

    def FunctionGradients ( self, state ):
        """Evaluate the function and its gradients."""
        state.DensitiesToTarget ( )
        energy = 0.0
        for closure in state.target.qcState.fockClosures: energy += ( closure ( ) )
        state.FockFromTarget ( )
        state.energy                 = energy
        state.numberOfFunctionCalls += 1

    def Initialize ( self, state ):
        """Initialization before iteration."""
        try:
            self.FunctionGradients ( state )
        except Exception as error:
            state.HandleError ( error )

    def Iterate ( self, target, log = logFile ):
        """Apply the algorithm to a target."""
        state = self.__class__._stateObject.FromTarget ( target, self.maximumHistory, self.useODA )
        self.LogStart     ( state, log = log )
        self.Initialize   ( state )
        self.LogIteration ( state )
        while ( self.Continue ( state ) ):
            self.Iteration    ( state )
            self.LogIteration ( state )
        self.LogStop ( state )
        return state.Finalize ( )

    def Iteration ( self, state ):
        """Perform an iteration."""
        try:
            self.ModifyFockMatrices ( state )
            self.MakeDensities      ( state )
            self.FunctionGradients  ( state )
        except Exception as error:
            state.HandleError ( error )
        state.numberOfIterations += 1

    def LogIteration ( self, state ):
        """Log an iteration."""
        if ( state.table is not None ) and ( state.numberOfIterations % self.logFrequency == 0 ):
            state.table.Entry ( "{:d}".format ( state.numberOfIterations ) )
            state.table.Entry ( "{:20.8f}".format ( state.energy         ) )
            state.table.Entry ( "{:20.8f}".format ( state.dEOld          ) )
            state.table.Entry ( "{:20.8f}".format ( state.rmsDifference  ) )
            state.table.Entry ( "{:20.8f}".format ( state.diisError      ) )
            if state.numberOfIterations == 0:
                state.table.Entry  ( "Init." )
                state.table.EndRow ( )
            elif state.diisOn:
                state.table.Entry ( "DIIS" )
                state.table.Entry ( "{:d}".format ( state.history ) )
            elif self.useODA:
                state.table.Entry ( "ODA" )
                state.table.Entry ( "{:12.4g}".format ( state.rcaMu ) )
            else:
                state.table.Entry ( "Damping" )
                state.table.Entry ( "{:12.4g}".format ( state.damp  ) )

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        if LogFileActive ( log ):
            state.log = log
            table = log.GetTable ( columns = [ 10, 20, 20, 20, 20, 8, 12 ] )
            table.Start   ( )
            table.Heading ( "Cycle"              )
            table.Heading ( "Energy"             )
            table.Heading ( "Energy Change"      )
            table.Heading ( "RMS Density Change" )
            table.Heading ( "Max. DIIS Error"    )
            table.Heading ( "Operation", columnSpan = 2 )
            state.table = table

    def LogStop ( self, state ):
        """Stop logging."""
        if state.table is not None:
            state.table.Stop ( )
            state.table = None
            if state.statusMessage is not None: state.log.Paragraph ( state.statusMessage )
            items = [ ( "Converged"     , "{:s}".format ( repr ( state.isConverged  ) ) ) ,
                      ( "Function Value", "{:g}".format ( state.energy                ) ) ,
                      ( "Iterations"    , "{:d}".format ( state.numberOfIterations    ) ) ,
                      ( "Function Calls", "{:d}".format ( state.numberOfFunctionCalls ) ) ]
            state.log.SummaryOfItems ( items, title = self.__class__._classLabel + " Report" )

    def MakeDensities ( self, state ):
        """Make the densities."""
        old           = state.workS
        scratch       = state.target.scratch
        rmsDifference = 0.0
        for ( d, f, o ) in state.currentFrame:
            o.MakeFromFock ( f, scratch, orthogonalizer = state.orthogonalizer )
            d.CopyTo ( old )
            d.MakeFromEigenSystem ( o.occupancyHandler.numberOccupied, o.occupancies, o.orbitals )
            old.Add ( d, scale = -1.0 )
            rmsDifference = max ( rmsDifference, old.RootMeanSquare ( ) )
        state.densitiesAreValid = True
        state.rmsDifference     = rmsDifference

    def ModifyFockMatrices ( self, state ):
        """Modify the Fock matrices for subsequent diagonalization."""
        dE        = state.energy - state.eOld
        diisError = 0.0
        if state.densitiesAreValid:
            if not self.useODA:
                ( state.damp, state.dEAverage ) = self.DavidsonDampingFactor ( state.numberOfIterations, dE, state.dEOld, state.dEAverage, state.damp )
            diisError = self.DIISCurrentFrame ( state )
            if state.numberOfIterations > 0:
                if self.useODA:
                    if   state.diisOn                : state.diisOn = ( diisError < self.rcaOnset  )
                    elif state.numberOfIterations > 1: state.diisOn = ( diisError < self.diisOnset )
                else:
                    if fabs ( dE ) < self.dampEnergyTolerance:
                        state.diisOn = state.diisOn or ( state.damp < self.dampTolerance )
                    else:
                        if state.diisOn: state.damp = self.dampTolerance
                        state.diisOn = False
                if   state.diisOn: self.DIISIterate ( state )
                elif self.useODA : state.energy = self.ODAIterate ( state, state.eOld, state.energy )
                else:              self.DavidsonDampingIterate ( state )
            if self.useODA: self.ODASave             ( state )
            else:           self.DavidsonDampingSave ( state )
        state.dEOld     = dE
        state.diisError = diisError
        state.eOld      = state.energy
        state.energy    = None

    def ODAIterate ( self, state, eOld, energy ):
        """Apply the ODA procedure."""
        # . Calculate the polynomial coefficients.
        A = 0.0
        C = 0.0
        for ( s, ( d, f, _ ) ) in enumerate ( state.currentFrame ):
            dH  = state.storeODA[s][0]
            fH  = state.storeODA[s][1]
            A  += ( d.TraceOfProduct ( f  ) - dH.TraceOfProduct ( f  ) )
            C  += ( d.TraceOfProduct ( fH ) - dH.TraceOfProduct ( fH ) )
        D  = eOld
        A += ( C + 2.0 * ( D - energy ) )
        B  = energy - A - C - D
        # . Find mu by minimizing the cubic polynomial in the range [0,1].
        # . Find the mu at the boundary with the smallest energy (either 0 or 1).
        fMin   = energy
        muBMin = 1.0  
        if eOld < energy:
            fMin   = eOld
            muBMin = 0.0
        # . Solve the polynomial.
        if A == 0.0:
            # . Line.
            if B == 0.0:
                mu = muBMin
            # . Quadratic.
            else:
                muT = - C / ( 2.0 * B )
                if ( B > 0.0 ) and ( muT > 0.0 ) and ( muT < 1.0 ): mu = muT
                else: mu = muBMin
        # . Cubic.
        else:
            fac = B * B - 3.0 * A * C
            # . Real turning points.
            if fac >= 0.0:
                muT1 = ( - B + sqrt ( fac ) ) / ( 3.0 * A )
                muT2 = ( - B - sqrt ( fac ) ) / ( 3.0 * A )
                if ( muT1 > 0.0 ) and ( muT1 < 1.0 ):
                    fac = ( ( A * muT1 + B ) * muT1 + C ) * muT1 + D
                    if fac < fMin:
                        fMin   = fac
                        muBMin = muT1
                if ( muT2 > 0.0 ) and ( muT2 < 1.0 ):
                    fac = ( ( A * muT2 + B ) * muT2 + C ) * muT2 + D
                    if fac < fMin:
                        fMin   = fac
                        muBMin = muT2
                mu = muBMin
            # . No real turning points.
            else: mu = muBMin
        # . Constraint mu.
        mu = max ( mu, self.minimumMu )
        # . Form the final matrices.
        nu = 1.0 - mu
        for ( s, ( d, f, _ ) ) in enumerate ( state.currentFrame ):
            d.Scale ( mu ) ; d.Add ( state.storeODA[s][0], scale = nu )
            f.Scale ( mu ) ; f.Add ( state.storeODA[s][1], scale = nu )
        # . Save the mu parameter.
        state.rcaMu = mu
        # . Reset energy (assuming the interpolation is a reasonable approximation).
        energy = ( ( A * mu + B ) * mu + C ) * mu + D
        return energy

    def ODASave ( self, state ):
        """Save the ODA data."""
        for ( s, ( d, f, _ ) ) in enumerate ( state.currentFrame ):
            d.CopyTo ( state.storeODA[s][0] )
            f.CopyTo ( state.storeODA[s][1] )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

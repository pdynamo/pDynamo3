"""Classes for performing quasi-Newton multidimensional minimization."""

import math

from  .MultiDimensionalMinimizer      import MultiDimensionalMinimizer      , \
                                             MultiDimensionalMinimizerState
from  .ObjectiveFunctionIteratorError import ObjectiveFunctionIteratorError
from ..Arrays                         import Array                          , \
                                             StorageType     
from ..LinearAlgebra                  import EigenPairs      

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QuasiNewtonMinimizerState ( MultiDimensionalMinimizerState ):
    """Quasi-Newton minimizer state."""

    _attributable = dict ( MultiDimensionalMinimizerState._attributable )
    _attributable.update ( { "eigenValues"  : None ,
                             "eigenVectors" : None ,
                             "fOld"         : None ,
                             "gOld"         : None ,
                             "hessian"      : None ,
                             "rfoMatrix"    : None ,
                             "rfoScale"     : None ,
                             "trustRadius"  : None ,
                             "w"            : None } )

    def SetUp ( self ):
        """Set up the state."""
        super ( QuasiNewtonMinimizerState, self ).SetUp ( )
        # . Allocation.
        n = self.numberOfVariables
        p = n + 1
        self.eigenValues  = Array.WithExtent  ( p    ) ; self.eigenValues.Set ( 0.0 )
        self.gOld         = Array.WithExtent  ( n    ) ; self.gOld.Set        ( 0.0 )
        self.w            = Array.WithExtent  ( n    ) ; self.w.Set           ( 0.0 )
        self.eigenVectors = Array.WithExtents ( p, p )
        self.rfoMatrix    = Array.WithExtent  ( p, storageType = StorageType.Symmetric ) ; self.rfoMatrix.Set ( 0.0 )
        # . Initialization.
        self.rfoScale = 1.0e+00 / math.sqrt ( float ( n ) ) # . The RFO scale factor.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QuasiNewtonMinimizer ( MultiDimensionalMinimizer ):
    """Quasi-Newton minimization."""

    _attributable = dict ( MultiDimensionalMinimizer._attributable )
    _classLabel   = "Quasi-Newton Minimizer"
    _stateObject  = QuasiNewtonMinimizerState
    _summarizable = dict ( MultiDimensionalMinimizer._summarizable )
    _attributable.update ( { "initialStep"          : 0.3e+00                  ,
                             "maximumTrust"         : 1.0e+00                  ,
                             "minimumTrust"         : 1.0e-02                  ,
                             "trustScaleDownBound"  : 0.25e+00                 ,
                             "trustScaleFactor"     : 2.0e+00                  ,
                             "trustScaleUpBound"    : 0.75e+00                 ,
                             "trustUpdateTolerance" : 1.0e-03                  } )
    _summarizable.update ( { "initialStep"          : "Initial Step"           ,
                             "maximumTrust"         : "Maximum Trust Radius"   ,
                             "minimumTrust"         : "Minimum Trust Radius"   ,
                             "trustScaleDownBound"  : "Trust Scale Down Bound" ,
                             "trustScaleFactor"     : "Trust Scale Factor"     ,
                             "trustScaleUpBound"    : "Trust Scale Up Bound"   ,
                             "trustUpdateTolerance" : "Trust Update Tolerance" } )

    def DetermineStep ( self, state ):
        """Get the step."""
        if state.numberOfIterations > 0:
            # . Calculate the predicted change in the function value.
            state.hessian.VectorMultiply ( state.d, state.w )
            state.w.Scale ( 0.5 )
            state.w.Add ( state.gOld )
            dFPredicted = state.d.Dot ( state.w )
            # . Calculate the change in the function value.
            dFActual = state.f - state.fOld
            # . Update the hessian.
            state.g.CopyTo ( state.w )
            state.w.Add ( state.gOld, scale = -1.0 )
            state.hessian.Update   ( state.d, state.w )
            # . Update the trust radius.
            self.UpdateTrustRadius ( state, dFActual, dFPredicted )
        # . Set up and solve the RFO equations.
        # . Needs to be more efficient.
        scale = 1.0 / math.sqrt ( state.rfoScale )
        for i in range ( state.numberOfVariables ):
            for j in range ( i+1 ):
                state.rfoMatrix[i,j] = state.hessian[i,j] * scale**2
            state.rfoMatrix[state.numberOfVariables,i] = state.g[i] * scale
        state.rfoMatrix[state.numberOfVariables,state.numberOfVariables] = 0.0
        EigenPairs ( state.rfoMatrix, state.eigenValues, state.eigenVectors )
        # . The step corresponds to the appropriate elements of the intermediately-normalized eigenVector of smallest eigenValue.
        scale /= state.eigenVectors[state.numberOfVariables,0]
        for i in range ( state.numberOfVariables ): state.d[i] = state.eigenVectors[i,0]
        state.d.Scale ( scale )
        # . Apply constraints.
        state.objectiveFunction.ApplyConstraintsToVector ( state.d )
        # . Check the step length.
        stepSize = state.d.Norm2 ( )
        if stepSize > state.trustRadius:
            state.d.Scale ( state.trustRadius / stepSize )
            stepSize = state.trustRadius
        # . Save various old values.
        state.fOld = state.f
        state.g.CopyTo ( state.gOld )
        # . Determine the new point.
        state.x.Add ( state.d )
        # . Finish up.
        state.stepType = "N"

    def Initialize ( self, state ):
        """Initialization before iteration."""
        super ( QuasiNewtonMinimizer, self ).Initialize ( state )
        # . Get the starting hessian.
        if not hasattr ( state.objectiveFunction, "StartingHessian" ): raise ObjectiveFunctionIteratorError ( "Objective function has no hessian evaluator." )
        state.hessian     = state.objectiveFunction.StartingHessian ( state.x )
        # . Other initialization.
        state.fOld        = state.f
        state.trustRadius = self.initialStep

    def UpdateTrustRadius ( self, state, dFActual, dFPredicted ):
        """Update the trust radius."""
        if ( math.fabs ( dFActual ) > self.trustUpdateTolerance ) and ( math.fabs ( dFPredicted ) > self.trustUpdateTolerance ):
            ratio = dFActual / dFPredicted
            if   ratio < self.trustScaleDownBound: state.trustRadius /=             self.trustScaleFactor
            elif ratio > self.trustScaleUpBound  : state.trustRadius *= math.sqrt ( self.trustScaleFactor )
        if state.trustRadius < self.minimumTrust: state.trustRadius = self.minimumTrust
        if state.trustRadius > self.maximumTrust: state.trustRadius = self.maximumTrust

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

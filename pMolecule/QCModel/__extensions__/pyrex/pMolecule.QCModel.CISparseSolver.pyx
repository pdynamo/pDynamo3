"""A aparse CI solver using the primme library."""

import math

from  pCore              import DataType
from  pScientific.Arrays import Array
from .QCModelError       import QCModelError

# . Eigenvectors are stored row-wise!

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CISparseSolver:
    """A CI sparse solver using the primme library."""

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

    def _Initialize ( self ):
        """Initialization."""
        # . Python objects.
        self.eigenValues    = None
        self.eigenVectors   = None
        self.iWork          = None
        self.matrix         = None
        self.preconditioner = None
        self.residualNorms  = None
        self.rWork          = None
        # . Primme.
        primme_initialize ( &(self.cObject) )
        self.cObject.eps                           = 1.0e-12
        self.cObject.maxMatvecs                    = 100000
        self.cObject.printLevel                    = 0                    
        self.cObject.correctionParams.precondition = 1

    def ApplyMatrix ( self, RealArray1D x not None ,
                            RealArray1D y not None ):
        """Apply the matrix to x in y."""
        CISparseSolver_ApplyMatrix ( RealArray1D_PointerToData ( x.cObject ) ,
                                     RealArray1D_PointerToData ( y.cObject ) ,
                                     NULL, &(self.cObject) )

    def ApplyPreconditioner ( self, RealArray1D x not None ,
                                    RealArray1D y not None ):
        """Apply the preconditioner to x in y."""
        CISparseSolver_ApplyPreconditioner ( RealArray1D_PointerToData ( x.cObject ) ,
                                             RealArray1D_PointerToData ( y.cObject ) ,
                                             NULL, &(self.cObject) )

    def CheckSolution ( self, referenceEigenvalues = None ):
        """Check the solution."""
        temporary  = Array.WithExtent ( self.cObject.n )
        deviation1 = 0.0
        deviation2 = 0.0
        for i in range ( self.cObject.numEvals ):
            row = self.eigenVectors[i,:]
            self.ApplyMatrix ( row, temporary )
            temporary.Add ( row, scale = -self.eigenValues[i] )
            deviation1 = max ( deviation1, temporary.AbsoluteMaximum ( ) )
            deviation2 = max ( deviation2, math.fabs ( row.Norm2 ( ) - 1.0 ) )
        results = { "Eigenvector Error"   : deviation1 ,
                    "Normalization Error" : deviation2 }
        if referenceEigenvalues is not None:
            deviation3 = 0.0
            for ( e, e0 ) in zip ( self.eigenValues, referenceEigenvalues ):
                deviation3 = max ( deviation3, math.fabs ( e - e0 ) )
            results["Eigenvalue Error"] = deviation3
        return results

    @classmethod
    def FromArrays ( selfClass, RealArray1D           eigenValues  not None ,
                                RealArray2D           eigenVectors not None ,
                                SparseSymmetricMatrix matrix       not None ,
                                RealArray1D           preconditioner        ):
        """Set up the solver from arrays."""
        cdef CISparseSolver self
        cdef int            n, v
        cdef CRealArray1D  *cPreconditioner = NULL
        ( v, n ) = eigenVectors.shape # . Row-wise storage.
        isOK     = ( n > 0 ) and ( v > 0 ) and ( v == len ( eigenValues ) ) and ( n == matrix.shape[0] ) and \
                   ( ( preconditioner is None ) or ( ( preconditioner is not None ) and ( n == len ( preconditioner ) ) ) )
        if not isOK: raise QCModelError ( "Error setting up CI sparse matrix solver." )
        self = selfClass ( )
        # . Primme data.
        # . Basic initialization.
        if preconditioner is not None: cPreconditioner = preconditioner.cObject
        self.cObject.n                   = n # . Matrix size.
        self.cObject.numEvals            = v # . Number of eigenvectors.
        self.cObject.matrixMatvec        = CISparseSolver_ApplyMatrix
        self.cObject.applyPreconditioner = CISparseSolver_ApplyPreconditioner
        self.cObject.matrix              = <void*> matrix.cObject
        self.cObject.preconditioner      = <void*> cPreconditioner
        # . Assignment.
        self.eigenValues    = eigenValues 
        self.eigenVectors   = eigenVectors
        self.matrix         = matrix
        self.preconditioner = preconditioner
        # . Allocate space.
        primme_set_method ( DYNAMIC, &(self.cObject) ) # . Needs to be done here after n, etc. set up.
        dprimme ( NULL, NULL, NULL, &(self.cObject) )
        self.residualNorms    = Array.WithExtent ( n )
        self.iWork            = Array.WithExtent ( self.cObject.intWorkSize  // sizeof ( CInteger ) , dataType = DataType.Integer )
        self.rWork            = Array.WithExtent ( self.cObject.realWorkSize // sizeof ( CReal    ) )
        self.cObject.intWork  = IntegerArray1D_PointerToData ( self.iWork.cObject )
        self.cObject.realWork = RealArray1D_PointerToData    ( self.rWork.cObject )
        # . Finish up.
        return self

    def SetOptions ( self, errorTolerance                     = None ,
                           maximumMatrixVectorMultiplications = None ,
                           printLevel                         = None ,
                           usePreconditioning                 = None ):
        """Set some options."""
        if errorTolerance                     is not None: self.cObject.eps        = errorTolerance
        if maximumMatrixVectorMultiplications is not None: self.cObject.maxMatvecs = maximumMatrixVectorMultiplications
        if printLevel                         is not None: self.cObject.printLevel = printLevel
        if usePreconditioning                 is not None: 
            if usePreconditioning: self.cObject.correctionParams.precondition = 1
            else:                  self.cObject.correctionParams.precondition = 0

    def Solve ( self ):
        """Solve for the eigenvalues and eigenvectors."""
        cdef int returnCode
        # . Simple case which gives an error in dprimme.
        if self.cObject.n == 1:
            self.eigenVectors.Set ( 1.0 )
            CISparseSolver_ApplyMatrix ( RealArray2D_PointerToData ( self.eigenVectors.cObject ) ,
                                         RealArray1D_PointerToData ( self.eigenValues.cObject  ) ,
                                         NULL, &(self.cObject) )
            return { "Converged"                     : True ,
                     "Converged Pairs"               : 1    ,
                     "Matrix-Vector Multiplications" : 1    ,
                     "Return Code"                   : 0    }
        # . Normal cases.
        else:
#            primme_display_params ( self.cObject )
            self.cObject.initSize = 0 # . Necessary for restart.
            returnCode = dprimme ( RealArray1D_PointerToData ( self.eigenValues.cObject   ) ,
                                   RealArray2D_PointerToData ( self.eigenVectors.cObject  ) ,
                                   RealArray1D_PointerToData ( self.residualNorms.cObject ) ,
                                   &(self.cObject) )
            return { "Converged"                     : ( returnCode == 0 ) and ( self.cObject.initSize == self.cObject.numEvals ) ,
                     "Converged Pairs"               : self.cObject.initSize         ,
                     "Matrix-Vector Multiplications" : self.cObject.stats.numMatvecs ,
                     "Return Code"                   : returnCode                    }

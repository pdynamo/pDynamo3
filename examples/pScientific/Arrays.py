"""Tests for real N-D arrays."""

import math, os, os.path

from Definitions        import outPath
from pCore              import CPUTime             , \
                               DeepClone           , \
                               logFile             , \
                               ShallowClone        , \
                               TestScriptExit_Fail , \
                               YAMLPickle          , \
                               YAMLUnpickle
from pScientific.Arrays import Array               , \
                               ArrayPrint          , \
                               ArrayPrint2D        , \
                               Flatten             , \
                               Reshape

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_AValue    = 10.0
_CValue    = 20.0
_Power     = 3.9
_Print     = True
_Tolerance = 1.0e-06
_VSet      = [ -13.0, -7.0, 11.0, 17.0 ]

# . Derived options.
_APower    = math.pow ( _AValue, _Power )

# . __getitem__ is used instead of [] so that indexes and slices can be predefined.
# . Indexes.
_Indexes   = [ [ 1, 8 ] ,
               [ ( 1, 1 ), ( 2, 4 ) ] ,
               [ ( 3, 0, 3 ), ( 5, 1, 2 ) ] ,
               [ ( 4, 2, 1, 1 ), ( 0, 0, 0, 2 ) ] ]

# . Shapes.
_Shapes    = [ ( 9, ), ( 3, 5 ), ( 6, 2, 4 ), ( 7, 3, 2, 3 ) ]

# . Reshaping (all good).
_Reshapes  = [ [ (  9,  1 ), (          3, 3 ), ( 1, 3, 3 ), ( 1, 9, 1, 1 ) ] ,
               [ (  5,  3 ), (       5, 1, 3 ) ] ,
               [ ( 24,  2 ), ( 3, 2, 2, 2, 2 ) ] ,
               [ (  7, 18 ), (       7, 9, 2 ) ] ]

# . Slices (with flattening flag and reshaping tests).
_Slices    = [ [ ( slice ( 1, 6, 2 ), True, [] ) ] ,
               [ ( ( 2, slice ( None, None, None ) ), True, [] ) ,
                 ( ( slice ( 0, 2, 1 ), slice ( 1, 2, 1 ) ), True, [] ) ] ,
               [ ( ( 5, slice ( 0, 1, 1 ), 2 ), True, [] ) ,
                 ( ( slice ( None, None, None ), 1, slice ( 2, 3, 1 ) ), True, [ ( ( 3, 2, 1 ), True ) ] ) ,
                 ( ( slice ( 0, 5, 2 ), slice ( None, None, None ), slice ( 1, 3, 1 ) ), False, [ ( ( 3, 1, 2, 2 ), True ) ] ) ] ,
               [ ( ( 3, 2, slice ( 0, 1, 1 ), 1 ), True, [] ) ,
                 ( ( slice ( 2, 5, 1 ), slice ( 0, 2, 2 ), 1, 2 ), True, [] ) ,
                 ( ( slice ( None, None, None ), 0, slice ( None, None, None ), slice ( None, None, None ) ), False, [] ) ,
                 ( ( slice ( 3, 7, 1 ), slice ( None, None, None ), slice ( 0, 2, 1 ), slice ( 0, 2, 2 ) ), True, [ ( ( 2, 2, 3, 2 ), True ), ( ( 12, 2 ), True ) ] ) ] ]

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK     = True
cpuTimer = CPUTime ( )

# . Paths.
outPath = os.path.join ( outPath, "arrays" )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Run the tests.
numberFailed = 0
for ( i, shape ) in enumerate ( _Shapes ):

    # . Allocation.
    A  = Array.WithShape ( shape )
    iA = A.iterator
    A.Set ( _AValue )
    if _Print: ArrayPrint ( A, title = "Array A with Shape {:s}".format ( repr ( list ( shape ) ) ), useLabelRange = True )
    if math.fabs ( iA.Sum ( ) - float ( A.size ) * _AValue ) > _Tolerance: numberFailed += 1

    # . Indexing (this is implicitly done by ArrayPrint anyway).
    for index in _Indexes[i]:
        value = A.__getitem__ ( index )
        if math.fabs ( value - _AValue ) > _Tolerance: numberFailed += 1
    
    # . Deep cloning.
    C  = DeepClone ( A )
    iC = C.iterator
    C.Set ( _CValue )
    if _Print: ArrayPrint ( C, title = "Array C (clone of A)", useLabelRange = True )
    if math.fabs ( iC.Sum ( ) - float ( C.size ) * _CValue ) > _Tolerance: numberFailed += 1

    # . Copying.
    A.CopyTo ( C )
    if _Print: ArrayPrint ( C, title = "Array C (after copy from A)", useLabelRange = True )
    if math.fabs ( iC.Sum ( ) - iA.Sum ( ) ) > _Tolerance: numberFailed += 1

    # . Adding.
    iC.Add ( iA, scale = -1.0 )
    if iC.AbsoluteMaximum ( ) > _Tolerance: numberFailed += 1

    # . Serialization.
    name     = "array_{:d}_{:d}_{:d}.yaml".format ( i, A.rank, A.size )
    yamlPath = os.path.join ( outPath, name )
    YAMLPickle ( yamlPath, A )
    newA     = YAMLUnpickle ( yamlPath )
    iNewA    = newA.iterator
    if _Print: ArrayPrint ( newA, title = "Array A (after unpickling)", useLabelRange = True )
    iNewA.Add ( iA, scale = -1.0 )
    if iNewA.AbsoluteMaximum ( ) > _Tolerance: numberFailed += 1

    # . Other operations.
    A.CopyTo ( C )
    iC.Power ( _Power )
    if _Print: ArrayPrint ( C, title = "Array A (after power)", itemFormat = "{:.1f}", useLabelRange = True )
    if math.fabs ( iC.Sum ( ) - float ( C.size ) * _APower ) > _Tolerance: numberFailed += 1

    # . Flattening and reshaping.
    F = Flatten ( A )
    if _Print: ArrayPrint ( F, title = "Array A (after flattening)", itemFormat = "{:.1f}", useLabelRange = True )
    if F is not A: # . Error if A = F.
        iF = F.iterator 
        iF.Add ( iA, scale = -1.0 )
        if iF.AbsoluteMaximum ( ) > _Tolerance: numberFailed += 1
    for newShape in _Reshapes[i]:
        R = Reshape ( A, newShape )
        if _Print: ArrayPrint ( R, title = "Array A (after reshaping)", itemFormat = "{:.1f}", useLabelRange = True )
        if R is not A: # . Error if A = F.
            iR = R.iterator
            iR.Add ( iA, scale = -1.0 )
            if iR.AbsoluteMaximum ( ) > _Tolerance: numberFailed += 1

    # . View operations.
    for ( v, ( index, canFlatten, reshapes ) ) in enumerate ( _Slices[i] ):

        # . Get the view.
        view  = A.__getitem__ ( index )
        iView = view.iterator
        value = _VSet[v]

        # . Setting.
        A.Set    ( _AValue )
        view.Set (   value )
        if _Print:
            ArrayPrint ( view, title = "Rank {:d} View (after set to {:.1f})".format ( view.rank, value ), useLabelRange = True )
            ArrayPrint ( A   , title = "Array A (after view set to {:.1f})".format   (            value ), useLabelRange = True )
        if math.fabs ( iA.Sum ( ) - float ( A.size - view.size ) * _AValue - float ( view.size ) * value ) > _Tolerance: numberFailed += 1

        # . Flattening.
        if view.rank > 1:
            try:
                vF = Flatten ( view )
                if _Print: ArrayPrint ( vF, title = "Flattened view", useLabelRange = True )
                if not canFlatten: numberFailed += 1
            except:
                if canFlatten: numberFailed += 1

        # . Reshaping.
        for ( newShape, canReshape ) in reshapes:
            try:
                vR = Reshape ( view, newShape )
                if _Print: ArrayPrint ( vR, title = "Reshaped view", useLabelRange = True )
                if not canReshape: numberFailed += 1
            except:
                if canReshape: numberFailed += 1

        # . Cloning.
        cD = DeepClone    ( view ) # . Independent of A and view.
        cS = ShallowClone ( view ) # . Dependent on A and equivalent to view.
        if _Print:
            ArrayPrint ( cD, title = "Deep Clone of Rank {:d} View".format    ( view.rank ), useLabelRange = True )
            ArrayPrint ( cS, title = "Shallow Clone of Rank {:d} View".format ( view.rank ), useLabelRange = True )
        if math.fabs ( cD.iterator.Sum ( ) - float ( cD.size ) * value ) > _Tolerance: numberFailed += 1
        if math.fabs ( cS.iterator.Sum ( ) - float ( cS.size ) * value ) > _Tolerance: numberFailed += 1

        cD.Set ( - 2.0 * _AValue )
        cS.Set ( -       _AValue ) # . Changes A and view as well.
        if _Print:
            ArrayPrint ( cD, title = "Deep Clone of Rank {:d} View (after set)".format    ( view.rank ), useLabelRange = True )
            ArrayPrint ( cS, title = "Shallow Clone of Rank {:d} View (after set)".format ( view.rank ), useLabelRange = True )
        if math.fabs ( iA.Sum ( ) - float ( A.size - 2 * cS.size ) * _AValue ) > _Tolerance: numberFailed += 1

        # . Serialization.
        name     = "view_{:d}_{:d}_{:d}_{:d}.yaml".format ( i, v, view.rank, view.size )
        yamlPath = os.path.join ( outPath, name )
        YAMLPickle ( yamlPath, view )
        newV     = YAMLUnpickle ( yamlPath ) # . newV independent of A (and view and cS) after unpickling.
        iNewV    = newV.iterator
        if _Print: ArrayPrint ( newV, title = "Rank {:d} View (after unpickling)".format ( newV.rank ), useLabelRange = True )
        iNewV.Add ( iView, scale = -1.0 )
        if iNewV.AbsoluteMaximum ( ) > _Tolerance: numberFailed += 1

# . Footer.
isOK = ( numberFailed <= 0 )
if isOK: logFile.Paragraph ( "All tests passed successfully." )
else:    logFile.Paragraph ( "{:d} of the tests failed."      )
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )

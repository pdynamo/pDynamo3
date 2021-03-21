"""Array printing."""

from  pCore         import Align         , \
                           DataType      , \
                           logFile       , \
                           LogFileActive
from .ArrayIterator import ArrayIterator

# . The condition function has arguments of ( value, indices ).

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Item formatting defaults.
_DefaultItemFormat  = "{!r}"
_DefaultItemsPerRow =    9
_DefaultItemWidth   =   10
_ItemFormat         = { DataType.Boolean : "{!r}" , DataType.Real : "{:18.8f}" , DataType.Integer : "{:12d}" }
_ItemsPerRow        = { DataType.Boolean :    9   , DataType.Real :     6      , DataType.Integer :     9    }
_ItemWidth          = { DataType.Boolean :   10   , DataType.Real :    19      , DataType.Integer :    10    }

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def ArrayPrint ( array, conditionFunction = None    ,
                        itemFormat        = None    ,
                        itemsPerRow       = None    ,
                        itemWidth         = None    ,
                        labels            = None    ,
                        log               = logFile ,
                        title             = None    ,
                        useLabelRange     = False   ):
    """Printing with an optional condition function and label range format."""
    ( isOK, itemFormat, itemsPerRow, itemWidth ) = _ProcessInputOptions ( array, itemFormat, itemsPerRow, itemWidth, log )
    itemsPerRow = min ( itemsPerRow, array.size )
    if isOK:
        # . Ensure that labels exist.
        ( labels, labelWidths ) = _CheckDimensionLabels ( array, labels )
        # . Get the iterator.
        iterator = ArrayIterator ( array, conditionFunction = conditionFunction, includeIndices = True )
        # . Find the table column format.
        if useLabelRange: columns = labelWidths + [ 2 ] + labelWidths + itemsPerRow * [ itemWidth ]
        else:             columns = itemsPerRow * ( labelWidths + [ itemWidth ] )
        # . Initialize the table.
        table = log.GetTable ( columns = columns )
        table.Start ( )
        if title is not None: table.Title ( title )
        # . Output the table.
        if useLabelRange:
            empty   = [ "" ] * ( array.rank + 1 )
            entries = []
            first   = None
            last    = empty
            for ( item, indices ) in iterator:
                iLabels = [ labels[i][j] for ( i, j ) in enumerate ( indices ) ]
                if first is None: first =            iLabels
                else:             last  = [ " -" ] + iLabels
                entries.append ( itemFormat.format ( item ) )
                if len ( entries ) == itemsPerRow:
                    for item in ( first + last + entries ): table.Entry ( item )
                    entries = []
                    first   = None
                    last    = empty
            if len ( entries ) > 0:
                for item in ( first + last + entries ): table.Entry ( item )
        else:
            for ( item, indices ) in iterator:
                for ( i, j ) in enumerate ( indices ):
                    table.Entry ( labels[i][j] )
                table.Entry ( itemFormat.format ( item ) )
        # . Finish up.
        table.Stop ( )

#-----------------------------------------------------------------------------------------------------------------------------------
def ArrayPrint2D ( array, columnLabels    = None    ,
                          columnSelection = None    ,
                          extraColumnData = None    ,
                          itemFormat      = None    ,
                          itemsPerRow     = None    ,
                          itemWidth       = None    ,
                          log             = logFile ,
                          rowLabels       = None    ,
                          rowSelection    = None    ,
                          title           = None    ):
    """Printing of 2-D arrays with optional extra column data and column and row selections.

    |extraColumnData| must be a list of 3-tuples that contain:
    (1) the tag to use for the data
    (2) the data to be output as strings
    (3) the align for the data strings
    """
    ( isOK, itemFormat, itemsPerRow, itemWidth ) = _ProcessInputOptions ( array, itemFormat, itemsPerRow, itemWidth, log )
    if isOK and ( array.rank == 2 ):
        # . Sort out the row and column selections.
        if columnSelection is None: columnSelection = range ( array.shape[1] )
        if rowSelection    is None: rowSelection    = range ( array.shape[0] )
        numberOfColumns = len ( columnSelection )
        itemsPerRow     = min ( itemsPerRow, numberOfColumns )
        # . Ensure that labels exist.
        ( labels, labelWidths ) = _CheckDimensionLabels ( array, [ rowLabels, columnLabels ] )
        # . Sort out the row labels.
        rowLabels  = labels     [0]
        indexWidth = labelWidths[0]
        # . Sort out the column labels.
        columnLabelSets = [ ( None, labels[1], "center" ) ]
        itemWidth       = max ( itemWidth, labelWidths[1] )
        if extraColumnData is not None:
            for ( tag, items, align ) in extraColumnData:
                indexWidth = max ( indexWidth, len ( tag ) )
                for item in items: itemWidth = max ( itemWidth, len ( item ) )
            columnLabelSets.extend ( extraColumnData )
        # . Tables.
        ( numberOfTables, lastItemsPerRow ) = divmod ( numberOfColumns, itemsPerRow )
        if lastItemsPerRow > 0: numberOfTables += 1
        for t in range ( numberOfTables ):
            if ( t == numberOfTables - 1 ) and ( lastItemsPerRow > 0 ): numberOfItems = lastItemsPerRow
            else:                                                       numberOfItems = itemsPerRow
            columnIndices = [ columnSelection[c] for c in range ( t * itemsPerRow, t * itemsPerRow + numberOfItems ) ]
            columns       = [ indexWidth ] + numberOfItems * [ itemWidth ]
            table         = log.GetTable ( columns = columns )
            table.Start ( )
            if title is not None: table.Title ( title )
            for ( tag, labels, align ) in columnLabelSets:
                table.Entry ( tag, align = Align.Left )
                for c in columnIndices: table.Entry ( labels[c], align = align )
            for r in rowSelection:
                table.Entry ( rowLabels[r], align = Align.Left )
                for c in columnIndices: table.Entry ( itemFormat.format ( array[r,c] ) )
            table.Stop ( )

#-----------------------------------------------------------------------------------------------------------------------------------
def _CheckDimensionLabels ( array, labels ):
    """Check the dimension labels."""
    # . Each dimension is validated separately.
    # . labels      - labels for each entry of each dimension (a list of lists).
    # . labelWidths - maximum width for the labels of each dimension.
    rank  = array.rank
    shape = array.shape
    # . Special 1-D case.
    if ( labels is not None ) and ( rank == 1 ) and ( len ( labels ) == shape[0] ):
        labels = [ labels ]
    else:
        isOK = [ False for i in range ( rank ) ]
        if ( labels is None ) or ( len ( labels ) != rank ):
            labels = [ [] for i in range ( rank ) ]
        else:
            for ( i, ( e, items ) ) in enumerate ( zip ( shape, labels ) ):
                isOK[i] = ( items is not None ) and ( len ( items ) == e )
        if not all ( isOK ):
            labelsD = [ repr ( i ) for i in range ( max ( shape ) ) ]
            for ( i, ( e, done ) ) in enumerate ( zip ( shape, isOK ) ):
                if not done: labels[i] = labelsD[0:e]
    # . Label widths.
    labelWidths = []
    for labelsD in labels:
        width = 0
        for label in labelsD: width = max ( width, len ( label ) )
        labelWidths.append ( width + 2 )
    return ( labels, labelWidths )

#-----------------------------------------------------------------------------------------------------------------------------------
def _ProcessInputOptions ( array, itemFormat, itemsPerRow, itemWidth, log ):
    """Process the input options."""
    isOK = ( array.size > 0 ) and LogFileActive ( log )
    if isOK:
        if itemFormat  is None: itemFormat  = _ItemFormat.get  ( array.block.dataType, _DefaultItemFormat  )
        if itemsPerRow is None: itemsPerRow = _ItemsPerRow.get ( array.block.dataType, _DefaultItemsPerRow )
        if itemWidth   is None: itemWidth   = _ItemWidth.get   ( array.block.dataType, _DefaultItemWidth   )
    return ( isOK, itemFormat, itemsPerRow, itemWidth )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

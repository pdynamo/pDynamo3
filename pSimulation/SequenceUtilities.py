"""Utilities for dealing with sequences assuming a PDB-like format."""

import string

from collections import defaultdict
from pBabel      import PDBModel          , \
                        PDBModelEntity
from pCore       import Clone             , \
                        logFile           , \
                        LogFileActive
from pMolecule   import Sequence          , \
                        SequenceComponent , \
                        SequenceEntity
from pScientific import PeriodicTable

# . Move eventually to sequence?

#===================================================================================================================================
# . Create a sequence consisting of one atom per component named after elemental type.
#===================================================================================================================================
def CreateElementSequence ( system, entityDescription = None, entityLabel = "E" ):
    """Create a sequence of a single entity using element symbols for atom names and components and one atom per component."""
    # . Initialization.
    minorSeparator = Sequence._attributable["fieldSeparator"]
    sequence       = Sequence ( )
    entity         = SequenceEntity.WithOptions ( label = entityLabel )
    sequence.AddChild ( entity )
    # . Process the atoms.
    for ( index, atom ) in enumerate ( system.atoms ):
        symbol    = PeriodicTable.Symbol ( atom.atomicNumber ).upper ( )
        component = SequenceComponent.WithOptions ( genericLabel = symbol, label = symbol + minorSeparator + repr ( index + 1 ) )
        entity.AddChild    ( component )
        component.AddChild ( atom      )
        atom.index = index
        atom.label = symbol
    # . Finish up.
    system.__dict__["_sequence"] = sequence

#===================================================================================================================================
# . Create a homogeneous isolate sequence consisting of a single entity.
#===================================================================================================================================
def CreateHomogeneousIsolateSequence ( system ,
                                       atomData              = [ ( 8, "O" ), ( 1, "H1" ), ( 1, "H2" ) ] ,
                                       componentGenericLabel = "HOH"   ,
                                       entityDescription     = "Water" ,
                                       entityLabel           = "W"     ,
                                       firstComponentNumber  = 1       ):
    """Create a sequence assuming a system composed of homogeneous isolates with water as the default."""
    # . Initialization.
    minorSeparator = Sequence._attributable["fieldSeparator"]
    sequence       = Sequence.WithDefaults ( )
    entity         = SequenceEntity.WithOptions ( label = entityLabel )
    sequence.AddChild ( entity )
    # . Loop over the atoms per isolate.
    atomIndex = 0
    for ( index, isolate ) in enumerate ( system.connectivity.isolateIndices ):
        if len ( isolate ) != len ( atomData ): raise ValueError ( "Incompatible atom data and isolate lengths." )
        component = SequenceComponent.WithOptions ( genericLabel = componentGenericLabel, label = componentGenericLabel + minorSeparator + repr ( index + 1 ) )
        entity.AddChild ( component )
        for ( ( atomicNumber, atomLabel ), i ) in zip ( atomData, isolate ):
            atom = system.atoms[i]
            if ( atomicNumber != atom.atomicNumber ): raise ValueError ( "Atomic number mismatch." )
            component.AddChild ( atom )
            atom.index = atomIndex
            atom.label = atomLabel
            atomIndex += 1
    # . Finish up.
    system.__dict__["_sequence"] = sequence

#===================================================================================================================================
# . Determine a unique entity label given a set of systems.
#===================================================================================================================================
def DetermineUniqueEntityLabel ( *arguments, **options ):
    """Determine a unique entity label."""
    # . Get existing labels.
    labels = set ( )
    for arg in arguments:
        sequence = getattr ( arg, "sequence", None )
        if sequence is not None:
            for entity in sequence.children:
                fields = sequence.ParseLabel ( entity.label, fields = 1 )
                labels.add ( fields[0].upper ( ) )
    # . Get free labels.
    freelabels = list ( set ( string.ascii_uppercase ).difference ( labels ) )
    if len ( freelabels ) <= 0: raise ValueError ( "There are no unique entity labels remaining." )
    # . Label specified and not taken so use it.
    label = options.get ( "label", None )
    if ( label is not None ) and ( label in freelabels ):
        unique = label
    # . Use the first free label.
    else:
        freelabels.sort ( )
        unique = freelabels.pop ( 0 )
    return unique

#===================================================================================================================================
# . Print the component frequencies or sequences of the entities in a PDB model or sequence.
#===================================================================================================================================
# . Private functions.
def _EntityFrequencies ( entity ):
    """Get entity component frequencies."""
    frequencies = defaultdict ( int )
    for component in entity.children:
        frequencies[component.genericLabel] += 1
    if ( entity.label is None ) or ( len ( entity.label ) <= 0 ): label = None
    else:                                                         label = entity.label
    results = []
    for key in sorted ( frequencies.keys ( ) ):
        results.extend ( [ key, "{:d}".format ( frequencies[key] ) ] )
    return ( label, results )

def _EntitySequence ( entity ):
    """Get entity sequence."""
    if ( entity.label is None ) or ( len ( entity.label ) <= 0 ): label = None
    else:                                                         label = entity.label
    return ( label, [ component.label for component in entity.children ] )

def PrintComponentData ( target, doFrequencies = True, log = logFile ):
    """Print the component frequencies or sequences for an entity or for the entities of a PDB model or sequence."""
    if LogFileActive ( log ):
        # . Options.
        if doFrequencies: dataTag = "Component Frequencies"
        else:             dataTag = "Component Sequence"
        if isinstance ( PDBModel, PDBModelEntity ): targetTag = "PDB Model"
        else:                                       targetTag = "Sequence"
        # . Gather data.
        data = []
        # . Models or sequences.
        if isinstance ( target, ( PDBModel, Sequence ) ):
            title = "{:s} {:s}".format ( targetTag, dataTag )
            data.append ( ( title, [] ) )
            for entity in target.children:
                if doFrequencies: ( label, items ) = _EntityFrequencies ( entity )
                else:             ( label, items ) = _EntitySequence    ( entity )
                if len ( items ) > 0:
                    if label is None: title = "Unnamed Entity"
                    else:             title = "Entity \"{:s}\"".format ( label )
                    data.append ( ( title, items ) )
        # . Entities.
        elif isinstance ( target, ( PDBModelEntity, SequenceEntity ) ):
            if doFrequencies: ( label, items ) = _EntityFrequencies ( target )
            else:             ( label, items ) = _EntitySequence    ( target )
            if len ( items ) > 0:
                if label is None:
                    title = "{:s} for Unnamed {:s} Entity".format ( dataTag, targetTag )
                else:
                    title = "{:s} for {:s} Entity \"{:s}\"".format ( dataTag, targetTag, label )
                data.append ( ( title, items ) )
        # . Argument error.
        else: raise TypeError ( "Invalid argument." )
        # . Title and item lengths and item numbers.
        iLength = 0
        iNumber = 0
        tLength = 0
        for ( title, items ) in data:
            tLength = max ( tLength, len ( title ) )
            iNumber = max ( iNumber, len ( items ) )
            if len ( items ) > 0: iLength = max ( iLength, max ( [ len ( item ) for item in items ] ) )
        iLength += 2
        tLength += 2
        # . Get columns.
        pageWidth = max ( tLength, min ( 80, iLength * iNumber ) ) # . 80 is default page width.
        nColumns  = ( pageWidth + iLength - 1 ) // iLength
        if doFrequencies and ( nColumns % 2 != 0 ): nColumns += 1 # . Must be even for frequencies (as pairs).
        # . Printing.
        table = log.GetTable ( columns = nColumns * [ iLength ] )
        table.Start ( )
        for ( r, ( title, items ) ) in enumerate ( data ):
            if r == 0: table.Title   ( title )
            else:      table.Section ( title )
            for item in items: table.Entry ( item )
        table.Stop ( )

#===================================================================================================================================
# . Renumber the components in the entities of a sequence.
#===================================================================================================================================
def RenumberEntityComponents ( system, entityLabels = None ):
    """Renumber the components in the entities of a sequence."""
    if system.sequence is not None:
        for entity in system.sequence.children:
            if ( entityLabels is None ) or ( entity.label in entityLabels ):
                for ( i, component ) in enumerate ( entity.children ):
                    fields = system.sequence.ParseLabel ( component.label )
                    if len ( fields ) < 1: fields.append ( "" )
                    fields[1] =  repr ( i + 1 )
                    component.label = system.sequence.MakeLabel ( *fields )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

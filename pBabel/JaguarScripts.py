"""Scripts for manipulating files from the Jaguar program."""

import math

from  pCore                  import Clone         , \
                                    logFile       , \
                                    LogFileActive
from  pScientific            import PeriodicTable
from  pScientific.Arrays     import Array         , \
                                    StorageType
from .JaguarInputFileReader  import JaguarInputFileReader
from .JaguarOutputFileReader import JaguarOutputFileReader

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The tolerance for bond order values.
_DEFAULTBONDORDERTOLERANCE = 0.1

# . The tolerance for charge deviations.
_DEFAULTCHARGETOLERANCE = 1.0e-3

# . The occupancy tolerance for including orbitals in the density matrix calculation.
_OCCUPANCYTOLERANCE = 1.0e-10

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def JaguarBondOrders ( inFile, outFile, bondOrdertolerance = _DEFAULTBONDORDERTOLERANCE, chargetolerance = _DEFAULTCHARGETOLERANCE, log = logFile, QCROSS = False ):
    """Calculate bond orders given Jaguar input and output files."""

    # . Parse the input file.
    inFile = JaguarInputFileReader.FromPath ( inFile )
    inFile.Parse ( )
    inFile.Summary ( log = logFile )

    # . Parse the output file.
    outFile = JaguarOutputFileReader.FromPath ( outFile )
    outFile.Parse   ( )
    outFile.Summary ( log = logFile )

    # . Get data from the files.
    # . Input.
    atomicNumbersi = getattr ( inFile,  "atomicNumbers", None )
    coordinates3   = getattr ( inFile,  "coordinates3",  None )
    orbitalSets    = getattr ( inFile,  "orbitalSets",   None )

    # . Output.
    atomicNumberso = getattr ( outFile, "atomicNumbers", None )
    nFunctions     = getattr ( outFile, "nFunctions",    None )
    overlap        = getattr ( outFile, "overlap",       None )

    # . Check for coordinates.
    if coordinates3 is None: raise ValueError ( "The coordinate data is missing from the input file." )

    # . Check the orbital data.
    QOK = ( orbitalSets is not None ) and ( len ( orbitalSets ) > 0 )
    if QOK:
        if   ( len ( orbitalSets ) == 1 ) and ( ""      in orbitalSets )                               : spinDensitiesRestricted = True
        elif ( len ( orbitalSets ) == 2 ) and ( "alpha" in orbitalSets ) and ( "beta" in orbitalSets ) : spinDensitiesRestricted = False
        else: QOK = False
    if not QOK: raise ValueError ( "Invalid orbital data on input file." )
    if spinDensitiesRestricted: nBasisi = inFile.orbitalSets[""]     [1]
    else:                       nBasisi = inFile.orbitalSets["alpha"][1]

    # . Check the overlap.
    if ( overlap is None ): raise ValueError ( "The overlap matrix is missing from the output file." )
    nBasiso = overlap.Dimension ( )

    # . Check the array giving the number of functions per atom.
#    print nFunctions, len ( nFunctions ), len ( atomicNumberso ), sum ( nFunctions ), nBasiso
    if ( nFunctions is None ) or ( len ( nFunctions ) != len ( atomicNumberso ) or ( sum ( nFunctions ) != nBasiso ) ): raise ValueError ( "Basis function data on the output file is missing or invalid." )

    # . Create the function index array.
    findices = []
    first = 0
    last  = 0
    for f in nFunctions:
        last  += f
        findices.append ( ( first, last ) )
        first += f

    # . Check for compatibility between the data.
    QOK = ( atomicNumbersi is not None ) and ( atomicNumberso is not None ) and ( len ( atomicNumbersi ) == len ( atomicNumberso ) ) and ( nBasisi == nBasiso )
    if QOK:
        for ( i, j ) in zip ( atomicNumbersi, atomicNumberso ):
            if i != j:
                QOK = False
                break
    if not QOK: raise ValueError ( "The systems on the input and output files are incompatible." )

    # . Set the keys for the calculation.
    if spinDensitiesRestricted: keys = [ "" ]
    else:               keys = [ "alpha", "beta" ]

    # . Get the densities multiplied by the overlap.
    ps = {}
    for key in keys:
        p      = _MakeDensity ( orbitalSets[key] )
        result = Array.WithExtents ( nBasisi, nBasisi )
        result.Set ( 999.0 )
        p.PostMultiply ( overlap, result )
        ps[key] = result
    # . Scale ps correctly for the spin-restricted case.
    if spinDensitiesRestricted: ps[""].Scale ( 2.0 )

    # . If cross terms are not required condense the ps matrices.
    if ( not QCROSS ) and ( not spinDensitiesRestricted ):
        tps  = ps.pop ( keys[0] )
        for key in keys[1:]: tps.Add ( ps[key] )
        ps   = { "": tps }
        keys = [ "" ]

    # . Get the bond-orders.
    bondOrders = {}
    for key1 in keys:
        for key2 in keys:
            bondOrders[ ( key1, key2 ) ] = _MakeBondOrders ( ps[key1], ps[key2], findices )

    # . Make the total bond-order if necessary.
    bokeys = list ( bondOrders.keys ( ) )
    if len ( bokeys ) > 1:
        bokeys.sort ( )
        tbo    = Clone ( bondOrders[bokeys[0]] )
        for key in bokeys[1:]: tbo.Add ( bondOrders[key] )
        tkey   = ( "", "" )
        bokeys.append ( tkey )
        bondOrders[tkey] = tbo

    # . Compute the electronic contribution to the Mulliken charges.
    qmulliken = Array.WithExtent ( len ( atomicNumbersi ) )
    qmulliken.Set ( 0.0 )
    for key in keys: _MakeElectronicMullikenCharges ( ps[key], findices, qmulliken )

    # . Determine atom valencies.
    free      = Array.WithExtent ( len ( atomicNumbersi ) ) ; free.Set      ( 0.0 )
    valencies = Array.WithExtent ( len ( atomicNumbersi ) ) ; valencies.Set ( 0.0 )
    tbo       = bondOrders[ ( "", "" ) ]
    for i in range ( len ( atomicNumbersi ) ):
        valencies[i] = ( - 2.0 * qmulliken[i] ) - tbo[i,i]
        totalbo = 0.0
        for j in range ( len ( atomicNumbersi ) ): totalbo += tbo[i,j]
        free[i] = ( - 2.0 * qmulliken[i] ) - totalbo

    # . Add in the core contributions to the Mulliken charges.
    for ( i, q ) in enumerate ( atomicNumbersi ): qmulliken[i] += float ( q )
    if outFile.QECP:
        for ( i, q ) in enumerate ( outFile.ecpElectrons ): qmulliken[i] -= float ( q )

    # . Output the results.
    if LogFileActive ( log ):

        # . Get the spin label.
        if spinDensitiesRestricted: spinlabel = "Spin Restricted"
        else:               spinlabel = "Spin Unrestricted"

        # . Create the atom names.
        atomnames = []
        for ( i, n ) in enumerate ( atomicNumbersi ):
            atomnames.append ( PeriodicTable.Symbol ( n, index = i + 1 ) )

        # . Atom data.
        columns = [ 10, 20, 20, 20 ]
        if not spinDensitiesRestricted: columns.append ( 20 )
        table   = log.GetTable ( columns = columns )
        table.Start ( )
        table.Title ( "Atom Data (" + spinlabel + ")" )
        table.Heading ( "Atom"            )
        table.Heading ( "Charge"          )
        table.Heading ( "Self Bond Order" )
        table.Heading ( "Valence"         )
        if not spinDensitiesRestricted: table.Heading ( "Free Valence" )
        for ( i, ni ) in enumerate ( atomnames ):
            table.Entry ( ni )
            table.Entry ( "{:.3f}".format ( qmulliken[i] ) )
            table.Entry ( "{:.3f}".format ( tbo[i,i]     ) )
            table.Entry ( "{:.3f}".format ( valencies[i] ) )
            if not spinDensitiesRestricted: table.Entry ( "{:.3f}".format ( free[i] ) )
        table.Stop ( )

        # . Bond orders.
        for key in bokeys:
            orders = bondOrders[key]
            table  = log.GetTable ( columns = [ 10, 10, 20, 20 ] )
            table.Start ( )
            if key == ( "", "" ): table.Title ( "Total Bond Orders" )
            else:                 table.Title ( key[0].title ( ) + "/" + key[1].title ( ) + " Bond Orders" )
            table.Heading ( "Atom 1"   )
            table.Heading ( "Atom 2"   )
            table.Heading ( "Order"    )
            table.Heading ( "Distance" )
            for ( i, ni ) in enumerate ( atomnames ):
                for ( j, nj ) in enumerate ( atomnames[0:i] ):
                    b = orders[i,j]
                    if math.fabs ( b ) > bondOrdertolerance:
                        table.Entry ( ni )
                        table.Entry ( nj )
                        table.Entry ( "{:.3f}".format ( b ) )
                        table.Entry ( "{:.3f}".format ( coordinates3.Distance ( i, j ) ) )
            table.Stop ( )

        # . Checks on the calculation.
        # . Free valence.
        if spinDensitiesRestricted:
            deviation = free.AbsoluteMaximum ( )
            if deviation > chargetolerance: log.Paragraph ( "Warning: the largest deviation between the free valence values and zero is {:.3f}.".format ( deviation ) )

        # . Total charge.
        deviation = math.fabs ( qmulliken.Sum ( ) - float ( inFile.charge ) )
        if deviation > chargetolerance: log.Paragraph ( "Warning: the total charge deviates from the formal charge by {:.3f}.".format ( deviation ) )

        # . Check for a valid reference set of Mulliken charges.
        qreference = getattr ( outFile, "qmulliken", None )
        if ( qreference is not None ) and ( len ( qreference ) == len ( atomicNumbersi ) ):
            qmulliken.Add ( qreference, scale = -1.0 )
            deviation = qmulliken.AbsoluteMaximum ( )
            if deviation > chargetolerance: log.Paragraph ( "Warning: the largest deviation between calculated and read Mulliken charges is {:.3f}.".format ( deviation ) )

    # . Finish up.
    return bondOrders

#===================================================================================================================================
# . Private functions.
# . These functions could be made more efficient.
#===================================================================================================================================
def _MakeBondOrders ( ps1, ps2, findices ):
    """Make bond orders."""
    natoms = len ( findices )
    orders = Array.WithExtent ( natoms, storageType = StorageType.Symmetric )
    orders.Set ( 0.0 )
    for ( i, ( ifirst, ilast ) ) in enumerate ( findices ):
        for ( j, ( jfirst, jlast ) ) in enumerate ( findices[0:i+1] ):
            b = 0.0
            for m in range ( ifirst, ilast ):
                for n in range ( jfirst, jlast ):
                    b += ps1[m,n] * ps2[n,m]
            orders[i,j] = b
    return orders

def _MakeDensity ( orbitalset ):
    """Make a density."""
    ( norbitals, nBasis, energies, occupancies, vectors ) = orbitalset
    p = Array.WithExtent ( nBasis, storageType = StorageType.Symmetric )
    p.Set ( 0.0 )
    for i in range ( norbitals ):
        o = occupancies[i]
        if math.fabs ( o ) > _OCCUPANCYTOLERANCE:
            for m in range ( nBasis ):
                for n in range ( m + 1 ):
                    p[m,n] += o * vectors[m,i] * vectors[n,i]
    return p

def _MakeElectronicMullikenCharges ( ps, findices, q ):
    """Determine an electronic contribution to the Mulliken charges."""
    for ( i, ( first, last ) ) in enumerate ( findices ):
        for j in range ( first, last ): q[i] -= ps[j,j]

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

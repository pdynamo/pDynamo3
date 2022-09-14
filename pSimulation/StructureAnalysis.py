"""Helper functions for analyzing structures."""

import math

from  pCore                 import Align , \
                                   logFile
from  pScientific.Arrays    import Array
from  pScientific.Geometry3 import SelfPairList_FromCoordinates3
from .CrystalUtilities      import CrystalExpandToP1 , \
                                   CrystalGetImageBondPairs

#===================================================================================================================================
# . Identify possible bonds in a structure.
#===================================================================================================================================
def IdentifyPossibleBonds ( system, atomRadii = None, bondSafetyFactor = 0.45, closeContactDistance = 0.80, excludeMMBonds = True, log = logFile ):
    """Identify possible bonds in a structure using geometrical criteria only.

    By default, only bonds that do not appear in the MM model are identified.
    Both active and inactive MM terms are used in the check.

    Very short bonds (close contacts) are flagged.
    """

    # . Initialization.
    closeContacts = None
    imageResults  = None
    possibleBonds = None

    # . Check for coordinates.
    coordinates3 = system.coordinates3
    if coordinates3 is not None:

        # . If there are undefined coordinates do nothing.
        if system.coordinates3.numberUndefined <= 0:

            # . Get some options.
            if atomRadii is None: atomRadii = Array.FromIterable ( [ atom.covalentRadius for atom in system.atoms ] )

            # . Get the MM bonds.
            mmbonds = None
            if excludeMMBonds and ( system.mmState is not None ):
                mmbonds = system.mmState.Get12Exclusions ( )
                if ( mmbonds is not None ) and ( len ( mmbonds ) <= 0 ): mmbonds = None

            # . Get close contacts (no matter whether MM or not) and possible bonds (excluding MM).
            closeContacts = SelfPairList_FromCoordinates3 ( coordinates3, safety = closeContactDistance )
            possibleBonds = SelfPairList_FromCoordinates3 ( coordinates3, exclusions = mmbonds, radii = atomRadii, safety = bondSafetyFactor )

            # . Output.
            # . Close contacts.
            if ( len ( closeContacts ) > 0 ):
                table = log.GetTable ( columns = 4 * [ 5, 5, 10, 5 ] )
                table.Start ( )
                table.Title ( "Close Contacts" )
                for ( i, j ) in closeContacts:
                    d = coordinates3.Distance ( i, j )
                    table.Entry ( "{:d}".format ( i ) )
                    table.Entry ( "{:d}".format ( j ) )
                    table.Entry ( "{:6.3f}".format ( d ) )
                    table.Entry ( " *" )
                table.Stop ( )

            # . Possible bonds.
            if ( len ( possibleBonds ) > 0 ):
                table = log.GetTable ( columns = 4 * [ 5, 5, 10, 5 ] )
                table.Start ( )
                table.Title ( "Possible Bonds (Excluding Close Contacts)" )
                nbonds = 0
                for ( i, j ) in possibleBonds:
                    d = coordinates3.Distance ( i, j )
                    if ( d > closeContactDistance ):
                        nbonds += 1
                        table.Entry ( "{:d}".format ( i ) )
                        table.Entry ( "{:d}".format ( j ) )
                        table.Entry ( "{:6.3f}".format ( d ) )
                        table.Entry ( "" )
                if nbonds == 0: table.Entry ( "There are no possible bonds which are not also close contacts.", align = Align.Center, columnSpan = 16 )
                table.Stop ( )

            # . The system has translational symmetry.
            if ( system.symmetry is not None ) and hasattr ( system.symmetry, "transformations" ):

                # . Get the pairlists.
                imageResults = CrystalGetImageBondPairs ( system, radii = atomRadii, safety = bondSafetyFactor )

                # . Output.
                if ( imageResults is not None ) and ( len ( imageResults ) > 0 ):
                    for ( t, a, b, c, pairs ) in imageResults:
                        table = log.GetTable ( columns = 4 * [ 5, 5, 10, 5 ] )
                        table.Start ( )
                        table.Title ( "Possible Bonds With Image {:d} {:d} {:d} {:d}".format ( t, a, b, c ) )
                        for ( i, j, d ) in pairs:
                            table.Entry ( "{:d}".format ( i ) )
                            table.Entry ( "{:d}".format ( j ) )
                            table.Entry ( "{:6.3f}".format ( d ) )
                            if ( d <= closeContactDistance ): table.Entry ( " *" )
                            else:                             table.Entry ( "" )
                        table.Stop ( )

    # . Finish up.
    return ( closeContacts, possibleBonds, imageResults )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

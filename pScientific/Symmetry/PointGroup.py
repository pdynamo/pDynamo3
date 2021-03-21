"""Point group classes and functions."""

import glob, math, os, os.path

from itertools       import combinations                  , \
                            combinations_with_replacement
from   pCore         import AttributableObject            , \
                            YAMLUnpickle
from  .SymmetryError import SymmetryError
from ..Arrays        import Array

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_CharacterMatchTolerance = 0.1 # . Dimensionless.
_UnassignedIRLabel       = "?"

#===================================================================================================================================
# . Utility class.
#===================================================================================================================================
class _IRIterator:
    """Iterator over combinations of IRs of a specific size and degeneracy."""

    def __init__ ( self, pointGroup, nIRs, degeneracy ):
        """Constructor."""
        self.combinations = self._MakeCombinations ( pointGroup.degeneracies, pointGroup.irreducibleRepresentations, nIRs, degeneracy )
        self.count        = 0

    def __iter__ ( self ): return self

    def __next__ ( self ):
        if self.count < len ( self.combinations ):
            items       = self.combinations[self.count]
            self.count += 1
            return items
        else:
            raise StopIteration

    def _MakeCombinations ( self, degeneracies, labels, nIRs, degeneracy ):
        """Make all combinations of indices - repeats allowed."""
        # . Each IR has a minimum degeneracy of 1, so the maximum IR possible is when (nIRs-1) IRs have a degeneracy of 1
        # . and the remaining IR has a degeneracy of ( degeneracy - ( nIRs - 1 ) ).
        upper   = degeneracy + 1 - nIRs
        allowed = [ r for r in range ( len ( degeneracies ) ) if degeneracies[r] <= upper ]
        if nIRs == 1:
            combinations = [ ( ( r, ), labels[r] ) for r in allowed ]
        else:
            combinations = []
            for indices in combinations_with_replacement ( allowed, nIRs ):
                if sum ( [ degeneracies[r] for r in indices ] ) == degeneracy:
                    combinations.append ( ( indices, "/".join ( sorted ( set ( [ labels[r] for r in indices ] ) ) ) ) )
        return combinations

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PointGroup ( AttributableObject ):
    """Point group class."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "characterMatchTolerance"     : _CharacterMatchTolerance ,
                             "characterSymmetryOperations" : None                     ,
                             "characterTable"              : None                     ,
                             "cInfinityRotations"          : None                     ,
                             "degeneracies"                : None                     ,
                             "irreducibleRepresentations"  : None                     ,
                             "label"                       : ""                       ,
                             "maximumDegeneracy"           : 0                        ,
                             "principalAxisOperation"      : None                     ,
                             "symmetryOperations"          : None                     } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( PointGroup, self )._CheckOptions ( )
        # . Find the maximum degeneracy.
        if self.characterTable is not None:
            e                      = self.characterSymmetryOperations.index ( "E" )
            self.degeneracies      = [ int ( round ( self.characterTable[c,e] ) ) for c in range ( len ( self.irreducibleRepresentations ) ) ]
            self.maximumDegeneracy = max ( self.degeneracies )
        # . Find the principal axis.
        if self.symmetryOperations is not None:
            if "CInfinity" in self.symmetryOperations: self.principalAxisOperation = "CInfinity"
            else:
                maximumN = 0
                for key in self.symmetryOperations.keys ( ):
                   if key.startswith ( "C" ): maximumN = max ( maximumN, int ( key[1:] ) )
                if maximumN > 0:
                    key = "C{:d}".format ( maximumN )
                    if self.symmetryOperations[key] == 1: self.principalAxisOperation = key

    def _IRIterator ( self, nIRs, degeneracy ):
        """An iterator over IRs."""
        return _IRIterator ( self, nIRs, degeneracy )

    def IdentifyIrreducibleRepresentationsX ( self, characters, itemValues, degeneracyTolerance, maximumIRs = 1 ):
        """Identify irreducible representations of items from their calculated characters."""
        # . Accidental degeneracy is handled by allowing combinations of IRs.
        # . Basic initialization.
        numberClasses  = len ( self.characterSymmetryOperations )
        numberItems    = characters.rows
        labels         = [ _UnassignedIRLabel for i in range ( numberItems ) ]
        irCharacters   = Array.WithExtent ( numberClasses )
        itemCharacters = Array.WithExtent ( numberClasses )
        # . Determine groups of possibly degenerate items.
        oldV   = itemValues[0]
        groups = [ [ 0 ] ]
        for i in range ( 1, numberItems ):
            v = itemValues[i]
            if math.fabs ( v - oldV ) <= degeneracyTolerance:
                groups[-1].append ( i )
            else:
                groups.append ( [ i ] )
            oldV = v
        # . Loop over groups of items.
        for group in groups:
            indices = list ( group )
            while len ( indices ) > 0:
                done = False
                # . The number of IRs to consider in the match in cases of accidental degeneracy.
                for nIRs in range ( 1, maximumIRs + 1 ):
                    # . Different combinations of indices starting with the largest.
                    for degeneracy in range ( min ( len ( indices ), nIRs * self.maximumDegeneracy ), ( nIRs-1 ), -1 ):
                        for combination in combinations ( indices, degeneracy ):
                            itemCharacters.Set ( 0.0 )
                            for itemIndex in combination:
                                itemCharacters.Add ( characters[itemIndex,:] )
                            # . All possible IR combinations with the same degeneracy as the combination.
                            for ( iRs, label ) in self._IRIterator ( nIRs, degeneracy ):
                                irCharacters.Set ( 0.0 )
                                for r in iRs:
                                    irCharacters.Add ( self.characterTable[r,:] )
                                matched = True
                                for c in range ( numberClasses ):
                                    if math.fabs ( irCharacters[c] - itemCharacters[c] ) > self.characterMatchTolerance:
                                        matched = False
                                        break
                                if matched:
                                    for itemIndex in combination:
                                        indices.remove ( itemIndex )
                                        labels[itemIndex] = label
                                    done = True
                                    break
                            if done: break
                        if done: break
                    if done: break
                # . No match found.
                if not done:
                    break
        #unassigned = [ i for i in range ( characters.rows ) if labels[i] == _UnassignedIRLabel ]
        # . Finish up.
        return labels

    def IdentifyIrreducibleRepresentations ( self, characters, itemValues, degeneracyTolerance, maximumIRs = 1 ):
        """Identify irreducible representations of items from their calculated characters."""
        # . Accidental degeneracy is handled by allowing combinations of IRs.
        # . Basic initialization.
        numberClasses  = len ( self.characterSymmetryOperations )
        numberItems    = characters.rows
        degeneracies   = self.degeneracies
        irLabels       = self.irreducibleRepresentations
        labels         = [ _UnassignedIRLabel for i in range ( numberItems ) ]
        irCharacters   = Array.WithExtent ( numberClasses )
        itemCharacters = Array.WithExtent ( numberClasses )
        # . Determine groups of possibly degenerate items.
        oldV   = itemValues[0]
        groups = [ [ 0 ] ]
        for i in range ( 1, numberItems ):
            v = itemValues[i]
            if math.fabs ( v - oldV ) <= degeneracyTolerance:
                groups[-1].append ( i )
            else:
                groups.append ( [ i ] )
            oldV = v
        # . Loop over groups of items.
        for group in groups:
            indices = list ( group )
            while len ( indices ) > 0:
                done = False
                # . The number of IRs to consider in the match in cases of accidental degeneracy.
                for nIRs in range ( 1, maximumIRs + 1 ):
                    # . Different combinations of indices starting with the largest.
                    for degeneracy in range ( min ( len ( indices ), nIRs * self.maximumDegeneracy ), ( nIRs-1 ), -1 ):
                        # . Find allowed IRs.
                        # . Each IR has a minimum degeneracy of 1, so the maximum IR possible is when (nIRs-1) IRs have a degeneracy of 1
                        # . and the remaining IR has a degeneracy of ( degeneracy - ( nIRs - 1 ) ).
                        allowed = [ r for r in range ( len ( degeneracies ) ) if degeneracies[r] <= ( degeneracy + 1 - nIRs ) ]
                        # . Loop over combinations of indices.
                        for combination in combinations ( indices, degeneracy ):
                            itemCharacters.Set ( 0.0 )
                            for itemIndex in combination:
                                itemCharacters.Add ( characters[itemIndex,:] )
                            # . All possible IR combinations with the same degeneracy as the combination.
                            for iRs in combinations_with_replacement ( allowed, nIRs ):
                                if sum ( [ degeneracies[r] for r in iRs ] ) == degeneracy:
                                    irCharacters.Set ( 0.0 )
                                    for r in iRs:
                                        irCharacters.Add ( self.characterTable[r,:] )
                                    matched = True
                                    for c in range ( numberClasses ):
                                        if math.fabs ( irCharacters[c] - itemCharacters[c] ) > self.characterMatchTolerance:
                                            matched = False
                                            break
                                    if matched:
                                        label = "/".join ( sorted ( set ( [ irLabels[r] for r in iRs ] ) ) )
                                        for itemIndex in combination:
                                            indices.remove ( itemIndex )
                                            labels[itemIndex] = label
                                        done = True
                                        break
                            if done: break
                        if done: break
                    if done: break
                # . No match found.
                if not done:
                    break
        #unassigned = [ i for i in range ( characters.rows ) if labels[i] == _UnassignedIRLabel ]
        # . Finish up.
        return labels

    def OperationKey ( self ):
        """Return the operation key."""
        if self.__dict__.get ( "_operationKey", None ) is None:
            keys = list ( self.symmetryOperations.keys ( ) )
            keys.sort ( )
            items = []
            for key in keys:
                n = self.symmetryOperations[key]
                if n == 1: items.append ( key )
                else:      items.append ( repr ( n ) + "*" + key )
            self._operationKey = " ".join ( items )
        return self._operationKey

    @classmethod
    def FromYAML ( selfClass, path ):
        """Constructor from a YAML file."""
        from math import cos, pi, sin
        # . Get basic data.
        data                        = YAMLUnpickle ( path )
        label                       = data["Label"]
        characterSymmetryOperations = data["Character Symmetry Operations"]
        cInfinityRotations          = data.get ( "CInfinity Rotations", [] )
        rows                        = data["Character Table"    ]
        symmetryOperations          = data["Symmetry Operations"]
        # . Process character table data.
        irreducibleRepresentations = list ( rows.keys ( ) )
        irreducibleRepresentations.sort ( )
        characterTable = Array.WithExtents ( len ( irreducibleRepresentations ), len ( characterSymmetryOperations ) )
        characterTable.Set ( 0.0 )
        for ( i, iR ) in enumerate ( irreducibleRepresentations ):
            row = rows[iR]
            for ( j, entry ) in enumerate ( row ):
                if isinstance ( entry, str ): value = eval ( entry )
                else:                         value = entry
                characterTable[i,j] = value
        # . Create the point group.
        self = selfClass.WithOptions ( label                       = label                       ,
                                       characterSymmetryOperations = characterSymmetryOperations ,
                                       characterTable              = characterTable              ,
                                       cInfinityRotations          = cInfinityRotations          ,
                                       irreducibleRepresentations  = irreducibleRepresentations  ,
                                       symmetryOperations          = symmetryOperations          )
        # . Finish up.
        return self

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def PointGroups_FromYAML ( inPath = None ):
    """Read point groups from a directory containing YAML files.

    Return as a dictionary with symmetry operation keys.
    """
    # . Get the path.
    if inPath is None:
        try:    inPath = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "pointGroups" )
        except: raise SymmetryError ( "Unable to find point group directory." )
    # . Initialization.
    pointGroups = {}
    # . Loop over all files.
    paths = glob.glob ( os.path.join ( inPath, "*.yaml" ) )
    for path in paths:
        if path.find ( "ReadMe" ) < 0:
            group = PointGroup.FromYAML ( path )
            key   = group.OperationKey ( )
            pointGroups[key] = group
    # . Finish up.
    return pointGroups

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

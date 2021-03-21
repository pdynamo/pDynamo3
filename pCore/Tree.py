"""Basic tree classes."""

from .AttributableObject import AttributableObject

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Label data.
_DefaultLabel = ""

# . Separators for fields and labels.
_FieldSeparator = "."
_LabelSeparator = ":"

# . The wildcard character.
_WildCard = "*"

#===================================================================================================================================
# . Base class for tree nodes.
#===================================================================================================================================
class TreeNode ( AttributableObject ):
    """Base class for tree nodes."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "label"  : None ,
                             "parent" : None ,
                             "_path"  : None ,
                             "_root"  : None } )

    def ClearPathInformation ( self ):
        """Clear path information."""
        self._path = None
        self._root = None

    def ReturnAncestorAtGivenLevel ( self, level ):
        """Return the ancestor item at a given level."""
        item = None
        if level >= 0:
            parent = self
            for i in range ( level ):
                parent = parent.parent
                if parent is None: break
            item = parent
        return item

    def ReturnFirstAncestorOfGivenClass ( self, itemClass ):
        """Return the first ancestor of a given class."""
        item   = None
        parent = self
        while parent is not None:
            if isinstance ( parent, itemClass ):
                item = parent
                break
            parent = parent.parent
        return item

    @property
    def path ( self ):
        """Find the path for the node."""
        if self._path is None:
            if self.root is self:
                self._path = self.label
            else:
                labels = []
                root   = self
                while root.parent is not None:
                    labels.append ( root.label )
                    root = root.parent
                labels.reverse ( )
                self._path = root.MakePath ( *labels )
        return self._path

    @property
    def root ( self ):
        """Find the root for the node."""
        if self._root is None:
            root = self
            while root.parent is not None:
                root = root.parent
            self._root = root
        return self._root

#===================================================================================================================================
# . Terminal or leaf nodes of a tree.
#===================================================================================================================================
class TreeLeafNode ( TreeNode ):
    """The leaf node of a tree."""
    pass

#===================================================================================================================================
# . Branch nodes of a tree.
#===================================================================================================================================
class TreeBranchNode ( TreeNode ):
    """The branch node of a tree."""

    _attributable = dict ( TreeNode._attributable )
    _attributable.update ( { "children"   : list ,
                             "childIndex" : dict } )

    def AddChild ( self, child, toFollow = None ):
        """Add a child to the node."""
        isOK = True
        # . Check if correct type?
        isOK = ( child.label is not None ) and ( child.label not in self.childIndex )
        if isOK:
            child.ClearPathInformation ( )
            child.parent = self
            self.childIndex[child.label] = child
            if toFollow is None: other = None
            else:                other = self.childIndex.get ( toFollow, None )
            if other is None:
                self.children.append ( child )
            else:
                index = self.children.index ( other ) + 1
                self.children.insert ( index, atom )
        return isOK

    def GetDescendantFromLabels ( self, *labels ):
        """Return a descendant given a sequence of labels."""
        descendant = self
        item       = None
        for label in labels:
            try:    item = descendant.childIndex.get ( label, None )
            except: item = None
            if item is None: break
            else: descendant = item
        return item

    def GetDescendantFromPath ( self, path ):
        """Return a descendant given a path."""
        labels = self.ParsePath ( path )
        item   = self.GetDescendantFromLabels ( *labels )
        return item

    def GetDescendantPaths ( self, level = -1 ):
        """Return a list of the paths of all descendants.

        Level < 0   return all leaf node paths
              = 0   return the path of self
              > 0   return all paths of descendants of the given level
        """
        paths = []
        if level == 0:
            paths.append ( self.path )
        else:
            if level > 0: level -= 1
            for child in self.children:
                if isinstance ( child, TreeLeafNode ):
                    if level <= 0: paths.append ( child.path )
                else:
                    paths.extend ( child.GetDescendantPaths ( level = level ) )
        return paths

    def GetDescendantsFromLabels ( self, *labels ):
        """Return the sequence of descendants corresponding to a sequence of labels."""
        descendant = self
        item       = None
        items      = []
        for label in labels:
            if descendant is None:
                items.append ( None )
            else:
                try:    item = descendant.childIndex.get ( label, None )
                except: item = None
                items.append ( item )
                descendant = item
        return items

#===================================================================================================================================
# . Root node of a tree.
#===================================================================================================================================
class TreeRootNode ( TreeBranchNode ):
    """The root node of a tree."""

    _attributable = dict ( TreeBranchNode._attributable )
    _attributable.update ( { "defaultLabel"   : _DefaultLabel   ,
                             "fieldSeparator" : _FieldSeparator ,
                             "labelSeparator" : _LabelSeparator ,
                             "wildCard"       : _WildCard       } )

    # . Should always strip whitespace here when making or parsing.

    def MakeLabel ( self, *fields ):
        """Make a label from fields."""
        n = len ( fields )
        for field in reversed ( fields ):
            if field == self.defaultLabel: n -= 1
            else:                          break
        return self.fieldSeparator.join ( fields[0:n] )

    def MakePath ( self, *labels ):
        """Make a path from labels."""
        return self.labelSeparator.join ( labels )

    def ParseLabel ( self, label, fields = None ):
        """Parse a label into a given number of fields."""
        tokens = label.split ( self.fieldSeparator )
        n = len ( tokens )
        for token in reversed ( tokens ):
            if token == self.defaultLabel: n -= 1
            else:                          break
        tokens = tokens[0:n]
        if fields is not None:
            n = len ( tokens )
            if   fields > n: tokens.extend ( ( fields - n ) * [ self.defaultLabel ] )
            elif fields < n: tokens = tokens[0:max(fields,0)]
        return tokens

    def ParsePath ( self, path, labels = None ):
        """Parse a path."""
        tokens = path.split ( self.labelSeparator )
        n = len ( tokens )
        for token in reversed ( tokens ):
            if token == self.defaultLabel: n -= 1
            else:                          break
        tokens = tokens[0:n]
        if labels is not None:
            n = len ( tokens )
            if   labels > n: tokens.extend ( ( labels - n ) * [ self.defaultLabel ] )
            elif labels < n: tokens = tokens[0:max(labels,0)]
        return tokens

    def ParsePathToFields ( self, path, fields = None ):
        """Parse a path into a given number of fields per label."""
        labels = self.ParsePath ( path )
        if fields is not None:
            tokens = []
            for ( n, label ) in zip ( fields, labels ):
                tokens.extend ( self.ParseLabel ( label, fields = n ) )
            if len ( fields ) > len ( labels ):
                for n in fields[len ( labels ):-1]:
                    tokens.extend ( n * [ self.defaultLabel ] )
        else: tokens = labels
        return tokens

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass

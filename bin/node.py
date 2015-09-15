#
# Copyright (C) 2015 INRA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Frederic Escudie - Plateforme bioinformatique Toulouse'
__copyright__ = 'Copyright (C) 2015 INRA'
__license__ = 'GNU General Public License'
__version__ = '0.1.0'
__email__ = 'frogs@toulouse.inra.fr'
__status__ = 'dev'


import json


class Node:
    """
    @summary: A node is an element of a tree structure. It is linked with its parent and its children and it is described by several metadata.
    """
    def __init__(self, name, parent_node=None, children_nodes=None, metadata=None):
        """
        @param name: [str] The node name.
        @param parent_node: [Node] The parent node.
        @param children_nodes: [list] List of children nodes.
        @param metadata: [dict] The metadata.
        """
        self.children = dict()
        if children_nodes is not None:
            for current_child in children_nodes:
                self.add_child( current_child )
        self.parent = parent_node
        self.name = name
        self.metadata = metadata if metadata is not None else dict()

    def __str__(self):
        """
        @summary: Returns a string representation of the node.
        @return: [str] The representation of the node.
        """
        node_str = self.name
        node_str += "\n\tParent=" + (self.parent.name if self.parent is not None else "-")
        node_str += "\n\tChilds=" + ", ".join(self.children.keys())
        metadata = list()
        for key in sorted(self.metadata.keys()):
            metadata.append(key + ":" + str(self.metadata[key]))
        node_str += "\n\tMetadata=" + ", ".join(metadata)
        return node_str

    def has_child(self, name=None):
        """
        @summary: Returns true if the node has the specified child (name is set) or at least one child (name is None).
        @param name: [str] the name of the searched child.
        @return: [bool]
        """
        if name is None:
            return len(self.children) > 0
        else:
            return self.children.has_key( name )

    def get_child(self, name):
        """
        @summary: Returns the specified child node.
        @param name: [str] the name of the searched child.
        @return: [Node] The child.
        """
        if not self.has_child( name ):
            raise Exception( self.name + " doesn't have child named '" + name + "'." )
        return self.children[name]

    def get_children(self):
        """
        @summary: Returns the children of the node.
        @return: [list] The list of children nodes.
        """
        return self.children.values()

    def get_parent(self):
        """
        @summary: Returns the parent node.
        @return: [Node] The parent node.
        """
        return self.parent

    def get_ancestors(self):
        """
        @summary: Returns the ancestors of the node.
        @return: [list] The list ancestors nodes. The nodes are in parent to child order.
        """
        ancestors = list()
        if self.parent is not None:
            ancestors.extend( self.parent.get_ancestors() )
            ancestors.extend( [self.parent] )
        return ancestors

    def get_depth(self):
        """
        @summary: Returns the depth of the node (= the branch length/ = the number of ancestors nodes).
        @return: [int] The depth of the node. The depth for the root element is 0.
        """
        if self.parent is None:
            return 0
        else:
            return( self.parent.get_depth() + 1 )

    def add_child(self, child):
        """
        @summary: Adds node as child.
        @param child: [Node] the node added.
        """
        if self.children.has_key( child.name ):
            raise Exception( "Duplicated child name '" + child.name + "' in node '" + self.name + "'." )
        child.parent = self
        self.children[child.name] = child

    def to_newick(self, distance_tag=None):
        """
        @summary: Returns the representation of the tree rooted by the node in newick format.
        @param distance_tag: [str] The metadata tag for the node distance (default: 'dist'). The distance is not necessary to use this method.
        @returns: [str] the newick representation of the tree.
        """
        if distance_tag is None: distance_tag = "dist"
        if not self.has_child():
            if self.metadata.has_key( distance_tag ):
                return '"' + self.name + '":' + str(self.metadata[distance_tag])
            else:
                return '"' + self.name + '"'
        else:
            children_newick = list()
            for child_name in self.children:
                child = self.children[child_name]
                children_newick.append( child.to_newick(distance_tag) )
            return '(' + ','.join(children_newick) + ')"' + self.name + '"'

    def to_extended_newick(self):
        """
        @summary: Returns the representation of the tree rooted by the node in extended newick format. In extended newick the distance tag is replaced by a json represetation of the metadata.
        @returns: [str] the extended newick representation of the tree.
        """
        if not self.has_child():
            if len(self.metadata.keys()) != 0:
                return '"' + self.name + '":' + json.dumps(self.metadata)
            else:
                return '"' + self.name + '"'
        else:
            children_newick = list()
            for child_name in self.children:
                child = self.children[child_name]
                children_newick.append( child.to_extended_newick() )
            if len(self.metadata.keys()) != 0:
                return '(' + ','.join(children_newick) + ')"' + self.name + '":' + json.dumps(self.metadata)
            else:
                return '(' + ','.join(children_newick) + ')"' + self.name + '"'
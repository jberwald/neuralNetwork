import random
import networkx as NX
from numpy import cumsum
try:
    from jb.neural_network.core import graph_object as go
except:
    # probably on the mac
    from neural_network.core import graph_object as go

def make_subgraphs( num_clusters, edim=5, idim=5 ):
    """
    Create a list of correctly-labeled graphs to form into a larger cluster.
    """
    subclust = go.network_of_clusters( num_clusters, edim, idim )
    return go.relabel_nodes( subclust )

def make_tree( graph_dict, num_ ):
    """
    From a dict of graphs, create a self-similar/binary tree by
    connecting subgraphs in appropriate order.

    key |--> ( successor nodes )

    Eg., { G0 : ( G1, G2 ),
             G1 : ( G3, G4 ),
             G2 : ( G5, G6 ) }

    We assume a binary tree structure.
    """
    Tree = NX.MultiDiGraph()

    print "levels", len( graph_dict )

    # keys in dict are parents at each level of the tree
    for parent in graph_dict:
        # grab children
        children = graph_dict[ parent ]
        # add all required edges to Tree from parent and children
        Tree.add_edges_from( parent.edges() )
        for g in children:
            Tree.add_edges_from( g.edges() )
        # add edges between children/parent
        for i in num_edges:
            children2parent( Tree, children, parent )
            parent2children( Tree, parent, children )
    return Tree

def choose_node( G ):
    """
    Random choice of node.
    """
    return random.choice( G.graph.nodes() )
        
def children2parent( G, children, parent ):
    """
    Create random directed eges from children to parent.
    """
    for c in children:
        u = choose_node( c )
        v = choose_node( parent )
        while u[1] == 'e' and v[1] == 'e':
            u = choose_node( c )
            v = choose_node( parent )
        G.add_edge( u, v )

def parent2children( G, parent, children ):
    """
    Create directed eges from parent to children. This is identical to
    above, except swtiches direction of edges. This is 'downstream'
    towards the leaves of the tree. Random choice of edge using choose_node() for nodes in the edge.
    """
    for c in children:
        u = choose_node( parent )
        v = choose_node( c )
        while u[1] == 'e' and v[1] == 'e':
            u = choose_node( parent )
            v = choose_node( c )
        G.add_edge( u, v )
        

if __name__ == "__main__":

    levels = 3
    num_edges = 2
    nc = sum( [ 2**i for i in range( levels ) ] )
    num_levels = cumsum( [2**i for i in range( levels-1 )] )
    
    subgraphs = make_subgraphs( nc )

    gdict = {}
    for i in range( num_levels[-1] ):
        print i
        print 2*i+1, 2*i +2
        gdict[ subgraphs[i] ] = ( subgraphs[2*i+1], subgraphs[2*i+2] )
    
    Tree = make_tree( gdict, num_edges )
                            

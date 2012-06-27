"""
Take existing graphs and shuffle EI and IE edges to created new
graphs.
"""
import networkx as nx
from random import *
from neural_network.core import graph_object as GO

class NewGraph(GO.Graph):

    def __init__(self, fname='graph.0.5.5.pkl', etype='ei'):
        GO.Graph.__init__(self)
        
        self.fname = fname
        self.etype = etype
        self.graph = self.read_graph(self.fname)
        
        self.ei = [e for e in self.graph.edges(data=True) if e[2]=='ei']
        self.ie = [e for e in self.graph.edges(data=True) if e[2]=='ie']

    def read_graph(self, fname):
        return nx.read_gpickle(fname)
    
    def make_edges(self, ebunch):
        """
        Place edges from ebunch into graph
        """
        for e in ebunch:
            self.graph.add_edge(e[0],e[1],self.etype)

    def repeat_edges(self, ebunch):
        """
        Check for edge repeats in edge list.

        Return True if there are repeat edges. OW, return False.
        """
        if len(ebunch) != len(set(ebunch)):
            return True
        else:
            return False

        
    def remove_edges(self, ebunch):
        """
        remove ebunch edges from graph
        """
        for e in ebunch:
            try:
                self.graph.remove_edge(e[0],e[1])
            except:
                # probably means we have a repeat edge. 
                continue

    def shuffle_edges(self):
        """
        extract type = {'ei', 'ie'} edges--predecessors and
        successors--and shuffle the predecessors.
        """
        if self.etype == 'ei':
            ebunch = self.ei
        elif self.etype == 'ie':
            ebunch = self.ie
        else:
            print "What type of edges are we shuffling?"
            exit(1)


        print "orig ebunch ", ebunch
        print
        
        self.remove_edges(ebunch)

        pre = [e[0] for e in ebunch]
        suc = [e[1] for e in ebunch]

        shuffle(pre)
        newedges = [(pre[i],suc[i],self.etype) for i in range(len(ebunch))]

        while self.repeat_edges(newedges):
            print "shuffling again..."
            shuffle(pre)
            newedges = [(pre[i],suc[i],self.etype) for i in range(len(ebunch))]
        

        print "newedges final ", newedges
        self.make_edges(newedges)

        if self.etype == 'ei':
            self.ei = [e for e in self.graph.edges(data=True) if e[2]=='ei']
        else:
            self.ie = [e for e in self.graph.edges(data=True) if e[2]=='ie']

        

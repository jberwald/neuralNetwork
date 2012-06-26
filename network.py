import numpy
import networkx as NX
import random


class Cell( object ):
    """
    Base class for Network composed of numerous Cells. 
    """
    def __init__( self, node_id, node_type, epsilon ):
        """
        node_id : unique integer identifying node
        node_type : excitatory or inhibitory
        epsilon: roughly equivalent to the oscillatory mode of the cell
        """
        # cell params and other essentials
	self.params = { 'alpha' : 4.0, 
		       'theta' : 0.1, 
		       'theta_x' : 0.1,  
		       'theta_I' : 0.1, 
		       'v_In' : 30., 
		       'v_Ex' : -30., 
		       'beta' : 0.1, 
		       'alpha_In' : 4., 
		       'beta_In' : 0.1,
		       'alpha_x' : 1.0, 
		       'beta_x' : 4.,
		       }                              
        self.node = node_id
        self.type = node_type
        self.epsilon = epsilon

	
    # def __repr__( self ):
    #     return str( (self.node, self.type) )

class Network( Cell ):
    """
    Inherits from Cell. Build a network from cells. Typically, will
    construct a network of Cell objects from a given NetworkX graph.

    A network without epsilons associated to each cell can be
    built. In order to run a neural network simulation, though, a set
    of epsilon values, as either a vector or path to a file, must be
    provided using the keyword arg 'epsilon_vec'.
    """
    def __init__( self, graph, random_couple=False, **args ):
        """
        """
        fargs = { 'g_InEx' : 0.4, 
                  'g_ExIn' : 0.4, 
                  'g_InIn' : 0.4, 
                  'min_rand': 0.3,
                  'max_rand': 0.5,
                  'epsilon_vec': None
                  }
        fargs.update( args )
        self.graph = graph
        self.DIM = len( graph )
        self.random_couple = random_couple
        # read in epsilon vector is necessary
        if hasattr( fargs['epsilon_vec'], '__iter__' ):
            self.epsvec = fargs['epsilon_vec']
        else:
            try:
                ev = numpy.loadtxt( fargs['epsilon_vec'] )
                self.epsvec = ev
            except IOError:
                raise
        # make the network
        self.make_network( **fargs )
        # now set some parameters from the graph
        self.EDIM = len( [u for u in self.network.nodes() if u[1]=='e'] )
        self.IDIM = len( [u for u in self.network.nodes() if u[1]=='i'] )
        self.DIM = 4*( self.EDIM+self.IDIM )
        self.NDIM = len( self.network )

    def __len__( self ):
        return len( self.network.nodes() )

    def make_network( self, **args ):
        """
        Construct a cell bunch of Cells and hook them together with
        appropriate coupling strengths.
        """
        nbunch = self.graph.nodes() 
        ebunch = self.graph.edges()

        # create empty network graph
        self.network = NX.DiGraph()
        # fill network nodes
        for i,u in enumerate( nbunch ):
            self.network.add_node( u, cell=Cell( u[0], u[1], self.epsvec[i] ) )
        # set up edge weights
        if not self.random_couple:
            weights = { 'ei': args['g_ExIn'],
                        'ie': args['g_InEx'],
                        'ii': args['g_InIn']
                        }
        for edge in ebunch:
            etype = edge[0][1] + edge[1][1]
            if edge in self.network.edges():
                continue
            if not self.random_couple:
                wt = weights[etype]
            else:
                wt = random.uniform( args['min_rand'], args['max_rand'] )
            self.network.add_edge( edge[0], edge[1],  type=etype, weight=wt )
  
    def coupling_matrices( self ):
        """
        Populate coupling matrix C with structure of coupling

        clusters : number of clusters connected by one or two edges.

        base_size  : number of nodes in each cluster
        
        interclust  : strength of connection between edges

        (Shift by 4 in all cases since U has partial zeros in
        exc. portion of x-column)
        """
        self.C_ei =\
            numpy.mat(numpy.zeros((self.DIM, self.DIM)))
        self.C_ie =\
            numpy.mat(numpy.zeros((self.DIM, self.DIM)))
        self.C_ii =\
            numpy.mat(numpy.zeros((self.DIM, self.DIM)))

        nbunch = self.network.nodes()
        nbunch.sort()

        ebunch = self.network.edges( data=True )
        for edge in ebunch:
            # predecessor: sends voltage to v
            u = edge[0] 
            # receiving node: thus, this is the row index
            v = edge[1] 
            # E->I
            if edge[2]['type'] == 'ei':
            #if u[1]=='e' and v[1]=='i':    
                s_coord = 2*self.NDIM + u[0] # at Inh. s(v[0])
                self.C_ei[self.EDIM+v[0], s_coord] = edge[2]['weight']
            # I->E
            elif edge[2]['type'] == 'ie':
                # the excitatory "s slots"
                # considering "e voltage" with coupling from "i cell"
                s_coord = 2*self.NDIM + self.EDIM + u[0]  # as Exc. s slot
                self.C_ie[v[0], s_coord] = edge[2]['weight']
            # I->I
            elif edge[2]['type'] == 'ii':
                # same, since EDIM just skips the excitatory "s slots"
                s_coord = 2*self.NDIM + self.EDIM + u[0] # at Inh. s slot
                # still have to skip the voltage "e slots" here
                self.C_ii[self.EDIM + v[0], s_coord] = edge[2]['weight']
            # self.C_ei = self.C_ei.T
            # self.C_ie = self.C_ie.T
            # self.C_ii = self.C_ii.T

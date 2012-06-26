import network
import networkx as NX
import numpy
import optparse, sys, shutil
import cPickle as pkl
try:
    from neural_network.C import neurons_file_cells as neurons
    from neural_network.core import timeseries as TS
except ImportError:
    print "Could not load neurons.so (or similar) module. We'll try to continue..." 
    pass

def make_input_dict( neural_network, extE=0.5, extI=0.3 ):
    """
    cell_input dictionary
    """
    nn = neural_network
    # grab universal params from single cell
    cell = nn.network.node[ (0,'e') ]['cell']
    # make external input vectors - constant values for each type of cell
    externE = extE * numpy.ones( nn.EDIM )
    externI = extI * numpy.ones( nn.IDIM )
    # populate cell_input
    cell_input = cell.params
    cell_input['epsilon'] = nn.epsvec
    cell_input['external_input'] = numpy.concatenate( (externE, externI) )
    return cell_input
              
def load_graph( gname ):
    """
    gname points to a graph stored (pickled) on disk.
    """
    return NX.read_gpickle( gname )

def load_connvec( fname ):
    """
    First try fname by itself. If that fails, then try appending
    the 'default' prefix.

    default prefix is '/data/jberwald/neurons/conns/'

    Returns dictionary keyed by 'g_XY' for the three types of
    connections.
    """
    try:
        with open( fname, 'r'  ) as fh:
            C = pkl.load( fh )
    except IOError:
        prefix = '/data/jberwald/neurons/conns/'
        with open( prefix + fname, 'r' ) as fh:
            C = pkl.load( fh )
    return C

def load_epsvec( fname ):
    """
    First try  fname by itself. If that fails, then try appending
    the 'default' prefix.

    default prefix is '/data/jberwald/neurons/epsilons/'

    Returns vector of epsilon values for a network.
    """
    try:
        E = numpy.loadtxt( fname )
    except IOError:
        prefix = '/data/jberwald/neurons/epsilons/'
        E = numpy.loadtxt( prefix + fname )
    return E

def rand_init_vector( nx, ni ):
    """
    Return intial conditions: [v0,...,vn,w0,...,wn, etc.],
    accounting for number of clusters (nc).
    """
    v = numpy.random.uniform(-0.5, 0.5, nx+ni)
    w = numpy.random.uniform(-0.4, 0.4, nx+ni)
    s = numpy.random.uniform(0, 1, nx+ni)
    xe = numpy.zeros(nx)
    x = numpy.random.uniform(0, 0.5, ni)
    initvec = numpy.concatenate( (v, w, s, xe, x) )
    return initvec

def run_simulation( neural_network, tfinal, cell_input, uinit, **args ):
    """
    Peel off necessary network attributes from params and pass to C
    module.
    
    This function calls neurons.solve_network(). The solution is
    stored in self.soln.
    
    Required args passed to neurons.solve_network() in this order:
    
    tfinal, 
    numpy.asarray(cpl_ei), 
    numpy.asarray(cpl_ie), 
    numpy.asarray(cpl_ii), 
    uinit,
    cell_input, 
    dt, 
    len(ex_nodes), 
    len(in_nodes)
    gamma,
    sigma,
    (tree -- bool)
    
    Returns python list of time series of each cell's equations, as
    well as times at each time step and average exctitatory
    voltage. Also returns final state of ODE solver for 'restart' of
    longer tests.
    """
    fargs = {'dt' : 1e-4,
             'gamma' : 1.,
             'sigma' : 0.,
             'tree': 0
             }
    fargs.update( args )
    
    nn = neural_network
    # Call network_fitness(), C extension
    soln = neurons.network_fitness( float( tfinal ), 
                                   numpy.asarray( nn.C_ei ),
                                   numpy.asarray( nn.C_ie ),
                                   numpy.asarray( nn.C_ii),
                                   uinit,
                                   cell_input, 
                                   fargs['dt'], 
                                   nn.EDIM,
                                   nn.IDIM,
                                   fargs['gamma'],
                                   fargs['sigma'],
                                   fargs['tree'] )
    return soln

if __name__ == "__main__":

    basedir_help = "Path to folders containing graphs, epsilon vectors, "\
                   "connection vectors, etc. [/data/jberwald/neurons]"
    tfinal_help = "Length of simulation, number of time steps. [1e4]"
    graph_name_help = "Name of graph file, instead of a number."
    epsvec_help = "Path to epsilon vector file."
    connvec_help = "Path to connection vector file."
    save_help = "Directory where we want to save results."
    sigma_help = "Gaussian noise level (variance=sigma, mean=0) [0.0]"
    gamma_help = "Amplitude of noise (with variance sigma) [1.0]"
    random_cpl_help = "Generate coupling strengths randomly. "\
                      "See min and max args below. [False]"
    min_rand_help = "Lower bound for random connection strength."
    max_rand_help = "Upper bound for random conection strength."
    write_help = "Write time series to disk while it is being solved. [False]"
    debug_help = "Toggle debugging messages. [False]"
    dryrun_help = "Set up the network and exit. Do not call the solver. [False]"

    
    parser = optparse.OptionParser()
    parser.usage = "python run_neurons.py [options]\n\n"\
                   "Example: python run_network.py -t 5000 --graph_name='graph.0.5.5.pkl' "\
                   "--save='/data/jberwald/neurons/test/' "

    parser.add_option("--basedir", "-b",
                      help=basedir_help,
                      type="string",
                      action="store",
                      dest="basedir",
                      default="/data/jberwald/neurons")
    parser.add_option("--tfinal", "-t",
                      help=tfinal_help,
                      type="int",
                      action="store",
                      dest="tfinal",
                      default=1e3)
    parser.add_option("--graph_name",
                      help=graph_name_help,
                      type="string",
                      action="store",
                      dest="graph_name",
                      default='graph.1.30.pkl')
    parser.add_option("--epsvec", "-e",
                      help=epsvec_help,
                      type="string",
                      action="store",
                      dest="epsvec",
                      default='eps.30.157')
    parser.add_option("--connvec", "-c",
                      help=connvec_help,
                      type="string",
                      action="store",
                      dest="connvec")
    parser.add_option("--save", "-s",
                      help=save_help,
                      type="string",
                      action="store",
                      dest="savedir",
                      default="/data/jberwald/neurons/test/")
    parser.add_option("--debug", "-d",
                      help=debug_help,
                      action="store_true",
                      dest="debug",
                      default=False)
    parser.add_option("--sigma",
                      help=sigma_help,
                      type="float",
                      action="store",
                      dest="sigma",
                      default=0.0)
    parser.add_option("--gamma",
                      help=gamma_help,
                      type="float",
                      action="store",
                      dest="gamma",
                      default=1.0)
    parser.add_option("--random_couple",
                      help=random_cpl_help,
                      action="store_true",
                      dest="rand_cpl",
                      default=False)
    parser.add_option( "--min_rand",
                       help=min_rand_help,
                       type="float",
                       action="store",
                       dest="min_rand",
                       default=0.3)
    parser.add_option( "--max_rand",
                       help=max_rand_help,
                       type="float",
                       action="store",
                       dest="max_rand",
                       default=0.5)    
    parser.add_option("--dryrun",
                      help=dryrun_help,
                      action="store_true",
                      dest="dryrun",
                      default=False)
    
    global options
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
    # get on with running a simulation
    else:
        # from neural_network.C import neurons_file_cells as neurons
        basedir = options.basedir #'/home/jberwald/data/jberwald/'

        save_dir = options.savedir
        if save_dir == None:
            raise ValueError, "Must provide directory to save results"

        if options.debug:
            if options.graph_name == '':
                options.graph_name = 'graph.1.30.pkl'
            print "tfinal", options.tfinal
            print "epsvec", options.epsvec
            print "graph", options.graph_name
            print "connvec", options.connvec

        # set up network
        G = load_graph( options.graph_name )
        if options.connvec:
            connvec = load_connvec( options.connvec )
        if not options.rand_cpl:
            epsvec = load_epsvec( options.epsvec )
            #epsvec = numpy.random.random( 30 )
        # initialize the network
        if options.connvec:
            if options.debug:
                print "connvec", connvec
            neuralNetwork = network.Network( G, epsilon_vec=epsvec,
                                             g_ExIn=connvec['g_ExIn'],
                                             g_InEx=connvec['g_InEx'],
                                             g_InIn=connvec['g_InIn'] )
        else:
            neuralNetwork = network.Network( G, epsilon_vec=epsvec,
                                             random_couple=options.rand_cpl,
                                             min_rand=options.min_rand,
                                             max_rand=options.max_rand )
        # create the coupling matrices
        neuralNetwork.coupling_matrices()
        # initial conditions (linspace, not random)
        initvec = rand_init_vector( neuralNetwork.EDIM, neuralNetwork.IDIM )

        if not options.dryrun:
            # set up simulation
            cell_input = make_input_dict( neuralNetwork )
            soln = run_simulation( neuralNetwork, options.tfinal, cell_input, initvec )
               
            # if necessary, move tmp time series to proper file
            TEMP = '/data/jberwald/neurons/tmp_data/neurons_tmp'
            save_vec = options.savedir + 'vec.out'
            #save_vec = '/home/jberwald/Dropbox/Projects/neuralNet/avgvec.out'
            print "save", save_vec
            shutil.move( TEMP, save_vec )

            # Then, use TimeSeries class to get equally-spaced bins...
            A = numpy.loadtxt( save_vec )
            tseries = TS.Timeseries( tvec=A[:,0], svec=A[:,1] )
            # bin data points to get evenly-spaced data
            tseries.bin_average()
            numpy.savetxt( options.savedir + 'avgvec.out',
                           tseries.avgvec )
            # save connections matrices if using random coupling strengths
            numpy.savetxt( options.savedir + 'C_ei', neuralNetwork.C_ei )
            numpy.savetxt( options.savedir + 'C_ie', neuralNetwork.C_ie )
            numpy.savetxt( options.savedir + 'C_ii', neuralNetwork.C_ii )
            
            # masked arrays for finding the mean
            M_ei = numpy.ma.masked_equal( neuralNetwork.C_ei, 0 )
            M_ie = numpy.ma.masked_equal( neuralNetwork.C_ie, 0 )
            M_ii = numpy.ma.masked_equal( neuralNetwork.C_ii, 0 )
            arrs = { 'mean_cpl_ei': M_ei,
                     'mean_cpl_ie': M_ie,
                     'mean_cpl_ii': M_ii
                     }
            # save the means
            for k in arrs:
                with open( options.savedir + k, 'w' ) as fh:
                    fh.write( str( arrs[k].mean() ) )
                                        
                                        


from pylab import *
import numpy
from neural_network.core import timeseries
from neural_network.netutils import utils

D = utils.Data_Handlers()

def noisy_vec(v, sigma=0.1):
    """
    Add Gaussian noise to each point of a vector
    """
    nv = []
    for x in v:
        # must be better way to do this
        nv.append(numpy.random.normal(x, sigma))  
    return nv

def plot_avgvec(vec, vdir, sigma):
    """
    """
    pname = 'ts.n'+str(sigma)+'.eps'
    outfile = vdir+pname

    ptitle = "Time series of voltage (binned average)"
    xl = "Time ('sec')"
    yl = "Average voltage"
    linestyle = 'b-'

    fig = figure()
    title(ptitle)
    xlabel(xl)
    ylabel(yl)
    plot(vec, linestyle)
    savefig(outfile)

    close(fig)

def run(vecdir, rseed=None, **args):
    """
    Read in vector, add noise to it, compute new dfa, plot noisy
    vector and dfa.

    vecdir : where original avgvec is located (without avgvolt... suffix)
    
    rseed  : for reproducibility of randomness. 
    """
    fargs = {'sigma': 0.1,
             'avgvec': 'avgvolt.sr.1.v.0.out'}

    for k in args.keys():
        fargs[k] = args[k]

    sigma = fargs['sigma']

    print "dir ", vecdir
    print "sigma ", sigma

    if not vecdir.endswith('/'): vecdir += '/'

    avname = fargs['avgvec']
    av = fromfile(vecdir + avname, sep='\n')

    # seed the generator
    if rseed == None:
        rseed = D.load_pickle(vecdir+'noise/rseed.pkl')
        
    numpy.random.seed(seed=rseed)
    nv= noisy_vec(av, sigma=sigma)

    if rseed == None:
        D.write_pkl(rseed, vecdir+'noise/rseed.pkl')

    ts = timeseries.Timeseries()
    ts.avgvec = nv
    ts.avgvoltloc = vecdir+'noise/avgvec.n'+str(sigma)+'.out'
    
    savetxt(ts.avgvoltloc, ts.avgvec, delimiter='\n')
    
    ts.compute_dfa()
    alpha = ts.alpha_min()
    ma = alpha[0]
    ts.linregress(xmin=1.3, xmax=ts.dfa[ma,0])
    ts.plot_dfa(num='.n'+str(sigma), opt=True, dir=vecdir+'noise/')

    plot_avgvec(ts.avgvec, vecdir+'noise/', sigma=sigma)
    
    
if __name__ == "__main__":

    d = '/data/jberwald/tests/fitness_tests/network_topology/'
    suffix = 'n10g0e1c1s0/'
    
    avgvecdir = d+suffix

    run(avgvecdir)

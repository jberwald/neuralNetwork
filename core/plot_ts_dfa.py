from pylab import *
import os
from neural_network.netutils import utils

D=utils.Data_Handlers()

def plot_ts_dfa(fdir):
    """
    fdir should be 'eps0.xxx....'
    """
    if not fdir.endswith('/'):
        fdir += '/'

    fdir += 'vec/'

    # DFA
    try:
        dfa = fromfile(fdir+'voltage.0.sr1.dfa', sep=' ')
    except:
        print "probably empty dfa file in ", fdir
        pass
        
    dfa = dfa.reshape((len(dfa)/2, 2))
    slope, intercept, r ,err = D.linregress(dfa[:,0], dfa[:,1])
    D.linreg_plot(dfa, slope, intercept, fdir+'dfaplot')

    # TS
    avgvolt = fromfile(fdir+'avgvolt.0.sr1.out', sep='\n')
    figure()
    title("Average voltage")
    xlabel("time")
    ylabel("voltage")
    plot(avgvolt)
    savefig(fdir+'ts.eps')

def plot_loop(fdir, start):
    """
    """
    if not fdir.endswith('/'):
        fdir += '/'
    d = os.listdir(fdir)
    print d
    for subd in d:
        if subd.startswith(start):
            print fdir+subd
            try:
                plot_ts_dfa(fdir+subd)
            except:
                print "probably empty dfa file in ", subd
                continue
    

"""
Make small utility to loop through all past sims in netowrk topy. folder and find optimal crossover. save the areas and min tuple, save dfaplot.opt.eps
"""
import os
from numpy import *
from neural_network.core import timeseries as T
from neural_network.netutils import utils

D = utils.Data_Handlers()

basepth = '/data/jberwald/tests/fitness_tests/network_topology/'

def read_dfa(dir, dfa_num=0):
    """
    Read and return avgvolt.*.dfa
    """
    if not dir.endswith('/'):
        dir+='/'
    av = fromfile(dir+'avgvolt.sr.1.v.'+str(dfa_num)+'.dfa', sep='\n')
    return av.reshape((len(av)/2,2))

def make_ts(dfa_dir):
    """
    """
    ts = T.Timeseries()
    ts.dfa = read_dfa(dfa_dir)
    return ts

def make_analysis_dir(dir):
    """
    Make the analysis/ directory to store results in
    """
    try:
        if not dir.endswith('/'): dir+='/'
        os.mkdir(dir+'analysis')
        return True
    except OSError:
        print "does analysis/ directory exists already"
        return False

def find_opt_regression(ts_obj, d):
    """
    For avgvolt.X.Y.dfa in d, find optimized linear regression
    crossovers.
    """
    TS = make_ts(d)
    make_analysis_dir(d)
    return TS.alpha_min()
    
def run(all=False, dfafile=None):

    if all:
        dpre = basepth+'n10g0e'
        
        for e in range(1,11):
            for c in range(1,11):
                print "e ", str(e)
                print "c ", str(c)
                
                thedir = dpre + str(e) + 'c'+str(c)+'s0/'
                TS = make_ts(thedir) # timeseries w/dfa attribute "filled"
                make_analysis_dir(thedir)
                try:
                    optLR = TS.alpha_min()
                except ValueError:
                    print 'Problems in alpha_min()!'
                    
                    D.write_pkl(optLR, thedir+'analysis/opt_linreg.dat')
                    min_dfax = TS.dfa[optLR[0],0]
                    TS.linregress(xmin=1.3,xmax=min_dfax)
                    TS.plot_dfa(opt=True, dir=thedir+'analysis/')
    else:
        thedir = dfafile
        TS = make_ts(thedir) # timeseries with dfa attribute "filled"
        D.write_pkl(optLR, thedir+'analysis/opt_linreg.dat')
        min_dfax = TS.dfa[optLR[0],0]
        TS.linregress(xmin=1.3,xmax=min_dfax)
        TS.plot_dfa(opt=True, dir=thedir+'analysis/')

if __name__ == "__main__":

    run()
        

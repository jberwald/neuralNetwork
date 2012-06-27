#!/usr/bin python
# ----------------------------------------------------------------------
# 
#   data_extract.py
# 
#   Jesse Berwald
# 
#   Version:    0.1  
#   Repository: .svn_repo/neuralNet (if it existed) 
# 
#   Opened:        Jan 2009
# 
# ----------------------------------------------------------------------
"""
Extract results from searches
"""
import os, pylab, numpy
import networkx as NX
from neural_network.netutils import utils
from neural_network.netutils import view_pkl as vp

IO = utils.IOutils()
DH = utils.Data_Handlers()

class Data:

    def __init__(self, svec=None, tvec=None, dvec=None, dir=None, inhib=None,
                 network_size=None, start_eps=None, 
                 final_eps=None, gsyns=None, single=True,
                 binscale=2):
        """
        svec   : voltage (state) vector
        tvec   : time vector corresponding to svec
        dvec   : dfa vector
        
        """
        self.svec = svec
        self.tvec = tvec
        self.dfa = dvec
        self.inhib = inhib
        #self.fft = None
        self.avgvec = None
        self.bins = None

        # info
        self.nsize = network_size
        self.startpt = start_eps
        self.final_eps = final_eps
        self.gsyns = gsyns
        self.single = single
        self.stop_code = ''
        self.search_iters = None


        # bin spacing
        self.step = 0.0
        self.binscale= binscale

        # complexity attributes
        self.slope = 0.0
        self.intercept = 0.0
        self.area = -1

        self.dir = dir

    def linregress(self):
        self.slope, self.intercept, r, sterr = D.linregress(self.dfa[:,0], 
                                                            self.dfa[:,1])

    def _dir_exists(self, dir, subdir):
        dlist = os.listdir(dir)

        print dlist
        print subdir[1:-1]

        if subdir[1:-1] in dlist:
            return 0
        else:
            os.mkdir(dir+subdir)
            return 0

    def plot_series(self, outfile=None, **args):
        """ 
        Plot binned time series.

        Options:

        vectype 

        :nobin:   use original time and voltage vector to plot
        :fft:     use fftfreq and fft vectors
        (neither => use avgvec and binned time vector)

        """
        fargs = {'title': None,
                 'usetex': True,
                 'vectype': 'avg',
                 'show': False,
                 'currdir': 'dir_',
                 'inhib': False}

        for k in args.keys():
            fargs[k] = args[k]

        if fargs['usetex']:
            params = {'axes.labelsize': 14,
                      'text.fontsize': 14,
                      'legend.fontsize': 10,
                      'font.family': 'serif',
                      'text.usetex': True}
            pylab.rcParams.update(params)

        v = self.svec
        t = self.tvec
        
        n1 = 0
        n2 = numpy.where(t < 20000)[0]

        vsnip = v[n1:n2[-1]:100]
        tsnip = t[n1:n2[-1]:100]

        title = "Time series of voltage (non-binned)\n"+\
            "eps: "+str(self.final_eps)+"\n"+\
            "g's: "+str(self.gsyns)
        xlabel = "Time ('sec')"
        ylabel = "Voltage"
        linestyle = 'b-'

        # trim off inner-most directory and replace with "dfaplots/"
        trim = self.dir.rfind('/s_')

        if self.single is False:
            # check whether save dir exist
            self._dir_exists(self.dir[:trim+1]+'plots', self.dir[trim:]) 
            outfile = self.dir[:trim+1]+\
                'plots'+self.dir[trim:]+'ts_plot'+'_'+fargs['currdir']
        else:
            outfile = self.dir+'ts_plot'
#         if not outfile.endswith('.eps'):
#             outfile += '.eps'

        print "plot name ", outfile

        fig = pylab.figure()
        pylab.title(title)
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)
        
#         pylab.text(2,1.8,"eps: "+str(self.final_eps))
#         pylab.text(2,1.6,"g's: "+str(self.gsyns)) 
        

        #pylab.plot(t, v, linestyle)
        pylab.plot(tsnip, vsnip, linestyle)
        
        if self.inhib is not None:
            nx, ny = self.inhib.shape
            for i in range(ny):
                pylab.plot(self.inhib[:,i])

        print "plotted"

        pylab.savefig(outfile)

        if fargs['show']:
            pylab.show()

        pylab.close(fig)

    def plot_dfa(self, outfile=None, **args):
        """
        Plot DFA curve with linear regression line using
        utils.Data_Handlers.linreg_plot()
        """
        fargs = {'show': False,
                 'sigma': 0,
                 'srate': 1,
                 'currdir': 'dir_'
                 }

        for k in args.keys():
            fargs[k] = args[k]

        self.linregress()
        
        # trim off end directory and replace with "dfaplots/"
        trim = self.dir.rfind('/s_')

        if self.single is False:
            # check whether save dir exist
            self._dir_exists(self.dir[:trim+1]+'plots', self.dir[trim:])  

            outfile = self.dir[:trim+1]+\
                'plots'+self.dir[trim:]+'dfa_plot'+'_'+fargs['currdir']
        else:
            outfile = self.dir[:trim+1]+'dfa_plot'+'_'+fargs['currdir']

        D.linreg_plot(self.dfa, self.slope, self.intercept, outfile, 
                      srate=fargs['srate'], sigma=fargs['sigma'], **args)

        if fargs['show']:
            pylab.show()

    def average(self, svec): 
        ## axis = 1 ?? 
        return numpy.average(svec)

    def bin_average(self):
        """ 
        Returns sampled timeseries obtained by binning with window
        size kT, where T = max(t_{i+1} - t_i), k=integer scaling of T,
        and taking the average of the resulting voltage vector over
        each bin.
        """
        tdiff = numpy.diff(self.tvec)
        self.steptup = (min(tdiff), max(tdiff))
        self.step = self.binscale * self.steptup[1]

        nbins = int(numpy.ceil(self.tvec[-1]/(self.step)))
        
        # this truncates the end of tvec
        endpt = self.step * (nbins-1)
        self.bins = numpy.linspace(0, endpt, nbins)
        
        self.avgvec = self.binavg()

    def binavg(self):
        indices = self.bins_to_indices()
        avglist = [numpy.average(self.svec[k]) for k in indices]      
        return numpy.asarray(avglist)

    def bins_to_indices(self):
        """ Return indices of tvec that fall into bins
        """
        n = len(self.bins)
        rbins = []

        leftend = 0
        for i in xrange(n-1):
            b1 = numpy.where(self.tvec[leftend:] >= self.bins[i])
            b2 = numpy.where(self.tvec[leftend:] <= self.bins[i+1])

            interval = numpy.intersect1d(leftend+b1[0],leftend+b2[0])
            leftend = interval[-1]
            rbins.append(interval)

        return rbins


#------------------------------------------

def extract_vectors(dir, type='state'):
    """
    Extract and return final dfa vector from search

    dir : upper level search directory,
    eg. /data1/jberwald/simplex_search/s_N
    """
    if dir.endswith('/'): pass
    else: dir += '/'
    
    d = os.listdir(dir)
    
    for f in d: 
        if f.startswith('eps0.'):
            dir += f + '/vec/'
            dv = os.listdir(dir)
            break

    if type is 'state':
        # just voltage vectors
        vec=[u for u in dv if u.startswith('voltage') and u.endswith('.out')]
        vec.sort()
    elif type is 'dfa':
        vec = [u for u in dv if u.startswith('voltage') and u.endswith('.dfa')]
        vec.sort()
    elif type is 'time':
        vec = [u for u in dv if u.startswith('tvec') and u.endswith('.out')]
        vec.sort()
    elif type is 'inhib':
        vec = [u for u in dv if u.startswith('tvec') and u.endswith('.out')]
        vec.sort()
    else:
        raise ValueError("Unrecognized type for vector")

    # grab first and last from sorted list
    v0 = vec[0]
    v1 = vec[-1]

    if type is 'state' or type is 'dfa' or type is 'inhib':
        if int(v1[8]) > int(v0[8]):
            # read in v1 (last time series
            v = IO.read_text_vec(dir+v1)        
        else:
            # last one put first in list, so read in v0
            v = IO.read_text_vec(dir+v0)
    else:
        # type = 'time'
        if int(v1[5]) > int(v0[5]):
            # read in v1
            v = IO.read_text_vec(dir+v1)        
        else:
            # last one put first in list, so read in v0
            v = IO.read_text_vec(dir+v0)

    if type is 'dfa':
        return reshape(v)
    else:
        return v

def extract_single_vec(fname, type='dfa'):
    v = IO.read_text_vec(fname)        
    if type is 'dfa':
        return reshape(v)
    else:
        return v
    
def reshape(vec):
    """
    reshape dfa vector into two-column array
    """
    n = len(vec)
    return vec.reshape((n/2,2))

def finalvals(dir):
    """ 
    grab final values from search, split into eps and g arrays
    """
    fname = dir + 'finalvals.simplex_search'
    if fname in dir:
        f = vp.load_pkl(fname)
        return f[0:2]  # final fitness and search space values
    else:
        return None

def _split_final_vals(fv, ncells):
    """
    parse final values and return 10**(eps) and gsyns
    """
    vals = fv[0]
    fit = fv[1]
    
    eps = vals[:ncells]
    g = vals[ncells:]
    return eps, g

def _start_eps(dir, sep='\n'):
    fname = dir + 'eps'
    return IO.fromfile(fname, sep=sep)

def network_size(eps):
    """
    Obtain network size from len(epsilon)
    """
    return len(eps)

def _extract_ivec(dir, ncells=2):
    if not dir.endswith('/'):
        dir+='/'
    d = os.listdir(dir+'vec/')
    for f in d: 
        if f.startswith('inhib'):
            v = numpy.fromfile(f, sep='\n')
            D.inhib = v.reshape((len(v)/ncells, ncells))
            break
        else: continue

def init_data(dir, info=False, single=True, avgvec=False, inhib=True):
    """
    For a single search results directory extract info
    """
    if not dir.endswith('/'):
        dir+='/'

    # vectors
    svec = extract_vectors(dir, type='state')
    tvec = extract_vectors(dir, type='time')
    dfa = extract_vectors(dir, type='dfa')

    # info
    if info is True:
        fvals = finalvals(dir)
        seps = _start_eps(dir)
        nsize = network_size(seps)

        feps, gsyns = _split_final_vals(fvals, nsize)
        if len(gsyns)==0:
            gsyns = None

        if inhib:
            ivec = _extract_ivec(dir)
        else:
            ivec = None

        data = Data(svec, tvec, dfa,
                    network_size=nsize,
                    start_eps=seps,
                    final_eps=10**(feps),
                    gsyns=gsyns,
                    single=single,
                    dir=dir)
        
    else:
        if inhib:
            ivec = _extract_ivec(dir)
        else:
            ivec = None
        data = Data(svec, tvec, dfa, dir, inhib=ivec)

    if avgvec:
        data.bin_average()

    return data

def do_directory(dir, start_at=None, dfaplot=False, tsplot=False):
    """
    Analyze all searches in given directory dir

    start_at  : int n -> start analyzing at directory s_n
    """
    if dir.endswith('/'):
        pass
    else:
        dir+='/'

    dlist = os.listdir(dir)
    
    if start_at is not None:
        # only analyze directories past start_at
        for f in dlist:
            if f.endswith('plots'):
                continue
            print "dir: ", dir+f
            sub_num = int((dir+f).lstrip(dir+'s_'))
            if sub_num < start_at:
                continue
            else:
                d = init_data(dir+f, info=True, single=False)
            if dfaplot:
                d.plot_dfa(currdir=f, eps=d.final_eps,
                           seps=d.startpt, g=d.gsyns)
            if tsplot:
                d.plot_series(currdir=f)
            
    else:
        for f in dlist:
            if f.endswith('dfaplots'):
                continue
            print "dir: ", dir+f
            d = init_data(dir+f, info=True, single=False)
            if dfaplot:
                d.plot_dfa(currdir=f, eps=d.final_eps, 
                           seps=d.startpt, g=d.gsyns)
            if tsplot:
                d.plot_series(currdir=f)

    print 'Done!'

def fitness(m, a=0, alpha=0.1, func=1):
    """
    """
    if func is 1:
        return (1-m)**2 + alpha*(a)
    elif func is 2:
        return 1./(1-m)
    elif func is 3:
        return (1-m)**2   



def analyze_slope_area(dir, alpha=0.1, fitfunc=1):
    """
    Analyze slope and area for dfa curve. Then filter through
    different fitness functions.

    fitfunc  : 1 == (1-m)**2 + scaled area
             : 2 == 1/(1-m)
             : 3 == (1-m)**2
    """
    dir += 'vec/'
    dv = extract_single_vec(dir+'voltage.0.sr1.dfa')
    lin = DH.linregress(dv[:,0], dv[:,1])
    m = lin[0]
    icept = lin[1]

    if fitfunc is 1:
        area = DH.linreg_dfa_area(dv, slope=m, intercept=icept)
        print "area ", area
        fit = fitness(m, area, alpha=alpha, func=fitfunc)
    elif fitfunc is 2:
        fit = fitness(m, func=fitfunc)
    elif fitfunc is 3:
        fit = fitness(m, func=fitfunc)

    print "fit ", fit

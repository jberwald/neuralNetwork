#!/usr/bin python
import numpy, pylab, time
from netutils import fileio, utils
#import fitfunc

#from matplotlib.cbook import report_memory

D = utils.Data_Handlers()

class Timeseries:
    """
    Object to hold methods specific to a timeseries.
    
    tvec -- time vector for alignment purpose and storage

    svec -- state vector to be analyzed.Numpy ndarray object. Can be
    1d or nd. If nd then we have not averaged the voltage components
    (mostly testing/plotting if nd).

    """
    def __init__(self, tvec=None, svec=None, inhib=None,
                 binscale=4, srate=1, testdir=None):

        # array objects
        self.tvec = tvec
        self.svec = svec
        self.inhib=inhib
        self.dfa = None
        self.fft = None
        self.avgvec = None
        self.bins = None
        self.spike_dfa = None
        self.srate = srate

        # bin spacing
        self.step = 0.0
        self.binscale = binscale
        
        # all parameters and network structure
        self.testdir = testdir
        
        # complexity attributes
        self.slope = 0.0
        self.intercept = 0.0
        self.spike_slope = 0.0
        self.spike_intercept = 0.0
        self.area = -1
        self.dfa_partition = []
        self.spike_dfa_partition = []
        self.finite_diff = numpy.empty([])
        self.spike_gaps = None
        self.spikes = None
        self.alpha = None
        
        # I/O stuff
        if self.testdir is None:
            print "Setting testdir to /data/jberwald/debug"
            self.testdir = "/data/jberwald/debug"
        self.voltloc  = ''
        self.avgvoltloc = ''
        self.spikeloc = ''

    def __del__(self):
        pass

    def save_vecs(self, vec='a', suffix='out', trunc=False):
        """
        save vectors of different types. safety first: try-except
        stuff.
        """
        if not self.testdir.endswith('/'): 
            self.testdir += '/'
        if not fileio.check_dir(self.testdir):
            fileio.mkdir(self.testdir)
        self._savevec(vectype=vec,
                      trunc=trunc,
                      suffix=suffix)

    def _savevec(self, vec=None, vectype='v', suffix=None, **args):
        """
        Optional: 

        threshold 
        -- label file with 'clip' if we thresholded data.

        savedir 
        -- default is neurons.vecstore

        vectype 
           t : time
           v : original voltage
           a : averaged/binned voltage
           i : inhibitory voltage (indiv. cells)
           g : time gaps vector

        """
        fargs = {'trunc': 0 }
        for k in args:
            fargs[k] = args[k]
        
        if not self.testdir.endswith('/'): self.testdir += '/'
        savedir = self.testdir

        dot = '.'
        if suffix == None:
            suf = 'out'
        else:
            suf = suffix
        if 'trunc' in args.keys():
            if args['trunc'] != -1:
                suf = dot.join( ['trunc', str(args['trunc']), suf] )
#                suf.insert(-1,'trunc_'+str(args['trunc']))
        if vectype is 't':
            fullsuf = dot.join(['tvec', suf])
            fname = savedir + fullsuf
            numpy.savetxt(fname, self.tvec[::self.srate], fmt='%1.15f')
        elif vectype is 'v':
            fullsuf = dot.join(['voltage', suf])
            fname = savedir + fullsuf
            numpy.savetxt(fname, self.svec[::self.srate], fmt='%1.15f')
            self.voltloc = fname
        elif vectype is 'a':
            fullsuf = dot.join(['avgvolt', suf])
            fname = savedir + fullsuf
            numpy.savetxt(fname, self.avgvec[::self.srate], fmt='%1.15f')
            self.avgvoltloc = fname
        elif vectype is 'i': 
            fullsuf = dot.join(['inhibvolt',suf])
            fname = savedir + fullsuf
            numpy.savetxt(fname, self.inhib[::self.srate], fmt='%1.15f')
        elif vectype is 'g':
            fullsuf = dot.join(['timegap',suf])
            fname = savedir + fullsuf
            numpy.savetxt(fname, self.spike_gaps[::self.srate], fmt='%1.15f')
            self.spikeloc = fname
        elif vectype is 'd':
            fullsuf = dot.join(['dfa',suf])
            fname = savedir + fullsuf

            print "path to dfa", fname
            
            numpy.savetxt(fname, self.dfa, fmt='%1.15f')
            self.dfaloc = fname
            
        else:
            if 'suffix' not in args.keys():
                raise KeyError, "Must provide suffix for _savevec()"
            fullsuf = dot.join([args['suffix'],'sr',str(self.srate),
                                str(self.neurons.direction),
                                str(self.neurons.counter),'out'])
            fname = savedir + fullsuf
            numpy.savetxt(fname, vec, fmt='%1.15f')

    #=================
    # Binning
    #=================
    def average(self, svec): 
        ## axis = 1 ?? 
        return numpy.average( svec )

    def bin_average(self, binscale=None):
        """ Returns sampled timeseries obtained by binning with window
        size kT, where T = max(t_{i+1} - t_i), k=integer scaling of T,
        and taking the average of the resulting voltage vector over
        each bin. This has the effect of making each step on the
        x-axis equal.
        """
        if binscale:
            self.binscale = binscale
        tdiff = numpy.diff(self.tvec)
        self.steptup = (min(tdiff), max(tdiff))
        self.step = self.binscale * self.steptup[1]
        # self.nbins = int(numpy.ceil( self.tvec[-1]/(self.step)) )
        self.nbins = int( numpy.ceil( len(self.tvec)/self.step ) )
        # this truncates the end of tvec
        self.endpt = int( self.tvec.max() )#self.step * (self.nbins-1)
        self.bins = numpy.linspace(0, self.endpt, self.nbins)
        self.avgvec = self.do_binavg()

    def do_binavg(self):
        """
        Returns the average voltage value over each evenly spaced bin.
        """
        indices = self.bins_to_indices()
        avglist = []
        if indices[1]==0:
            # This checks if there are multiple bins at the beginning
            # that do not contain any tvec value.
            z = numpy.where(indices==0)[0]
            a = z[-1]
            indices = indices[a:]
            self.bins = self.bins[a:] #need to trim tvec as well
        # number of bins 
        n = len(indices)

        for k in xrange(n-1):
            #avglist = [numpy.average(self.svec[k] for k in indices)]
            if indices[k]==indices[k+1]:
                continue
            avglist.append( numpy.average(self.svec[indices[k]:indices[k+1]]) )
            # if k > 3020:
            #     print "bin", self.svec[indices[k]:indices[k+1]]
            #avgarr[k] = numpy.average(self.svec[indices[k]:indices[k+1]])
        return numpy.ascontiguousarray( avglist )

    def bins_to_indices(self):
        """
        Return indices of tvec that fall into bins
        """
        n = len(self.bins)
        dbins = numpy.digitize( self.tvec, self.bins )
        h = numpy.histogram(dbins, bins=numpy.arange( n ))        
        return numpy.cumsum(h[0])

    #=============
    # DFA stats
    #=============
    def calc_crossover(self):
        
        x = self.dfa[:,0]
        y = self.dfa[:,1]

        linreg = D.linregress(x,y)

        # store in output list
        slope = linreg[0]
        intercept = linreg[1]
        liny = D.make_line(slope, x, intercept)
        indices = D.dfa_relative_linreg(y, liny)

        # keep just the indices around the crossover
        indtrunc = D.find_hump(indices)
        indices = indtrunc

        # find index where dfa curve and lin. reg. curve are greatest
        # distance apart (approx.)
        mdindex = D.find_max_dist_index(x[indices], 
                                        y[indices], 
                                        slope,
                                        intercept)
        try:
            crossover_index = indices[0] + mdindex
        except:
            crossover_index = 0
        self.crossover = self.dfa[crossover_index,0]

    def partition_dfa(self, dfa, xmin, xmax):
        """
        """
        u = numpy.where(dfa[:,0] > xmin)[0]
        v = numpy.where(dfa[:,0] < xmax)[0]
        m = numpy.intersect1d(u,v)
        return (dfa[:m[0]], dfa[m], dfa[m[-1]:-1])

    def triple_linregress(self, xmin, xmax, spikes=False):
        """
        """
        if spikes: 
            dfa = self.spike_dfa
        else: 
            dfa = self.dfa

        dfa_parts = self.partition_dfa(dfa, xmin, xmax)
        for d in dfa_parts:
            m, icept, r, err = D.linregress(d[:,0],d[:,1])
            if spikes:
                self.spike_dfa_partition.append((m, icept, d))
            else:
                self.dfa_partition.append((m, icept, d))

    def linregress(self, xmin=0, xmax=-1, spikes=False):
        """
        Perform linear regression on dfa vector. 

        Optional args used to find linear regression on segment of DFA
        curve. xmin, xmax expected to be exponents (eg, xmin =>
        10^xmin window size). Call triple_linregress with these args.
        """
        if spikes:
            dfa = self.spike_dfa
        else:
            dfa = self.dfa

        if xmin != 0 and xmax != -1:
            # set attributes for three separate linear regression
            # lines
            self.triple_linregress(xmin, xmax, spikes=spikes) 
        
        # we always need this info if doing spikes
        elif spikes:
            self.spike_slope, self.spike_intercept, r, sterr = \
                D.linregress(dfa[:,0], dfa[:,1])
        else:
            self.slope, self.intercept, r, sterr =\
                D.linregress(dfa[:,0], dfa[:,1])

    def set_alpha_range(self, xmin, xmax):
        """
        Set range of indices on DFA x-axis to choose alpha from
        """
        if xmax == -1:
            xmax = self.dfa[:,0][-1]
        d0 = numpy.where(self.dfa[:,0] >= xmin)[0]
        d1 = numpy.where(self.dfa[:,0] <= xmax)[0]
        return numpy.intersect1d(d0,d1)

    def alpha_area(self, lines, xind, alpha):
        """
        Find area between DFA curve and lines. Assume single parameter
        alpha, giving ony two regions to sum area over. Return sum of
        abs(area).
        """
        xt0 = lines[0]
        xt1 = lines[1]
        a0 = numpy.sum(numpy.abs(xt0 - self.dfa[:,1][xind[0]:alpha]))
        a1 = numpy.sum(numpy.abs(xt1 - self.dfa[:,1][alpha:xind[-1]]))
        return a0+a1
        
    def alpha_min(self, xmin=2.0, xmax=-1):
        """
        Method to call to calculate minimizing alpha parameter.
        """
        area = {}

        # list of indices
        alphaRange = self.set_alpha_range(xmin, xmax)
        aleft = alphaRange[0]
        aright = alphaRange[-1]
        for alpha in alphaRange[2:-2]: #start away from endpt
            reg1 = D.linregress(self.dfa[:,0][aleft:alpha], 
                                self.dfa[:,1][aleft:alpha])
            reg2 = D.linregress(self.dfa[:,0][alpha:aright], 
                                self.dfa[:,1][alpha:aright])
            lr1 = D.make_line(reg1[0], 
                              self.dfa[:,0][aleft:alpha],
                              reg1[1])
            lr2 = D.make_line(reg2[0], 
                              self.dfa[:,0][alpha:aright],
                              reg2[1])
            area[alpha] = self.alpha_area((lr1,lr2), alphaRange, alpha)
        
        areatup = (D.get_key(area, min(area.values())), area)
        self.alpha = (self.dfa[areatup[0],0].item(), areatup[0])        
        return areatup
        
    #==================
    # Spike analysis
    #==================            
    def time_gaps(self, high_t):
        """
        Calculate time gaps between spikes. high_t is array of indices
        where svec is above threshold == (svec.mean + window); see
        threshold_spikes().
        """
        self.spikeIndices = self.spike_detector(high_t)
        tvec = self.tvec[high_t]
        self.interspike = numpy.diff(tvec[self.spikeIndices])

        print 'interspike'
        print len(self.interspike)

        #return numpy.diff(self.spikes)

    def spike_detector(self, thresh):
        """
        Detect basic change in derivative of vec using numpy.diff. To
        detect spikes we look for local maxima ==> where diff of a
        sign vector == 2.

        over_thresh -- indices cooresponding to svec values >
        threshold

        returns     -- indices of spikes
        """
        print 'svec '
        print len(self.svec)

        pylab.plot(self.svec)

        print numpy.diff(self.svec[thresh])
        dvec = numpy.diff(self.svec[thresh])
        print "dvec ", dvec
        print len(dvec)
        
        sig_dvec = numpy.sign(dvec) # pos/neg slopes
        print "sig ", sig_dvec[:100]
        print len(sig_dvec)

        extrema = numpy.diff(sig_dvec) # local extrema

        print "extrema ", extrema[:100]
        print len(extrema)

        # indices of local maxima
        vmax = numpy.where(extrema == 2)[0] 
        print 'vmax ', vmax
        print len(vmax)
        return vmax

    def threshold_spikes(self, buf=0.2):
        """
        Look for "spikes" above threshold: Threshold the true voltage
        timeseries by mean+buf. Pass this thresholded time
        series to spike_detector()
        """
        vmean = self.svec.mean()
        high_thresh = numpy.where(self.svec > (vmean+buf))[0]

        print "thresh ", high_thresh

        self.spike_gaps = self.time_gaps(high_thresh)
        
    #================
    # Fnite diff
    #================
    def compute_finite_diff(self, vec):
        """ 
        Compute first-order finite difference derivative.
        """
        self.finite_diff = fitfunc.finite_diff(vec)

    #================
    # Sample Entropy
    #================
    def sampen(self, m=2, r=0.2):
        """
        m == maximum epoch length
        r == tolerance
        """
        outfile = self.testdir+'sampen.out'
        fitfun.sampen(self.avgvoltloc, outfile=outfile, m=m, r=r)

    def sampen_l2norm(self, fname, start=1):
        """
        Calculate L2 norm of sample entropy function from <start> to
        last entry.
        """
        samp = self.extract_sampen(fname=fname, start=start)
        return fitfunc.L2_norm()
        

    def extract_sampen(self,fname='sampen.out', start=1):
        return fitfunc.extract_sampen(self.testdir+fname,start=start)
        

     

    #================
    # Plotting
    #================
    def plot_dfa(self, outfile=None, **args):
        """
        Plot DFA curve with linear regression line using
        utils.Data_Handlers.linreg_plot()
        """
        fargs = {'save': True,
                 'show': False,
                 'eps': None,
                 'srate': 1,
                 }
        
        for k in args.keys():
            fargs[k] = args[k]

        if self.slope == 0.0:
            self.linregress()
        
        if outfile is None:
            outfile = self.testdir+'dfaplot_'
            if 'spikes' in args.keys():
                if args['spikes']: outfile += 'spikes_'
        else:
            if 'dir' in args.keys():
                if not args['dir'].endswith('/'): args['dir']+='/'
                outfile = args['dir'] + 'dfa_plot'
                if 'spikes' in args.keys():
                    if args['spikes']: outfile += 'spikes_'
            else:
                if fargs['save']:
                    print "Provide file to save dfa plot to"
                    return -1

        if 'direction' in args.keys():
            outfile += str(args['direction'])+'.'
        if 'node' in args.keys():
            outfile += str(args['node'])+'.'
        if 'num' in args.keys():
            outfile += str(args['num'])
        if 'trunc' in args.keys():
            outfile += '.trunc_'+str(args['trunc'])
        if 'opt' in args.keys():
            outfile += '.opt'

        # if self.neurons is None:
        #     # using Timeseries standalone
        #     srate = sigma = 0
        #     show = fargs['show']
        # else:
        #     srate = self.neurons.srate
        #     sigma = self.neurons.sigma
        #     show = self.neurons.showplots

        if 'spikes' in args.keys():
            if args['spikes']: 
                dfavec = self.spike_dfa
                slope = self.spike_slope
                intercept = self.spike_intercept
                part = self.spike_dfa_partition
            else: 
                dfavec = self.dfa
                slope = self.slope
                intercept = self.intercept
                part = self.dfa_partition
        else: 
            dfavec = self.dfa
            dfavec = self.dfa
            slope = self.slope
            intercept = self.intercept
            part = self.dfa_partition

        D.linreg_plot(dfavec, slope, intercept, outfile, 
                      srate=fargs['srate'],
                      tex=True,
                      sigma=fargs['sigma'],
                      dfa_parts=part, 
                      save=fargs['save'],
                      eps=fargs['eps'])

        if fargs['show']:
            pylab.show()

    def plot_series(self, outfile=None, **args):
        """ 
        Plot DFA curve with a linear regression line using

        Options:

        vectype 

        :nobin:   use original time and voltage vector to plot
        :fft:     use fftfreq and fft vectors
        (neither => use avgvec and binned time vector)

        """
        fargs = {'title': None,
                 'usetex': True,
                 'vectype': 'avg',
                 'dfaout': None,
                 'num': 0,
                 'direction': 0}

        for k in args.keys():
            fargs[k] = args[k]

        if fargs['usetex']:
            params = {'axes.labelsize': 14,
                      'text.fontsize': 14,
                      'legend.fontsize': 10,
                      'font.family': 'serif',
                      'text.usetex': True}
            pylab.rcParams.update(params)

        if fargs['vectype'] is 'nobin':
            # plot original vector against time
            v = self.svec
            t = self.tvec
            title = "Time series of non-binned voltage"
            xlabel = "Time ('sec')"
            ylabel = "Average voltage"
            linestyle = 'b-'

        elif fargs['vectype'] is 'fft':
            v = self.fft
            t = self.fftfreq
            title = "Fourier Transform:\n"+\
                " with external input "+str(self.Network.extE)
            xlabel = 'Frequency'
            ylabel = 'Amplitude'
            linestyle = 'r-'

        elif fargs['vectype'] is 'dfa':
            self.plot_dfa(outfile=fargs['dfaout'])   
            return None

        else:
            v = self.avgvec
            t = self.bins[1:-1]
            title = "Time series of voltage (binned average)"
            xlabel = "Time ('sec')"
            ylabel = "Average voltage"
            linestyle = 'b-'

        if outfile is None:
            if fargs['vectype'] is 'fft':
                outfile = self.testdir + 'fft.eps'

            elif fargs['vectype'] is 'nobin':
                outfile = self.testdir + 'ts_nobin.eps'
            else:
                outfile = self.testdir + 'ts_binned.'+\
                    str(fargs['direction'])
                if 'node' in args.keys():
                    outfile += '.'+str(fargs['node'])
                if 'num' in args.keys():
                    outfile += '.'+str(fargs['num'])
                if 'trunc' in args.keys():
                    outfile += '.trunc_'+str(fargs['trunc'])
        if not outfile.endswith('.eps'):
            outfile += '.eps'

        fig = pylab.figure()
        pylab.title(title)
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)
        pylab.plot(t, v, linestyle)
        pylab.savefig(outfile)

        if self.neurons.showplots:
            pylab.show()

        #print "memory usage ", report_memory()

        fig.clf()
        pylab.close(fig)


    def compute_dfa(self, vectype='avg', trunc=False, minbox=4, max_boxes=10):
        """
        Compute the DFA for a Timeseries object. 

        max_boxes ==> npts/max_boxes (dfa default npts/4)
        minbox ==> min window size
        """
        n = len(self.avgvec)/max_boxes
        
        fname_out = D.run_dfa(self.avgvoltloc, minbox=minbox, maxbox=n)
        dfaout = numpy.fromfile(fname_out, sep=' ')
        
        nx = len(dfaout)
        rdfa = dfaout.reshape((nx/2,2))

        # for testing purposes
        if trunc is not False:
            inds = numpy.where(rdfa[:,0] < trunc)[0]
            maxind = inds[-1]
        else:
            maxind = nx-1
        self.dfa = rdfa[0:maxind]

    def compute_time_gap_dfa(self, trunc=False, minbox=4, max_boxes=10):
        """
        Compute the DFA for time gaps between spikes for a Timeseries
        object.

        max_boxes ==> npts/max_boxes (dfa default npts/4)
        minbox ==> min window size
        """
        if self.spike_gaps == None:
            # self.voltage_threshold()
            print "need to compute voltage_threshold first!"
            exit(1)

        n = len(self.spike_gaps)/max_boxes
        
        fname_out = D.run_dfa(self.spikeloc, minbox=minbox, maxbox=n)
        dfaout = numpy.fromfile(fname_out, sep=' ')
        
        nx = len(dfaout)
        rdfa = dfaout.reshape((nx/2,2))

        # for testing purposes
        if trunc is not False:
            inds = numpy.where(rdfa[:,0] < trunc)[0]
            maxind = inds[-1]
        else:
            maxind = nx-1
        self.spike_dfa = rdfa[0:maxind]

    def compute_fft(self, vectype='avg'):
        """
        Compute FFT of Timeseries object. 
        """
        fftfreq = numpy.fft.fftfreq(len(self.bins), d=self.step)
        self.fft = numpy.abs(numpy.fft.rfft(self.avgvec))

        n = len(self.fft)
        self.fftfreq = fftfreq[:n]
 
        

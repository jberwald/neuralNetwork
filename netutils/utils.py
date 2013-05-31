#!/usr/bin python
# ----------------------------------------------------------------------
# 
#   utils.py
# 
#   Jesse Berwald
# 
#   Version:       
#   Repository: .svn_repo/neuralNet   
# 
#   Opened:        Thu 20 Dec 2007
# 
# ----------------------------------------------------------------------
"""
Modules containing utility functions for data manipulation, IO
operations, statistical analysis, plotting, etc.


Classes:

Data_Handlers()   : Plotting and some of the array manipulation needed
to plot (eg., slicing array to retain only excitatory cell data).

IOutils()         : class containing methods to assist in input/output and
pickling of graphs and data.

Stats()           : Statistical analysis functions.

Network_Anaylyzer : Object designed to handle the analysis of neural
networks
"""
import csv, gc, time, os
import pylab, numpy, math
import cPickle as pkl
import subprocess as sp
import tempfile

#-----------------------------------------------
#
# Statistical Functions
#
#-----------------------------------------------
class Stats:
    """
    Functions required for linear regression. Does not require SciPy
    (so more usable on 64-bit machines since some fortran-linked
    integration utilities don't actually work on 64-bit architecture).
    Source code from scipy.stats.stats
    """

    def __init__(self):
        pass

    def linregress(self, *args):
        """ 
        Calculates a regression line on two arrays, x and y, corresponding 
        to x,y pairs.  
        If a single 2D array is passed, linregress finds dim with 2 levels 
        and splits data into x,y pairs along that dim. 
        
        Returns: slope, intercept, r, two-tailed prob, stderr-of-the-estimate 
        
        Source code from scipy.stats.stats. 
        """ 
        TINY = 1.0e-20 
        if len(args) == 1:  # more than 1D array? 
            args = numpy.asarray(args[0]) 
            if len(args) == 2: 
                x = args[0] 
                y = args[1] 
            else: 
                x = args[:,0] 
                y = args[:,1] 
        else: 
            x = numpy.asarray(args[0]) 
            y = numpy.asarray(args[1]) 
            n = len(x) 
            xmean = numpy.mean(x,None) 
            ymean = numpy.mean(y,None) 
            xm,ym = x-xmean, y-ymean 
            r_num = numpy.add.reduce(xm*ym) 
            r_den = math.sqrt(self.ss(xm)*self.ss(ym)) 
        if r_den == 0.0: 
            r = 0.0 
        else: 
            r = r_num / r_den 
            if (r > 1.0): r = 1.0 # from numerical error 
        slope = r_num / self.ss(xm) 
        intercept = ymean - slope*xmean 
        sampstd = self.samplestd(y)

        # why did this throw math domain errors ?!?!
        #sterrest = math.sqrt(1-r*r)*self.samplestd(y) 
        sterrest = math.sqrt(1-r*r)*sampstd
        return slope, intercept, r, sterrest 

    def ss(self, a):
        """
        Sum of squares. Assume one axis only.
        """
        #a, axis = _chk_asarray(a, axis) 
        return numpy.sum(a*a) 

    def betai(self, a, b, x):
        """
        """
        x = numpy.asarray(x) 
        x = numpy.where(x < 1.0, x, 1.0)  # if x > 1 then return 1.0 
        return special.betainc(a, b, x) 

    def samplestd(self, a):
        """
        assume only one axis
        """
        return math.sqrt(self.samplevar(a)) 

    def samplevar(self, a, axis=0):
        """
        scipy.stats samplevar func
        """
        mn = numpy.expand_dims(numpy.mean(a, axis), axis) 
        deviations = a - mn 
        n = a.shape[axis] 
        svar = self.ss(deviations) / float(n) 
        return svar 

    def fft_linregress(self, time_series, xstart=10, xstop=1e5, xstep=100):
        """
        Calculate the slope of linear regression of the loglog plot of the fft
        """
        fourier = numpy.fft.rfft(time_series - numpy.average(time_series))
        xaxis = numpy.arange(xstart, xstop, xstep)
        (slope, intercept, r, stderr)=\
            self.linregress(numpy.log10(xaxis),numpy.log10(\
                numpy.abs(fourier[xstart: xstop: xstep])))
        return slope, stderr

    def decimal_place(self, x):
        a = math.log10(x)
        if a < 1:
            return abs(math.floor(a))
        elif a > 1:
            return math.ceil(a)
        

class IOutils:
    """
    Input/Output utilities.
    """
    def __init__(self):
        pass

    def getopts(self, argv):
        """put argv into a dictionary"""
        opts = {}
        while argv:
            if argv[0][0] == '-':
                opts[argv[0]] = argv[1]
                argv = argv[2:]
            else:
                argv = argv[1:]
        return opts

    def write_data(self, dvec, fname):
        """
        Write vector to newline-separated file

        dvec    -- data vector
        fname   -- string object containing file name
        """
        if type(fname) == tuple:
            # assume fname is mkstemp
            fd = fname[0]
            fpath = fname[1]
            fh = os.fdopen(fd, 'wb')
        else:
            fh = open(fname,'wb')
        writer = csv.writer(fh, delimiter='\n')
        writer.writerow(dvec)
        fh.close()
        return fname

    def append_data(self, data, fname):
        """
        file fname assumed to exist already. append str(data) to end
        of file.
        """
        try:
            fh = open(fname, 'a')
        except IOError:
            print "error opening file ", fname
        fh.write(str(data)+'\n')
        fh.close()        

    def write_value(self, value, fname):
        """
        Write a network value to file fname
        """
        fh = open(fname,'wb')
        fh.write(str(value))
        fh.close()

    def read_text_vec(self, fname, ncols=1, sep=' '):
        """ 
        Read a text file and return a vector. If ncols = 2 assumed to be
        dfa data, so reshape.
        """
        fh = open(fname, 'r')
        vec = numpy.fromfile(fh, sep=sep)
        fh.close()
        if ncols == 2:
            n = len(vec)
            vec = vec.reshape((n/2,2))
        return vec

    def call_oct_stats(self):
        fh = open("oct_stats.m",'rb')
        sp.call("/usr/bin/octave", stdin=fh)

    def fromfile(self, fname, sep='\n'):
        return numpy.fromfile(fname, sep=sep)

    def dfa_data_to_array(self, fname):
        """
        Returns 2d array of dfa data for plotting and analysis.
        """
        fh = open(fname, 'rb')
        data = numpy.fromfile(fh, sep=' ')
        fh.close()
        n = data.size
        rdata = data.reshape((n/2,2)) 
        return rdata

    def get_key(self, d, value):
        """
        Return list of keys in dict d that pair with value.
        """
        return [item[0] for item in d.items() if item[1]==value]

    #--------------------------------------------
    #
    # Pickle and YAML I/O
    #
    #---------------------------------------------
    def pickle_graph(self, G, path):
        """
        Write a graph in pickle format to fname.

        G     -- NetworkX graph object
        fname -- file name in str format
        Returns 0
        """
        f = open(path, 'wb')
        pkl.dump(G,f)
        f.close()
        return 0

    def write_pkl(self, obj, path):
        """
        Write an object in pickle format to fname.

        obj     -- object to pickle
        fname -- file name in str format
        Returns 0
        """
        with open(path, 'wb') as fh:
            pkl.dump(obj, fh)

    def load_pickle_graph(self, path):
        """
        Load a previously pickled graph with networkx utils
        """
        try:
            import networkx as NX
        except ImportError:
            raise
        G = NX.read_gpickle(path)
        return G
    
    def load_pickle(self, path):
        """
        Load stored pickle file
        """
        fh = open(path, 'r')
        data = pkl.load(fh)
        return data

    def write_yaml(self, G, path, mode='r'):
        """
        Write object G to file in YAML format
        """
        try:
            import yaml
            fh = open(path, mode=mode)
            yaml.dump(G,fh)
        except ImportError:
            print "Cannot find yaml, using cPickle"
            path = path[:-4]+'pkl'
            fh = open(path, mode=mode)
            pkl.dump(G,fh)
        fh.close()
        return 0

    def load_yaml_graph(self, path):
        """
        Load a yaml file into a graph object
        """
        try:
            import yaml
        except ImportError:
            raise ImportError, "Cannot find yaml package"
        f = open(path, 'rb')
        G = yaml.load(f)
        return G

    def node_count(self, G):
        """
        Return len(exc. nodes), len(inh. nodes)
        """
        enodes = [u for u in G if u[1]=='e']
        inodes = [u for u in G if u[1]=='i']
        return (len(enodes), len(inodes))


#------------------------------------------------
#
# Plotting and array manipulation needed to plot
#
#------------------------------------------------
class Data_Handlers(IOutils, Stats):
    """
    Data handlers. 
    """
    def __init__(self):
        pass

    def fft_freq_shift(self, data, high_freq, spacing=0.005):
        """
        Create frequency bin for fft.  Truncate at 1/2 length (for
        rfft), and slice off initial frequencies, making sure to at
        least start at bin 1 to avoid log(0).

        Default spacing is dt in ODE step.
 
        Returns fftfreq (truncated by 1/2) and shift 
        """
        freq = numpy.fft.fftfreq(len(data),d=spacing)
        len_freq = len(freq)

        # slice off highest frequency (fluctuation 
        freq_trunc = freq[0: len_freq/2]

        # take highest index for truncation
        shift = numpy.max(numpy.where(freq_trunc < high_freq))

        if shift == 0:
            shift = 1  # make sure we don't try to eval log(0)
        
        return freq_trunc, shift


    def run_dfa(self, fname, fodir=None, minbox=None, maxbox=None):
        """
        Call dfa.c. Courtesy CK Peng, et al. at PhysioNet.org
        """
        fin = file(fname,'rb')
        if fname[0:4] == '/tmp' or fname[0:3]=='tmp':
            foname = fname+".dfa"
            fout = file(foname,'wb')
        else:
            foname = fname[0:-4]+".dfa"
            fout = file(foname,'wb')
        
        cmd = ['dfa']
        if minbox is not None:
            cmd.append( '-l' )
            cmd.append( str(int(minbox)) )
        if maxbox is not None:
            cmd.append( '-u' )
            cmd.append( str(int(maxbox)) )
        sp.call( cmd, stdin=fin, stdout=fout)
        fin.close()
        fout.close()
        return foname


    def fft_srate(self, data, spacing, srate=1):
        """
        FFT of data sampled at srate. 

        spacing is the average steps length from 'h' vector of GSL ODE
        solver.
        """
        freq = numpy.fft.fftfreq(len(data), d=spacing)
        fx = freq[0: len(freq)/2]
        fy = abs(numpy.fft.rfft(data))
        return fx, fy[1:]

    #-----------------------------------------
    #
    # Plotting
    #
    #-----------------------------------------
    def plot_errorbars(self, xt, yt, err, fname, **args):
        """
        Assumed to be plotting error in y direction.
        """
        fargs = {'title': "Variance about DFA mean",
                 'state': 'v',
                 'node': 0 }

        s = fargs['title']+'; IC changes in '+fargs['state']+\
            ', node '+str(fargs['node'])

        fig = pylab.figure()
        pylab.title(s)
        pylab.errorbar(xt, yt, err, fmt='.', ecolor='r')
        pylab.xlabel('Log of window size')
        pylab.ylabel('Log of DFA')
        pylab.savefig(fname)
        pylab.close(fig)


    def scatter_plot(self, xt, yt, title='scatter plot',
                     xlabel='',ylabel='', fname='scatter'):
        """
        Scatter plot (xt,yt)
        """
        pylab.figure()
        pylab.title(title)
        pylab.xlabel(xlabel)
        pylab.ylabel(ylabel)
        pylab.scatter(xt,yt)
        pylab.savefig(fname+'.eps')

    def scatter_dfas(self, dfalist, fname):
        """
        Takes as sequence object containing nx2 arrays and scatter
        plots them on one canvas.

        dfalist   : list of dfa vectors to scatter
        fname     : where to save it (full path)
        """
        d=dfalist
        n = len(d)
        fig = pylab.figure()
        
        for i in range(n):
            pylab.scatter(d[i][:,0],d[i][:,1], s=5)
            
        pylab.xlabel("Log of window size")
        pylab.ylabel("Log of DFA")

        if not fname.endswith(".eps"): fname+='.eps'
        pylab.savefig(fname)
        pylab.close(fig)
            

    def plot_datafile(self, fname, sep='\n', **args):
        """
        Simple utility function to read in 2D data and plot it.
        """
        fargs = {'color': 'b',
                 'markers': '',
                 'width': 1.5 }
        for k in args:
            fargs[k] = args[k]
        
        data = numpy.fromfile( fname, sep=sep )
        data.resize( (len(data)/2,2) )
        nx, ny = data.shape
        if nx > ny:
            data = data.T
        ## color and markers for data. Eg., 'bo' plot arg
        linestyle = fargs['color'] + fargs['markers']
        pylab.plot( data[0], data[1], linestyle, lw=fargs['width'] )
        pylab.show()

    def plot_loglog_lin_regress(self, xdata, ydata, 
                                slope, intercept, **args):
                       
        """
        Plot the loglog of the fft of data and the linear regression
        line to the fft.

        keywords:
        
            -- save, should we save the plot 
            
            -- path, required if save=True. expects relative path from
               working directory, eg. 'temp/fft_plots'
            
            -- scatter, plot FFT spectrum as scatter plot (not valid
               for lin.regression plot)
        """
        fargs = {'srate': 1,
                 'save': True,
                 'path': 'temp/',
                 'shift': 1 }
        
        for k in args.keys():
            fargs[k] = args[k]
        
        
        lendata = len(xdata)
        y = numpy.log10(
            intercept + slope*xdata[fargs['shift']::fargs['srate']])

        print "y ", y

        print "ydata ", ydata

        pylab.figure()
        pylab.clf()
        pylab.title("Linear Regression of log-log of FFT")
        pylab.plot(numpy.log10(xdata[fargs['shift']::fargs['srate']]),
                   numpy.log10(numpy.abs(
                    ydata[fargs['shift']::fargs['srate']])),
                   'b-', label="Excit. fft")
        pylab.plot(numpy.log10(xdata[fargs['shift']::fargs['srate']]), y, 
                   'r-', label="lin. reg. line, slope="+str(round(slope,4)))
        pylab.grid(True)
        pylab.legend(loc=3)
        if fargs['save']:
            path = fargs['path']
            pylab.savefig(path+"loglog_linreg."+str(fargs['srate'])+".eps",
                          dpi=80)

    def linreg_plot(self, data, slope, intercept, fname, srate=1, sigma=0,
                    tex=False, dfa_parts=None, **args):
        """
        data   -- 2D dfa data
        """
        fargs = {'g': None,
                 'eps': None,
                 'seps': None,
                 'fdir': 'dir_',
                 'cells': 10 }

        for k in args.keys():
            fargs[k] = args[k]

        params = {'axes.labelsize': 14,
                  'text.fontsize': 14,
                  'legend.fontsize': 10,
                  'font.family': 'serif',
                  'text.usetex': tex }
          #'figure.figsize': fig_size}          'xtick.labelsize': 10,
          #'ytick.labelsize': 10, 
        pylab.rcParams.update(params)

        # create figure handle 
        fig = pylab.figure()
        
        Ptitle = "DFA"
            
#         if srate == 1:
#             pass
#         else:
#             Ptitle += ", srate "+str(srate)
        if sigma == 0:
            pass
        else:
            Ptitle += ", noise="+str(sigma)
        
        pylab.title(Ptitle)
        pylab.plot(data[:,0], data[:,1], linewidth=1.5,
                       label="DFA of "+str(fargs['cells'])+" cell network")


        print dfa_parts
        
        if type(dfa_parts) == list:
            d = dfa_parts
            color = ['g','r','m']
            for i in range(len(d)):
                ## d[i][0]=slope; d[i][1]=intercept, d[i][2]=dfa
                slope = d[i][0]
                intercept = d[i][1]
                dfa = d[i][2]
                y = slope*dfa[:,0] + intercept
                pylab.plot(dfa[:,0], y,
                           color=color[i],
                           linewidth=1.5,
                           label="linear reg. slope="+str(round(slope,6)))

        else:
            y = slope*data[:,0] + intercept
            pylab.plot(data[:,0], y, 
                       label="linear regression, slope="+str(round(slope,8)))
        if fargs['eps'] is not None:
            epstext = r'$\epsilon = $'+str(fargs['eps'])
            pylab.text(1, -1.3, epstext, fontsize=9, ha='left')
        if fargs['seps'] is not None:
            pylab.text(1.5, -1.55, 'start eps: '+str(fargs['seps']))
        if fargs['g'] is not None:
            pylab.text(1.5, -1.8, 'g_syns: '+str(fargs['g']))

        if not tex:
            pylab.xlabel('window size n')
            pylab.ylabel('DFA(n)')
        elif tex:
            pylab.xlabel(r'$\log(w)$')
            pylab.ylabel(r'$F(w)$')

        pylab.legend(loc=2)

        try:
            if fargs['save'] == True:
                pylab.savefig(fname+'.pdf')
                fig.clf()
                pylab.close(fig)
            else:
                pass
        except KeyError:
            pass
        finally:
            # same old default behavior
            pylab.savefig(fname+'.pdf')
            fig.clf()
            pylab.close(fig)
        

    def linreg_crossover_plot(self, data, slope1, intercept1, 
                              slope2, intercept2, 
                              indices1, indices2,
                              fname='/', srate=1, **args):
        """
        data   -- 1D time series
        """
        fargs = {'tex': True,
                 'sigma': 0,
                 'srate': 1,
                 'do_single': False}

        for k in args.keys():
            fargs[k] = args[k]

        sigma = fargs['sigma']
        srate = fargs['srate']

        if fargs['tex']:
            params = {'axes.labelsize': 14,
                      'text.fontsize': 14,
                      'legend.fontsize': 10,
                      'font.family': 'serif',
                      'text.usetex': True }
            pylab.rcParams.update(params)

        dx = data[:,0]
        dy = data[:,1]

        y1 = self.make_line(slope1, dx[indices1], intercept1)
        y2 = self.make_line(slope2, dx[indices2], intercept2)
        dx = data[:,0]
        dy = data[:,1]
        # plot shit
        pylab.figure()
        pylab.clf()
        pylab.title("DFA plot: sample rate "+str(srate)+", sigma "+str(sigma))
        pylab.plot(dx, dy, linewidth=1.5, label="DFA")
        pylab.plot(dx[indices1], y1, 
                   label="linear regression, slope="+str(round(slope1, 4)))
        pylab.plot(dx[indices2], y2, 
                   label="linear regression, slope="+str(round(slope2, 4)))
        pylab.xlabel('Log of window size')
        pylab.ylabel('DFA')
        #pylab.legend(loc=2)
        pylab.legend(loc=0)
        pylab.savefig(fname+'.eps')

    #-----------------------------------
    # Data Analysis tools
    #-----------------------------------
    def linreg_dfa_area(self, dfacurve, lr=None, slope=None, intercept=None):
        """
        Calculate area between dfa curve and linear regression line
        between x1 and x2, where x1 = index where linear reg. line
        goes from above to below dfa curve and x2 is the minimum index
        where the linear reg. line comes back out from under the dfa
        curve.

        Must provide either lr, linear regression line as 2D array, or
        slope and intercept with which to build line from dfacurve
        x-axis. If lr is passed in, it must have been built from the
        dfacurve prior to calling linreg_dfa_area.
        
        Returns   the    area   in   above    region   calculated   as
        sum((i+1 - i)*(DFA[i:i+1] - lr[i:i+1])**2) over indexes in region.
        """
        if lr is None:
            if slope is None or intercept is None:
                raise ValueError(
                    "If no linreg line, must provide slope and intercept!")
            else:
                lr = self.make_line(slope, dfacurve[:,0], intercept)
        else:
            if len(dfacurve[:,0]) != len(lr):
                raise ValueError("len of dfacurve and lr are not equal")

        if self.is_ragged(dfacurve, lr):
            return 0
        else:
            left_end, right_end = self.find_endpts(dfacurve, lr)
        
        if abs(left_end - right_end) == 0:
            # calc area only to the right of crossing
            right_end = len(dfacurve)-1
            # or return 0 ??

        area = self.dfa_linreg_area(dfacurve, lr, left_end, right_end)
        #print "area ", area
        return area

    def dfa_linreg_area(self, dfa, lr, lpt, rpt):
        """
        Find L2 area between dfa and lr between under "hump".
        """
        return numpy.sum((dfa[lpt:rpt,1] - lr[lpt:rpt])**2)

    def find_above_below(self, dfa, lr):
        """
        Find where dfa curve is above/below linreg. line. 

        Returns two sets of indices
        """
        below = numpy.where(dfa[:,1] < lr)[0]
        above = numpy.where(dfa[:,1] > lr)[0]
        return below, above

    def is_ragged(self, dfa, lr):
        """
        If dfa and lr cross many times, flag this so that we can approximate
        area as zero.
        """
        below, above = self.find_above_below(dfa, lr)
        crossings = self.crossings(above, below+1)
        if len(crossings) > 4:
            return True
        else:
            return False

    def crossings(self, above, below):
        return [x for x in above if x in below] 

    def find_endpts(self, dfa, lr):
        """
        """
        below, above = self.find_above_below(dfa, lr)
        crossings = self.crossings(above, below+1)
        if len(crossings) > 4:
            # if many crossings probably ragged
            # should have already passed this above
            return below[0]+1, below[1]+1
        else:
            # hopefully more regular curve
            # return indices where curves cross
            return crossings[0], crossings[-1]

    def find_endpt(self, dfacurve, lr):
        """
        Find last reasonable point to calculate area with before dfa
        curve becomes jumpy.
        """

        dmax = dfacurve[-1,0]

        if dmax > 3.5:
            trunc = numpy.where(dfacurve[:,0] < 3.5)[0]
            x2 = trunc[-1]
        else:
            trunc = numpy.where(dfacurve[:,0] < 3.0)[0]
            x2 = trunc[-1]

        return x2

    def crossover_point(self, data, cr=None, shift=None):
        """
        Simplistic way to set crossover point of dfa data, which is
        assumed to be a 2d array of (box size n)x(DFA(n)).

        Returns indices data[:,0] to the left and right of cr.

        TODO - Find a true crossover point
        """
        rx = data[:,0]
        ry = data[:,1]

        print 'data ', data
        print rx
        print ry

        if cr == None:
            w = abs(rx[0] - rx[-1])/2 + shift
        else: 
            # if we have a specific crossover pt
            w = cr
        indl = numpy.where(rx < rx[0] + w)
        indr = numpy.where(rx >= rx[0] + w)
        return indl, indr

    def plot_data(self, tdata, ndata, grid=False):
        """
        Basic plot with no labels
        """
        pylab.clf()
        if grid == True:
            pylab.grid(True)
        pylab.plot(tdata,ndata)

    def draw_graph(self, G, layout="spring"):
        """
        Draw graph, coloring Ex and In nodes different colors.
        """
        ex_nodes = [u for u in G.nodes() if u[1]=="e"]
        in_nodes = [u for u in G.nodes() if u[1]=="i"]
        if layout != "spring":
            if layout == "circular":
                pos = NX.circular_layout(G)
            else:
                print "Unknown layout. Using spring_layout"
                layout = "spring"
        elif layout == "spring":
            pos = NX.spring_layout(G)
        NX.draw_networkx_nodes(G,pos,nodelist=ex_nodes,node_color='r')
        NX.draw_networkx_nodes(G,pos,nodelist=in_nodes,node_color='b')
        NX.draw_networkx_edges(G,pos)
        n = len(G.nodes())
        pylab.savefig("neural_network_"+str(n)+".png")
        #pylab.show()

    def state_avg(self, data, **args):
        """ 
        data -- state array with shape: (states , num_eqns)

        num_cells -- assumed (for now) that these are excitatory cells.

        Returns average voltage of excitatory cells (more code needed
        for inh.)x
        """
        fargs = {'num_cells': 10,
                 'num_exc': 5 }

        for k in args.keys():
            fargs[k] = args[k]

        eqns = 4*fargs['num_cells']
        N = len(data) / (eqns)
        
        rdata = data.reshape(N, eqns)

        return self.avg_cells(rdata)

    def avg_cells(self, data, **args):
        """
        data is state array of shape (steps, num_eqns)
        """
        fargs = {'num_exc': 5}

        edata = data[:,0:fargs['num_exc']]
        avg_vec = (edata.sum(axis=1))/(fargs['num_exc'])
        avg = numpy.average(avg_vec)
        eavg = avg_vec - avg
        return eavg

    def make_line(self, m, x, b):
        """ x is an array of indices. m=slope, b=intercept"""
        return m * x + b

    def dfa_relative_linreg(self, dy, ly):
        """
        Indices where dfa is above/below lin. regression line. This
        assumes that dfa cruve is not too "wiggly". Not sure how to
        handle wiggly case....
        """
        dfa_above = numpy.where(dy > ly)[0]
        dfa_below = numpy.where(dy <= ly)[0]

        if 0 in dfa_below:
            # then dfa curve starts below linreg
            # want indices where dy *above* ly
            return dfa_above
        elif 0 in dfa_above:
            # then dfa curve starts above linreg
            # want indices where dy *below* ly
            return dfa_below

    def find_hump(self, indices):
        """
        Find the first place the linear regression line crosses the
        DFA
        """
        diff = numpy.absolute(indices[1:]-indices[:-1])
        stop = numpy.where(diff != 1)[0]

        if len(stop) == 0:
            return indices
        elif stop[0] < float(len(indices))/3:
            if len(stop) > 1:
                return indices[0:stop[1]]
            else: 
                return indices[0:stop[0]]
        else:
            return indices[0:stop[0]]

    def find_right_endpt(self, noise, srate):
        """
        Grab the right number from the dictionary.

        This is slow with the file handling....

        return 
        """
        fh = open("scutoff.pkl", 'rb')
        d = pkl.load(fh)
        fh.close()

        # hack for sigma 0-45
        noise = 2*noise

        cutoffs = d[int(noise)]
        w = numpy.where(cutoffs == float(srate))[0]
        return cutoffs[w[0],1]

    def find_linreg_perp_intersect(self, pt, m, bl):
        """
        pt  -- point on DFA curve
        m   -- slope of linear regression line
        bl  -- intercept of linear regression line

        Given point pt, find x so that for the two lines below
        intersect:

        [y - y2 = -(1/m)(x - x2)] == [y - y1 = m(x - x1)]

        where x2=pt[0], y2=pt[1] are on the DFA curve.
        """
        x = (1 / (m + 1/m))*(pt[0]/m + pt[1] - bl)
        y = m * x + bl
        z = (x,y)
        return self.l2norm((x-pt[0], y-pt[1]))

    def l2norm(self, v):
        return pylab.l2norm(v)

    def find_max_dist_index(self, dx, dy, m, intercept):
        """
        Find index where dy and ly curves are max distance apart.
        """
        lindist = []
        nx = len(dx)
        ny = len(dy)
        if ny != ny:
            # don't know when this would ever happen, but...
            raise ValueError, "Array lengths don't match!"
        for i in xrange(nx):
            pt = [dx[i], dy[i]]
            lindist.append(self.find_linreg_perp_intersect(pt, m, intercept))
        linarr = numpy.asarray(lindist)
#        dist = [self.l2norm(v) for v in linarr]
        try:
            M = max(linarr)
        except:
            print "dist ", linarr, " dfa ", dx, dy
            return 0
            
        w = numpy.where(linarr == M)[0]
        return w[0]

    def timeseries2dfa(self, svec=None, **args):
        """
        Typical usage: pass in vector svec and either find the average
        of excitatory cells or pass it along to dfa analysis if svec
        is already averaged.

        If svec is not passed a value, then *must* pass timeseries2dfa
        the key/value pair path='path/to/state_array'

        Returns 2d dfa vector 
        """
        fargs = {'num_cells': 5,
                 'noise': -1,
                 'do_avg': True,
                 'tex': True,
                 'srate': 1,
                 'do_single': False,
                 'dfa': None,
                 'path': 'temp/'}

        for k in args.keys():
            fargs[k] = args[k]

        fpath = fargs['path']
        srate = fargs['srate']

        if svec != None:
            # if svec is the full state array, we must convert it to
            # an average of excitatory cells first
            if fargs['do_avg'] and svec != None:
                eavg = self.state_avg(svec)
            else:
                eavg = svec
            fname = self.write_data(eavg[::srate], 
                                    fpath + 'excAvg_srate'+str(srate)+\
                                        '.dat')
        else:
            # exc. avg file already exists
            #fname = fpath + 'excAvg_srate'+str(srate)+'.dat'
            s = self.read_state_array(fpath)
            eavg = self.state_avg(s, fargs['num_cells'])
            fname = self.write_data(eavg[::srate], 
                                    fpath + 'excAvg_srate'+str(srate)+\
                                        '.dat')
        print 'fname ', fname

        dfa_fname = self.run_dfa(fname)
        return self.dfa_data_to_array(dfa_fname)

    def calc_crossover_at_srate(self, fpath, srate, svec=None, **args):
        """
        1) read in state.out array which resides at fpath/state.out
        2) calculate excitatory average (time series)
        3) write 2) to file *at given srate* 
        4) run DFA function on file from 3)
        4') read in array written in 4)
        5) calculate single linear regression line to data returned in 4)
        6) determine where likely crossover point is on DFA curve
        7) calculate two linear regression curves from result in 6),
        return slopes m1, m2
        8) record m1-m2 in array:

        max noise    ...    max noise/max noise
        .                          .
        .                          .
        .                          .
        min noise(0) ...    max sample rate/min noise


        args:  num_cells  -- 5 excitatory cells for now
               noise      -- for cutting off extraneous data at end of data
                             (see find_right_endpt())
               do_avg     -- True for neural network,
                             False for single time series (eg. annie data)

        Note: Assumed that averaging over exc. cells. state_avg() is
        inadequate to handle other options yet.

        return list of slopes and intercepts:

        [m1-m2, m1, m2, m (single lin. reg slope),
        intercept1, intercept2, intercept (for m), 
        dfa_data]
        """
        fargs = {'num_cells': 5,
                 'noise': -1,
                 'do_avg': True,
                 'tex': True,
                 'srate': 1,
                 'do_single': False,
                 'dfa_data': None,
                 'plot': False}

        for k in args.keys():
            fargs[k] = args[k]

        if fargs['dfa_data'] == None:
            rdata = self.timeseries2dfa(svec=svec, fpath=fpath)
        else:
            rdata = fargs['dfa_data']
            
        # DFA data
        x = rdata[:,0] # x values for both dfa and linreg line(s)
        y = rdata[:,1]

        if fargs['noise'] != -1:
            right_endpt = self.find_right_endpt(fargs['noise'], 
                                                srate)
            rt = numpy.where(x <= right_endpt)[0]
            rtstop = rt[-1]
        else:
            rtstop = len(x)

        # single linear reg. to dfa curve
        linreg = self.linregress(x[0:rtstop], y[0:rtstop])

        # store in output list
        slope = linreg[0]
        intercept = linreg[1]
        liny = self.make_line(slope, x, intercept)
        
        indices = self.dfa_relative_linreg(y, liny)

        # keep just the indices around the crossover
        indtrunc = self.find_hump(indices)
        indices = indtrunc

        # find index where dfa curve and lin. reg. curve are greatest
        # distance apart (approx.)
        mdindex = self.find_max_dist_index(x[indices], 
                                           y[indices], 
                                           slope,
                                           intercept)
        crossover_index = indices[0] + mdindex

        if fargs['noise'] != -1:
            right_endpt = self.find_right_endpt(fargs['noise'], 
                                                srate)
            rt = numpy.where(x <= right_endpt)[0]
            rtstop = rt[-1]
        else:
            rtstop = len(x)

        # left and right sets of indices
        indl = numpy.arange(0, crossover_index)
        indr = numpy.arange(crossover_index, rtstop)

        slope1, intercept1, r1, sterrest1 = \
            self.linregress(x[indl],y[indl])
        slope2, intercept2, r2, sterrest2 = \
            self.linregress(x[indr],y[indr])
        
        #print "slopes ", slope1, slope2, " with srate ", srate

        # plot stuff
        if fargs['plot']:
            plotname = fpath + "dfa_srate"+str(srate)+"_cutoff.pdf"
            pylab.figure()
            pylab.title("DFA curve with linear regression line")
            pylab.plot(x,y, label="DFA curve")
            pylab.plot(x,liny, label="slope "+str(round(slope,4)))
            pylab.xlabel("Log of window size")
            pylab.ylabel("log DFA")
            pylab.legend(loc=2)
            pylab.savefig(plotname)
    
            # linear regression lines to DFA curve
            if fargs['do_single']:
                crossplot = fpath +\
                    "DFA_co.single_srate"+str(srate)+"."+str(fargs['noise'])
            else:
                crossplot = fpath +\
                    "DFA_co_srate"+str(srate)+"."+str(fargs['noise'])
            self.linreg_crossover_plot(rdata, 
                                       slope1, intercept1, 
                                       slope2, intercept2, 
                                       indl, indr,
                                       fname=crossplot, 
                                       srate=srate, 
                                       tex=fargs['tex'])

        return (slope1 - slope2, slope1, slope2, slope, 
                intercept1, intercept2, intercept, rdata)

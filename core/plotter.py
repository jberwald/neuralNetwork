#!/usr/bin python
# ----------------------------------------------------------------------
# 
#   plotter.py
# 
#   Jesse Berwald
# 
#   Version:    0.1  
# 
#   Opened:     Mar 19, 2009
# 
# ----------------------------------------------------------------------
"""
Class(es) to be inherited by wrappers in other modules. We create N
matplotlib figure objects. These are saved then cleared as needed for
better memory management.
"""

class Plot:
    """
    """
    def __init__(self, nfigs=10, gui=False):
        if gui:
            import pylab
        else:
            import matplotlib
            matplotlib.use('Agg') # non-gui plotting 
            import pylab

    def clfigs(self):
        """ Clear all figs """
        for fig in figlist:
            fig.clf()

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
                 'fdir': 'dir_'}

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
            
        if srate == 1:
            pass
        else:
            Ptitle += ", srate "+str(srate)
        if sigma == 0:
            pass
        else:
            Ptitle += ", noise="+str(sigma)
        
        pylab.title(Ptitle)
        pylab.plot(data[:,0], data[:,1], linewidth=1.5,
                       label="DFA")
        if type(dfa_parts) == list:
            d = dfa_parts
            color = ['g','r','m']
            for i in range(len(d)):
                # d[i][0]=slope; d[i][1]=intercept, d[i][2]=dfa
                slope = d[i][0]
                intercept = d[i][1]
                dfa = d[i][2]
                y = slope*dfa[:,0] + intercept
                pylab.plot(dfa[:,0], y, color=color[i],linewidth=1.5,label="linear reg. slope="+str(round(slope,6)))

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
            pylab.xlabel('Log of window size')
            pylab.ylabel('DFA')

        pylab.legend(loc=2)
        
        if 'save' in fargs.keys():
            if fargs['save'] == True:
                pylab.savefig(fname+'.eps', dpi=80)
                fig.clf()
                pylab.close(fig)
            else:
                pass
        else:
            # same old default behavior
            pylab.savefig(fname+'.eps', dpi=80)
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

    def plot_sampen(self, arr, norm=None, outfile='sampenplot.eps', **args):
        """
        arr  : array extracted from sample entropy file.
        
        """
        fargs = {'show': False,
                 'tex': True }

        for k in args.keys():
            fargs[k] = args[k]

        if fargs['tex']:
            params = {'axes.labelsize': 14,
                      'text.fontsize': 14,
                      'legend.fontsize': 10,
                      'font.family': 'serif',
                      'text.usetex': True }
            pylab.rcParams.update(params)
        
        if not outfile.endswith('.eps') or not outfile.endswith('.pdf'):
            outfile+='.pdf'

        ptitle = 'Sample entropy'
        xs = r'$m$'
        ys = r'$SE(m)$'

        fig = pylab.figure()
        pylab.title(ptitle)
        nx = pylab.arange(1,len(arr)+1)
        #pylab.fill(nx, arr, 'b', alpha=0.5)
        pylab.plot(nx, arr, 'bo-')
        pylab.xlabel(xs)
        pylab.ylabel(ys)
        pylab.grid(True)
        if norm is not None:
            n = str(round(norm,4))
            t = r'$||SE|| = $'+n
            ycoord = (arr[0]+arr[1])/2
            pylab.text(4,ycoord, t, fontsize=16)
        if fargs['show']:
            pylab.show()
        fig.savefig(outfile)
        pylab.close(fig)

    def plot_mse(self, arr, window, r, outfile='mse.eps', norm=None, **args):
        """
        Take extracted values from mse computations and plot
        them. (See mse -h for explanation of options.)

        window:      m value 
        r:           tolerance (r value)
        """
        fargs = {'show': False,
                 'tex': True,
                 'keep_up': False }

        for k in args.keys():
            fargs[k] = args[k]

        if fargs['tex']:
            params = {'axes.labelsize': 14,
                      'text.fontsize': 14,
                      'legend.fontsize': 10,
                      'font.family': 'serif',
                      'text.usetex': True }
            pylab.rcParams.update(params)
        
        if not outfile.endswith('.eps') or not outfile.endswith('.pdf'):
            outfile+='.pdf'

        ptitle = 'Multiscale entropy'
        xs = r'$m$'
        ys = r'$MSE(m)$'

        print 'outfile ', outfile

        fig = pylab.figure()
        pylab.title(ptitle)
        #nx = pylab.arange(1,len(arr)+1)
        #pylab.fill(nx, arr, 'b', alpha=0.5)
        #pylab.plot(nx, arr, 'bo-')
        pylab.fill_between(arr[:,0], arr[:,1])
        pylab.xlabel(xs)
        pylab.ylabel(ys)
        pylab.grid(True)
        if norm is not None:
            n = str(round(norm,4))
            t = r'$||MSE|| = $'+n
            ycoord = (arr[0,1]+arr[1,1])/2
            pylab.text(4,ycoord, t, fontsize=16)
        if fargs['show']:
            pylab.show()
        fig.savefig(outfile)
        if not fargs['keep_up']:
            pylab.close(fig)

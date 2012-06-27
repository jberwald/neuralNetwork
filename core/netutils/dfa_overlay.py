import pylab, os
from neural_network.netutils import utils

D = utils.Data_Handlers()

def linreg_line(data):
    """
    data  --  2d array
    """
    return D.linregress(data[:,0],data[:,1])

def linreg_plot(data1, data2, fname, srate=1, sigma=0,
                tex=True, **args):
    """
    data   -- 1D time series
    """
    fargs = {'g': None,
             'eps': None,
             'seps': None,
             'fdir': 'dir_'}

    for k in args.keys():
        fargs[k] = args[k]


    if len(data1.shape) == 1:
        data1 = data1.reshape((len(data1)/2,2))
    if len(data2.shape) == 1:
        data2 = data2.reshape((len(data2)/2,2))

    m1, intercept1, r1, err1 = linreg_line(data1)
    m2, intercept2, r2, err2 = linreg_line(data2)
    

    params = {'axes.labelsize': 14,
              'text.fontsize': 14,
              'legend.fontsize': 10,
              'font.family': 'serif',
              'text.usetex': tex }
      #'figure.figsize': fig_size}          'xtick.labelsize': 10,
      #'ytick.labelsize': 10, 
    pylab.rcParams.update(params)

    y1 = m1*data1[:,0] + intercept1
    y2 = m2*data1[:,0] + intercept2

    # plot shit
    pylab.figure()
    pylab.clf()

    Ptitle = "DFA plot"

    if srate == 1:
        pass
    else:
        Ptitle += ", srate "+str(srate)
    if sigma == 0:
        pass
    else:
        Ptitle += ", noise="+str(sigma)

    pylab.title(Ptitle)

    #dfa curves
    pylab.plot(data1[:,0], data1[:,1], linewidth=1.2, label="DFA1")
    pylab.plot(data1[:,0], data2[:,1], linewidth=1.2, label="DFA2")
               
    # linreg lines
    pylab.plot(data1[:,0], y1, 
               label="linear regression (1), slope="+str(round(m1,8)))
    pylab.plot(data1[:,0], y2,
               label="linear regression (2), slope="+str(round(m2,8)))
    if fargs['eps'] is not None:
        pylab.text(1.5, -1.3, 'eps       :'+str(fargs['eps']))
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
    pylab.show()
 #   pylab.savefig(fname+'.'+str(round(slope,8))+'.eps', dpi=80)

def plot_noisy_overlay(plist, fdir):
    """
    Create plot with objects from plist.

    plist  : list of tuples, (dfa array, noise level)
    """
    params = {'axes.labelsize': 14,
              'text.fontsize': 14,
              'legend.fontsize': 10,
              'font.family': 'serif',
              'text.usetex': True }
    pylab.rcParams.update(params)

    fig = pylab.figure()
    ax = fig.gca()

    ax.set_title('Overlay of DFA plots for noisy timeseries')

    for p in plist:
        ax.plot(p[0][:,0], p[0][:,1], label='sigma='+str(p[1]))
    
    ax.set_xlabel('log of window size')
    ax.set_ylabel('DFA')
    pylab.legend(loc=2, shadow=True)
    pylab.savefig(fdir+'dfa_noisy_overlay.eps')
    pylab.show()

def noisy_dfa_vecs(fdir):
    """
    Return list of tuples of dfa vectors with associated sigma (noise
    level).
    """
    if not fdir.endswith('noise/'): 
        raise ValueError("Probably not pulling from noise directory")
    if not fdir.endswith('/'): fdir += '/'

    dlist = os.listdir(fdir)
    dlist.sort()
    dfalist = []
    
    for d in dlist:
        if d.endswith('.dfa'):
            v = pylab.fromfile(fdir+d, sep='\n')
            sigma = extract_noise_level(d)
            dfalist.append((v.reshape((len(v)/2,2)), sigma))
    return dfalist
    
def extract_noise_level(fname):
    """
    Extract nX.Y in 'avgvec.nX.Y.dfa'
    """
    # extract filename from full path
    slash = fname.rfind('/')
    avgvec = fname[slash+1:]
    s = avgvec.split('.')

    # strip 'n'; if sigma=1, don't add '.dfa'
    if s[2] == 'dfa':
        sigma = '.'.join((s[1][1:],'0'))
    else:
        sigma = '.'.join((s[1][1:],s[2]))

    return sigma

def noisy_overlay(fdir, avgvec='avgvolt.sr.1.v.0.dfa', **args):
    """
    Overlay on one plot a sequence of noisy dfa plots.

    fdir  : top level of simulation directory (eg. .../n10g0e1c1s0/)
    """
    noisedir = fdir+'noise/'
    ndvecs = noisy_dfa_vecs(noisedir)

    av = pylab.fromfile(fdir+avgvec,sep='\n')
    avdfa = av.reshape((len(av)/2,2))
    
    ndvecs.append((avdfa, '0.0'))

    plot_noisy_overlay(ndvecs, fdir)

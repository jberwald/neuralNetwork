"""
Make subplot with two figures. One dfa the other finite differences of
that dfa curve.
"""
# import matplotlib
# matplotlib.use('TKAgg')
from matplotlib.font_manager import FontProperties
from pylab import *
from neural_network.netutils import fileio
from neural_network.core import fitfunc, timeseries

TS = timeseries.Timeseries(None,None)
FP = FontProperties()

def make_plots(numfiles):

    n = numfiles

    # read in dfa data
    dfa = [fromfile('avgvolt.sr.1.0.'+str(i)+'.dfa',sep='\n') for i in range(n)]

    #reshape dfa data
    for i in range(n):
        dfa[i] = dfa[i].reshape((len(dfa[i])/2,2))

    # save dir
    if not fileio.check_dir('fin_diff'):
        fileio.mkdir('fin_diff')

    # plot data
    fig = figure()
    subplots_adjust(hspace=0.4, wspace=0.4)
    FP.set_size('xx-small')
    for i in range(0,n,2):
        clf()
        
        print "plot ", i

        suptitle("Comparison of Finite Differences and DFA curves - simulations "+str(i)+" and "+str(i+1))

        TS.dfa_partition = []

        TS.dfa = d = dfa[i]
        TS.linregress(xmin=1.5, xmax=3.3)

        subplot(221)
#        d = dfa[i]
        title('Finite difference of DFA', fontsize='small')
        xlabel('log of window size', fontsize='x-small')
        ylabel('Finite difference', fontsize='x-small')
        fd = fitfunc.finite_diff(d)
        plot(d[:-1,0],fd)

        subplot(222)
        dpart = TS.dfa_partition
        color = ['g','r','m']
        title('DFA with linear regressions', fontsize='small')
        xlabel('log of window size', fontsize='x-small')
        ylabel('DFA', fontsize='x-small')
        plot(d[:,0], d[:,1], linewidth=1.5,
                   label="DFA")
        for j in range(len(dpart)):
            # d[i][0]=slope; d[i][1]=intercept, d[i][2]=dfa
            slope = dpart[j][0]
            intercept = dpart[j][1]
            dx = dpart[j][2]
            y = slope*dx[:,0] + intercept
            plot(dx[:,0], y, color=color[j],linewidth=1,label="linear reg. slope="+str(round(slope,6)))
        legend(prop=FP,shadow=True, loc=4)

        #next one
        TS.dfa_partition = []

        TS.dfa = d = dfa[i+1]
        TS.linregress(xmin=1.5, xmax=3.3)

        subplot(223)
        #title('Finite difference of DFA')
        xlabel('log of window size', fontsize='x-small')
        ylabel('Finite difference', fontsize='x-small')
        fd = fitfunc.finite_diff(d)
        plot(d[:-1,0],fd)

        subplot(224)
        dpart = TS.dfa_partition
        color = ['g','r','m']

        plot(d[:,0], d[:,1], linewidth=1.5,
                   label="DFA")
        for j in range(len(dpart)):
            # d[i][0]=slope; d[i][1]=intercept, d[i][2]=dfa
            slope = dpart[j][0]
            intercept = dpart[j][1]
            dx = dpart[j][2]
            y = slope*dx[:,0] + intercept
            plot(dx[:,0], y, color=color[j],linewidth=1,label="linear reg. slope="+str(round(slope,6)))
        legend(prop=FP,shadow=True, loc=4)
            
        savefig('fin_diff/fin_diff_dfa.'+str(i)+'_'+str(i+1)+'.eps')


def plot_dfa_fin_diff(tsobj, num=0):
    """
    Plot a single instance of comparison. 
    """
    fig = figure()
    subplots_adjust(hspace=0.4, wspace=0.4)
    FP.set_size('xx-small')

    suptitle("Comparison of Finite Difference and DFA curve")

    d = tsobj.dfa

    subplot(221)
    #d = dfa[i]
    title('Finite difference of DFA', fontsize='small')
    xlabel('log of window size', fontsize='x-small')
    ylabel('Finite difference', fontsize='x-small')
    fd = fitfunc.finite_diff(d)
    plot(d[:-1,0],fd)

    subplot(222)
    dpart = tsobj.dfa_partition

    color = ['g','r','m']
    title('DFA with linear regressions', fontsize='small')
    xlabel('log of window size', fontsize='x-small')
    ylabel('DFA', fontsize='x-small')
    plot(d[:,0], d[:,1], linewidth=1.5,
               label="DFA")
    for j in range(len(dpart)):
        # d[i][0]=slope; d[i][1]=intercept, d[i][2]=dfa
        slope = dpart[j][0]
        intercept = dpart[j][1]
        dx = dpart[j][2]
        y = slope*dx[:,0] + intercept
        plot(dx[:,0], y, color=color[j],linewidth=1,label="linear reg. slope="+str(round(slope,6)))
    legend(prop=FP,shadow=True, loc=4)

    savefig(tsobj.testdir+'findiff_dfa.'+str(num)+'.eps')
    close(fig)

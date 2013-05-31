import numpy, pylab
from neural_network.netutils import utils 

D = utils.Data_Handlers()

def plot(fname, dir):

    dfavec = numpy.fromfile(fname, sep=' ')

    rd = len(dfavec)
    rdfa = dfavec.reshape((rd/2,2))

    slope, intercept, r, sterr = D.linregress(rdfa[:,0], rdfa[:,1])

    dfaname = dir + 'dfa_fitness.eps'

    D.linreg_plot(rdfa, slope, intercept, dfaname, srate=10, sigma=0.0)
    
    pylab.show()
    

